!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Partition F and homogenize P/dPdF
!--------------------------------------------------------------------------------------------------
submodule(homogenization) mechanical


  interface

    module subroutine mechanical_pass_init
    end subroutine mechanical_pass_init

    module subroutine mechanical_isostrain_init
    end subroutine mechanical_isostrain_init

    module subroutine mechanical_RGC_init(num_homogMech)
      class(tNode), pointer, intent(in) :: &
        num_homogMech                                                                               !< pointer to mechanical homogenization numerics data
    end subroutine mechanical_RGC_init


    module subroutine mechanical_isostrain_partitionDeformation(F,avgF)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
    end subroutine mechanical_isostrain_partitionDeformation

    module subroutine mechanical_RGC_partitionDeformation(F,avgF,instance,of)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
      integer,                          intent(in)  :: &
        instance, &
        of
    end subroutine mechanical_RGC_partitionDeformation


    module subroutine mechanical_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
      real(pReal),   dimension (3,3),       intent(out) :: avgP                                     !< average stress at material point
      real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                               !< average stiffness at material point

      real(pReal),   dimension (:,:,:),     intent(in)  :: P                                        !< partitioned stresses
      real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                     !< partitioned stiffnesses
      integer,                              intent(in)  :: instance
    end subroutine mechanical_isostrain_averageStressAndItsTangent

    module subroutine mechanical_RGC_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,instance)
      real(pReal),   dimension (3,3),       intent(out) :: avgP                                     !< average stress at material point
      real(pReal),   dimension (3,3,3,3),   intent(out) :: dAvgPdAvgF                               !< average stiffness at material point

      real(pReal),   dimension (:,:,:),     intent(in)  :: P                                        !< partitioned stresses
      real(pReal),   dimension (:,:,:,:,:), intent(in)  :: dPdF                                     !< partitioned stiffnesses
      integer,                              intent(in)  :: instance
    end subroutine mechanical_RGC_averageStressAndItsTangent


    module function mechanical_RGC_updateState(P,F,avgF,dt,dPdF,ip,el) result(doneAndHappy)
      logical, dimension(2) :: doneAndHappy
      real(pReal), dimension(:,:,:),     intent(in)    :: &
        P,&                                                                                         !< partitioned stresses
        F                                                                                           !< partitioned deformation gradients
      real(pReal), dimension(:,:,:,:,:), intent(in) :: dPdF                                         !< partitioned stiffnesses
      real(pReal), dimension(3,3),       intent(in) :: avgF                                         !< average F
      real(pReal),                       intent(in) :: dt                                           !< time increment
      integer,                           intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
    end function mechanical_RGC_updateState


    module subroutine mechanical_RGC_results(instance,group)
      integer,          intent(in) :: instance                                                      !< homogenization instance
      character(len=*), intent(in) :: group                                                         !< group name in HDF5 file
    end subroutine mechanical_RGC_results

  end interface

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_init(num_homog)

  class(tNode), pointer, intent(in) :: &
    num_homog

  class(tNode), pointer :: &
    num_homogMech

  print'(/,a)', ' <<<+-  homogenization:mechanical init  -+>>>'

  allocate(homogenization_dPdF(3,3,3,3,discretization_nIPs*discretization_Nelems), source=0.0_pReal)
  homogenization_F0 = spread(math_I3,3,discretization_nIPs*discretization_Nelems)                   ! initialize to identity
  homogenization_F = homogenization_F0                                                              ! initialize to identity
  allocate(homogenization_P(3,3,discretization_nIPs*discretization_Nelems),        source=0.0_pReal)

  num_homogMech => num_homog%get('mech',defaultVal=emptyDict)
  if (any(homogenization_type == HOMOGENIZATION_NONE_ID))      call mechanical_pass_init
  if (any(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)) call mechanical_isostrain_init
  if (any(homogenization_type == HOMOGENIZATION_RGC_ID))       call mechanical_RGC_init(num_homogMech)

end subroutine mechanical_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition F onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_partition(subF,ip,el)

  real(pReal), intent(in), dimension(3,3) :: &
    subF
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point
    el                                                                                              !< element number

  integer :: co
  real(pReal), dimension (3,3,homogenization_Nconstituents(material_homogenizationAt(el))) :: Fs


  chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))

    case (HOMOGENIZATION_NONE_ID) chosenHomogenization
      Fs(1:3,1:3,1) = subF

    case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
      call mechanical_isostrain_partitionDeformation(Fs,subF)

    case (HOMOGENIZATION_RGC_ID) chosenHomogenization
      call mechanical_RGC_partitionDeformation(Fs,subF,ip,el)

  end select chosenHomogenization

  do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
    call phase_mechanical_setF(Fs(1:3,1:3,co),co,ip,el)
  enddo


end subroutine mechanical_partition


!--------------------------------------------------------------------------------------------------
!> @brief Average P and dPdF from the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_homogenize(dt,ip,el)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
       ip, &                                                                                        !< integration point
       el                                                                                           !< element number

  integer :: co,ce
  real(pReal) :: dPdFs(3,3,3,3,homogenization_Nconstituents(material_homogenizationAt(el)))
  real(pReal) :: Ps(3,3,homogenization_Nconstituents(material_homogenizationAt(el)))


  ce = (el-1)* discretization_nIPs + ip
  chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))

    case (HOMOGENIZATION_NONE_ID) chosenHomogenization
        homogenization_P(1:3,1:3,ce)            = phase_mechanical_getP(1,ip,el)
        homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = phase_mechanical_dPdF(dt,1,ip,el)

    case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
      do co = 1, homogenization_Nconstituents(material_homogenizationAt(el))
        dPdFs(:,:,:,:,co) = phase_mechanical_dPdF(dt,co,ip,el)
        Ps(:,:,co) = phase_mechanical_getP(co,ip,el)
      enddo
      call mechanical_isostrain_averageStressAndItsTangent(&
        homogenization_P(1:3,1:3,ce), &
        homogenization_dPdF(1:3,1:3,1:3,1:3,ce),&
        Ps,dPdFs, &
        homogenization_typeInstance(material_homogenizationAt(el)))

    case (HOMOGENIZATION_RGC_ID) chosenHomogenization
      do co = 1, homogenization_Nconstituents(material_homogenizationAt(el))
        dPdFs(:,:,:,:,co) = phase_mechanical_dPdF(dt,co,ip,el)
        Ps(:,:,co) = phase_mechanical_getP(co,ip,el)
      enddo
      call mechanical_RGC_averageStressAndItsTangent(&
        homogenization_P(1:3,1:3,ce), &
        homogenization_dPdF(1:3,1:3,1:3,1:3,ce),&
        Ps,dPdFs, &
        homogenization_typeInstance(material_homogenizationAt(el)))

  end select chosenHomogenization

end subroutine mechanical_homogenize


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
module function mechanical_updateState(subdt,subF,ip,el) result(doneAndHappy)

  real(pReal), intent(in) :: &
    subdt                                                                                           !< current time step
  real(pReal), intent(in), dimension(3,3) :: &
    subF
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point
    el                                                                                              !< element number
  logical, dimension(2) :: doneAndHappy

  integer :: co
  real(pReal) :: dPdFs(3,3,3,3,homogenization_Nconstituents(material_homogenizationAt(el)))
  real(pReal) :: Fs(3,3,homogenization_Nconstituents(material_homogenizationAt(el)))
  real(pReal) :: Ps(3,3,homogenization_Nconstituents(material_homogenizationAt(el)))


  if (homogenization_type(material_homogenizationAt(el)) == HOMOGENIZATION_RGC_ID) then
      do co = 1, homogenization_Nconstituents(material_homogenizationAt(el))
        dPdFs(:,:,:,:,co) = phase_mechanical_dPdF(subdt,co,ip,el)
        Fs(:,:,co)        = phase_mechanical_getF(co,ip,el)
        Ps(:,:,co)        = phase_mechanical_getP(co,ip,el)
      enddo
      doneAndHappy = mechanical_RGC_updateState(Ps,Fs,subF,subdt,dPdFs,ip,el)
  else
    doneAndHappy = .true.
  endif

end function mechanical_updateState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to file.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_results(group_base,h)

  character(len=*), intent(in) :: group_base
  integer, intent(in)          :: h

  character(len=:), allocatable :: group

  group = trim(group_base)//'/mech'
  call results_closeGroup(results_addGroup(group))

  select case(homogenization_type(h))

    case(HOMOGENIZATION_rgc_ID)
      call mechanical_RGC_results(homogenization_typeInstance(h),group)

  end select

  !temp = reshape(homogenization_F,[3,3,discretization_nIPs*discretization_Nelems])
  !call results_writeDataset(group,temp,'F',&
  !                          'deformation gradient','1')
  !temp = reshape(homogenization_P,[3,3,discretization_nIPs*discretization_Nelems])
  !call results_writeDataset(group,temp,'P',&
  !                          '1st Piola-Kirchhoff stress','Pa')

end subroutine mechanical_results


end submodule mechanical
