!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Partition F and homogenize P/dPdF
!--------------------------------------------------------------------------------------------------
submodule(homogenization) mechanical


  interface

    module subroutine pass_init
    end subroutine pass_init

    module subroutine isostrain_init
    end subroutine isostrain_init

    module subroutine RGC_init(num_homogMech)
      class(tNode), pointer, intent(in) :: &
        num_homogMech                                                                               !< pointer to mechanical homogenization numerics data
    end subroutine RGC_init


    module subroutine isostrain_partitionDeformation(F,avgF)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
    end subroutine isostrain_partitionDeformation

    module subroutine RGC_partitionDeformation(F,avgF,ce)
      real(pReal),   dimension (:,:,:), intent(out) :: F                                            !< partitioned deformation gradient
      real(pReal),   dimension (3,3),   intent(in)  :: avgF                                         !< average deformation gradient at material point
      integer,                          intent(in)  :: &
        ce
    end subroutine RGC_partitionDeformation


    module function RGC_updateState(P,F,avgF,dt,dPdF,ce) result(doneAndHappy)
      logical, dimension(2) :: doneAndHappy
      real(pReal), dimension(:,:,:),     intent(in)    :: &
        P,&                                                                                         !< partitioned stresses
        F                                                                                           !< partitioned deformation gradients
      real(pReal), dimension(:,:,:,:,:), intent(in) :: dPdF                                         !< partitioned stiffnesses
      real(pReal), dimension(3,3),       intent(in) :: avgF                                         !< average F
      real(pReal),                       intent(in) :: dt                                           !< time increment
      integer,                           intent(in) :: &
        ce                                                                                          !< cell
    end function RGC_updateState


    module subroutine RGC_results(ho,group)
      integer,          intent(in) :: ho                                                            !< homogenization type
      character(len=*), intent(in) :: group                                                         !< group name in HDF5 file
    end subroutine RGC_results

  end interface

  integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable :: &
    homogenization_type                                                                             !< type of each homogenization

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

  call material_parseHomogenization2()

  allocate(homogenization_dPdF(3,3,3,3,discretization_nIPs*discretization_Nelems), source=0.0_pReal)
  homogenization_F0 = spread(math_I3,3,discretization_nIPs*discretization_Nelems)                   ! initialize to identity
  homogenization_F = homogenization_F0                                                              ! initialize to identity
  allocate(homogenization_P(3,3,discretization_nIPs*discretization_Nelems),        source=0.0_pReal)

  num_homogMech => num_homog%get('mech',defaultVal=emptyDict)
  if (any(homogenization_type == HOMOGENIZATION_NONE_ID))      call pass_init
  if (any(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)) call isostrain_init
  if (any(homogenization_type == HOMOGENIZATION_RGC_ID))       call RGC_init(num_homogMech)

end subroutine mechanical_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition F onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_partition(subF,ce)

  real(pReal), intent(in), dimension(3,3) :: &
    subF
  integer,     intent(in) :: &
    ce

  integer :: co
  real(pReal), dimension (3,3,homogenization_Nconstituents(material_homogenizationID(ce))) :: Fs


  chosenHomogenization: select case(homogenization_type(material_homogenizationID(ce)))

    case (HOMOGENIZATION_NONE_ID) chosenHomogenization
      Fs(1:3,1:3,1) = subF

    case (HOMOGENIZATION_ISOSTRAIN_ID) chosenHomogenization
      call isostrain_partitionDeformation(Fs,subF)

    case (HOMOGENIZATION_RGC_ID) chosenHomogenization
      call RGC_partitionDeformation(Fs,subF,ce)

  end select chosenHomogenization

  do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
    call phase_set_F(Fs(1:3,1:3,co),co,ce)
  enddo


end subroutine mechanical_partition


!--------------------------------------------------------------------------------------------------
!> @brief Average P and dPdF from the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_homogenize(dt,ce)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: ce

  integer :: co


  homogenization_P(1:3,1:3,ce)            = phase_P(1,ce)
  homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = phase_mechanical_dPdF(dt,1,ce)
  do co = 2, homogenization_Nconstituents(material_homogenizationID(ce))
    homogenization_P(1:3,1:3,ce)            = homogenization_P(1:3,1:3,ce) &
                                            + phase_P(co,ce)
    homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = homogenization_dPdF(1:3,1:3,1:3,1:3,ce) &
                                            + phase_mechanical_dPdF(dt,co,ce)
  enddo

  homogenization_P(1:3,1:3,ce)            = homogenization_P(1:3,1:3,ce) &
                                          / real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)
  homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = homogenization_dPdF(1:3,1:3,1:3,1:3,ce) &
                                          / real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end subroutine mechanical_homogenize


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
module function mechanical_updateState(subdt,subF,ce) result(doneAndHappy)

  real(pReal), intent(in) :: &
    subdt                                                                                           !< current time step
  real(pReal), intent(in), dimension(3,3) :: &
    subF
  integer,     intent(in) :: &
    ce
  logical, dimension(2) :: doneAndHappy

  integer :: co
  real(pReal) :: dPdFs(3,3,3,3,homogenization_Nconstituents(material_homogenizationID(ce)))
  real(pReal) :: Fs(3,3,homogenization_Nconstituents(material_homogenizationID(ce)))
  real(pReal) :: Ps(3,3,homogenization_Nconstituents(material_homogenizationID(ce)))


  if (homogenization_type(material_homogenizationID(ce)) == HOMOGENIZATION_RGC_ID) then
      do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
        dPdFs(:,:,:,:,co) = phase_mechanical_dPdF(subdt,co,ce)
        Fs(:,:,co)        = phase_F(co,ce)
        Ps(:,:,co)        = phase_P(co,ce)
      enddo
      doneAndHappy = RGC_updateState(Ps,Fs,subF,subdt,dPdFs,ce)
  else
    doneAndHappy = .true.
  endif

end function mechanical_updateState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to file.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_results(group_base,ho)

  character(len=*), intent(in) :: group_base
  integer, intent(in)          :: ho

  character(len=:), allocatable :: group

  group = trim(group_base)//'/mechanical'
  call results_closeGroup(results_addGroup(group))

  select case(homogenization_type(ho))

    case(HOMOGENIZATION_rgc_ID)
      call RGC_results(ho,group)

  end select

  !temp = reshape(homogenization_F,[3,3,discretization_nIPs*discretization_Nelems])
  !call results_writeDataset(group,temp,'F',&
  !                          'deformation gradient','1')
  !temp = reshape(homogenization_P,[3,3,discretization_nIPs*discretization_Nelems])
  !call results_writeDataset(group,temp,'P',&
  !                          '1st Piola-Kirchhoff stress','Pa')

end subroutine mechanical_results


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization2()

  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogMech

  integer :: h

  material_homogenization => config_material%get('homogenization')

  allocate(homogenization_type(size(material_name_homogenization)), source=HOMOGENIZATION_undefined_ID)

  do h=1, size(material_name_homogenization)
    homog => material_homogenization%get(h)
    homogMech => homog%get('mechanical')
    select case (homogMech%get_asString('type'))
      case('pass')
        homogenization_type(h) = HOMOGENIZATION_NONE_ID
      case('isostrain')
        homogenization_type(h) = HOMOGENIZATION_ISOSTRAIN_ID
      case('RGC')
        homogenization_type(h) = HOMOGENIZATION_RGC_ID
      case default
        call IO_error(500,ext_msg=homogMech%get_asString('type'))
    end select
  enddo

end subroutine material_parseHomogenization2


end submodule mechanical
