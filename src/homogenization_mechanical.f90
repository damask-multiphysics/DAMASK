!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Partition F and homogenize P/dPdF
!--------------------------------------------------------------------------------------------------
submodule(homogenization) mechanical


  interface

    module subroutine pass_init()
    end subroutine pass_init

    module subroutine isostrain_init()
    end subroutine isostrain_init

    module subroutine RGC_init()
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

  type :: tOutput                                                                                   !< requested output (per phase)
    character(len=pStringLen), allocatable, dimension(:) :: &
      label
  end type tOutput
  type(tOutput), allocatable, dimension(:) :: output_mechanical

  enum, bind(c); enumerator :: &
    MECHANICAL_UNDEFINED_ID, &
    MECHANICAL_PASS_ID, &
    MECHANICAL_ISOSTRAIN_ID, &
    MECHANICAL_RGC_ID
  end enum
  integer(kind(MECHANICAL_UNDEFINED_ID)), dimension(:),   allocatable :: &
    mechanical_type                                                                             !< type of each homogenization

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_init()

  print'(/,1x,a)', '<<<+-  homogenization:mechanical init  -+>>>'

  call parseMechanical()

  allocate(homogenization_dPdF(3,3,3,3,discretization_Ncells), source=0.0_pReal)
  homogenization_F0 = spread(math_I3,3,discretization_Ncells)
  homogenization_F = homogenization_F0
  allocate(homogenization_P(3,3,discretization_Ncells),source=0.0_pReal)

  if (any(mechanical_type == MECHANICAL_PASS_ID))      call pass_init()
  if (any(mechanical_type == MECHANICAL_ISOSTRAIN_ID)) call isostrain_init()
  if (any(mechanical_type == MECHANICAL_RGC_ID))       call RGC_init()

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


  chosenHomogenization: select case(mechanical_type(material_homogenizationID(ce)))

    case (MECHANICAL_PASS_ID) chosenHomogenization
      Fs(1:3,1:3,1) = subF

    case (MECHANICAL_ISOSTRAIN_ID) chosenHomogenization
      call isostrain_partitionDeformation(Fs,subF)

    case (MECHANICAL_RGC_ID) chosenHomogenization
      call RGC_partitionDeformation(Fs,subF,ce)

  end select chosenHomogenization

  do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
    call phase_set_F(Fs(1:3,1:3,co),co,ce)
  end do


end subroutine mechanical_partition


!--------------------------------------------------------------------------------------------------
!> @brief Average P and dPdF from the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_homogenize(Delta_t,ce)

  real(pReal), intent(in) :: Delta_t
  integer, intent(in) :: ce

  integer :: co


  homogenization_P(1:3,1:3,ce)            = phase_P(1,ce)*material_v(1,ce)
  homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = phase_mechanical_dPdF(Delta_t,1,ce)*material_v(1,ce)
  do co = 2, homogenization_Nconstituents(material_homogenizationID(ce))
    homogenization_P(1:3,1:3,ce)            = homogenization_P(1:3,1:3,ce) &
                                            + phase_P(co,ce)*material_v(co,ce)
    homogenization_dPdF(1:3,1:3,1:3,1:3,ce) = homogenization_dPdF(1:3,1:3,1:3,1:3,ce) &
                                            + phase_mechanical_dPdF(Delta_t,co,ce)*material_v(co,ce)
  end do

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


  if (mechanical_type(material_homogenizationID(ce)) == MECHANICAL_RGC_ID) then
      do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
        dPdFs(:,:,:,:,co) = phase_mechanical_dPdF(subdt,co,ce)
        Fs(:,:,co)        = phase_F(co,ce)
        Ps(:,:,co)        = phase_P(co,ce)
      end do
      doneAndHappy = RGC_updateState(Ps,Fs,subF,subdt,dPdFs,ce)
  else
    doneAndHappy = .true.
  end if

end function mechanical_updateState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to file.
!--------------------------------------------------------------------------------------------------
module subroutine mechanical_results(group_base,ho)

  character(len=*), intent(in) :: group_base
  integer, intent(in)          :: ho

  integer :: ou
  character(len=:), allocatable :: group


  group = trim(group_base)//'/mechanical'
  call results_closeGroup(results_addGroup(group))

  select case(mechanical_type(ho))

    case(MECHANICAL_RGC_ID)
      call RGC_results(ho,group)

  end select

  do ou = 1, size(output_mechanical(1)%label)

    select case (output_mechanical(ho)%label(ou))
      case('F')
        call results_writeDataset(reshape(homogenization_F,[3,3,discretization_nCells]),group,'F', &
                                  'deformation gradient','1')
      case('P')
        call results_writeDataset(reshape(homogenization_P,[3,3,discretization_nCells]),group,'P', &
                                  'first Piola-Kirchhoff stress','Pa')
    end select
  end do

end subroutine mechanical_results


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine parseMechanical()

  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    mechanical

  integer :: ho


  material_homogenization => config_material%get('homogenization')

  allocate(mechanical_type(size(material_name_homogenization)), source=MECHANICAL_UNDEFINED_ID)
  allocate(output_mechanical(size(material_name_homogenization)))

  do ho=1, size(material_name_homogenization)
    homog => material_homogenization%get(ho)
    mechanical => homog%get('mechanical')
#if defined(__GFORTRAN__)
    output_mechanical(ho)%label = output_as1dString(mechanical)
#else
    output_mechanical(ho)%label = mechanical%get_as1dString('output',defaultVal=emptyStringArray)
#endif
    select case (mechanical%get_asString('type'))
      case('pass')
        mechanical_type(ho) = MECHANICAL_PASS_ID
      case('isostrain')
        mechanical_type(ho) = MECHANICAL_ISOSTRAIN_ID
      case('RGC')
        mechanical_type(ho) = MECHANICAL_RGC_ID
      case default
        call IO_error(500,ext_msg=mechanical%get_asString('type'))
    end select
  end do

end subroutine parseMechanical


end submodule mechanical
