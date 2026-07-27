!--------------------------------------------------------------------------------------------------
!> @author Tilak Raj Pant, Indian Institute of Science
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
submodule(homogenization) electrical

  type :: tDataContainer
    real(pREAL), dimension(:,:), allocatable :: E
    real(pREAL), dimension(:,:), allocatable :: J
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: homog_electrical_current

  type :: tParameters
    character(len=pSTRLEN), allocatable, dimension(:) :: output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param

contains


!--------------------------------------------------------------------------------------------------
!> @brief Allocate storage and read output specifications.
!--------------------------------------------------------------------------------------------------
module subroutine electrical_init()

  type(tDict), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationElectrical
  integer :: ho, Nmembers


  print'(/,1x,a)', '<<<+-  homogenization:electrical init  -+>>>'

  configHomogenizations => config_material%get_dict('homogenization')
  allocate(param  (size(configHomogenizations)))
  allocate(homog_electrical_current(size(configHomogenizations)))

  do ho = 1, size(configHomogenizations)
    Nmembers = count(material_ID_homogenization == ho)
    allocate(homog_electrical_current(ho)%E(3,Nmembers), source=0.0_pREAL)
    allocate(homog_electrical_current(ho)%J(3,Nmembers), source=0.0_pREAL)

    configHomogenization => configHomogenizations%get_dict(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('electrical')) then
        configHomogenizationElectrical => &
          configHomogenization%get_dict('electrical')
#if defined(__GFORTRAN__)
        prm%output = output_as1dStr(configHomogenizationElectrical)
#else
        prm%output = configHomogenizationElectrical%get_as1dStr('output', &
                       defaultVal=emptyStrArray)
#endif
      else
        prm%output = emptyStrArray
      end if
    end associate
  end do

end subroutine electrical_init


!--------------------------------------------------------------------------------------------------
!> @brief Check whether electrical homogenization is active.
!--------------------------------------------------------------------------------------------------
module function homogenization_electrical_active() result(active)

  logical :: active


  active = any(electrical_active(:))

end function homogenization_electrical_active


!--------------------------------------------------------------------------------------------------
!> @brief Homogenized electrical conductivity via rule of mixtures.
!--------------------------------------------------------------------------------------------------
module function homogenization_sigma(ce) result(sigma)

  integer, intent(in)         :: ce
  real(pREAL), dimension(3,3) :: sigma
  integer :: co


  sigma = phase_sigma(1,ce) * material_v(1,ce)
  do co = 2, homogenization_Nconstituents(material_ID_homogenization(ce))
    sigma = sigma + phase_sigma(co,ce) * material_v(co,ce)
  end do

end function homogenization_sigma

!--------------------------------------------------------------------------------------------------
!> @brief Partition macro electric field onto the individual phase constituents of a cell.
!--------------------------------------------------------------------------------------------------
module subroutine electrical_partition(E,ce)

  real(pREAL), intent(in), dimension(3) :: E
  integer, intent(in) :: ce
  integer :: co

  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    call phase_electrical_setField(E,co,ce)
  end do

end subroutine electrical_partition

!--------------------------------------------------------------------------------------------------
!> @brief Compute J for one cell and store in J.
!--------------------------------------------------------------------------------------------------
module subroutine electrical_homogenize(E,ce)

  real(pREAL), intent(in), dimension(3) :: E
  integer, intent(in) :: ce
  integer :: co, ho, en
  real(pREAL), dimension(3) :: J

  ho = material_ID_homogenization(ce)
  en = material_entry_homogenization(ce)

  homog_electrical_current(ho)%E(:,en) = E

  J = matmul(phase_sigma(1,ce), E) * material_v(1,ce)

  do co = 2, homogenization_Nconstituents(ho)
    J = J + matmul(phase_sigma(co,ce), E) * material_v(co,ce)
  end do

  homog_electrical_current(ho)%J(:,en) = J

end subroutine electrical_homogenize

!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine electrical_result(ho, group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(ho))
    do o = 1, size(prm%output)
      select case(trim(prm%output(o)))
        case('E')
          call result_writeDataset(homog_electrical_current(ho)%E, group, 'E', 'electric field',  'V/m')
        case('J')
          call result_writeDataset(homog_electrical_current(ho)%J, group, 'J', 'current density', 'A/m^2')
      end select
    end do
  end associate

end subroutine electrical_result

end submodule electrical
