!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) thermal

  use lattice

  interface

    module subroutine pass_init
    end subroutine pass_init

    module subroutine isotemperature_init
    end subroutine isotemperature_init

  end interface

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: T, dot_T
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_init()

  class(tNode), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationThermal
  integer :: ho


  print'(/,a)', ' <<<+-  homogenization:thermal init  -+>>>'


  configHomogenizations => config_material%get('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%T(count(material_homogenizationID==ho)), source=300.0_pReal)
    allocate(current(ho)%dot_T(count(material_homogenizationID==ho)), source=0.0_pReal)
    configHomogenization => configHomogenizations%get(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('thermal')) then
        configHomogenizationThermal => configHomogenization%get('thermal')
#if defined (__GFORTRAN__)
        prm%output = output_as1dString(configHomogenizationThermal)
#else
        prm%output = configHomogenizationThermal%get_as1dString('output',defaultVal=emptyStringArray)
#endif
      else
        prm%output = emptyStringArray
      endif
    end associate
  enddo

  call pass_init()

end subroutine thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_partition(ce)

  integer,     intent(in) :: ce

  real(pReal) :: T, dot_T
  integer :: co


  T     = current(material_homogenizationID(ce))%T(material_homogenizationEntry(ce))
  dot_T = current(material_homogenizationID(ce))%dot_T(material_homogenizationEntry(ce))
  do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
    call phase_thermal_setField(T,dot_T,co,ce)
  enddo

end subroutine thermal_partition


!--------------------------------------------------------------------------------------------------
!> @brief Homogenized thermal viscosity.
!--------------------------------------------------------------------------------------------------
module function homogenization_mu_T(ce) result(mu)

  integer, intent(in) :: ce
  real(pReal) :: mu

  integer :: co


  mu = phase_mu_T(1,ce)
  do co = 2, homogenization_Nconstituents(material_homogenizationID(ce))
    mu = mu + phase_mu_T(co,ce)
  enddo

  mu = mu / real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end function homogenization_mu_T


!--------------------------------------------------------------------------------------------------
!> @brief Homogenized thermal conductivity in reference configuration.
!--------------------------------------------------------------------------------------------------
module function homogenization_K_T(ce) result(K)

  integer, intent(in) :: ce
  real(pReal), dimension(3,3) :: K

  integer :: co


  K = phase_K_T(1,ce)
  do co = 2, homogenization_Nconstituents(material_homogenizationID(ce))
    K = K + phase_K_T(co,ce)
  enddo

  K = K / real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end function homogenization_K_T


!--------------------------------------------------------------------------------------------------
!> @brief Homogenized heat generation rate.
!--------------------------------------------------------------------------------------------------
module function homogenization_f_T(ce) result(f)

  integer, intent(in) :: ce
  real(pReal) :: f

  integer :: co


  f = phase_f_T(material_phaseID(1,ce),material_phaseEntry(1,ce))
  do co = 2, homogenization_Nconstituents(material_homogenizationID(ce))
    f = f + phase_f_T(material_phaseID(co,ce),material_phaseEntry(co,ce))
  enddo

  f = f/real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end function homogenization_f_T


!--------------------------------------------------------------------------------------------------
!> @brief Set thermal field and its rate (T and dot_T).
!--------------------------------------------------------------------------------------------------
module subroutine homogenization_thermal_setField(T,dot_T, ce)

  integer, intent(in) :: ce
  real(pReal),   intent(in) :: T, dot_T


  current(material_homogenizationID(ce))%T(material_homogenizationEntry(ce)) = T
  current(material_homogenizationID(ce))%dot_T(material_homogenizationEntry(ce)) = dot_T


end subroutine homogenization_thermal_setField


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine thermal_results(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(ho))
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case('T')
          call results_writeDataset(current(ho)%T,group,'T','temperature','K')
      end select
    enddo outputsLoop
  end associate

end subroutine thermal_results

end submodule thermal
