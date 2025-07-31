submodule(homogenization) chemical

  interface

    module subroutine pass_init
    end subroutine pass_init

  end interface

  type :: tDataContainer
    real(pREAL), dimension(:,:), allocatable :: comp                                                ! Average composition
    real(pREAL), dimension(:,:), allocatable :: mu                                                  ! Diffusion potential
  end type tDataContainer

!  integer, dimension(:), allocatable :: &
!    homogenization_chemical_Ncomponents
!  integer :: &
!    homogenization_chemical_maxNcomponents
!
  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters), dimension(:), allocatable :: &
    param


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine chemical_init()

  class(tDict), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationChemical
  integer :: ho, Nmembers


  print'(/,1x,a)', '<<<+-  homogenization:chemical init  -+>>>'


  configHomogenizations => config_material%get_dict('homogenization')
  allocate(param(size(configHomogenizations)))
  allocate(current(size(configHomogenizations)))
  allocate(homogenization_chemical_Ncomponents(size(configHomogenizations)))

  do ho = 1, size(configHomogenizations)
    Nmembers = count(material_ID_homogenization == ho)
    configHomogenization => configHomogenizations%get_dict(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('chemical')) then
        configHomogenizationChemical => configHomogenization%get_dict('chemical')
        homogenization_chemical_Ncomponents(ho) = configHomogenizationChemical%get_asInt('N_components',defaultVal=0) ! should probably not be defined in homogenization?
#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(configHomogenizationChemical)
#else
        prm%output = configHomogenizationChemical%get_as1dStr('output',defaultVal=emptyStrArray)
#endif
        allocate(current(ho)%comp(homogenization_chemical_Ncomponents(ho), Nmembers),source=0.0_pREAL)
        allocate(current(ho)%mu(homogenization_chemical_Ncomponents(ho)-1, Nmembers),source=0.0_pREAL)

      end if
    end associate
  end do

  homogenization_chemical_maxNcomponents = maxval(homogenization_chemical_Ncomponents)

  call pass_init()

end subroutine chemical_init


!--------------------------------------------------------------------------------------------------
!> @brief Check if chemical homogemization description is present in the configuration file
!--------------------------------------------------------------------------------------------------
module function homogenization_chemical_active() result(active)

  logical :: active

  active = any(chemical_active(:))

end function homogenization_chemical_active


!--------------------------------------------------------------------------------------------------
!> @brief Partition chemical composition onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine chemical_partition(Delta_t, ce)

  real(pREAL), intent(in) :: Delta_t
  integer,     intent(in) :: ce

  real(pREAL), dimension(:),allocatable :: comp, mu
  integer :: co

  mu    = current(material_ID_homogenization(ce))%mu(:,material_entry_homogenization(ce))
  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    comp =  phase_calculate_composition(mu,co,ce)   ! need to calculate composition for each constituent for given diffusion potential
    !call phase_chemical_setField(comp,Delta_t,co,ce)
  end do

end subroutine chemical_partition

!--------------------------------------------------------------------------------------------------
!> @brief Calculate average composition at material point.
!--------------------------------------------------------------------------------------------------
module function homogenization_composition(mu,Delta_t,ce) result(comp)

  real(pREAL), dimension(:), intent(in) :: mu
  integer,     intent(in) :: ce
  real(pREAL), intent(in) :: Delta_t

  real(pREAL), dimension(:), allocatable :: comp
  integer :: co, N_components


  N_components = homogenization_chemical_Ncomponents(material_ID_homogenization(ce))
  allocate(comp(N_components),source = 0.0_pREAL)
  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    comp = comp + phase_calculate_composition(mu,co,ce)
  end do

  comp = comp/real(homogenization_Nconstituents(material_ID_homogenization(ce)),pREAL)

end function homogenization_composition


!--------------------------------------------------------------------------------------------------
!> @brief Calculate composition tangent at material point.
!--------------------------------------------------------------------------------------------------
module function homogenization_compositionTangent(mu,Delta_t,ce) result(comp_tangent)

  real(pREAL), dimension(:), intent(in) :: mu
  integer,     intent(in) :: ce
  real(pREAL), intent(in) :: Delta_t

  real(pREAL), dimension(:,:),allocatable :: comp_tangent
  integer :: co, N_components

  N_components = homogenization_chemical_Ncomponents(material_ID_homogenization(ce))
  allocate(comp_tangent(N_components-1,N_components-1),source = 0.0_pREAL)

  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    comp_tangent = comp_tangent + phase_compositionTangent(mu,co,ce)
  end do

  comp_tangent = comp_tangent/real(homogenization_Nconstituents(material_ID_homogenization(ce)),pREAL)

end function homogenization_compositionTangent


!--------------------------------------------------------------------------------------------------
!> @brief Retrieve mobility at each material point.
!--------------------------------------------------------------------------------------------------
module function homogenization_mobility(ce) result(mobility)

  integer, intent(in) :: ce
  real(pREAL), dimension(:,:), allocatable :: mobility

  integer :: co,N_components

  N_components = homogenization_chemical_Ncomponents(material_ID_homogenization(ce))
  allocate(mobility(N_components-1,N_components-1),source = 0.0_pREAL)
  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    mobility = mobility + phase_get_mobility(co,ce)
  end do

  mobility = mobility/homogenization_Nconstituents(material_ID_homogenization(ce))

end function homogenization_mobility


!--------------------------------------------------------------------------------------------------
!> @brief Set chemical field (mu)
!--------------------------------------------------------------------------------------------------
module subroutine homogenization_chemical_setField(mu, comp, Delta_t, ce)

  real(pREAL), dimension(:),  intent(in) :: mu
  real(pREAL), dimension(:),  intent(in) :: comp
  real(pREAL), intent(in) :: Delta_t
  integer, intent(in) :: ce


  current(material_ID_homogenization(ce))%mu(:,material_entry_homogenization(ce)) = mu
  current(material_ID_homogenization(ce))%comp(:,material_entry_homogenization(ce)) = comp  !NOTE: Average composition

  call chemical_partition(Delta_t,ce)

end subroutine homogenization_chemical_setField


!--------------------------------------------------------------------------------------------------
!> @brief Writes result to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine chemical_result(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o, ou

  associate(prm => param(ho))
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case('comp')
          do ou = 1, homogenization_chemical_Ncomponents(ho)
              call result_writeDataset(current(ho)%comp(ou,:),group,trim(material_name_species(ou)), &
                                       'concentration of '//trim(material_name_species(ou)),'mole fraction')
          end do
        case('mu')
          do ou = 1, homogenization_chemical_Ncomponents(ho)-1
              call result_writeDataset(current(ho)%mu(ou,:),group,'mu_'//trim(material_name_species(ou)), &
                                       'total potential of '//trim(material_name_species(ou)),'J')
          end do


      end select
    end do outputsLoop
  end associate

end subroutine chemical_result


end submodule chemical



