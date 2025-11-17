!----------------------------------------------------------------------------------------------------
!> @brief internal microstructure state for all chemical sources and kinematics constitutive models
!----------------------------------------------------------------------------------------------------
submodule(phase) chemical

  type :: tChemicalParameters
    real(pREAL) ::                 V_m = 0.0_pREAL                                                  !< molar volume
  end type tChemicalParameters

  integer(kind(UNDEFINED)),  dimension(:), allocatable :: &
    chemical_energy

  type :: tDataContainer             ! ?? not very telling name. Better: "fieldQuantities" ??
    real(pREAL), dimension(:,:), allocatable :: C, dot_C, C0
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable  :: current         ! ?? not very telling name. Better: "field" ?? MD: current(ho)%T(en) reads quite good

  type(tChemicalParameters), dimension(:), allocatable :: param


  interface


    module function quadEnergy_init() result(myChemicalEnergy)
      logical, dimension(:), allocatable :: mychemicalEnergy
    end function quadEnergy_init

    module function quadEnergy_composition(mu_chemical,ph,en) result(comp)
      real(pREAL), dimension(:), intent(in) :: mu_chemical
      integer, intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(:),allocatable :: comp
    end function quadEnergy_composition

    module function quadEnergy_compositionTangent(mu_chemical,ph,en) result(comp_tangent)
      real(pREAL), dimension(:), intent(in) :: mu_chemical
      integer, intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(:,:),allocatable :: comp_tangent
    end function quadEnergy_compositionTangent

    module function quadEnergy_mobility(ph,en) result(mobility)
      integer, intent(in) :: ph, en
      real(pREAL), dimension(:,:),allocatable :: mobility
    end function quadEnergy_mobility

    module subroutine quadEnergy_results(ph,comp,group)
      integer,          intent(in) :: ph
      real(pREAL), dimension(:,:),intent(in) :: comp
      character(len=*), intent(in) :: group
    end subroutine quadEnergy_results


 end interface

contains

!----------------------------------------------------------------------------------------------
!< @brief initializes chemical energy and sources/kinematics (eventually)
!----------------------------------------------------------------------------------------------
module subroutine chemical_init(phases)

  type(tDict), pointer :: &
    phases

  type(tDict), pointer :: &
    phase, &
    chemical, &
    components, &
    component

  integer :: &
    ph, &
    Nmembers


  print'(/,a)', ' <<<+-  phase:chemical init  -+>>>'

  allocate(current(size(phases)))
  allocate(chemical_energy(size(phases)),source=UNDEFINED)
  allocate(param(size(phases)))

  phases => config_material%get_dict('phase')
  do ph = 1, size(phases)
    phase => phases%get_dict(ph)
    chemical => phase%get_dict('chemical',defaultVal=emptyDict)
    param(ph)%V_m = chemical%get_asReal('V_m',defaultVal=1.0_pREAL)
  end do

  !initialize chemical energy model
  where(quadEnergy_init())       chemical_energy = CHEMICAL_QUADENERGY

  do ph = 1,size(phases)
    Nmembers = count(material_ID_phase == ph)
    if (chemical_energy(ph) == UNDEFINED) then
      allocate(current(ph)%C(1,Nmembers),source=0.0_pREAL)
      allocate(current(ph)%dot_C(1,Nmembers),source=0.0_pREAL)
      allocate(current(ph)%C0(1,Nmembers),source=0.0_pREAL)
    end if
  end do


end subroutine chemical_init


!----------------------------------------------------------------------------------------------
!< @brief Calculates composition for constituent.
!----------------------------------------------------------------------------------------------
module function phase_calculate_composition(mu,co,ce) result(conc)

  real(pREAL), intent(in), dimension(:) :: mu
  integer, intent(in) :: co, ce
  real(pREAL), dimension(:), allocatable :: conc
  real(pREAL), dimension(:), allocatable :: mu_chemical

  integer :: &
    ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  mu_chemical = mu
  chemicalEnergyType: select case (chemical_energy(ph))

    case (CHEMICAL_QUADENERGY)
      conc = quadEnergy_composition(mu_chemical,ph,en)

  end select chemicalEnergyType

end function phase_calculate_composition


!----------------------------------------------------------------------------------------------
!< @brief Retrieves composition for constituent.
!----------------------------------------------------------------------------------------------
module function phase_get_mobility(co,ce) result(mobility)

  integer, intent(in) :: co, ce
  real(pREAL), dimension(:,:),allocatable :: mobility

  integer :: &
    ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  chemicalEnergyType: select case (chemical_energy(ph))

    case (CHEMICAL_QUADENERGY)
      mobility = quadEnergy_mobility(ph,en)

  end select chemicalEnergyType

end function phase_get_mobility


!----------------------------------------------------------------------------------------------
!< @brief Retrieves composition tangent for constituent.
!----------------------------------------------------------------------------------------------
module function phase_compositionTangent(mu,co,ce) result(comp_tangent)

  real(pREAL), dimension(:), intent(in) :: mu
  integer, intent(in) :: co, ce
  real(pREAL), dimension(:,:),allocatable :: comp_tangent
  real(pREAL), dimension(:),  allocatable :: mu_chemical

  integer :: &
    ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  mu_chemical = mu
  chemicalEnergyType: select case (chemical_energy(ph))

    case (CHEMICAL_QUADENERGY)
      comp_tangent = quadEnergy_compositionTangent(mu_chemical,ph,en)

  end select chemicalEnergyType

end function phase_compositionTangent


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
!--------------------------------------------------------------------------------------------------
module subroutine chemical_forward()

  integer :: ph

  do ph = 1, size(current)
    current(ph)%C0 = current(ph)%C
  end do

end subroutine chemical_forward


!----------------------------------------------------------------------------------------------
!< @brief writes chemical results to HDF5 output file
!----------------------------------------------------------------------------------------------
module subroutine chemical_result(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph


  if (chemical_energy(ph) /= UNDEFINED) &
    call result_closeGroup(result_addGroup(group//'chemical'))

  chemicalEnergyType: select case (chemical_energy(ph))

    case (CHEMICAL_QUADENERGY)
      call quadEnergy_results(ph,current(ph)%C,group//'chemical/')

  end select chemicalEnergyType

end subroutine chemical_result


!--------------------------------------------------------------------------------------------------
!> @brief checks if a chemical energy module is active or not
!--------------------------------------------------------------------------------------------------
function chemical_active(chemical_label)  result(active_chemical)

  character(len=*), intent(in)       :: chemical_label                                              !< type of chemical energy model
  logical, dimension(:), allocatable :: active_chemical

  type(tDict), pointer :: &
    phases, &
    phase, &
    chemical
  integer :: ph

  phases => config_material%get_dict('phase')
  allocate(active_chemical(size(phases)), source = .false. )
  do ph = 1, size(phases)
    phase => phases%get_dict(ph)
    chemical => phase%get_dict('chemical',defaultVal=emptyDict)
    if (chemical%get_asStr('type',defaultVal='none') == chemical_label) active_chemical(ph) = .true.
  end do

end function chemical_active


end submodule chemical
