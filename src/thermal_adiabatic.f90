!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for adiabatic temperature evolution
!--------------------------------------------------------------------------------------------------
module thermal_adiabatic
  use prec
  use config
  use material
  use results
  use constitutive
  use YAML_types
  use crystallite
  use lattice

  implicit none
  private

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters
  
  type(tparameters),             dimension(:), allocatable :: &
    param
 
  public :: &
    thermal_adiabatic_init, &
    thermal_adiabatic_updateState, &
    thermal_adiabatic_getSourceAndItsTangent, &
    thermal_adiabatic_getSpecificHeat, &
    thermal_adiabatic_getMassDensity, &
    thermal_adiabatic_results
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_init
 
  integer :: maxNinstance,h,NofMyHomog
  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogThermal

  print'(/,a)', ' <<<+-  thermal_adiabatic init  -+>>>'; flush(6)
  
  maxNinstance = count(thermal_type == THERMAL_adiabatic_ID)
  if (maxNinstance == 0) return
  
  allocate(param(maxNinstance))
  
  material_homogenization => config_material%get('homogenization')
  do h = 1, size(material_name_homogenization)
    if (thermal_type(h) /= THERMAL_adiabatic_ID) cycle
    homog => material_homogenization%get(h)
    homogThermal => homog%get('thermal')
 
    associate(prm => param(thermal_typeInstance(h)))

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(homogThermal)
#else
    prm%output = homogThermal%get_asStrings('output',defaultVal=emptyStringArray)
#endif

    NofMyHomog=count(material_homogenizationAt==h)
    thermalState(h)%sizeState = 1
    allocate(thermalState(h)%state0   (1,NofMyHomog), source=thermal_initialT(h))
    allocate(thermalState(h)%subState0(1,NofMyHomog), source=thermal_initialT(h))
    allocate(thermalState(h)%state    (1,NofMyHomog), source=thermal_initialT(h))
 
    thermalMapping(h)%p => material_homogenizationMemberAt
    deallocate(temperature(h)%p)
    temperature(h)%p => thermalState(h)%state(1,:)
    deallocate(temperatureRate(h)%p)
    allocate  (temperatureRate(h)%p(NofMyHomog), source=0.0_pReal)
  
    end associate
  enddo
 
end subroutine thermal_adiabatic_init


!--------------------------------------------------------------------------------------------------
!> @brief  calculates adiabatic change in temperature based on local heat generation model  
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_updateState(subdt, ip, el)
  
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    subdt

  logical, dimension(2)  :: &
    thermal_adiabatic_updateState
  integer :: &
    homog, &
    offset
  real(pReal) :: &
    T, Tdot, dTdot_dT  
 
  homog  = material_homogenizationAt(el)
  offset = material_homogenizationMemberAt(ip,el)
  
  T = thermalState(homog)%subState0(1,offset)
  call thermal_adiabatic_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
  T = T + subdt*Tdot/(thermal_adiabatic_getSpecificHeat(ip,el)*thermal_adiabatic_getMassDensity(ip,el))
  
  thermal_adiabatic_updateState = [     abs(T - thermalState(homog)%state(1,offset)) &
                                     <= 1.0e-2_pReal &
                                   .or. abs(T - thermalState(homog)%state(1,offset)) &
                                     <= 1.0e-6_pReal*abs(thermalState(homog)%state(1,offset)), &
                                   .true.]
 
  temperature    (homog)%p(thermalMapping(homog)%p(ip,el)) = T  
  temperatureRate(homog)%p(thermalMapping(homog)%p(ip,el)) = &
    (thermalState(homog)%state(1,offset) - thermalState(homog)%subState0(1,offset))/(subdt+tiny(0.0_pReal))
  
end function thermal_adiabatic_updateState


!--------------------------------------------------------------------------------------------------
!> @brief returns heat generation rate
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
 
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), intent(in) :: &
    T
  real(pReal), intent(out) :: &
    Tdot, dTdot_dT
  integer :: &
    homog
 
  Tdot = 0.0_pReal
  dTdot_dT = 0.0_pReal
 
  homog  = material_homogenizationAt(el)
  call constitutive_thermal_getRateAndItsTangents(TDot, dTDot_dT, T, crystallite_S, crystallite_Lp, ip, el)

  Tdot = Tdot/real(homogenization_Nconstituent(homog),pReal)
  dTdot_dT = dTdot_dT/real(homogenization_Nconstituent(homog),pReal)
 
end subroutine thermal_adiabatic_getSourceAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized specific heat capacity
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_getSpecificHeat(ip,el)
 
  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  
  real(pReal) :: &
    thermal_adiabatic_getSpecificHeat
  integer :: &
    grain
   
  thermal_adiabatic_getSpecificHeat = 0.0_pReal
  
  do grain = 1, homogenization_Nconstituent(material_homogenizationAt(el))
    thermal_adiabatic_getSpecificHeat = thermal_adiabatic_getSpecificHeat & 
                                      + lattice_c_p(material_phaseAt(grain,el))
  enddo
 
  thermal_adiabatic_getSpecificHeat = thermal_adiabatic_getSpecificHeat &
                                    / real(homogenization_Nconstituent(material_homogenizationAt(el)),pReal)
  
end function thermal_adiabatic_getSpecificHeat
 
 
!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized mass density
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_getMassDensity(ip,el)
    
  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal) :: &
    thermal_adiabatic_getMassDensity
  integer :: &
    grain
   
  thermal_adiabatic_getMassDensity = 0.0_pReal
 
  do grain = 1, homogenization_Nconstituent(material_homogenizationAt(el))
    thermal_adiabatic_getMassDensity = thermal_adiabatic_getMassDensity &
                                     + lattice_rho(material_phaseAt(grain,el))
  enddo
 
  thermal_adiabatic_getMassDensity = thermal_adiabatic_getMassDensity &
                                   / real(homogenization_Nconstituent(material_homogenizationAt(el)),pReal)
 
end function thermal_adiabatic_getMassDensity


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group

  integer :: o
  
  associate(prm => param(damage_typeInstance(homog)))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('T') 
        call results_writeDataset(group,temperature(homog)%p,'T',&
                                  'temperature','K')
    end select
  enddo outputsLoop
  end associate

end subroutine thermal_adiabatic_results

end module thermal_adiabatic
