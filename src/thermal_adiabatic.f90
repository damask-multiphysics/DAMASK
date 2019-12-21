!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for adiabatic temperature evolution
!--------------------------------------------------------------------------------------------------
module thermal_adiabatic
  use prec
  use config
  use numerics
  use material
  use results
  use source_thermal_dissipation
  use source_thermal_externalheat
  use crystallite
  use lattice

  implicit none
  private

  enum, bind(c) 
    enumerator :: undefined_ID, &
                  temperature_ID
  end enum
  
  type :: tParameters
    integer(kind(undefined_ID)), dimension(:), allocatable :: &
      outputID
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
 
  integer :: maxNinstance,o,h,NofMyHomog   
  character(len=pStringLen), dimension(:), allocatable :: outputs
 
  write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_ADIABATIC_label//' init  -+>>>'; flush(6)
  
  maxNinstance = count(thermal_type == THERMAL_adiabatic_ID)
  if (maxNinstance == 0) return
  
  allocate(param(maxNinstance))
  
  do h = 1, size(thermal_type)
    if (thermal_type(h) /= THERMAL_adiabatic_ID) cycle
    associate(prm => param(thermal_typeInstance(h)),config => config_homogenization(h))
              
    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))

    do o=1, size(outputs)
      select case(outputs(o))
        case('temperature')
          prm%outputID = [prm%outputID, temperature_ID]
      end select
    enddo

    NofMyHomog=count(material_homogenizationAt==h)
    thermalState(h)%sizeState = 1
    allocate(thermalState(h)%state0   (1,NofMyHomog), source=thermal_initialT(h))
    allocate(thermalState(h)%subState0(1,NofMyHomog), source=thermal_initialT(h))
    allocate(thermalState(h)%state    (1,NofMyHomog), source=thermal_initialT(h))
 
    nullify(thermalMapping(h)%p)
    thermalMapping(h)%p => mappingHomogenization(1,:,:)
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
  offset = mappingHomogenization(1,ip,el)
  
  T = thermalState(homog)%subState0(1,offset)
  call thermal_adiabatic_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
  T = T + subdt*Tdot/(thermal_adiabatic_getSpecificHeat(ip,el)*thermal_adiabatic_getMassDensity(ip,el))
  
  thermal_adiabatic_updateState = [     abs(T - thermalState(homog)%state(1,offset)) &
                                     <= err_thermal_tolAbs &
                                   .or. abs(T - thermalState(homog)%state(1,offset)) &
                                     <= err_thermal_tolRel*abs(thermalState(homog)%state(1,offset)), &
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
 
  real(pReal) :: &
    my_Tdot, my_dTdot_dT
  integer :: &
    phase, &
    homog, &
    instance, &
    grain, &
    source, &
    constituent
    
  homog  = material_homogenizationAt(el)
  instance = thermal_typeInstance(homog)
   
  Tdot = 0.0_pReal
  dTdot_dT = 0.0_pReal
  do grain = 1, homogenization_Ngrains(homog)
    phase = material_phaseAt(grain,el)
    constituent = material_phasememberAt(grain,ip,el)
    do source = 1, phase_Nsources(phase)
      select case(phase_source(source,phase))                                                   
        case (SOURCE_thermal_dissipation_ID)
         call source_thermal_dissipation_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                              crystallite_S(1:3,1:3,grain,ip,el), &
                                                              crystallite_Lp(1:3,1:3,grain,ip,el), &
                                                              phase)
 
        case (SOURCE_thermal_externalheat_ID)
         call source_thermal_externalheat_getRateAndItsTangent(my_Tdot, my_dTdot_dT, &
                                                               phase, constituent)
 
        case default
         my_Tdot = 0.0_pReal
         my_dTdot_dT = 0.0_pReal
      end select
      Tdot = Tdot + my_Tdot
      dTdot_dT = dTdot_dT + my_dTdot_dT
    enddo  
  enddo
  
  Tdot = Tdot/real(homogenization_Ngrains(homog),pReal)
  dTdot_dT = dTdot_dT/real(homogenization_Ngrains(homog),pReal)
 
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
  
   
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    thermal_adiabatic_getSpecificHeat = thermal_adiabatic_getSpecificHeat + &
     lattice_specificHeat(material_phaseAt(grain,el))
  enddo
 
  thermal_adiabatic_getSpecificHeat = &
    thermal_adiabatic_getSpecificHeat/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
  
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
 
   
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    thermal_adiabatic_getMassDensity = thermal_adiabatic_getMassDensity + &
     lattice_massDensity(material_phaseAt(grain,el))
  enddo
 
  thermal_adiabatic_getMassDensity = &
    thermal_adiabatic_getMassDensity/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end function thermal_adiabatic_getMassDensity


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group
  integer :: o
  
  associate(prm => param(damage_typeInstance(homog)))

  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
    
      case (temperature_ID)
        call results_writeDataset(group,temperature(homog)%p,'T',&
                                  'temperature','K')
    end select
  enddo outputsLoop
  end associate

end subroutine thermal_adiabatic_results

end module thermal_adiabatic
