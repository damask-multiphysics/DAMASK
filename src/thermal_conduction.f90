!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for temperature evolution from heat conduction
!--------------------------------------------------------------------------------------------------
module thermal_conduction
  use prec
  use material
  use config
  use lattice
  use results
  use crystallite
  use source_thermal_dissipation
  use source_thermal_externalheat
 
  implicit none
  private

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters
  
  type(tparameters),             dimension(:), allocatable :: &
    param 
 
  public :: &
    thermal_conduction_init, &
    thermal_conduction_getSourceAndItsTangent, &
    thermal_conduction_getConductivity33, &
    thermal_conduction_getSpecificHeat, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_putTemperatureAndItsRate, &
    thermal_conduction_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_init

  
  integer :: maxNinstance,NofMyHomog,h
 
  write(6,'(/,a)') ' <<<+-  thermal_'//THERMAL_CONDUCTION_label//' init  -+>>>'; flush(6)
  
  maxNinstance = count(thermal_type == THERMAL_conduction_ID)
  if (maxNinstance == 0) return
  
  allocate(param(maxNinstance))
  
  do h = 1, size(thermal_type)
    if (thermal_type(h) /= THERMAL_conduction_ID) cycle
    associate(prm => param(thermal_typeInstance(h)),config => config_homogenization(h))
              
    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)
  
    NofMyHomog=count(material_homogenizationAt==h)
    thermalState(h)%sizeState = 0
    allocate(thermalState(h)%state0   (0,NofMyHomog))
    allocate(thermalState(h)%subState0(0,NofMyHomog))
    allocate(thermalState(h)%state    (0,NofMyHomog))
 
    thermalMapping(h)%p => material_homogenizationMemberAt
    deallocate(temperature    (h)%p)
    allocate  (temperature    (h)%p(NofMyHomog), source=thermal_initialT(h))
    deallocate(temperatureRate(h)%p)
    allocate  (temperatureRate(h)%p(NofMyHomog), source=0.0_pReal)
    
    end associate
  enddo
 
end subroutine thermal_conduction_init


!--------------------------------------------------------------------------------------------------
!> @brief returns heat generation rate
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_getSourceAndItsTangent(Tdot, dTdot_dT, T, ip, el)
 
  integer, intent(in) :: &
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
    offset, &
    instance, &
    grain, &
    source, &
    constituent
    
  homog  = material_homogenizationAt(el)
  offset = material_homogenizationMemberAt(ip,el)
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
 
end subroutine thermal_conduction_getSourceAndItsTangent
 

!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized thermal conductivity in reference configuration
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getConductivity33(ip,el)
  
  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), dimension(3,3) :: &
    thermal_conduction_getConductivity33
  integer :: &
    grain
    
   
  thermal_conduction_getConductivity33 = 0.0_pReal
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    thermal_conduction_getConductivity33 = thermal_conduction_getConductivity33 + &
     crystallite_push33ToRef(grain,ip,el,lattice_thermalConductivity33(:,:,material_phaseAt(grain,el)))
  enddo
 
  thermal_conduction_getConductivity33 = thermal_conduction_getConductivity33 &
                                       / real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end function thermal_conduction_getConductivity33


!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized specific heat capacity
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getSpecificHeat(ip,el)
  
  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal) :: &
    thermal_conduction_getSpecificHeat
  integer :: &
    grain
   
  thermal_conduction_getSpecificHeat = 0.0_pReal
  
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    thermal_conduction_getSpecificHeat = thermal_conduction_getSpecificHeat &
                                       + lattice_specificHeat(material_phaseAt(grain,el))
  enddo
 
  thermal_conduction_getSpecificHeat = thermal_conduction_getSpecificHeat &
                                     / real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end function thermal_conduction_getSpecificHeat


!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized mass density
!--------------------------------------------------------------------------------------------------
function thermal_conduction_getMassDensity(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal) :: &
    thermal_conduction_getMassDensity
  integer :: &
    grain
   
  thermal_conduction_getMassDensity = 0.0_pReal
  
   
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    thermal_conduction_getMassDensity = thermal_conduction_getMassDensity &
                                      + lattice_massDensity(material_phaseAt(grain,el))
  enddo
 
  thermal_conduction_getMassDensity = thermal_conduction_getMassDensity &
                                    / real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end function thermal_conduction_getMassDensity


!--------------------------------------------------------------------------------------------------
!> @brief updates thermal state with solution from heat conduction PDE
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_putTemperatureAndItsRate(T,Tdot,ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    T, &
    Tdot
  integer :: &
    homog, &
    offset  
  
  homog  = material_homogenizationAt(el)
  offset = thermalMapping(homog)%p(ip,el)
  temperature    (homog)%p(offset) = T
  temperatureRate(homog)%p(offset) = Tdot

end subroutine thermal_conduction_putTemperatureAndItsRate
 

!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group

  integer :: o
  
  associate(prm => param(damage_typeInstance(homog)))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('temperature')                                                                           ! ToDo: should be 'T'
        call results_writeDataset(group,temperature(homog)%p,'T',&
                                  'temperature','K')
    end select
  enddo outputsLoop
  end associate

end subroutine thermal_conduction_results

end module thermal_conduction
