!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for adiabatic temperature evolution
!--------------------------------------------------------------------------------------------------
module thermal_adiabatic
  use prec
  use config
  use numerics
  use material
  use source_thermal_dissipation
  use source_thermal_externalheat
  use crystallite
  use lattice

  implicit none
  private
 
  integer,                       dimension(:,:),  allocatable, target, public :: &
    thermal_adiabatic_sizePostResult                                                                !< size of each post result output
  character(len=64),             dimension(:,:),  allocatable, target, public :: &
    thermal_adiabatic_output                                                                        !< name of each post result output
    
  integer,                       dimension(:),    allocatable, target, public :: &
    thermal_adiabatic_Noutput                                                                       !< number of outputs per instance of this thermal model 
 
  enum, bind(c) 
    enumerator :: undefined_ID, &
                  temperature_ID
  end enum
  integer(kind(undefined_ID)),   dimension(:,:),  allocatable :: & 
    thermal_adiabatic_outputID                                                                      !< ID of each post result output
 
 
  public :: &
    thermal_adiabatic_init, &
    thermal_adiabatic_updateState, &
    thermal_adiabatic_getSourceAndItsTangent, &
    thermal_adiabatic_getSpecificHeat, &
    thermal_adiabatic_getMassDensity, &
    thermal_adiabatic_postResults
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_adiabatic_init
 
  integer :: maxNinstance,section,instance,i,sizeState,NofMyHomog   
  character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  character(len=65536), dimension(:), allocatable :: outputs
 
  write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_ADIABATIC_label//' init  -+>>>'
  
  maxNinstance = count(thermal_type == THERMAL_adiabatic_ID)
  if (maxNinstance == 0) return
  
  allocate(thermal_adiabatic_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0)
  allocate(thermal_adiabatic_output         (maxval(homogenization_Noutput),maxNinstance))
           thermal_adiabatic_output = ''
  allocate(thermal_adiabatic_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
  allocate(thermal_adiabatic_Noutput        (maxNinstance),                               source=0) 
 
  
  initializeInstances: do section = 1, size(thermal_type)
    if (thermal_type(section) /= THERMAL_adiabatic_ID) cycle
    NofMyHomog=count(material_homogenizationAt==section)
    instance = thermal_typeInstance(section)
    outputs = config_homogenization(section)%getStrings('(output)',defaultVal=emptyStringArray)
    do i=1, size(outputs)
      select case(outputs(i))
        case('temperature')
              thermal_adiabatic_Noutput(instance) = thermal_adiabatic_Noutput(instance) + 1
              thermal_adiabatic_outputID(thermal_adiabatic_Noutput(instance),instance) = temperature_ID
              thermal_adiabatic_output(thermal_adiabatic_Noutput(instance),instance) = outputs(i)
              thermal_adiabatic_sizePostResult(thermal_adiabatic_Noutput(instance),instance) = 1
      end select
    enddo
 
 ! allocate state arrays
    sizeState = 1
    thermalState(section)%sizeState = sizeState
    thermalState(section)%sizePostResults = sum(thermal_adiabatic_sizePostResult(:,instance))
    allocate(thermalState(section)%state0   (sizeState,NofMyHomog), source=thermal_initialT(section))
    allocate(thermalState(section)%subState0(sizeState,NofMyHomog), source=thermal_initialT(section))
    allocate(thermalState(section)%state    (sizeState,NofMyHomog), source=thermal_initialT(section))
 
    nullify(thermalMapping(section)%p)
    thermalMapping(section)%p => mappingHomogenization(1,:,:)
    deallocate(temperature(section)%p)
    temperature(section)%p => thermalState(section)%state(1,:)
    deallocate(temperatureRate(section)%p)
    allocate  (temperatureRate(section)%p(NofMyHomog), source=0.0_pReal)
      
  enddo initializeInstances
 
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
    phase = phaseAt(grain,ip,el)
    constituent = phasememberAt(grain,ip,el)
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
     lattice_specificHeat(material_phase(grain,ip,el))
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
     lattice_massDensity(material_phase(grain,ip,el))
  enddo
 
  thermal_adiabatic_getMassDensity = &
    thermal_adiabatic_getMassDensity/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
end function thermal_adiabatic_getMassDensity


!--------------------------------------------------------------------------------------------------
!> @brief return array of thermal results
!--------------------------------------------------------------------------------------------------
function thermal_adiabatic_postResults(homog,instance,of) result(postResults)
 
  integer, intent(in) :: &
    homog, &
    instance, &
    of
 
  real(pReal), dimension(sum(thermal_adiabatic_sizePostResult(:,instance))) :: &
    postResults
 
  integer :: &
    o, c
 
  c = 0
 
  do o = 1,thermal_adiabatic_Noutput(instance)
     select case(thermal_adiabatic_outputID(o,instance))
  
       case (temperature_ID)
         postResults(c+1) = temperature(homog)%p(of)
         c = c + 1
     end select
  enddo
 
end function thermal_adiabatic_postResults

end module thermal_adiabatic
