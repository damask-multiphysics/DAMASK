!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for temperature evolution from heat conduction
!--------------------------------------------------------------------------------------------------
module thermal_conduction
  use prec
  use material
  use config
  use lattice
  use crystallite
  use source_thermal_dissipation
  use source_thermal_externalheat
 
  implicit none
  private
 
  integer,                     dimension(:,:), allocatable, target, public :: &
    thermal_conduction_sizePostResult                                                               !< size of each post result output
  character(len=64),           dimension(:,:), allocatable, target, public :: &
    thermal_conduction_output                                                                       !< name of each post result output
    
  integer,                     dimension(:),   allocatable, target, public :: &
    thermal_conduction_Noutput                                                                      !< number of outputs per instance of this damage 
 
  enum, bind(c) 
    enumerator :: undefined_ID, &
                  temperature_ID
  end enum
  integer(kind(undefined_ID)), dimension(:,:),  allocatable,          private :: & 
    thermal_conduction_outputID                                                                     !< ID of each post result output
 
 
  public :: &
    thermal_conduction_init, &
    thermal_conduction_getSourceAndItsTangent, &
    thermal_conduction_getConductivity33, &
    thermal_conduction_getSpecificHeat, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_putTemperatureAndItsRate, &
    thermal_conduction_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_conduction_init

  
  integer :: maxNinstance,section,instance,i
  integer :: sizeState
  integer :: NofMyHomog   
  character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  character(len=65536), dimension(:), allocatable :: outputs
 
  write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_CONDUCTION_label//' init  -+>>>'
  
  maxNinstance = count(thermal_type == THERMAL_conduction_ID)
  if (maxNinstance == 0) return
  
  allocate(thermal_conduction_sizePostResult (maxval(homogenization_Noutput),maxNinstance),source=0)
  allocate(thermal_conduction_output         (maxval(homogenization_Noutput),maxNinstance))
           thermal_conduction_output = ''
  allocate(thermal_conduction_outputID       (maxval(homogenization_Noutput),maxNinstance),source=undefined_ID)
  allocate(thermal_conduction_Noutput        (maxNinstance),                               source=0) 
 
  
  initializeInstances: do section = 1, size(thermal_type)
    if (thermal_type(section) /= THERMAL_conduction_ID) cycle
    NofMyHomog=count(material_homogenizationAt==section)
    instance = thermal_typeInstance(section)
    outputs = config_homogenization(section)%getStrings('(output)',defaultVal=emptyStringArray)
    do i=1, size(outputs)
      select case(outputs(i))
        case('temperature')
              thermal_conduction_Noutput(instance) = thermal_conduction_Noutput(instance) + 1
              thermal_conduction_outputID(thermal_conduction_Noutput(instance),instance) = temperature_ID
              thermal_conduction_output(thermal_conduction_Noutput(instance),instance) = outputs(i)
              thermal_conduction_sizePostResult(thermal_conduction_Noutput(instance),instance) = 1
      end select
    enddo
 
 
 ! allocate state arrays
    sizeState = 0
    thermalState(section)%sizeState = sizeState
    thermalState(section)%sizePostResults = sum(thermal_conduction_sizePostResult(:,instance))
    allocate(thermalState(section)%state0   (sizeState,NofMyHomog))
    allocate(thermalState(section)%subState0(sizeState,NofMyHomog))
    allocate(thermalState(section)%state    (sizeState,NofMyHomog))
 
    nullify(thermalMapping(section)%p)
    thermalMapping(section)%p => mappingHomogenization(1,:,:)
    deallocate(temperature    (section)%p)
    allocate  (temperature    (section)%p(NofMyHomog), source=thermal_initialT(section))
    deallocate(temperatureRate(section)%p)
    allocate  (temperatureRate(section)%p(NofMyHomog), source=0.0_pReal)
      
  enddo initializeInstances
 
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
  offset = mappingHomogenization(1,ip,el)
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
     crystallite_push33ToRef(grain,ip,el,lattice_thermalConductivity33(:,:,material_phase(grain,ip,el)))
  enddo
 
  thermal_conduction_getConductivity33 = &
    thermal_conduction_getConductivity33/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
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
    thermal_conduction_getSpecificHeat = thermal_conduction_getSpecificHeat + &
     lattice_specificHeat(material_phase(grain,ip,el))
  enddo
 
  thermal_conduction_getSpecificHeat = &
    thermal_conduction_getSpecificHeat/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
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
                                      + lattice_massDensity(material_phase(grain,ip,el))
  enddo
 
  thermal_conduction_getMassDensity = &
    thermal_conduction_getMassDensity/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
 
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
!> @brief return array of thermal results
!--------------------------------------------------------------------------------------------------
function thermal_conduction_postResults(homog,instance,of) result(postResults)
 
  integer,              intent(in) :: &
    homog, &
    instance, &
    of
 
  real(pReal), dimension(sum(thermal_conduction_sizePostResult(:,instance))) :: &
    postResults
 
  integer :: &
    o, c
 
  c = 0
  do o = 1,thermal_conduction_Noutput(instance)
     select case(thermal_conduction_outputID(o,instance))
  
       case (temperature_ID)
         postResults(c+1) = temperature(homog)%p(of)
         c = c + 1
     end select
  enddo
 
end function thermal_conduction_postResults

end module thermal_conduction
