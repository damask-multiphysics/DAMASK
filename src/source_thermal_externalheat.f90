!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!> @brief material subroutine for variable heat source
!--------------------------------------------------------------------------------------------------
module source_thermal_externalheat
  use prec, only: &
    pReal

  implicit none
  private
  integer,           dimension(:),   allocatable,         public, protected :: &
    source_thermal_externalheat_offset, &                                                                 !< which source is my current thermal dissipation mechanism?
    source_thermal_externalheat_instance                                                                  !< instance of thermal dissipation source mechanism

  integer,           dimension(:,:), allocatable, target, public :: &
    source_thermal_externalheat_sizePostResult                                                            !< size of each post result output

  character(len=64), dimension(:,:), allocatable, target, public :: &
    source_thermal_externalheat_output                                                                    !< name of each post result output
   
  integer,                       dimension(:),   allocatable, target, public :: &
    source_thermal_externalheat_Noutput                                                                   !< number of outputs per instance of this source 

  type, private :: tParameters                                                                       !< container type for internal constitutive parameters
    real(pReal), dimension(:), allocatable :: &
      time, &
      heat_rate
    integer :: &
     nIntervals
  end type tParameters

  type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_thermal_externalheat_init, &
    source_thermal_externalheat_dotState, &
    source_thermal_externalheat_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_init
  use debug, only: &
    debug_level,&
    debug_constitutive,&
    debug_levelBasic
  use material, only: &
    material_allocateSourceState, &
    material_phase, &
    phase_source, &
    phase_Nsources, &
    phase_Noutput, &
    SOURCE_thermal_externalheat_label, &
    SOURCE_thermal_externalheat_ID
  use config, only: &
    config_phase, &
    material_Nphase
 
 
  integer :: maxNinstance,instance,source,sourceOffset,NofMyPhase,p   
 
  write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_thermal_externalheat_label//' init  -+>>>'
 
  
  maxNinstance = count(phase_source == SOURCE_thermal_externalheat_ID)
  if (maxNinstance == 0) return
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
  
  allocate(source_thermal_externalheat_offset(material_Nphase), source=0)
  allocate(source_thermal_externalheat_instance(material_Nphase), source=0)
  
  do p = 1, material_Nphase
    source_thermal_externalheat_instance(p) = count(phase_source(:,1:p) == SOURCE_thermal_externalheat_ID)
    do source = 1, phase_Nsources(p)
      if (phase_source(source,p) == SOURCE_thermal_externalheat_ID) &
        source_thermal_externalheat_offset(p) = source
    enddo    
  enddo
    
  allocate(source_thermal_externalheat_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0)
  allocate(source_thermal_externalheat_output  (maxval(phase_Noutput),maxNinstance))
           source_thermal_externalheat_output = ''
  allocate(source_thermal_externalheat_Noutput(maxNinstance),                             source=0) 
 
  allocate(param(maxNinstance))
 
  do p=1, size(config_phase)
    if (all(phase_source(:,p) /= SOURCE_thermal_externalheat_ID)) cycle
    instance = source_thermal_externalheat_instance(p)
    sourceOffset = source_thermal_externalheat_offset(p)
    NofMyPhase=count(material_phase==p)
 
    param(instance)%time       = config_phase(p)%getFloats('externalheat_time')
    param(instance)%nIntervals = size(param(instance)%time) - 1
    
    
    param(instance)%heat_rate = config_phase(p)%getFloats('externalheat_rate',requiredSize = size(param(instance)%time))
    
    call material_allocateSourceState(p,sourceOffset,NofMyPhase,1,1,0)
    
  enddo

end subroutine source_thermal_externalheat_init


!--------------------------------------------------------------------------------------------------
!> @brief rate of change of state
!> @details state only contains current time to linearly interpolate given heat powers
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_dotState(phase, of)
  use material, only: &
    sourceState
 
  integer, intent(in) :: &
    phase, &
    of
  integer :: &
    sourceOffset
 
  sourceOffset = source_thermal_externalheat_offset(phase)
  
  sourceState(phase)%p(sourceOffset)%dotState(1,of) = 1.0_pReal                                     ! state is current time

end subroutine source_thermal_externalheat_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns local heat generation rate 
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_getRateAndItsTangent(TDot, dTDot_dT, phase, of)
  use material, only: &
    sourceState
 
  integer, intent(in) :: &
    phase, &
    of
  real(pReal),  intent(out) :: &
    TDot, &
    dTDot_dT
  integer :: &
    instance, sourceOffset, interval
  real(pReal) :: &
    frac_time   
 
  instance = source_thermal_externalheat_instance(phase)
  sourceOffset = source_thermal_externalheat_offset(phase)
 
  do interval = 1, param(instance)%nIntervals                                                       ! scan through all rate segments
    frac_time = (sourceState(phase)%p(sourceOffset)%state(1,of) - &
                 param(instance)%time(interval)) / &
                (param(instance)%time(interval+1) - &
                 param(instance)%time(interval))                                                    ! fractional time within segment
    if (     (frac_time <  0.0_pReal .and. interval == 1) &
        .or. (frac_time >= 1.0_pReal .and. interval == param(instance)%nIntervals) &
        .or. (frac_time >= 0.0_pReal .and. frac_time < 1.0_pReal) ) &
      TDot = param(instance)%heat_rate(interval  ) * (1.0_pReal - frac_time) + &
             param(instance)%heat_rate(interval+1) * frac_time                                      ! interpolate heat rate between segment boundaries...
                                                                                                    ! ...or extrapolate if outside of bounds
  enddo
  dTDot_dT = 0.0
 
end subroutine source_thermal_externalheat_getRateAndItsTangent
 
end module source_thermal_externalheat
