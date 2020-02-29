!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!> @brief material subroutine for variable heat source
!--------------------------------------------------------------------------------------------------
module source_thermal_externalheat
  use prec
  use debug
  use discretization
  use material
  use config

  implicit none
  private

  integer,           dimension(:),   allocatable :: &
    source_thermal_externalheat_offset, &                                                           !< which source is my current thermal dissipation mechanism?
    source_thermal_externalheat_instance                                                            !< instance of thermal dissipation source mechanism

  type :: tParameters                                                                              !< container type for internal constitutive parameters
    real(pReal), dimension(:), allocatable :: &
      time, &
      heat_rate
    integer :: &
     nIntervals
  end type tParameters

  type(tParameters), dimension(:), allocatable  :: param                                            !< containers of constitutive parameters (len Ninstance)


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

  integer :: Ninstance,sourceOffset,NofMyPhase,p

  write(6,'(/,a)') ' <<<+-  source_'//SOURCE_thermal_externalheat_label//' init  -+>>>'; flush(6)

  Ninstance = count(phase_source == SOURCE_thermal_externalheat_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(source_thermal_externalheat_offset  (size(config_phase)), source=0)
  allocate(source_thermal_externalheat_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    source_thermal_externalheat_instance(p) = count(phase_source(:,1:p) == SOURCE_thermal_externalheat_ID)
    do sourceOffset = 1, phase_Nsources(p)
      if (phase_source(sourceOffset,p) == SOURCE_thermal_externalheat_ID) then
        source_thermal_externalheat_offset(p) = sourceOffset
        exit 
      endif
    enddo

    if (all(phase_source(:,p) /= SOURCE_thermal_externalheat_ID)) cycle
    associate(prm => param(source_thermal_externalheat_instance(p)), &
              config => config_phase(p))

    prm%time       = config%getFloats('externalheat_time')
    prm%nIntervals = size(prm%time) - 1

    prm%heat_rate = config%getFloats('externalheat_rate',requiredSize = size(prm%time))

    NofMyPhase = count(material_phaseAt==p) * discretization_nIP
    call material_allocateSourceState(p,sourceOffset,NofMyPhase,1,1,0)

    end associate
  enddo

end subroutine source_thermal_externalheat_init


!--------------------------------------------------------------------------------------------------
!> @brief rate of change of state
!> @details state only contains current time to linearly interpolate given heat powers
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_dotState(phase, of)

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

  integer, intent(in) :: &
    phase, &
    of
  real(pReal),  intent(out) :: &
    TDot, &
    dTDot_dT

  integer :: &
    sourceOffset, interval
  real(pReal) :: &
    frac_time

  sourceOffset = source_thermal_externalheat_offset(phase)

  associate(prm => param(source_thermal_externalheat_instance(phase)))
  do interval = 1, prm%nIntervals                                                       ! scan through all rate segments
    frac_time = (sourceState(phase)%p(sourceOffset)%state(1,of) - prm%time(interval)) &
              / (prm%time(interval+1) - prm%time(interval))                             ! fractional time within segment
    if (     (frac_time <  0.0_pReal .and. interval == 1) &
        .or. (frac_time >= 1.0_pReal .and. interval == prm%nIntervals) &
        .or. (frac_time >= 0.0_pReal .and. frac_time < 1.0_pReal) ) &
      TDot = prm%heat_rate(interval  ) * (1.0_pReal - frac_time) + &
             prm%heat_rate(interval+1) * frac_time                                      ! interpolate heat rate between segment boundaries...
                                                                                        ! ...or extrapolate if outside of bounds
  enddo
  dTDot_dT = 0.0
  end associate

end subroutine source_thermal_externalheat_getRateAndItsTangent

end module source_thermal_externalheat
