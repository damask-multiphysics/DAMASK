!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoDuctile
  use prec
  use debug
  use IO
  use discretization
  use material
  use config
  use results

  implicit none
  private

  integer,                       dimension(:),           allocatable :: &
    source_damage_isoDuctile_offset, &                                                              !< which source is my current damage mechanism?
    source_damage_isoDuctile_instance                                                               !< instance of damage source mechanism

  type, private :: tParameters                                                                      !< container type for internal constitutive parameters
    real(pReal) :: &
      critPlasticStrain, &
      N, &
      aTol
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable, private :: param                                    !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_damage_isoDuctile_init, &
    source_damage_isoDuctile_dotState, &
    source_damage_isoDuctile_getRateAndItsTangent, &
    source_damage_isoDuctile_Results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_init

  integer :: Ninstance,sourceOffset,NofMyPhase,p
  character(len=pStringLen) :: extmsg = ''

  write(6,'(/,a)') ' <<<+-  source_'//SOURCE_DAMAGE_ISODUCTILE_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_source == SOURCE_DAMAGE_ISODUCTILE_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(source_damage_isoDuctile_offset  (size(config_phase)), source=0)
  allocate(source_damage_isoDuctile_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    source_damage_isoDuctile_instance(p) = count(phase_source(:,1:p) == SOURCE_DAMAGE_ISODUCTILE_ID)
    do sourceOffset = 1, phase_Nsources(p)
      if (phase_source(sourceOffset,p) == SOURCE_DAMAGE_ISODUCTILE_ID) then
        source_damage_isoDuctile_offset(p) = sourceOffset
        exit
      endif
    enddo

    if (all(phase_source(:,p) /= SOURCE_DAMAGE_ISODUCTILE_ID)) cycle
    associate(prm => param(source_damage_isoDuctile_instance(p)), &
              config => config_phase(p))

    prm%aTol              = config%getFloat('isoductile_atol',defaultVal = 1.0e-3_pReal)
    prm%N                 = config%getFloat('isoductile_ratesensitivity')
    prm%critPlasticStrain = config%getFloat('isoductile_criticalplasticstrain')

    ! sanity checks
    if (prm%aTol              <  0.0_pReal) extmsg = trim(extmsg)//' isoductile_atol'
    if (prm%N                 <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_ratesensitivity'
    if (prm%critPlasticStrain <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_criticalplasticstrain'

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') &
      call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ISODUCTILE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

    NofMyPhase=count(material_phaseAt==p) * discretization_nIP
    call material_allocateSourceState(p,sourceOffset,NofMyPhase,1,1,0)
    sourceState(p)%p(sourceOffset)%aTolState=prm%aTol

    end associate
  enddo

end subroutine source_damage_isoDuctile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_dotState(ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

  integer :: &
    phase, &
    constituent, &
    sourceOffset, &
    damageOffset, &
    homog

  phase = material_phaseAt(ipc,el)
  constituent = material_phasememberAt(ipc,ip,el)
  sourceOffset = source_damage_isoDuctile_offset(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)

  associate(prm => param(source_damage_isoDuctile_instance(phase)))
  sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
    sum(plasticState(phase)%slipRate(:,constituent))/(damage(homog)%p(damageOffset)**prm%N)/prm%critPlasticStrain ! ToDo: abs for slip rate?
  end associate

end subroutine source_damage_isoDuctile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

  integer, intent(in) :: &
    phase, &
    constituent
  real(pReal),  intent(in) :: &
    phi
  real(pReal),  intent(out) :: &
    localphiDot, &
    dLocalphiDot_dPhi

  integer :: &
    sourceOffset

  sourceOffset = source_damage_isoDuctile_offset(phase)

  dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine source_damage_isoDuctile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_isoDuctile_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_isoDuctile_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('isoductile_drivingforce')
        call results_writeDataset(group,stt,'tbd','driving force','tbd')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_isoDuctile_results

end module source_damage_isoDuctile
