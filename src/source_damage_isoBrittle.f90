!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoBrittle
  use prec
  use debug
  use IO
  use math
  use discretization
  use material
  use config
  use results

  implicit none
  private

  integer,                       dimension(:),           allocatable :: &
    source_damage_isoBrittle_offset, &
    source_damage_isoBrittle_instance

  type, private :: tParameters                                                                      !< container type for internal constitutive parameters
    real(pReal) :: &
      critStrainEnergy, &
      N, &
      aTol
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_damage_isoBrittle_init, &
    source_damage_isoBrittle_deltaState, &
    source_damage_isoBrittle_getRateAndItsTangent, &
    source_damage_isoBrittle_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_init

  integer :: Ninstance,sourceOffset,NofMyPhase,p
  character(len=pStringLen) :: extmsg = ''

  write(6,'(/,a)') ' <<<+-  source_'//SOURCE_DAMAGE_ISOBRITTLE_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_source == SOURCE_DAMAGE_ISOBRITTLE_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(source_damage_isoBrittle_offset  (size(config_phase)), source=0)
  allocate(source_damage_isoBrittle_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    source_damage_isoBrittle_instance(p) = count(phase_source(:,1:p) == SOURCE_DAMAGE_ISOBRITTLE_ID)
    do sourceOffset = 1, phase_Nsources(p)
      if (phase_source(sourceOffset,p) == SOURCE_DAMAGE_ISOBRITTLE_ID) then
        source_damage_isoBrittle_offset(p) = sourceOffset
        exit
      endif
    enddo

    if (all(phase_source(:,p) /= SOURCE_DAMAGE_ISOBRITTLE_ID)) cycle
    associate(prm => param(source_damage_isoBrittle_instance(p)), &
              config => config_phase(p))

    prm%aTol             = config%getFloat('isobrittle_atol',defaultVal = 1.0e-3_pReal)
    prm%N                = config%getFloat('isobrittle_n')
    prm%critStrainEnergy = config%getFloat('isobrittle_criticalstrainenergy')

    ! sanity checks
    if (prm%aTol             <  0.0_pReal) extmsg = trim(extmsg)//' isobrittle_atol'
    if (prm%N                <= 0.0_pReal) extmsg = trim(extmsg)//' isobrittle_n'
    if (prm%critStrainEnergy <= 0.0_pReal) extmsg = trim(extmsg)//' isobrittle_criticalstrainenergy'

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') &
      call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ISOBRITTLE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

    NofMyPhase = count(material_phaseAt==p) * discretization_nIP
    call material_allocateSourceState(p,sourceOffset,NofMyPhase,1,1,1)
    sourceState(p)%p(sourceOffset)%aTolState=prm%aTol

    end associate
  enddo

end subroutine source_damage_isoBrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_deltaState(C, Fe, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal),  intent(in), dimension(3,3) :: &
    Fe
  real(pReal),  intent(in), dimension(6,6) :: &
    C

  integer :: &
    phase, &
    constituent, &
    sourceOffset
  real(pReal), dimension(6) :: &
    strain
  real(pReal) :: &
    strainenergy

  phase = material_phaseAt(ipc,el)                                                                  !< phase ID at ipc,ip,el
  constituent = material_phasememberAt(ipc,ip,el)                                                   !< state array offset for phase ID at ipc,ip,el
  sourceOffset = source_damage_isoBrittle_offset(phase)

  strain = 0.5_pReal*math_sym33to6(matmul(transpose(Fe),Fe)-math_I3)

  associate(prm => param(source_damage_isoBrittle_instance(phase)))
  strainenergy = 2.0_pReal*sum(strain*matmul(C,strain))/prm%critStrainEnergy
  ! ToDo: check strainenergy = 2.0_pReal*dot_product(strain,matmul(C,strain))/param(instance)%critStrainEnergy

  if (strainenergy > sourceState(phase)%p(sourceOffset)%subState0(1,constituent)) then
    sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
      strainenergy - sourceState(phase)%p(sourceOffset)%state(1,constituent)
  else
    sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
      sourceState(phase)%p(sourceOffset)%subState0(1,constituent) - &
      sourceState(phase)%p(sourceOffset)%state(1,constituent)
  endif
  end associate

end subroutine source_damage_isoBrittle_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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

  sourceOffset = source_damage_isoBrittle_offset(phase)

  associate(prm => param(source_damage_isoBrittle_instance(phase)))
  localphiDot = (1.0_pReal - phi)**(prm%n - 1.0_pReal) &
              - phi*sourceState(phase)%p(sourceOffset)%state(1,constituent)
  dLocalphiDot_dPhi = - (prm%n - 1.0_pReal)* (1.0_pReal - phi)**max(0.0_pReal,prm%n - 2.0_pReal) &
                      - sourceState(phase)%p(sourceOffset)%state(1,constituent)
  end associate

end subroutine source_damage_isoBrittle_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_isoBrittle_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_isoBrittle_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('isobrittle_drivingforce')
        call results_writeDataset(group,stt,'tbd','driving force','tbd')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_isoBrittle_results

end module source_damage_isoBrittle
