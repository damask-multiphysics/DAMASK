!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_anisoBrittle
  use prec
  use debug
  use IO
  use math
  use discretization
  use material
  use config
  use lattice
  use results

  implicit none
  private

  integer,                       dimension(:),           allocatable :: &
    source_damage_anisoBrittle_offset, &                                                            !< which source is my current source mechanism?
    source_damage_anisoBrittle_instance                                                             !< instance of source mechanism

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      sdot_0, &
      n
    real(pReal), dimension(:), allocatable :: &
      critDisp, &
      critLoad
    real(pReal), dimension(:,:,:,:), allocatable :: &
      cleavage_systems
    integer :: &
      sum_N_cl
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_damage_anisoBrittle_init, &
    source_damage_anisoBrittle_dotState, &
    source_damage_anisobrittle_getRateAndItsTangent, &
    source_damage_anisoBrittle_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoBrittle_init

  integer :: Ninstance,sourceOffset,NipcMyPhase,p
  integer, dimension(:), allocatable :: N_cl
  character(len=pStringLen) :: extmsg = ''

  write(6,'(/,a)') ' <<<+-  source_'//SOURCE_DAMAGE_ANISOBRITTLE_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_source == SOURCE_DAMAGE_ANISOBRITTLE_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(source_damage_anisoBrittle_offset  (size(config_phase)), source=0)
  allocate(source_damage_anisoBrittle_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    source_damage_anisoBrittle_instance(p) = count(phase_source(:,1:p) == SOURCE_DAMAGE_ANISOBRITTLE_ID)
    do sourceOffset = 1, phase_Nsources(p)
      if (phase_source(sourceOffset,p) == SOURCE_DAMAGE_ANISOBRITTLE_ID) then
        source_damage_anisoBrittle_offset(p) = sourceOffset
        exit
      endif
    enddo

    if (all(phase_source(:,p) /= SOURCE_DAMAGE_ANISOBRITTLE_ID)) cycle
    associate(prm => param(source_damage_anisoBrittle_instance(p)), &
              config => config_phase(p))

    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

    N_cl = config%getInts('ncleavage',defaultVal=emptyIntArray)
    prm%sum_N_cl = sum(abs(N_cl))

    prm%n         = config%getFloat('anisobrittle_ratesensitivity')
    prm%sdot_0    = config%getFloat('anisobrittle_sdot0')

    prm%critDisp  = config%getFloats('anisobrittle_criticaldisplacement',requiredSize=size(N_cl))
    prm%critLoad  = config%getFloats('anisobrittle_criticalload',        requiredSize=size(N_cl))

    prm%cleavage_systems = lattice_SchmidMatrix_cleavage(N_cl,config%getString('lattice_structure'),&
                                                         config%getFloat('c/a',defaultVal=0.0_pReal))

    ! expand: family => system
    prm%critDisp = math_expand(prm%critDisp,N_cl)
    prm%critLoad = math_expand(prm%critLoad,N_cl)

    ! sanity checks
    if (prm%n            <= 0.0_pReal)  extmsg = trim(extmsg)//' anisobrittle_n'
    if (prm%sdot_0       <= 0.0_pReal)  extmsg = trim(extmsg)//' anisobrittle_sdot0'
    if (any(prm%critLoad <  0.0_pReal)) extmsg = trim(extmsg)//' anisobrittle_critLoad'
    if (any(prm%critDisp <  0.0_pReal)) extmsg = trim(extmsg)//' anisobrittle_critDisp'

    NipcMyPhase = count(material_phaseAt==p) * discretization_nIP
    call material_allocateState(sourceState(p)%p(sourceOffset),NipcMyPhase,1,1,0)
    sourceState(p)%p(sourceOffset)%atol = config%getFloat('anisobrittle_atol',defaultVal=1.0e-3_pReal)
    if(any(sourceState(p)%p(sourceOffset)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisobrittle_atol'

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ANISOBRITTLE_LABEL//')')

enddo

end subroutine source_damage_anisoBrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoBrittle_dotState(S, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal),  intent(in), dimension(3,3) :: &
    S

  integer :: &
    phase, &
    constituent, &
    sourceOffset, &
    damageOffset, &
    homog, &
    i
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit

  phase = material_phaseAt(ipc,el)
  constituent = material_phasememberAt(ipc,ip,el)
  sourceOffset = source_damage_anisoBrittle_offset(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)

  associate(prm => param(source_damage_anisoBrittle_instance(phase)))
  sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = 0.0_pReal
  do i = 1, prm%sum_N_cl
    traction_d    = math_tensordot(S,prm%cleavage_systems(1:3,1:3,1,i))
    traction_t    = math_tensordot(S,prm%cleavage_systems(1:3,1:3,2,i))
    traction_n    = math_tensordot(S,prm%cleavage_systems(1:3,1:3,3,i))

    traction_crit = prm%critLoad(i)*damage(homog)%p(damageOffset)**2.0_pReal

    sourceState(phase)%p(sourceOffset)%dotState(1,constituent) &
    = sourceState(phase)%p(sourceOffset)%dotState(1,constituent) &
    + prm%sdot_0 / prm%critDisp(i) &
      * ((max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**prm%n + &
         (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**prm%n + &
         (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**prm%n)
  enddo
  end associate

end subroutine source_damage_anisoBrittle_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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

  sourceOffset = source_damage_anisoBrittle_offset(phase)

  dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine source_damage_anisoBrittle_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine source_damage_anisoBrittle_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_anisoBrittle_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_anisoBrittle_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('anisobrittle_drivingforce')
        call results_writeDataset(group,stt,'tbd','driving force','tbd')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_anisoBrittle_results

end module source_damage_anisoBrittle
