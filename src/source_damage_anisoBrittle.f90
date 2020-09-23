!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule (constitutive:constitutive_damage) source_damage_anisoBrittle

  integer,                       dimension(:),           allocatable :: &
    source_damage_anisoBrittle_offset, &                                                            !< which source is my current source mechanism?
    source_damage_anisoBrittle_instance                                                             !< instance of source mechanism

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      sdot_0, &                                                                                     !< opening rate of cleavage planes
      n                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:), allocatable :: &
      critDisp, &                                                                                   !< critical displacement 
      critLoad                                                                                      !< critical load
    real(pReal), dimension(:,:,:,:), allocatable :: &
      cleavage_systems
    integer :: &
      sum_N_cl                                                                                      !< total number of cleavage planes
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function source_damage_anisoBrittle_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length  
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: Ninstance,sourceOffset,NipcMyPhase,p
  integer, dimension(:), allocatable :: N_cl
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  source_damage_anisoBrittle init  -+>>>'

  mySources = source_active('damage_anisoBrittle',source_length)
  Ninstance = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstance; flush(IO_STDOUT)
  if(Ninstance == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstance))
  allocate(source_damage_anisoBrittle_offset  (phases%length), source=0)
  allocate(source_damage_anisoBrittle_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p) 
    if(any(mySources(:,p))) source_damage_anisoBrittle_instance(p) = count(mySources(:,1:p))
    if(count(mySources(:,p)) == 0) cycle
    sources => phase%get('source')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_damage_anisoBrittle_offset(p) = sourceOffset
        associate(prm  => param(source_damage_anisoBrittle_instance(p)))
        src => sources%get(sourceOffset) 
  
        N_cl = src%get_asInts('N_cl',defaultVal=emptyIntArray)
        prm%sum_N_cl = sum(abs(N_cl))
  
        prm%n         = src%get_asFloat('q')
        prm%sdot_0    = src%get_asFloat('dot_o')
  
        prm%critDisp  = src%get_asFloats('s_crit',  requiredSize=size(N_cl))
        prm%critLoad  = src%get_asFloats('g_crit',  requiredSize=size(N_cl))
  
        prm%cleavage_systems = lattice_SchmidMatrix_cleavage(N_cl,phase%get_asString('lattice'),&
                                                             phase%get_asFloat('c/a',defaultVal=0.0_pReal))
  
        ! expand: family => system
        prm%critDisp = math_expand(prm%critDisp,N_cl)
        prm%critLoad = math_expand(prm%critLoad,N_cl)

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif
 
          ! sanity checks
        if (prm%n            <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (prm%sdot_0       <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_o'
        if (any(prm%critLoad <  0.0_pReal)) extmsg = trim(extmsg)//' g_crit'
        if (any(prm%critDisp <  0.0_pReal)) extmsg = trim(extmsg)//' s_crit'

        NipcMyPhase = count(material_phaseAt==p) * discretization_nIP
        call constitutive_allocateState(sourceState(p)%p(sourceOffset),NipcMyPhase,1,1,0)
        sourceState(p)%p(sourceOffset)%atol = src%get_asFloat('anisobrittle_atol',defaultVal=1.0e-3_pReal)
        if(any(sourceState(p)%p(sourceOffset)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisobrittle_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoBrittle)')
      endif
    enddo
  enddo

end function source_damage_anisoBrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_anisoBrittle_dotState(S, ipc, ip, el)

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
module subroutine source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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
module subroutine source_damage_anisoBrittle_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_anisoBrittle_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_anisoBrittle_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('f_phi')
        call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_anisoBrittle_results

end submodule source_damage_anisoBrittle
