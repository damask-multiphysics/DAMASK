!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(constitutive:constitutive_damage)  source_damage_anisoDuctile

  integer,                       dimension(:),           allocatable :: &
    source_damage_anisoDuctile_offset, &                                                            !< which source is my current damage mechanism?
    source_damage_anisoDuctile_instance                                                             !< instance of damage source mechanism

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      n                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:), allocatable :: &
      critPlasticStrain                                                                             !< critical plastic strain per slip system
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function source_damage_anisoDuctile_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length  
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    pl, &
    sources, &
    src
  integer :: Ninstance,sourceOffset,NipcMyPhase,p
  integer, dimension(:), allocatable :: N_sl
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  source_damage_anisoDuctile init  -+>>>'

  mySources = source_active('damage_anisoDuctile',source_length)
  Ninstance = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstance; flush(IO_STDOUT)
  if(Ninstance == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstance))
  allocate(source_damage_anisoDuctile_offset  (phases%length), source=0)
  allocate(source_damage_anisoDuctile_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p) 
    if(any(mySources(:,p))) source_damage_anisoDuctile_instance(p) = count(mySources(:,1:p))
    if(count(mySources(:,p)) == 0) cycle
    sources => phase%get('source')
    pl => phase%get('plasticity')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_damage_anisoDuctile_offset(p) = sourceOffset
        associate(prm  => param(source_damage_anisoDuctile_instance(p)))
        src => sources%get(sourceOffset) 

        N_sl = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
        prm%n                 = src%get_asFloat('q')
        prm%critPlasticStrain = src%get_asFloats('gamma_crit',requiredSize=size(N_sl))

        ! expand: family => system
        prm%critPlasticStrain = math_expand(prm%critPlasticStrain,N_sl)

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif
 
        ! sanity checks
        if (prm%n                     <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (any(prm%critPlasticStrain <  0.0_pReal)) extmsg = trim(extmsg)//' gamma_crit'

        NipcMyPhase=count(material_phaseAt==p) * discretization_nIP
        call constitutive_allocateState(sourceState(p)%p(sourceOffset),NipcMyPhase,1,1,0)
        sourceState(p)%p(sourceOffset)%atol = src%get_asFloat('anisoDuctile_atol',defaultVal=1.0e-3_pReal)
        if(any(sourceState(p)%p(sourceOffset)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisoductile_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoDuctile)')
      endif
    enddo
  enddo


end function source_damage_anisoDuctile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_anisoDuctile_dotState(ipc, ip, el)

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
  sourceOffset = source_damage_anisoDuctile_offset(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)

  associate(prm => param(source_damage_anisoDuctile_instance(phase)))
  sourceState(phase)%p(sourceOffset)%dotState(1,constituent) &
    = sum(plasticState(phase)%slipRate(:,constituent)/(damage(homog)%p(damageOffset)**prm%n)/prm%critPlasticStrain)
  end associate

end subroutine source_damage_anisoDuctile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_anisoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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

  sourceOffset = source_damage_anisoDuctile_offset(phase)

  dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine source_damage_anisoDuctile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_anisoDuctile_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_anisoDuctile_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_anisoDuctile_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('f_phi')
        call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_anisoDuctile_results

end submodule source_damage_anisoDuctile
