!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damagee) anisoductile

  integer,                       dimension(:),           allocatable :: &
    source_damage_anisoDuctile_offset, &                                                            !< which source is my current damage mechanism?
    source_damage_anisoDuctile_instance                                                             !< instance of damage source mechanism

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      q                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:), allocatable :: &
      gamma_crit                                                                                    !< critical plastic strain per slip system
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function anisoductile_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl, &
    sources, &
    src
  integer :: Ninstances,sourceOffset,Nconstituents,p
  integer, dimension(:), allocatable :: N_sl
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  phase:damage:anisoductile init  -+>>>'

  mySources = source_active('damage_anisoDuctile',source_length)
  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstances))
  allocate(source_damage_anisoDuctile_offset  (phases%length), source=0)
  allocate(source_damage_anisoDuctile_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p)
    if(any(mySources(:,p))) source_damage_anisoDuctile_instance(p) = count(mySources(:,1:p))
    if(count(mySources(:,p)) == 0) cycle
    mech  => phase%get('mechanics')
    pl    => mech%get('plasticity')
    sources => phase%get('source')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_damage_anisoDuctile_offset(p) = sourceOffset
        associate(prm  => param(source_damage_anisoDuctile_instance(p)))
        src => sources%get(sourceOffset)

        N_sl           = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
        prm%q          = src%get_asFloat('q')
        prm%gamma_crit = src%get_asFloats('gamma_crit',requiredSize=size(N_sl))

        ! expand: family => system
        prm%gamma_crit = math_expand(prm%gamma_crit,N_sl)

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif

        ! sanity checks
        if (prm%q              <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (any(prm%gamma_crit <  0.0_pReal)) extmsg = trim(extmsg)//' gamma_crit'

        Nconstituents=count(material_phaseAt==p) * discretization_nIPs
        call phase_allocateState(damageState(p)%p(sourceOffset),Nconstituents,1,1,0)
        damageState(p)%p(sourceOffset)%atol = src%get_asFloat('anisoDuctile_atol',defaultVal=1.0e-3_pReal)
        if(any(damageState(p)%p(sourceOffset)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisoductile_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoDuctile)')
      endif
    enddo
  enddo


end function anisoductile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine anisoductile_dotState(co, ip, el)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

  integer :: &
    ph, &
    me, &
    sourceOffset, &
    damageOffset, &
    homog

  ph = material_phaseAt(co,el)
  me = material_phasememberAt(co,ip,el)
  sourceOffset = source_damage_anisoDuctile_offset(ph)
  homog = material_homogenizationAt(el)
  damageOffset = material_homogenizationMemberAt(ip,el)

  associate(prm => param(source_damage_anisoDuctile_instance(ph)))
  damageState(ph)%p(sourceOffset)%dotState(1,me) &
    = sum(plasticState(ph)%slipRate(:,me)/(damagestate_h(homog)%state(1,damageOffset)**prm%q)/prm%gamma_crit)
  end associate

end subroutine anisoductile_dotState


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

  dLocalphiDot_dPhi = -damageState(phase)%p(sourceOffset)%state(1,constituent)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine source_damage_anisoDuctile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine anisoductile_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_anisoDuctile_instance(phase)), &
            stt => damageState(phase)%p(source_damage_anisoDuctile_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('f_phi')
        call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
    end select
  enddo outputsLoop
  end associate

end subroutine anisoductile_results

end submodule anisoductile
