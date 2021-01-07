!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating isotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule (constitutive:constitutive_damage) source_damage_isoDuctile

  integer,                       dimension(:),           allocatable :: &
    source_damage_isoDuctile_offset, &                                                              !< which source is my current damage mechanism?
    source_damage_isoDuctile_instance                                                               !< instance of damage source mechanism

  type:: tParameters                                                                                !< container type for internal constitutive parameters
    real(pReal) :: &
      gamma_crit, &                                                                                 !< critical plastic strain
      q
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function source_damage_isoDuctile_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length  
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: Ninstances,sourceOffset,Nconstituents,p
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  source_damage_isoDuctile init  -+>>>'

  mySources = source_active('damage_isoDuctile',source_length)
  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstances))
  allocate(source_damage_isoDuctile_offset  (phases%length), source=0)
  allocate(source_damage_isoDuctile_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p) 
    if(count(mySources(:,p)) == 0) cycle
    if(any(mySources(:,p))) source_damage_isoDuctile_instance(p) = count(mySources(:,1:p))
    sources => phase%get('source')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_damage_isoDuctile_offset(p) = sourceOffset
        associate(prm  => param(source_damage_isoDuctile_instance(p)))
        src => sources%get(sourceOffset) 

        prm%q          = src%get_asFloat('q')
        prm%gamma_crit = src%get_asFloat('gamma_crit')

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif
 
        ! sanity checks
        if (prm%q          <= 0.0_pReal) extmsg = trim(extmsg)//' q'
        if (prm%gamma_crit <= 0.0_pReal) extmsg = trim(extmsg)//' gamma_crit'

        Nconstituents=count(material_phaseAt==p) * discretization_nIPs
        call constitutive_allocateState(sourceState(p)%p(sourceOffset),Nconstituents,1,1,0)
        sourceState(p)%p(sourceOffset)%atol = src%get_asFloat('isoDuctile_atol',defaultVal=1.0e-3_pReal)
        if(any(sourceState(p)%p(sourceOffset)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' isoductile_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_isoDuctile)')
      endif
    enddo
  enddo


end function source_damage_isoDuctile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_isoDuctile_dotState(co, ip, el)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

  integer :: &
    phase, &
    constituent, &
    sourceOffset, &
    damageOffset, &
    homog

  phase = material_phaseAt(co,el)
  constituent = material_phasememberAt(co,ip,el)
  sourceOffset = source_damage_isoDuctile_offset(phase)
  homog = material_homogenizationAt(el)
  damageOffset = material_homogenizationMemberAt(ip,el)

  associate(prm => param(source_damage_isoDuctile_instance(phase)))
  sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
    sum(plasticState(phase)%slipRate(:,constituent))/(damage(homog)%p(damageOffset)**prm%q)/prm%gamma_crit
  end associate

end subroutine source_damage_isoDuctile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

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
module subroutine source_damage_isoDuctile_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(source_damage_isoDuctile_instance(phase)), &
            stt => sourceState(phase)%p(source_damage_isoDuctile_offset(phase))%state)
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case ('f_phi')
        call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
    end select
  enddo outputsLoop
  end associate

end subroutine source_damage_isoDuctile_results

end submodule source_damage_isoDuctile
