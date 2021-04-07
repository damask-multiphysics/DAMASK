!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damage) anisoductile

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      q                                                                                             !< damage rate sensitivity
    real(pReal), dimension(:), allocatable :: &
      gamma_crit                                                                                    !< critical plastic strain per slip system
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function anisoductile_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl, &
    sources, &
    src
  integer :: Ninstances,Nmembers,ph
  integer, dimension(:), allocatable :: N_sl
  character(len=pStringLen) :: extmsg = ''


  mySources = source_active('anisoductile')
  if(count(mySources) == 0) return

  print'(/,a)', ' <<<+-  phase:damage:anisoductile init  -+>>>'
  print'(a,i0)', ' # phases: ',count(mySources); flush(IO_STDOUT)


  phases => config_material%get('phase')
  allocate(param(phases%length))

  do ph = 1, phases%length
    if(mySources(ph)) then
      phase => phases%get(ph)
      mech  => phase%get('mechanical')
      pl    => mech%get('plastic')
      sources => phase%get('damage')


        associate(prm  => param(ph))
        src => sources%get(1)

        N_sl           = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
        prm%q          = src%get_asFloat('q')
        prm%gamma_crit = src%get_as1dFloat('gamma_crit',requiredSize=size(N_sl))

        ! expand: family => system
        prm%gamma_crit = math_expand(prm%gamma_crit,N_sl)

#if defined (__GFORTRAN__)
        prm%output = output_as1dString(src)
#else
        prm%output = src%get_as1dString('output',defaultVal=emptyStringArray)
#endif

        ! sanity checks
        if (prm%q              <= 0.0_pReal)  extmsg = trim(extmsg)//' q'
        if (any(prm%gamma_crit <  0.0_pReal)) extmsg = trim(extmsg)//' gamma_crit'

        Nmembers=count(material_phaseID==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,1,0)
        damageState(ph)%atol = src%get_asFloat('anisoDuctile_atol',defaultVal=1.0e-3_pReal)
        if(any(damageState(ph)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' anisoductile_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_anisoDuctile)')
      endif

  enddo


end function anisoductile_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine anisoductile_dotState(ph,me)

  integer, intent(in) :: &
    ph, &
    me


  associate(prm => param(ph))
    damageState(ph)%dotState(1,me) = sum(plasticState(ph)%slipRate(:,me)/(damage_phi(ph,me)**prm%q)/prm%gamma_crit)
  end associate

end subroutine anisoductile_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
module subroutine anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, ph,me)

  integer, intent(in) :: &
    ph, &
    me
  real(pReal),  intent(in) :: &
    phi
  real(pReal),  intent(out) :: &
    localphiDot, &
    dLocalphiDot_dPhi


  dLocalphiDot_dPhi = -damageState(ph)%state(1,me)

  localphiDot = 1.0_pReal &
              + dLocalphiDot_dPhi*phi

end subroutine anisoductile_getRateAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine anisoductile_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state)
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('f_phi')
          call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
      end select
    enddo outputsLoop
  end associate

end subroutine anisoductile_results

end submodule anisoductile
