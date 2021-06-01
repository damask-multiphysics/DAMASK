!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damage) isobrittle

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      W_crit                                                                                        !< critical elastic strain energy
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstances)

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function isobrittle_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: Nmembers,ph
  character(len=pStringLen) :: extmsg = ''


  mySources = source_active('isobrittle')
  if(count(mySources) == 0) return

  print'(/,a)', ' <<<+-  phase:damage:isobrittle init  -+>>>'
  print'(a,i0)', ' # phases: ',count(mySources); flush(IO_STDOUT)


  phases => config_material%get('phase')
  allocate(param(phases%length))

  do ph = 1, phases%length
    if(mySources(ph)) then
    phase => phases%get(ph)
    sources => phase%get('damage')

        associate(prm  => param(ph))
        src => sources%get(1)

        prm%W_crit = src%get_asFloat('W_crit')

#if defined (__GFORTRAN__)
        prm%output = output_as1dString(src)
#else
        prm%output = src%get_as1dString('output',defaultVal=emptyStringArray)
#endif

        ! sanity checks
        if (prm%W_crit <= 0.0_pReal) extmsg = trim(extmsg)//' W_crit'

        Nmembers = count(material_phaseID==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,1,1)
        damageState(ph)%atol = src%get_asFloat('isobrittle_atol',defaultVal=1.0e-3_pReal)
        if(any(damageState(ph)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' isobrittle_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_isobrittle)')
      endif

  enddo


end function isobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine isobrittle_deltaState(C, Fe, ph,me)

  integer, intent(in) :: ph,me
  real(pReal),  intent(in), dimension(3,3) :: &
    Fe
  real(pReal),  intent(in), dimension(6,6) :: &
    C

  real(pReal), dimension(6) :: &
    strain
  real(pReal) :: &
    strainenergy


  strain = 0.5_pReal*math_sym33to6(matmul(transpose(Fe),Fe)-math_I3)

  associate(prm => param(ph))
    strainenergy = 2.0_pReal*sum(strain*matmul(C,strain))/prm%W_crit
    ! ToDo: check strainenergy = 2.0_pReal*dot_product(strain,matmul(C,strain))/prm%W_crit

    damageState(ph)%deltaState(1,me) = merge(strainenergy - damageState(ph)%state(1,me), &
                                             damageState(ph)%subState0(1,me) - damageState(ph)%state(1,me), &
                                             strainenergy > damageState(ph)%subState0(1,me))
  end associate

end subroutine isobrittle_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine isobrittle_results(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), &
            stt => damageState(phase)%state)
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('f_phi')
          call results_writeDataset(stt,group,trim(prm%output(o)),'driving force','J/m³')
      end select
    enddo outputsLoop
  end associate

end subroutine isobrittle_results

end submodule isobrittle
