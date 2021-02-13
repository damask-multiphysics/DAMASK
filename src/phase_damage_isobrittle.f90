!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damagee) isobrittle

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
module function isobrittle_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: Ninstances,sourceOffset,Nconstituents,p
  character(len=pStringLen) :: extmsg = ''

  print'(/,a)', ' <<<+-  phase:damage:isobrittle init  -+>>>'

  mySources = source_active('isobrittle',source_length)
  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(phases%length))

  do p = 1, phases%length
    phase => phases%get(p)
    if(count(mySources(:,p)) == 0) cycle
    sources => phase%get('damage')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        associate(prm  => param(p))
        src => sources%get(sourceOffset)

        prm%W_crit = src%get_asFloat('W_crit')

#if defined (__GFORTRAN__)
        prm%output = output_asStrings(src)
#else
        prm%output = src%get_asStrings('output',defaultVal=emptyStringArray)
#endif

        ! sanity checks
        if (prm%W_crit <= 0.0_pReal) extmsg = trim(extmsg)//' W_crit'

        Nconstituents = count(material_phaseAt==p) * discretization_nIPs
        call phase_allocateState(damageState(p),Nconstituents,1,1,1)
        damageState(p)%atol = src%get_asFloat('isoBrittle_atol',defaultVal=1.0e-3_pReal)
        if(any(damageState(p)%atol < 0.0_pReal)) extmsg = trim(extmsg)//' isobrittle_atol'

        end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
        if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_isoBrittle)')
      endif
    enddo
  enddo


end function isobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_isoBrittle_deltaState(C, Fe, ph,me)

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

    if (strainenergy > damageState(ph)%subState0(1,me)) then
      damageState(ph)%deltaState(1,me) = strainenergy - damageState(ph)%state(1,me)
    else
      damageState(ph)%deltaState(1,me) = damageState(ph)%subState0(1,me) - damageState(ph)%state(1,me)
    endif
  end associate

end subroutine source_damage_isoBrittle_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
module subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

  integer, intent(in) :: &
    phase, &
    constituent
  real(pReal),  intent(in) :: &
    phi
  real(pReal),  intent(out) :: &
    localphiDot, &
    dLocalphiDot_dPhi


  associate(prm => param(phase))
    localphiDot = 1.0_pReal &
                - phi*damageState(phase)%state(1,constituent)
    dLocalphiDot_dPhi = - damageState(phase)%state(1,constituent)
  end associate

end subroutine source_damage_isoBrittle_getRateAndItsTangent


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
          call results_writeDataset(group,stt,trim(prm%output(o)),'driving force','J/m³')
      end select
    enddo outputsLoop
  end associate

end subroutine isobrittle_results

end submodule isobrittle
