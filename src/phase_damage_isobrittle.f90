!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:damage) isobrittle

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      W_crit                                                                                        !< critical elastic strain energy
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tIsobrittleState
    real(pREAL), pointer, dimension(:) :: &                                                         !< vectors along Nmembers
      r_W                                                                                           !< ratio between actual and critical strain energy density
  end type tIsobrittleState

  type(tParameters),      allocatable, dimension(:) :: param                                        !< containers of constitutive parameters (len Ninstances)
  type(tIsobrittleState), allocatable, dimension(:) :: &
    deltaState, &
    state

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function isobrittle_init() result(mySources)

  logical, dimension(:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    src
  integer :: Nmembers,ph
  character(len=:), allocatable :: &
    refs, &
    extmsg


  mySources = source_active('isobrittle')
  if (count(mySources) == 0) return

  print'(/,1x,a)', '<<<+-  phase:damage:isobrittle init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(deltaState(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (mySources(ph)) then
      phase => phases%get_dict(ph)
      src => phase%get_dict('damage')

      associate(prm => param(ph), dlt => deltaState(ph), stt => state(ph))

        prm%W_crit = src%get_asReal('G_crit')/src%get_asReal('l_c')

        print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
        refs = config_listReferences(src,indent=3)
        if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(src)
#else
        prm%output = src%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

        ! sanity checks
        if (prm%W_crit <= 0.0_pREAL) extmsg = trim(extmsg)//' W_crit'

        Nmembers = count(material_ID_phase==ph)
        call phase_allocateState(damageState(ph),Nmembers,1,0,1)
        damageState(ph)%atol = src%get_asReal('atol_phi',defaultVal=1.0e-9_pREAL)
        if (any(damageState(ph)%atol < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_phi'

        stt%r_W => damageState(ph)%state(1,:)
        dlt%r_W => damageState(ph)%deltaState(1,:)

      end associate


      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(damage_isobrittle)')
    end if

  end do


end function isobrittle_init


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
module subroutine isobrittle_deltaState(C, Fe, ph,en)

  integer, intent(in) :: ph,en
  real(pREAL),  intent(in), dimension(3,3) :: &
    Fe
  real(pREAL),  intent(in), dimension(6,6) :: &
    C

  real(pREAL), dimension(6) :: &
    epsilon
  real(pREAL) :: &
    r_W


  epsilon = math_33toVoigt6_strain(0.5_pREAL*(matmul(transpose(Fe),Fe)-math_I3))

  associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph))

    r_W = (0.5_pREAL*dot_product(epsilon,matmul(C,epsilon)))/prm%W_crit
    dlt%r_W(en) = merge(r_W - stt%r_W(en), 0.0_pREAL, r_W > stt%r_W(en))

  end associate

end subroutine isobrittle_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine isobrittle_result(phase,group)

  integer,          intent(in) :: phase
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(phase), stt => damageState(phase)%state) ! point to state and output r_W (is scalar, not 1D vector)

    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case ('r_W')
          call result_writeDataset(stt,group,trim(prm%output(o)),'ratio between actual and critical strain energy density','-')
      end select
    end do outputsLoop

  end associate

end subroutine isobrittle_result

end submodule isobrittle
