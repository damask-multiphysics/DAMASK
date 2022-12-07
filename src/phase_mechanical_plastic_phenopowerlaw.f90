!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) phenopowerlaw

  type :: tParameters
    real(pReal) :: &
      dot_gamma_0_sl = 1.0_pReal, &                                                                 !< reference shear strain rate for slip
      dot_gamma_0_tw = 1.0_pReal, &                                                                 !< reference shear strain rate for twin
      n_sl           = 1.0_pReal, &                                                                 !< stress exponent for slip
      n_tw           = 1.0_pReal, &                                                                 !< stress exponent for twin
      f_sat_sl_tw    = 1.0_pReal, &                                                                 !< push-up factor for slip saturation due to twinning
      c_1            = 1.0_pReal, &
      c_2            = 1.0_pReal, &
      c_3            = 1.0_pReal, &
      c_4            = 1.0_pReal, &
      h_0_sl_sl      = 1.0_pReal, &                                                                 !< reference hardening slip - slip
      h_0_tw_sl      = 1.0_pReal, &                                                                 !< reference hardening twin - slip
      h_0_tw_tw      = 1.0_pReal, &                                                                 !< reference hardening twin - twin
      a_sl           = 1.0_pReal
    real(pReal),               allocatable, dimension(:) :: &
      xi_inf_sl, &                                                                                  !< maximum critical shear stress for slip
      h_int, &                                                                                      !< per family hardening activity (optional)
      gamma_char                                                                                    !< characteristic shear for twins
    real(pReal),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      h_sl_tw, &                                                                                    !< slip resistance from twin activity
      h_tw_sl, &                                                                                    !< twin resistance from slip activity
      h_tw_tw                                                                                       !< twin resistance from twin activity
    real(pReal),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw                                                                                      !< total number of active twin systems
    logical :: &
      nonSchmidActive = .false.
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl, &
      systems_tw
  end type tParameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      xi_sl, &
      xi_tw, &
      gamma_sl, &
      gamma_tw
  end type tIndexDotState

  type :: tPhenopowerlawState
    real(pReal), pointer, dimension(:,:) :: &
      xi_sl, &
      xi_tw, &
      gamma_sl, &
      gamma_tw
  end type tPhenopowerlawState

!--------------------------------------------------------------------------------------------------
! containers for parameters, dot state index,  and state
  type(tParameters),         allocatable, dimension(:) :: param
  type(tIndexDotState),      allocatable, dimension(:) :: indexDotState
  type(tPhenopowerlawState), allocatable, dimension(:) :: state

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_phenopowerlaw_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl, N_tw
  real(pReal), dimension(:), allocatable :: &
    xi_0_sl, &                                                                                      !< initial critical shear stress for slip
    xi_0_tw, &                                                                                      !< initial critical shear stress for twin
    a                                                                                               !< non-Schmid coefficients
  character(len=pStringLen) :: &
    extmsg = ''
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('phenopowerlaw')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:phenopowerlaw init  -+>>>'
  print'(/,a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), stt => state(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech => phase%get_dict('mechanical')
    pl => mech%get_dict('plastic')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = lattice_labels_slip(N_sl,phase_lattice(ph))
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal=emptyRealArray)
        if (size(a) > 0) prm%nonSchmidActive = .true.
        prm%P_nS_pos = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%P_nS_neg = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%P_nS_pos = prm%P_sl
        prm%P_nS_neg = prm%P_sl
      end if
      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,pl%get_as1dFloat('h_sl-sl'),phase_lattice(ph))

      xi_0_sl             = pl%get_as1dFloat('xi_0_sl',   requiredSize=size(N_sl))
      prm%xi_inf_sl       = pl%get_as1dFloat('xi_inf_sl', requiredSize=size(N_sl))
      prm%h_int           = pl%get_as1dFloat('h_int',     requiredSize=size(N_sl), &
                                            defaultVal=[(0.0_pReal,i=1,size(N_sl))])

      prm%dot_gamma_0_sl  = pl%get_asFloat('dot_gamma_0_sl')
      prm%n_sl            = pl%get_asFloat('n_sl')
      prm%a_sl            = pl%get_asFloat('a_sl')
      prm%h_0_sl_sl       = pl%get_asFloat('h_0_sl-sl')

      ! expand: family => system
      xi_0_sl             = math_expand(xi_0_sl,      N_sl)
      prm%xi_inf_sl       = math_expand(prm%xi_inf_sl,N_sl)
      prm%h_int           = math_expand(prm%h_int,    N_sl)

      ! sanity checks
      if (    prm%dot_gamma_0_sl  <= 0.0_pReal)      extmsg = trim(extmsg)//' dot_gamma_0_sl'
      if (    prm%a_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' a_sl'
      if (    prm%n_sl            <= 0.0_pReal)      extmsg = trim(extmsg)//' n_sl'
      if (any(xi_0_sl             <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_0_sl'
      if (any(prm%xi_inf_sl       <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_inf_sl'

    else slipActive
      xi_0_sl = emptyRealArray
      allocate(prm%xi_inf_sl,prm%h_int,source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    N_tw         = pl%get_as1dInt('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%systems_tw = lattice_labels_twin(N_tw,phase_lattice(ph))
      prm%P_tw     = lattice_SchmidMatrix_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%h_tw_tw  = lattice_interaction_TwinByTwin(N_tw,pl%get_as1dFloat('h_tw-tw'),phase_lattice(ph))
      prm%gamma_char = lattice_characteristicShear_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))

      xi_0_tw             = pl%get_as1dFloat('xi_0_tw',requiredSize=size(N_tw))

      prm%c_1             = pl%get_asFloat('c_1',defaultVal=0.0_pReal)
      prm%c_2             = pl%get_asFloat('c_2',defaultVal=1.0_pReal)
      prm%c_3             = pl%get_asFloat('c_3',defaultVal=0.0_pReal)
      prm%c_4             = pl%get_asFloat('c_4',defaultVal=0.0_pReal)
      prm%dot_gamma_0_tw  = pl%get_asFloat('dot_gamma_0_tw')
      prm%n_tw            = pl%get_asFloat('n_tw')
      prm%f_sat_sl_tw     = pl%get_asFloat('f_sat_sl-tw')
      prm%h_0_tw_tw       = pl%get_asFloat('h_0_tw-tw')

      ! expand: family => system
      xi_0_tw       = math_expand(xi_0_tw,N_tw)

      ! sanity checks
      if (prm%dot_gamma_0_tw <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_gamma_0_tw'
      if (prm%n_tw           <= 0.0_pReal)  extmsg = trim(extmsg)//' n_tw'

    else twinActive
      xi_0_tw = emptyRealArray
      allocate(prm%gamma_char,source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
    end if twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h_0_tw_sl  = pl%get_asFloat('h_0_tw-sl')
      prm%h_sl_tw    = lattice_interaction_SlipByTwin(N_sl,N_tw,pl%get_as1dFloat('h_sl-tw'), &
                                                      phase_lattice(ph))
      prm%h_tw_sl    = lattice_interaction_TwinBySlip(N_tw,N_sl,pl%get_as1dFloat('h_tw-sl'), &
                                                      phase_lattice(ph))
    else slipAndTwinActive
      allocate(prm%h_sl_tw(prm%sum_N_sl,prm%sum_N_tw))                                              ! at least one dimension is 0
      allocate(prm%h_tw_sl(prm%sum_N_tw,prm%sum_N_sl))                                              ! at least one dimension is 0
      prm%h_0_tw_sl = 0.0_pReal
    end if slipAndTwinActive

!--------------------------------------------------------------------------------------------------
!  output pararameters

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%xi_sl = [startIndex,endIndex]
    stt%xi_sl => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_sl =  spread(xi_0_sl, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%xi_tw = [startIndex,endIndex]
    stt%xi_tw => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_tw =  spread(xi_0_tw, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%gamma_tw = [startIndex,endIndex]
    stt%gamma_tw => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!> @details asummes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
pure module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,               intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos,dot_gamma_sl_neg, &
    ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw,ddot_gamma_dtau_tw

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph))

  call kinetics_sl(Mp,ph,en,dot_gamma_sl_pos,dot_gamma_sl_neg,ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg)
  slipSystems: do i = 1, prm%sum_N_sl
    Lp = Lp + (dot_gamma_sl_pos(i)+dot_gamma_sl_neg(i))*prm%P_sl(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_sl_pos(i) * prm%P_sl(k,l,i) * prm%P_nS_pos(m,n,i) &
                       + ddot_gamma_dtau_sl_neg(i) * prm%P_sl(k,l,i) * prm%P_nS_neg(m,n,i)
  end do slipSystems

  call kinetics_tw(Mp,ph,en,dot_gamma_tw,ddot_gamma_dtau_tw)
  twinSystems: do i = 1, prm%sum_N_tw
    Lp = Lp + dot_gamma_tw(i)*prm%P_tw(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_tw(i)*prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
  end do twinSystems

  end associate

end subroutine phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function phenopowerlaw_dotState(Mp,ph,en) result(dotState)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pReal) :: &
    xi_sl_sat_offset,&
    sumF
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos,dot_gamma_sl_neg, &
    left_SlipSlip

  associate(prm => param(ph), stt => state(ph), &
            dot_xi_sl => dotState(indexDotState(ph)%xi_sl(1):indexDotState(ph)%xi_sl(2)), &
            dot_xi_tw => dotState(indexDotState(ph)%xi_tw(1):indexDotState(ph)%xi_tw(2)), &
            dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)), &
            dot_gamma_tw => dotState(indexDotState(ph)%gamma_tw(1):indexDotState(ph)%gamma_tw(2)))

    call kinetics_sl(Mp,ph,en, dot_gamma_sl_pos,dot_gamma_sl_neg)
    dot_gamma_sl = abs(dot_gamma_sl_pos+dot_gamma_sl_neg)
    call kinetics_tw(Mp,ph,en, dot_gamma_tw)
    sumF = sum(stt%gamma_tw(:,en)/prm%gamma_char)

    xi_sl_sat_offset = prm%f_sat_sl_tw*sqrt(sumF)
    left_SlipSlip = sign(abs(1.0_pReal-stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))**prm%a_sl, &
                             1.0_pReal-stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))

    dot_xi_sl = prm%h_0_sl_sl * (1.0_pReal + prm%c_1 * sumF**prm%c_2) * (1.0_pReal + prm%h_int) &
                * left_SlipSlip * matmul(prm%h_sl_sl,dot_gamma_sl) &
              + matmul(prm%h_sl_tw,dot_gamma_tw)

    dot_xi_tw = prm%h_0_tw_sl * sum(stt%gamma_sl(:,en))**prm%c_3 &
                * matmul(prm%h_tw_sl,dot_gamma_sl) &
              + prm%h_0_tw_tw * sumF**prm%c_4 * matmul(prm%h_tw_tw,dot_gamma_tw)

  end associate

end function phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_phenopowerlaw_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('xi_sl')
          call results_writeDataset(stt%xi_sl,group,trim(prm%output(ou)), &
                                    'resistance against plastic slip','Pa',prm%systems_sl)
        case('gamma_sl')
          call results_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                    'plastic shear','1',prm%systems_sl)

        case('xi_tw')
          call results_writeDataset(stt%xi_tw,group,trim(prm%output(ou)), &
                                    'resistance against twinning','Pa',prm%systems_tw)
        case('gamma_tw')
          call results_writeDataset(stt%gamma_tw,group,trim(prm%output(ou)), &
                                    'twinning shear','1',prm%systems_tw)

      end select

    end do

  end associate

end subroutine plastic_phenopowerlaw_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_sl(Mp,ph,en, &
                            dot_gamma_sl_pos,dot_gamma_sl_neg,ddot_gamma_dtau_sl_pos,ddot_gamma_dtau_sl_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl_pos, &
    dot_gamma_sl_neg
  real(pReal),                  intent(out), optional, dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau_sl_pos, &
    ddot_gamma_dtau_sl_neg

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau_sl_pos, &
    tau_sl_neg
  integer :: i

  associate(prm => param(ph), stt => state(ph))

    do i = 1, prm%sum_N_sl
      tau_sl_pos(i) =       math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,i))
      tau_sl_neg(i) = merge(math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,i)), &
                            0.0_pReal, prm%nonSchmidActive)
    end do

    where(dNeq0(tau_sl_pos))
      dot_gamma_sl_pos = prm%dot_gamma_0_sl * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &     ! 1/2 if non-Schmid active
                       * sign(abs(tau_sl_pos/stt%xi_sl(:,en))**prm%n_sl,  tau_sl_pos)
    else where
      dot_gamma_sl_pos = 0.0_pReal
    end where

    where(dNeq0(tau_sl_neg))
      dot_gamma_sl_neg = prm%dot_gamma_0_sl * 0.5_pReal &                                           ! only used if non-Schmid active, always 1/2
                       * sign(abs(tau_sl_neg/stt%xi_sl(:,en))**prm%n_sl,  tau_sl_neg)
    else where
      dot_gamma_sl_neg = 0.0_pReal
    end where

    if (present(ddot_gamma_dtau_sl_pos)) then
      where(dNeq0(dot_gamma_sl_pos))
        ddot_gamma_dtau_sl_pos = dot_gamma_sl_pos*prm%n_sl/tau_sl_pos
      else where
        ddot_gamma_dtau_sl_pos = 0.0_pReal
      end where
    end if
    if (present(ddot_gamma_dtau_sl_neg)) then
      where(dNeq0(dot_gamma_sl_neg))
        ddot_gamma_dtau_sl_neg = dot_gamma_sl_neg*prm%n_sl/tau_sl_neg
      else where
        ddot_gamma_dtau_sl_neg = 0.0_pReal
      end where
    end if

  end associate

end subroutine kinetics_sl


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress. Twinning is assumed to take place only in untwinned volume.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_tw(Mp,ph,en,&
                            dot_gamma_tw,ddot_gamma_dtau_tw)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal), dimension(param(ph)%sum_N_tw), intent(out) :: &
    dot_gamma_tw
  real(pReal), dimension(param(ph)%sum_N_tw), intent(out), optional :: &
    ddot_gamma_dtau_tw

  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    tau_tw
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    tau_tw = [(math_tensordot(Mp,prm%P_tw(1:3,1:3,i)),i=1,prm%sum_N_tw)]

    where(tau_tw > 0.0_pReal)
      dot_gamma_tw = (1.0_pReal-sum(stt%gamma_tw(:,en)/prm%gamma_char)) &                           ! only twin in untwinned volume fraction
                   * prm%dot_gamma_0_tw*(abs(tau_tw)/stt%xi_tw(:,en))**prm%n_tw
    else where
      dot_gamma_tw = 0.0_pReal
    end where

    if (present(ddot_gamma_dtau_tw)) then
      where(dNeq0(dot_gamma_tw))
        ddot_gamma_dtau_tw = dot_gamma_tw*prm%n_tw/tau_tw
      else where
        ddot_gamma_dtau_tw = 0.0_pReal
      end where
    end if

  end associate

end subroutine kinetics_tw

end submodule phenopowerlaw
