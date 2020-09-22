!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
submodule(constitutive:constitutive_plastic) plastic_phenopowerlaw

  type :: tParameters
    real(pReal) :: &
      gdot0_slip  = 1.0_pReal, &                                                                    !< reference shear strain rate for slip
      gdot0_twin  = 1.0_pReal, &                                                                    !< reference shear strain rate for twin
      n_slip      = 1.0_pReal, &                                                                    !< stress exponent for slip
      n_twin      = 1.0_pReal, &                                                                    !< stress exponent for twin
      spr         = 1.0_pReal, &                                                                    !< push-up factor for slip saturation due to twinning
      c_1         = 1.0_pReal, &
      c_2         = 1.0_pReal, &
      c_3         = 1.0_pReal, &
      c_4         = 1.0_pReal, &
      h0_SlipSlip = 1.0_pReal, &                                                                    !< reference hardening slip - slip
      h0_TwinSlip = 1.0_pReal, &                                                                    !< reference hardening twin - slip
      h0_TwinTwin = 1.0_pReal, &                                                                    !< reference hardening twin - twin
      a_slip      = 1.0_pReal
    real(pReal),               allocatable, dimension(:) :: &
      xi_slip_sat, &                                                                                !< maximum critical shear stress for slip
      H_int, &                                                                                      !< per family hardening activity (optional)
      gamma_twin_char                                                                               !< characteristic shear for twins
    real(pReal),               allocatable, dimension(:,:) :: &
      interaction_SlipSlip, &                                                                       !< slip resistance from slip activity
      interaction_SlipTwin, &                                                                       !< slip resistance from twin activity
      interaction_TwinSlip, &                                                                       !< twin resistance from slip activity
      interaction_TwinTwin                                                                          !< twin resistance from twin activity
    real(pReal),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw                                                                                      !< total number of active twin systems
    logical :: &
      nonSchmidActive = .false.
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tPhenopowerlawState
    real(pReal), pointer, dimension(:,:) :: &
      xi_slip, &
      xi_twin, &
      gamma_slip, &
      gamma_twin
  end type tPhenopowerlawState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),         allocatable, dimension(:) :: param
  type(tPhenopowerlawState), allocatable, dimension(:) :: &
    dotState, &
    state

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_phenopowerlaw_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    Ninstance, &
    p, i, &
    NipcMyPhase, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl, N_tw
  real(pReal), dimension(:), allocatable :: &
    xi_slip_0, &                                                                                    !< initial critical shear stress for slip
    xi_twin_0, &                                                                                    !< initial critical shear stress for twin
    a                                                                                               !< non-Schmid coefficients
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    pl

  print'(/,a)', ' <<<+-  plastic_phenopowerlaw init  -+>>>'

  myPlasticity = plastic_active('phenopowerlaw')
  Ninstance = count(myPlasticity)
  print'(a,i2)', ' # instances: ',Ninstance; flush(IO_STDOUT)
  if(Ninstance == 0) return

  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))

  phases => config_material%get('phase')
  i = 0
  do p = 1, phases%length
    phase => phases%get(p)

    if(.not. myPlasticity(p)) cycle
    i = i + 1
    associate(prm => param(i), &
              dot => dotState(i), &
              stt => state(i))
    pl  => phase%get('plasticity')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                           phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(phase%get_asString('lattice') == 'bcc') then
        a = pl%get_asFloats('nonSchmid_coefficients',defaultVal=emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos  = prm%P_sl
        prm%nonSchmid_neg  = prm%P_sl
      endif
      prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(N_sl, &
                                                                pl%get_asFloats('h_sl_sl'), &
                                                                phase%get_asString('lattice'))

      xi_slip_0       = pl%get_asFloats('xi_0_sl',   requiredSize=size(N_sl))
      prm%xi_slip_sat = pl%get_asFloats('xi_inf_sl', requiredSize=size(N_sl))
      prm%H_int       = pl%get_asFloats('h_int',     requiredSize=size(N_sl), &
                                         defaultVal=[(0.0_pReal,i=1,size(N_sl))])

      prm%gdot0_slip  = pl%get_asFloat('dot_gamma_0_sl')
      prm%n_slip      = pl%get_asFloat('n_sl')
      prm%a_slip      = pl%get_asFloat('a_sl')
      prm%h0_SlipSlip = pl%get_asFloat('h_0_sl_sl')

      ! expand: family => system
      xi_slip_0       = math_expand(xi_slip_0,      N_sl)
      prm%xi_slip_sat = math_expand(prm%xi_slip_sat,N_sl)
      prm%H_int       = math_expand(prm%H_int,      N_sl)

      ! sanity checks
      if (    prm%gdot0_slip  <= 0.0_pReal)      extmsg = trim(extmsg)//' dot_gamma_0_sl'
      if (    prm%a_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' a_sl'
      if (    prm%n_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' n_sl'
      if (any(xi_slip_0       <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_0_sl'
      if (any(prm%xi_slip_sat <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_inf_sl'

    else slipActive
      xi_slip_0 = emptyRealArray
      allocate(prm%xi_slip_sat,prm%H_int,source=emptyRealArray)
      allocate(prm%interaction_SlipSlip(0,0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    N_tw         = pl%get_asInts('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%P_tw                 = lattice_SchmidMatrix_twin(N_tw,phase%get_asString('lattice'),&
                                                           phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%interaction_TwinTwin = lattice_interaction_TwinByTwin(N_tw,&
                                                                pl%get_asFloats('h_tw_tw'), &
                                                                phase%get_asString('lattice'))
      prm%gamma_twin_char      = lattice_characteristicShear_twin(N_tw,phase%get_asString('lattice'),&
                                                                  phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      xi_twin_0       = pl%get_asFloats('xi_0_tw',requiredSize=size(N_tw))

      prm%c_1         = pl%get_asFloat('c_1',defaultVal=0.0_pReal)
      prm%c_2         = pl%get_asFloat('c_2',defaultVal=1.0_pReal)
      prm%c_3         = pl%get_asFloat('c_3',defaultVal=0.0_pReal)
      prm%c_4         = pl%get_asFloat('c_4',defaultVal=0.0_pReal)
      prm%gdot0_twin  = pl%get_asFloat('dot_gamma_0_tw')
      prm%n_twin      = pl%get_asFloat('n_tw')
      prm%spr         = pl%get_asFloat('f_sl_sat_tw')
      prm%h0_TwinTwin = pl%get_asFloat('h_0_tw_tw')

      ! expand: family => system
      xi_twin_0       = math_expand(xi_twin_0,N_tw)

      ! sanity checks
      if (prm%gdot0_twin <= 0.0_pReal)  extmsg = trim(extmsg)//' dot_gamma_0_tw'
      if (prm%n_twin     <= 0.0_pReal)  extmsg = trim(extmsg)//' n_tw'

    else twinActive
      xi_twin_0 = emptyRealArray
      allocate(prm%gamma_twin_char,source=emptyRealArray)
      allocate(prm%interaction_TwinTwin(0,0))
    endif twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h0_TwinSlip          = pl%get_asFloat('h_0_tw_sl')
      prm%interaction_SlipTwin = lattice_interaction_SlipByTwin(N_sl,N_tw,&
                                                                pl%get_asFloats('h_sl_tw'), &
                                                                phase%get_asString('lattice'))
      prm%interaction_TwinSlip = lattice_interaction_TwinBySlip(N_tw,N_sl,&
                                                                pl%get_asFloats('h_tw_sl'), &
                                                                phase%get_asString('lattice'))
    else slipAndTwinActive
      allocate(prm%interaction_SlipTwin(prm%sum_N_sl,prm%sum_N_tw))                                 ! at least one dimension is 0
      allocate(prm%interaction_TwinSlip(prm%sum_N_tw,prm%sum_N_sl))                                 ! at least one dimension is 0
      prm%h0_TwinSlip = 0.0_pReal
    endif slipAndTwinActive

!--------------------------------------------------------------------------------------------------
!  output pararameters

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(pl)
#else
    prm%output = pl%get_asStrings('output',defaultVal=emptyStringArray)
#endif

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
    sizeDotState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw
    sizeState = sizeDotState


    call constitutive_allocateState(plasticState(p),NipcMyPhase,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%xi_slip => plasticState(p)%state   (startIndex:endIndex,:)
    stt%xi_slip =  spread(xi_slip_0, 2, NipcMyPhase)
    dot%xi_slip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%xi_twin => plasticState(p)%state   (startIndex:endIndex,:)
    stt%xi_twin =  spread(xi_twin_0, 2, NipcMyPhase)
    dot%xi_twin => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_slip => plasticState(p)%state   (startIndex:endIndex,:)
    dot%gamma_slip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'
    ! global alias
    plasticState(p)%slipRate => plasticState(p)%dotState(startIndex:endIndex,:)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%gamma_twin => plasticState(p)%state   (startIndex:endIndex,:)
    dot%gamma_twin => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(phenopowerlaw)')

  enddo

end function plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!> @details asummes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
pure module subroutine plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)

  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,               intent(in) :: &
    instance, &
    of

  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_slip_pos,gdot_slip_neg, &
    dgdot_dtauslip_pos,dgdot_dtauslip_neg
  real(pReal), dimension(param(instance)%sum_N_tw) :: &
    gdot_twin,dgdot_dtautwin

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(instance))

  call kinetics_slip(Mp,instance,of,gdot_slip_pos,gdot_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg)
  slipSystems: do i = 1, prm%sum_N_sl
    Lp = Lp + (gdot_slip_pos(i)+gdot_slip_neg(i))*prm%P_sl(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtauslip_pos(i) * prm%P_sl(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + dgdot_dtauslip_neg(i) * prm%P_sl(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo slipSystems

  call kinetics_twin(Mp,instance,of,gdot_twin,dgdot_dtautwin)
  twinSystems: do i = 1, prm%sum_N_tw
    Lp = Lp + gdot_twin(i)*prm%P_tw(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtautwin(i)*prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
  enddo twinSystems

  end associate

end subroutine plastic_phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_phenopowerlaw_dotState(Mp,instance,of)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal) :: &
    c_SlipSlip,c_TwinSlip,c_TwinTwin, &
    xi_slip_sat_offset,&
    sumGamma,sumF
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    left_SlipSlip,right_SlipSlip, &
    gdot_slip_pos,gdot_slip_neg

  associate(prm => param(instance), stt => state(instance), dot => dotState(instance))

  sumGamma = sum(stt%gamma_slip(:,of))
  sumF     = sum(stt%gamma_twin(:,of)/prm%gamma_twin_char)

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
  c_SlipSlip = prm%h0_slipslip * (1.0_pReal + prm%c_1*sumF** prm%c_2)
  c_TwinSlip = prm%h0_TwinSlip * sumGamma**prm%c_3
  c_TwinTwin = prm%h0_TwinTwin * sumF**prm%c_4

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors
  left_SlipSlip  = 1.0_pReal + prm%H_int
  xi_slip_sat_offset = prm%spr*sqrt(sumF)
  right_SlipSlip = abs(1.0_pReal-stt%xi_slip(:,of) / (prm%xi_slip_sat+xi_slip_sat_offset)) **prm%a_slip &
                 * sign(1.0_pReal,1.0_pReal-stt%xi_slip(:,of) / (prm%xi_slip_sat+xi_slip_sat_offset))

!--------------------------------------------------------------------------------------------------
! shear rates
  call kinetics_slip(Mp,instance,of,gdot_slip_pos,gdot_slip_neg)
  dot%gamma_slip(:,of) = abs(gdot_slip_pos+gdot_slip_neg)
  call kinetics_twin(Mp,instance,of,dot%gamma_twin(:,of))

!--------------------------------------------------------------------------------------------------
! hardening
  dot%xi_slip(:,of) = c_SlipSlip * left_SlipSlip * &
                      matmul(prm%interaction_SlipSlip,dot%gamma_slip(:,of)*right_SlipSlip) &
                    + matmul(prm%interaction_SlipTwin,dot%gamma_twin(:,of))

  dot%xi_twin(:,of) = c_TwinSlip * matmul(prm%interaction_TwinSlip,dot%gamma_slip(:,of)) &
                    + c_TwinTwin * matmul(prm%interaction_TwinTwin,dot%gamma_twin(:,of))
  end associate

end subroutine plastic_phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_phenopowerlaw_results(instance,group)

  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))

      case('xi_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%xi_slip,   trim(prm%output(o)), &
                                                     'resistance against plastic slip','Pa')
      case('gamma_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma_slip,trim(prm%output(o)), &
                                                     'plastic shear','1')

      case('xi_tw')
        if(prm%sum_N_tw>0) call results_writeDataset(group,stt%xi_twin,   trim(prm%output(o)), &
                                                     'resistance against twinning','Pa')
      case('gamma_tw')
        if(prm%sum_N_tw>0) call results_writeDataset(group,stt%gamma_twin,trim(prm%output(o)), &
                                                     'twinning shear','1')

    end select
  enddo outputsLoop
  end associate

end subroutine plastic_phenopowerlaw_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_slip(Mp,instance,of, &
                              gdot_slip_pos,gdot_slip_neg,dgdot_dtau_slip_pos,dgdot_dtau_slip_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal),                  intent(out), dimension(param(instance)%sum_N_sl) :: &
    gdot_slip_pos, &
    gdot_slip_neg
  real(pReal),                  intent(out), optional, dimension(param(instance)%sum_N_sl) :: &
    dgdot_dtau_slip_pos, &
    dgdot_dtau_slip_neg

  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    tau_slip_pos, &
    tau_slip_neg
  integer :: i

  associate(prm => param(instance), stt => state(instance))

  do i = 1, prm%sum_N_sl
    tau_slip_pos(i) =       math_tensordot(Mp,prm%nonSchmid_pos(1:3,1:3,i))
    tau_slip_neg(i) = merge(math_tensordot(Mp,prm%nonSchmid_neg(1:3,1:3,i)), &
                            0.0_pReal, prm%nonSchmidActive)
  enddo

  where(dNeq0(tau_slip_pos))
    gdot_slip_pos = prm%gdot0_slip * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &              ! 1/2 if non-Schmid active
                  * sign(abs(tau_slip_pos/stt%xi_slip(:,of))**prm%n_slip,  tau_slip_pos)
  else where
    gdot_slip_pos = 0.0_pReal
  end where

  where(dNeq0(tau_slip_neg))
    gdot_slip_neg = prm%gdot0_slip * 0.5_pReal &                                                    ! only used if non-Schmid active, always 1/2
                  * sign(abs(tau_slip_neg/stt%xi_slip(:,of))**prm%n_slip,  tau_slip_neg)
  else where
    gdot_slip_neg = 0.0_pReal
  end where

  if (present(dgdot_dtau_slip_pos)) then
    where(dNeq0(gdot_slip_pos))
      dgdot_dtau_slip_pos = gdot_slip_pos*prm%n_slip/tau_slip_pos
    else where
      dgdot_dtau_slip_pos = 0.0_pReal
    end where
  endif
  if (present(dgdot_dtau_slip_neg)) then
    where(dNeq0(gdot_slip_neg))
      dgdot_dtau_slip_neg = gdot_slip_neg*prm%n_slip/tau_slip_neg
    else where
      dgdot_dtau_slip_neg = 0.0_pReal
    end where
  endif
  end associate

end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress. Twinning is assumed to take place only in untwinned volume.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(Mp,instance,of,&
                              gdot_twin,dgdot_dtau_twin)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal), dimension(param(instance)%sum_N_tw), intent(out) :: &
    gdot_twin
  real(pReal), dimension(param(instance)%sum_N_tw), intent(out), optional :: &
    dgdot_dtau_twin

  real(pReal), dimension(param(instance)%sum_N_tw) :: &
    tau_twin
  integer :: i

  associate(prm => param(instance), stt => state(instance))

  do i = 1, prm%sum_N_tw
    tau_twin(i)  = math_tensordot(Mp,prm%P_tw(1:3,1:3,i))
  enddo

  where(tau_twin > 0.0_pReal)
    gdot_twin = (1.0_pReal-sum(stt%gamma_twin(:,of)/prm%gamma_twin_char)) &                         ! only twin in untwinned volume fraction
              * prm%gdot0_twin*(abs(tau_twin)/stt%xi_twin(:,of))**prm%n_twin
  else where
    gdot_twin = 0.0_pReal
  end where

  if (present(dgdot_dtau_twin)) then
    where(dNeq0(gdot_twin))
      dgdot_dtau_twin = gdot_twin*prm%n_twin/tau_twin
    else where
      dgdot_dtau_twin = 0.0_pReal
    end where
  endif

  end associate

end subroutine kinetics_twin

end submodule plastic_phenopowerlaw
