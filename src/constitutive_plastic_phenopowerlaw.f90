!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
submodule(constitutive) plastic_phenopowerlaw

  type :: tParameters
    real(pReal) :: &
      gdot0_slip, &                                                                                 !< reference shear strain rate for slip
      gdot0_twin, &                                                                                 !< reference shear strain rate for twin
      n_slip, &                                                                                     !< stress exponent for slip
      n_twin, &                                                                                     !< stress exponent for twin
      spr, &                                                                                        !< push-up factor for slip saturation due to twinning
      c_1, &
      c_2, &
      c_3, &
      c_4, &
      h0_SlipSlip, &                                                                                !< reference hardening slip - slip
      h0_TwinSlip, &                                                                                !< reference hardening twin - slip
      h0_TwinTwin, &                                                                                !< reference hardening twin - twin
      a_slip
    real(pReal),               allocatable, dimension(:) :: &
      xi_slip_0, &                                                                                  !< initial critical shear stress for slip
      xi_twin_0, &                                                                                  !< initial critical shear stress for twin
      xi_slip_sat, &                                                                                !< maximum critical shear stress for slip
      nonSchmidCoeff, &
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
module subroutine plastic_phenopowerlaw_init

  integer :: &
    Ninstance, &
    p, i, &
    NipcMyPhase, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer, dimension(:), allocatable :: &
    N_sl, N_tw
  character(len=pStringLen) :: &
    extmsg = ''

  write(6,'(/,a)') ' <<<+-  plastic_'//PLASTICITY_PHENOPOWERLAW_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))

  do p = 1, size(phase_plasticity)
    if (phase_plasticity(p) /= PLASTICITY_PHENOPOWERLAW_ID) cycle
    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)), &
              config => config_phase(p))

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
    prm%c_1            = config%getFloat('twin_c',defaultVal=0.0_pReal)
    prm%c_2            = config%getFloat('twin_b',defaultVal=1.0_pReal)
    prm%c_3            = config%getFloat('twin_e',defaultVal=0.0_pReal)
    prm%c_4            = config%getFloat('twin_d',defaultVal=0.0_pReal)

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = config%getInts('nslip',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(N_sl)
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,config%getString('lattice_structure'),&
                                           config%getFloat('c/a',defaultVal=0.0_pReal))

      if(trim(config%getString('lattice_structure')) == 'bcc') then
        prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                              defaultVal = emptyRealArray)
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(N_sl,prm%nonSchmidCoeff,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(N_sl,prm%nonSchmidCoeff,-1)
      else
        allocate(prm%nonSchmidCoeff(0))
        prm%nonSchmid_pos  = prm%P_sl
        prm%nonSchmid_neg  = prm%P_sl
      endif
      prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(N_sl, &
                                                                config%getFloats('interaction_slipslip'), &
                                                                config%getString('lattice_structure'))

      prm%xi_slip_0   = config%getFloats('tau0_slip',   requiredSize=size(N_sl))
      prm%xi_slip_sat = config%getFloats('tausat_slip', requiredSize=size(N_sl))
      prm%H_int       = config%getFloats('h_int',       requiredSize=size(N_sl), &
                                         defaultVal=[(0.0_pReal,i=1,size(N_sl))])

      prm%gdot0_slip  = config%getFloat('gdot0_slip')
      prm%n_slip      = config%getFloat('n_slip')
      prm%a_slip      = config%getFloat('a_slip')
      prm%h0_SlipSlip = config%getFloat('h0_slipslip')

      ! expand: family => system
      prm%xi_slip_0   = math_expand(prm%xi_slip_0,  N_sl)
      prm%xi_slip_sat = math_expand(prm%xi_slip_sat,N_sl)
      prm%H_int       = math_expand(prm%H_int,      N_sl)

      ! sanity checks
      if (    prm%gdot0_slip  <= 0.0_pReal)      extmsg = trim(extmsg)//' gdot0_slip'
      if (    prm%a_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' a_slip'
      if (    prm%n_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' n_slip'
      if (any(prm%xi_slip_0   <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_slip_0'
      if (any(prm%xi_slip_sat <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_slip_sat'
    else slipActive
      allocate(prm%interaction_SlipSlip(0,0))
      allocate(prm%xi_slip_0(0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    N_tw         = config%getInts('ntwin', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(N_tw)
    twinActive: if (prm%sum_N_tw > 0) then
      prm%P_tw                 = lattice_SchmidMatrix_twin(N_tw,config%getString('lattice_structure'),&
                                                           config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%interaction_TwinTwin = lattice_interaction_TwinByTwin(N_tw,&
                                                                config%getFloats('interaction_twintwin'), &
                                                                config%getString('lattice_structure'))
      prm%gamma_twin_char      = lattice_characteristicShear_twin(N_tw,config%getString('lattice_structure'),&
                                                                  config%getFloat('c/a'))

      prm%xi_twin_0            = config%getFloats('tau0_twin',requiredSize=size(N_tw))

      prm%gdot0_twin           = config%getFloat('gdot0_twin')
      prm%n_twin               = config%getFloat('n_twin')
      prm%spr                  = config%getFloat('s_pr')
      prm%h0_TwinTwin          = config%getFloat('h0_twintwin')

      ! expand: family => system
      prm%xi_twin_0   = math_expand(prm%xi_twin_0,  N_tw)

      ! sanity checks
      if (prm%gdot0_twin <= 0.0_pReal)  extmsg = trim(extmsg)//' gdot0_twin'
      if (prm%n_twin     <= 0.0_pReal)  extmsg = trim(extmsg)//' n_twin'
    else twinActive
      allocate(prm%interaction_TwinTwin(0,0))
      allocate(prm%xi_twin_0(0))
      allocate(prm%gamma_twin_char(0))
    endif twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h0_TwinSlip          = config%getFloat('h0_twinslip')
      prm%interaction_SlipTwin = lattice_interaction_SlipByTwin(N_sl,N_tw,&
                                                                config%getFloats('interaction_sliptwin'), &
                                                                config%getString('lattice_structure'))
      prm%interaction_TwinSlip = lattice_interaction_TwinBySlip(N_tw,N_sl,&
                                                                config%getFloats('interaction_twinslip'), &
                                                                config%getString('lattice_structure'))
    else slipAndTwinActive
      allocate(prm%interaction_SlipTwin(prm%sum_N_sl,prm%sum_N_tw))                                 ! at least one dimension is 0
      allocate(prm%interaction_TwinSlip(prm%sum_N_tw,prm%sum_N_sl))                                 ! at least one dimension is 0
      prm%h0_TwinSlip = 0.0_pReal
    endif slipAndTwinActive

!--------------------------------------------------------------------------------------------------
!  output pararameters
    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
    sizeDotState = size(['tau_slip  ','gamma_slip']) * prm%sum_N_sl &
                 + size(['tau_twin  ','gamma_twin']) * prm%sum_N_tw
    sizeState = sizeDotState

    call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%xi_slip => plasticState(p)%state   (startIndex:endIndex,:)
    stt%xi_slip = spread(prm%xi_slip_0, 2, NipcMyPhase)
    dot%xi_slip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_resistance',defaultVal=1.0_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex)<=0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%xi_twin => plasticState(p)%state   (startIndex:endIndex,:)
    stt%xi_twin = spread(prm%xi_twin_0, 2, NipcMyPhase)
    dot%xi_twin => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_resistance',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_slip => plasticState(p)%state   (startIndex:endIndex,:)
    dot%gamma_slip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex)<=0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma_slip'
    ! global alias
    plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    stt%gamma_twin => plasticState(p)%state   (startIndex:endIndex,:)
    dot%gamma_twin => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex)<=0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_PHENOPOWERLAW_LABEL//')')

  enddo

end subroutine plastic_phenopowerlaw_init


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

      case('resistance_slip')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%xi_slip,   'xi_sl', &
                                                     'resistance against plastic slip','Pa')
      case('accumulatedshear_slip')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma_slip,'gamma_sl', &
                                                     'plastic shear','1')

      case('resistance_twin')
        if(prm%sum_N_tw>0) call results_writeDataset(group,stt%xi_twin,   'xi_tw', &
                                                     'resistance against twinning','Pa')
      case('accumulatedshear_twin')
        if(prm%sum_N_tw>0) call results_writeDataset(group,stt%gamma_twin,'gamma_tw', &
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
  logical :: nonSchmidActive

  associate(prm => param(instance), stt => state(instance))

  nonSchmidActive = size(prm%nonSchmidCoeff) > 0

  do i = 1, prm%sum_N_sl
    tau_slip_pos(i) =       math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))
    tau_slip_neg(i) = merge(math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,i)), &
                            0.0_pReal, nonSchmidActive)
  enddo

  where(dNeq0(tau_slip_pos))
    gdot_slip_pos = prm%gdot0_slip * merge(0.5_pReal,1.0_pReal, nonSchmidActive) &                  ! 1/2 if non-Schmid active
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
    tau_twin(i)  = math_mul33xx33(Mp,prm%P_tw(1:3,1:3,i))
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
