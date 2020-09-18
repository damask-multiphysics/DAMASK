!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystal plasticity model for bcc metals, especially Tungsten
!--------------------------------------------------------------------------------------------------
submodule(constitutive:constitutive_plastic) plastic_disloTungsten

  real(pReal), parameter :: &
    kB = 1.38e-23_pReal                                                                             !< Boltzmann constant in J/Kelvin

  type :: tParameters
    real(pReal) :: &
      D    = 1.0_pReal, &                                                                           !< grain size
      mu   = 1.0_pReal, &                                                                           !< equivalent shear modulus
      D_0  = 1.0_pReal, &                                                                           !< prefactor for self-diffusion coefficient
      Q_cl = 1.0_pReal                                                                              !< activation energy for dislocation climb
    real(pReal),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< magnitude of burgers vector [m]
      D_a, &
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations
      atomicVolume, &                                                                               !< factor to calculate atomic volume
      tau_0, &                                                                                      !< Peierls stress
      !* mobility law parameters
      delta_F, &                                                                                    !< activation energy for glide [J]
      v0, &                                                                                         !< dislocation velocity prefactor [m/s]
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      B, &                                                                                          !< friction coefficient
      kink_height, &                                                                                !< height of the kink pair
      w, &                                                                                          !< width of the kink pair
      omega                                                                                         !< attempt frequency for kink pair nucleation
    real(pReal),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      forestProjection
    real(pReal),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl                                                                                      !< total number of active slip system
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
    logical :: &
      dipoleFormation                                                                               !< flag indicating consideration of dipole formation
  end type                                                                                          !< container type for internal constitutive parameters

  type :: tDisloTungstenState
    real(pReal), dimension(:,:), pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tDisloTungstenState

  type :: tDisloTungstendependentState
    real(pReal), dimension(:,:), allocatable :: &
      Lambda_sl, &
      threshold_stress
  end type tDisloTungstendependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),              allocatable, dimension(:) :: param
  type(tDisloTungstenState),          allocatable, dimension(:) :: &
    dotState, &
    state
  type(tDisloTungstendependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_disloTungsten_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    Ninstance, &
    p, i, &
    NipcMyPhase, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,    dimension(:), allocatable :: &
    N_sl
  real(pReal),dimension(:), allocatable :: &
    rho_mob_0, &                                                                                    !< initial dislocation density
    rho_dip_0, &                                                                                    !< initial dipole density
    a                                                                                               !< non-Schmid coefficients
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    pl

  print'(/,a)', ' <<<+-  plastic_dislotungsten init  -+>>>'

  myPlasticity = plastic_active('disloTungsten')
  Ninstance = count(myPlasticity)
  print'(a,i2)', ' # instances: ',Ninstance; flush(6)
  if(Ninstance == 0) return
  
  print*, 'Cereceda et al., International Journal of Plasticity 78:242–256, 2016'
  print*, 'https://dx.doi.org/10.1016/j.ijplas.2015.09.002'

  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))
  allocate(dependentState(Ninstance))

  phases => config_material%get('phase')
  i = 0
  do p = 1, phases%length
    phase => phases%get(p)

    if(.not. myPlasticity(p)) cycle
    i = i + 1
    associate(prm => param(i), &
              dot => dotState(i), &
              stt => state(i), &
              dst => dependentState(i))
    pl  => phase%get('plasticity')

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(pl)
#else
    prm%output = pl%get_asStrings('output',defaultVal=emptyStringArray)
#endif

    ! This data is read in already in lattice
    prm%mu = lattice_mu(p)

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                           phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(trim(phase%get_asString('lattice')) == 'bcc') then
        a = pl%get_asFloats('nonSchmid_coefficients',defaultVal = emptyRealArray)
        prm%nonSchmid_pos = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos = prm%P_sl
        prm%nonSchmid_neg = prm%P_sl
      endif

      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,pl%get_asFloats('h_sl_sl'), &
                                                   phase%get_asString('lattice'))
      prm%forestProjection = lattice_forestProjection_edge(N_sl,phase%get_asString('lattice'),&
                                                           phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%forestProjection = transpose(prm%forestProjection)

      rho_mob_0       = pl%get_asFloats('rho_mob_0',     requiredSize=size(N_sl))
      rho_dip_0       = pl%get_asFloats('rho_dip_0',     requiredSize=size(N_sl))
      prm%v0          = pl%get_asFloats('v_0',           requiredSize=size(N_sl))
      prm%b_sl        = pl%get_asFloats('b_sl',          requiredSize=size(N_sl))
      prm%delta_F     = pl%get_asFloats('Q_s',           requiredSize=size(N_sl))

      prm%i_sl        = pl%get_asFloats('i_sl',          requiredSize=size(N_sl))
      prm%tau_0       = pl%get_asFloats('tau_peierls',   requiredSize=size(N_sl))
      prm%p           = pl%get_asFloats('p_sl',          requiredSize=size(N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(N_sl))])
      prm%q           = pl%get_asFloats('q_sl',          requiredSize=size(N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(N_sl))])
      prm%kink_height = pl%get_asFloats('h',             requiredSize=size(N_sl))
      prm%w           = pl%get_asFloats('w',             requiredSize=size(N_sl))
      prm%omega       = pl%get_asFloats('omega',         requiredSize=size(N_sl))
      prm%B           = pl%get_asFloats('B',             requiredSize=size(N_sl))

      prm%D               = pl%get_asFloat('D')
      prm%D_0             = pl%get_asFloat('D_0')
      prm%Q_cl            = pl%get_asFloat('Q_cl')
      prm%atomicVolume    = pl%get_asFloat('f_at')       * prm%b_sl**3.0_pReal
      prm%D_a             = pl%get_asFloat('D_a')        * prm%b_sl

      prm%dipoleformation = pl%get_asBool('dipole_formation_factor', defaultVal = .true.)

      ! expand: family => system
      rho_mob_0          = math_expand(rho_mob_0,          N_sl)
      rho_dip_0          = math_expand(rho_dip_0,          N_sl)
      prm%q              = math_expand(prm%q,              N_sl)
      prm%p              = math_expand(prm%p,              N_sl)
      prm%delta_F        = math_expand(prm%delta_F,        N_sl)
      prm%b_sl           = math_expand(prm%b_sl,           N_sl)
      prm%kink_height    = math_expand(prm%kink_height,    N_sl)
      prm%w              = math_expand(prm%w,              N_sl)
      prm%omega          = math_expand(prm%omega,          N_sl)
      prm%tau_0          = math_expand(prm%tau_0,          N_sl)
      prm%v0             = math_expand(prm%v0,             N_sl)
      prm%B              = math_expand(prm%B,              N_sl)
      prm%i_sl           = math_expand(prm%i_sl,           N_sl)
      prm%atomicVolume   = math_expand(prm%atomicVolume,   N_sl)
      prm%D_a            = math_expand(prm%D_a,            N_sl)

      ! sanity checks
      if (    prm%D_0          <= 0.0_pReal)  extmsg = trim(extmsg)//' D_0'
      if (    prm%Q_cl         <= 0.0_pReal)  extmsg = trim(extmsg)//' Q_cl'
      if (any(rho_mob_0        <  0.0_pReal)) extmsg = trim(extmsg)//' rho_mob_0'
      if (any(rho_dip_0        <  0.0_pReal)) extmsg = trim(extmsg)//' rho_dip_0'
      if (any(prm%v0           <  0.0_pReal)) extmsg = trim(extmsg)//' v_0'
      if (any(prm%b_sl         <= 0.0_pReal)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%delta_F      <= 0.0_pReal)) extmsg = trim(extmsg)//' Q_s'
      if (any(prm%tau_0        <  0.0_pReal)) extmsg = trim(extmsg)//' tau_peierls'
      if (any(prm%D_a          <= 0.0_pReal)) extmsg = trim(extmsg)//' D_a or b_sl'
      if (any(prm%atomicVolume <= 0.0_pReal)) extmsg = trim(extmsg)//' f_at or b_sl'

    else slipActive
      rho_mob_0= emptyRealArray; rho_dip_0 = emptyRealArray
      allocate(prm%b_sl,prm%D_a,prm%i_sl,prm%atomicVolume,prm%tau_0, &
               prm%delta_F,prm%v0,prm%p,prm%q,prm%B,prm%kink_height,prm%w,prm%omega, &
               source = emptyRealArray)
      allocate(prm%forestProjection(0,0))
      allocate(prm%h_sl_sl         (0,0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl
    sizeState = sizeDotState

    call constitutive_allocateState(plasticState(p),NipcMyPhase,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%rho_mob => plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_mob =  spread(rho_mob_0,2,NipcMyPhase)
    dot%rho_mob => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_rho',defaultVal=1.0_pReal)
    if (any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%rho_dip => plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_dip =  spread(rho_dip_0,2,NipcMyPhase)
    dot%rho_dip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_rho',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_sl => plasticState(p)%state(startIndex:endIndex,:)
    dot%gamma_sl => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = 1.0e-2_pReal
    ! global alias
    plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)

    allocate(dst%Lambda_sl(prm%sum_N_sl,NipcMyPhase),         source=0.0_pReal)
    allocate(dst%threshold_stress(prm%sum_N_sl,NipcMyPhase),  source=0.0_pReal)

    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(disloTungsten)')

  enddo

end function plastic_disloTungsten_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine plastic_disloTungsten_LpAndItsTangent(Lp,dLp_dMp, &
                                                         Mp,T,instance,of)
  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                 intent(in) :: &
    T                                                                                               !< temperature
  integer,                     intent(in) :: &
    instance, &
    of

  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    dot_gamma_pos,dot_gamma_neg, &
    ddot_gamma_dtau_pos,ddot_gamma_dtau_neg

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(instance))

  call kinetics(Mp,T,instance,of,dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg)
  do i = 1, prm%sum_N_sl
    Lp = Lp + (dot_gamma_pos(i)+dot_gamma_neg(i))*prm%P_sl(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_pos(i) * prm%P_sl(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + ddot_gamma_dtau_neg(i) * prm%P_sl(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo

  end associate

end subroutine plastic_disloTungsten_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloTungsten_dotState(Mp,T,instance,of)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal) :: &
    VacancyDiffusion
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos, gdot_neg,&
    tau_pos,&
    tau_neg, &
    v_cl, &
    dot_rho_dip_formation, &
    dot_rho_dip_climb, &
    dip_distance

  associate(prm => param(instance), stt => state(instance),dot => dotState(instance), dst => dependentState(instance))

  call kinetics(Mp,T,instance,of,&
                gdot_pos,gdot_neg, &
                tau_pos_out = tau_pos,tau_neg_out = tau_neg)

  dot%gamma_sl(:,of) = (gdot_pos+gdot_neg)                                                          ! ToDo: needs to be abs
  VacancyDiffusion = prm%D_0*exp(-prm%Q_cl/(kB*T))

  where(dEq0(tau_pos))                                                                              ! ToDo: use avg of pos and neg
    dot_rho_dip_formation = 0.0_pReal
    dot_rho_dip_climb     = 0.0_pReal
  else where
    dip_distance = math_clip(3.0_pReal*prm%mu*prm%b_sl/(16.0_pReal*PI*abs(tau_pos)), &
                             prm%D_a, &                                                             ! lower limit
                             dst%Lambda_sl(:,of))                                                   ! upper limit
    dot_rho_dip_formation = merge(2.0_pReal*dip_distance* stt%rho_mob(:,of)*abs(dot%gamma_sl(:,of))/prm%b_sl, & ! ToDo: ignore region of spontaneous annihilation
                                  0.0_pReal, &
                                  prm%dipoleformation)
    v_cl = (3.0_pReal*prm%mu*VacancyDiffusion*prm%atomicVolume/(2.0_pReal*pi*kB*T)) &
                  * (1.0_pReal/(dip_distance+prm%D_a))
    dot_rho_dip_climb = (4.0_pReal*v_cl*stt%rho_dip(:,of))/(dip_distance-prm%D_a)                   ! ToDo: Discuss with Franz: Stress dependency?
  end where

  dot%rho_mob(:,of) = abs(dot%gamma_sl(:,of))/(prm%b_sl*dst%Lambda_sl(:,of)) &                      ! multiplication
                    - dot_rho_dip_formation &
                    - (2.0_pReal*prm%D_a)/prm%b_sl*stt%rho_mob(:,of)*abs(dot%gamma_sl(:,of))        ! Spontaneous annihilation of 2 single edge dislocations
  dot%rho_dip(:,of) = dot_rho_dip_formation &
                    - (2.0_pReal*prm%D_a)/prm%b_sl*stt%rho_dip(:,of)*abs(dot%gamma_sl(:,of)) &      ! Spontaneous annihilation of a single edge dislocation with a dipole constituent
                    - dot_rho_dip_climb

  end associate

end subroutine plastic_disloTungsten_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloTungsten_dependentState(instance,of)

  integer,      intent(in) :: &
    instance, &
    of

  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    dislocationSpacing

  associate(prm => param(instance), stt => state(instance),dst => dependentState(instance))

  dislocationSpacing = sqrt(matmul(prm%forestProjection,stt%rho_mob(:,of)+stt%rho_dip(:,of)))
  dst%threshold_stress(:,of) = prm%mu*prm%b_sl &
                             * sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,of)+stt%rho_dip(:,of)))

  dst%Lambda_sl(:,of) = prm%D/(1.0_pReal+prm%D*dislocationSpacing/prm%i_sl)

  end associate

end subroutine plastic_disloTungsten_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloTungsten_results(instance,group)

  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('rho_mob')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_mob,trim(prm%output(o)), &
                                                     'mobile dislocation density','1/m²')
      case('rho_dip')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_dip,trim(prm%output(o)), &
                                                     'dislocation dipole density''1/m²')
      case('gamma_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma_sl,trim(prm%output(o)), &
                                                     'plastic shear','1')
      case('Lambda_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(group,dst%Lambda_sl,trim(prm%output(o)), &
                                                     'mean free path for slip','m')
      case('tau_pass')
        if(prm%sum_N_sl>0) call results_writeDataset(group,dst%threshold_stress,trim(prm%output(o)), &
                                                     'threshold stress for slip','Pa')
    end select
  enddo outputsLoop
  end associate

end subroutine plastic_disloTungsten_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,T,instance,of, &
                 dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg,tau_pos_out,tau_neg_out)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal),                  intent(out), dimension(param(instance)%sum_N_sl) :: &
    dot_gamma_pos, &
    dot_gamma_neg
  real(pReal),                  intent(out), optional, dimension(param(instance)%sum_N_sl) :: &
    ddot_gamma_dtau_pos, &
    ddot_gamma_dtau_neg, &
    tau_pos_out, &
    tau_neg_out
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    StressRatio, &
    StressRatio_p,StressRatio_pminus1, &
    dvel, vel, &
    tau_pos,tau_neg, &
    t_n, t_k, dtk,dtn, &
    needsGoodName                                                                                   ! ToDo: @Karo: any idea?
  integer :: j

  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

  do j = 1, prm%sum_N_sl
    tau_pos(j) = math_tensordot(Mp,prm%nonSchmid_pos(1:3,1:3,j))
    tau_neg(j) = math_tensordot(Mp,prm%nonSchmid_neg(1:3,1:3,j))
  enddo


  if (present(tau_pos_out)) tau_pos_out = tau_pos
  if (present(tau_neg_out)) tau_neg_out = tau_neg

  associate(BoltzmannRatio  => prm%delta_F/(kB*T), &
            dot_gamma_0     => stt%rho_mob(:,of)*prm%b_sl*prm%v0, &
            effectiveLength => dst%Lambda_sl(:,of) - prm%w)

  significantPositiveTau: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
    StressRatio = (abs(tau_pos)-dst%threshold_stress(:,of))/prm%tau_0
    StressRatio_p       = StressRatio** prm%p
    StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
    needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)

    t_n = prm%b_sl/(needsGoodName*prm%omega*effectiveLength)
    t_k = effectiveLength * prm%B /(2.0_pReal*prm%b_sl*tau_pos)

    vel = prm%kink_height/(t_n + t_k)

    dot_gamma_pos = dot_gamma_0 * sign(vel,tau_pos) * 0.5_pReal
  else where significantPositiveTau
    dot_gamma_pos = 0.0_pReal
  end where significantPositiveTau

  if (present(ddot_gamma_dtau_pos)) then
    significantPositiveTau2: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
      dtn = -1.0_pReal * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pReal-StressRatio_p)**(prm%q - 1.0_pReal) &
          * (StressRatio)**(prm%p - 1.0_pReal) / prm%tau_0
      dtk = -1.0_pReal * t_k / tau_pos

      dvel = -1.0_pReal * prm%kink_height * (dtk + dtn) / (t_n + t_k)**2.0_pReal

      ddot_gamma_dtau_pos = dot_gamma_0 * dvel* 0.5_pReal
    else where significantPositiveTau2
      ddot_gamma_dtau_pos = 0.0_pReal
    end where significantPositiveTau2
  endif

  significantNegativeTau: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
    StressRatio = (abs(tau_neg)-dst%threshold_stress(:,of))/prm%tau_0
    StressRatio_p       = StressRatio** prm%p
    StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
    needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)

    t_n = prm%b_sl/(needsGoodName*prm%omega*effectiveLength)
    t_k = effectiveLength * prm%B /(2.0_pReal*prm%b_sl*tau_pos)

    vel = prm%kink_height/(t_n + t_k)

    dot_gamma_neg = dot_gamma_0 * sign(vel,tau_neg) * 0.5_pReal
  else where significantNegativeTau
    dot_gamma_neg = 0.0_pReal
  end where significantNegativeTau

  if (present(ddot_gamma_dtau_neg)) then
    significantNegativeTau2: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
      dtn = -1.0_pReal * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pReal-StressRatio_p)**(prm%q - 1.0_pReal) &
          * (StressRatio)**(prm%p - 1.0_pReal) / prm%tau_0
      dtk = -1.0_pReal * t_k / tau_neg

      dvel = -1.0_pReal * prm%kink_height * (dtk + dtn) / (t_n + t_k)**2.0_pReal

      ddot_gamma_dtau_neg = dot_gamma_0 * dvel * 0.5_pReal
    else where significantNegativeTau2
      ddot_gamma_dtau_neg = 0.0_pReal
    end where significantNegativeTau2
  end if

  end associate
  end associate

end subroutine kinetics

end submodule plastic_disloTungsten
