!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystal plasticity model for bcc metals, especially Tungsten
!--------------------------------------------------------------------------------------------------
submodule(constitutive) plastic_disloUCLA

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
      atomicVolume, &
      tau_0, &
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

  type :: tDisloUCLAState
    real(pReal), dimension(:,:), pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tDisloUCLAState

  type :: tDisloUCLAdependentState
    real(pReal), dimension(:,:), allocatable :: &
      Lambda_sl, &
      threshold_stress
  end type tDisloUCLAdependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),              allocatable, dimension(:) :: param
  type(tDisloUCLAState),          allocatable, dimension(:) :: &
    dotState, &
    state
  type(tDisloUCLAdependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloUCLA_init

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

  write(6,'(/,a)') ' <<<+-  plastic_'//PLASTICITY_DISLOUCLA_LABEL//' init  -+>>>'; flush(6)

  write(6,'(/,a)') ' Cereceda et al., International Journal of Plasticity 78:242–256, 2016'
  write(6,'(a)')   ' https://dx.doi.org/10.1016/j.ijplas.2015.09.002'

  Ninstance = count(phase_plasticity == PLASTICITY_DISLOUCLA_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))
  allocate(dependentState(Ninstance))

  do p = 1, size(phase_plasticity)
    if (phase_plasticity(p) /= PLASTICITY_DISLOUCLA_ID) cycle
    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)), &
              dst => dependentState(phase_plasticityInstance(p)), &
              config => config_phase(p))

    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

    ! This data is read in already in lattice
    prm%mu = lattice_mu(p)

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = config%getInts('nslip',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = lattice_SchmidMatrix_slip(N_sl,config%getString('lattice_structure'),&
                                           config%getFloat('c/a',defaultVal=0.0_pReal))

      if(trim(config%getString('lattice_structure')) == 'bcc') then
        a = config%getFloats('nonschmid_coefficients',defaultVal = emptyRealArray)
        prm%nonSchmid_pos = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos = prm%P_sl
        prm%nonSchmid_neg = prm%P_sl
      endif

      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,config%getFloats('interaction_slipslip'), &
                                                   config%getString('lattice_structure'))
      prm%forestProjection = lattice_forestProjection_edge(N_sl,config%getString('lattice_structure'),&
                                                           config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%forestProjection = transpose(prm%forestProjection)

      rho_mob_0       = config%getFloats('rhoedge0',       requiredSize=size(N_sl))
      rho_dip_0       = config%getFloats('rhoedgedip0',    requiredSize=size(N_sl))
      prm%v0          = config%getFloats('v0',             requiredSize=size(N_sl))
      prm%b_sl        = config%getFloats('slipburgers',    requiredSize=size(N_sl))
      prm%delta_F     = config%getFloats('qedge',          requiredSize=size(N_sl))

      prm%i_sl        = config%getFloats('clambdaslip',    requiredSize=size(N_sl))
      prm%tau_0       = config%getFloats('tau_peierls',    requiredSize=size(N_sl))
      prm%p           = config%getFloats('p_slip',         requiredSize=size(N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(N_sl))])
      prm%q           = config%getFloats('q_slip',         requiredSize=size(N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(N_sl))])
      prm%kink_height = config%getFloats('kink_height',    requiredSize=size(N_sl))
      prm%w           = config%getFloats('kink_width',     requiredSize=size(N_sl))
      prm%omega       = config%getFloats('omega',          requiredSize=size(N_sl))
      prm%B           = config%getFloats('friction_coeff', requiredSize=size(N_sl))

      prm%D               = config%getFloat('grainsize')
      prm%D_0             = config%getFloat('d0')
      prm%Q_cl            = config%getFloat('qsd')
      prm%atomicVolume    = config%getFloat('catomicvolume')       * prm%b_sl**3.0_pReal
      prm%D_a             = config%getFloat('cedgedipmindistance') * prm%b_sl
      prm%dipoleformation = config%getFloat('dipoleformationfactor') > 0.0_pReal                    !should be on by default, ToDo: change to /key/-type key

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
      if (any(rho_mob_0        <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedge0'
      if (any(rho_dip_0        <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedgedip0'
      if (any(prm%v0           <  0.0_pReal)) extmsg = trim(extmsg)//' v0'
      if (any(prm%b_sl         <= 0.0_pReal)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%delta_F      <= 0.0_pReal)) extmsg = trim(extmsg)//' qedge'
      if (any(prm%tau_0        <  0.0_pReal)) extmsg = trim(extmsg)//' tau_0'
      if (any(prm%D_a          <= 0.0_pReal)) extmsg = trim(extmsg)//' cedgedipmindistance or b_sl'
      if (any(prm%atomicVolume <= 0.0_pReal)) extmsg = trim(extmsg)//' catomicvolume or b_sl'

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

    call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%rho_mob => plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_mob =  spread(rho_mob_0,2,NipcMyPhase)
    dot%rho_mob => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_rho',defaultVal=1.0_pReal)
    if (any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%rho_dip => plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_dip =  spread(rho_dip_0,2,NipcMyPhase)
    dot%rho_dip => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_rho',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_sl => plasticState(p)%state(startIndex:endIndex,:)
    dot%gamma_sl => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = 1.0e6_pReal                                         ! ARRG
    ! global alias
    plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)

    allocate(dst%Lambda_sl(prm%sum_N_sl,NipcMyPhase),         source=0.0_pReal)
    allocate(dst%threshold_stress(prm%sum_N_sl,NipcMyPhase),  source=0.0_pReal)

    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_DISLOUCLA_LABEL//')')

  enddo

end subroutine plastic_disloUCLA_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp, &
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

end subroutine plastic_disloUCLA_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloUCLA_dotState(Mp,T,instance,of)

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

end subroutine plastic_disloUCLA_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloUCLA_dependentState(instance,of)

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

end subroutine plastic_disloUCLA_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_disloUCLA_results(instance,group)

  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('edge_density')                                                                          ! ToDo: should be rho_mob
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_mob,'rho_mob',&
                                                     'mobile dislocation density','1/m²')
      case('dipole_density')                                                                        ! ToDo: should be rho_dip
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_dip,'rho_dip',&
                                                     'dislocation dipole density''1/m²')
      case('shear_rate_slip')                                                                       ! should be gamma
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma_sl,'dot_gamma_sl',&            ! this is not dot!!
                                                     'plastic shear','1')
      case('mfp_slip')                                                                              !ToDo: should be Lambda
        if(prm%sum_N_sl>0) call results_writeDataset(group,dst%Lambda_sl,'Lambda_sl',&
                                                     'mean free path for slip','m')
      case('threshold_stress_slip')                                                                 !ToDo: should be tau_pass
        if(prm%sum_N_sl>0) call results_writeDataset(group,dst%threshold_stress,'tau_pass',&
                                                     'threshold stress for slip','Pa')
    end select
  enddo outputsLoop
  end associate

end subroutine plastic_disloUCLA_results


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
    tau_pos(j) = math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,j))
    tau_neg(j) = math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,j))
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

end submodule plastic_disloUCLA
