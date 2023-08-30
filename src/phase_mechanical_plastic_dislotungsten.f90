!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystal plasticity model for bcc metals, especially Tungsten
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) dislotungsten

  type :: tParameters
    real(pREAL) :: &
      D    = 1.0_pREAL, &                                                                           !< grain size
      D_0  = 1.0_pREAL, &                                                                           !< prefactor for self-diffusion coefficient
      Q_cl = 1.0_pREAL                                                                              !< activation energy for dislocation climb
    real(pREAL),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< magnitude of Burgers vector [m]
      d_caron, &                                                                                    !< distance of spontaneous annhihilation
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations
      f_at, &                                                                                       !< factor to calculate atomic volume
      tau_Peierls, &                                                                                !< Peierls stress
      !* mobility law parameters
      Q_s, &                                                                                        !< activation energy for glide [J]
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      B, &                                                                                          !< friction coefficient
      h, &                                                                                          !< height of the kink pair
      w, &                                                                                          !< width of the kink pair
      omega                                                                                         !< attempt frequency for kink pair nucleation
    real(pREAL),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      forestProjection
    real(pREAL),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl                                                                                      !< total number of active slip system
    character(len=:),          allocatable               :: &
      isotropic_bound
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
    logical :: &
      dipoleFormation                                                                               !< flag indicating consideration of dipole formation
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl
  end type tParameters                                                                              !< container type for internal constitutive parameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tIndexDotState

  type :: tDislotungstenState
    real(pREAL), dimension(:,:), pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tDislotungstenState

  type :: tDislotungstenDependentState
    real(pREAL), dimension(:,:), allocatable :: &
      Lambda_sl, &
      tau_pass
  end type tDislotungstenDependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),                  allocatable, dimension(:) :: param
  type(tIndexDotState),               allocatable, dimension(:) :: indexDotState
  type(tDisloTungstenState),          allocatable, dimension(:) :: state
  type(tDisloTungstenDependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislotungsten_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,    dimension(:), allocatable :: &
    N_sl
  real(pREAL),dimension(:), allocatable :: &
    f_edge, &                                                                                       !< edge character fraction of total dislocation density
    rho_mob_0, &                                                                                    !< initial dislocation density
    rho_dip_0, &                                                                                    !< initial dipole density
    a                                                                                               !< non-Schmid coefficients
  character(len=:), allocatable :: &
    refs, &
    extmsg
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislotungsten')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:dislotungsten init  -+>>>'

  print'(/,1x,a)', 'D. Cereceda et al., International Journal of Plasticity 78:242–256, 2016'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2015.09.002'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  allocate(dependentState(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), &
              stt => state(ph), dst => dependentState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech  => phase%get_dict('mechanical')
    pl  => mech%get_dict('plastic')

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(pl,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(pl)
#else
    prm%output = pl%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

    prm%isotropic_bound = pl%get_asStr('isotropic_bound',defaultVal='isostrain')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))
      prm%P_sl = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        a = pl%get_as1dReal('a_nonSchmid',defaultVal = emptyRealArray)
        prm%P_nS_pos = crystal_nonSchmidMatrix(N_sl,a,+1)
        prm%P_nS_neg = crystal_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%P_nS_pos = prm%P_sl
        prm%P_nS_neg = prm%P_sl
      end if

      prm%dipoleformation = .not. pl%get_asBool('no_dipole_formation', defaultVal=.false.)

      prm%D    = pl%get_asReal('D')
      prm%D_0  = pl%get_asReal('D_0')
      prm%Q_cl = pl%get_asReal('Q_cl')

      f_edge          = math_expand(pl%get_as1dReal('f_edge',      requiredSize=size(N_sl), &
                                                    defaultVal=[(0.5_pREAL, i=1,size(N_sl))]),N_sl)
      rho_mob_0       = math_expand(pl%get_as1dReal('rho_mob_0',   requiredSize=size(N_sl)),N_sl)
      rho_dip_0       = math_expand(pl%get_as1dReal('rho_dip_0',   requiredSize=size(N_sl)),N_sl)
      prm%b_sl        = math_expand(pl%get_as1dReal('b_sl',        requiredSize=size(N_sl)),N_sl)
      prm%Q_s         = math_expand(pl%get_as1dReal('Q_s',         requiredSize=size(N_sl)),N_sl)
      prm%i_sl        = math_expand(pl%get_as1dReal('i_sl',        requiredSize=size(N_sl)),N_sl)
      prm%tau_Peierls = math_expand(pl%get_as1dReal('tau_Peierls', requiredSize=size(N_sl)),N_sl)
      prm%p           = math_expand(pl%get_as1dReal('p_sl',        requiredSize=size(N_sl)),N_sl)
      prm%q           = math_expand(pl%get_as1dReal('q_sl',        requiredSize=size(N_sl)),N_sl)
      prm%h           = math_expand(pl%get_as1dReal('h',           requiredSize=size(N_sl)),N_sl)
      prm%w           = math_expand(pl%get_as1dReal('w',           requiredSize=size(N_sl)),N_sl)
      prm%omega       = math_expand(pl%get_as1dReal('omega',       requiredSize=size(N_sl)),N_sl)
      prm%B           = math_expand(pl%get_as1dReal('B',           requiredSize=size(N_sl)),N_sl)
      prm%d_caron     = prm%b_sl *  pl%get_asReal('D_a')
      prm%f_at        = prm%b_sl**3*pl%get_asReal('f_at')

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'), &
                                                   phase_lattice(ph))

      prm%forestProjection = spread(          f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_edge (N_sl,phase_lattice(ph),phase_cOverA(ph)) &
                           + spread(1.0_pREAL-f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_screw(N_sl,phase_lattice(ph),phase_cOverA(ph))

      ! sanity checks
      if (    prm%D_0          <  0.0_pREAL)  extmsg = trim(extmsg)//' D_0'
      if (    prm%Q_cl         <= 0.0_pREAL)  extmsg = trim(extmsg)//' Q_cl'
      if (any(rho_mob_0        <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_mob_0'
      if (any(rho_dip_0        <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_dip_0'
      if (any(prm%b_sl         <= 0.0_pREAL)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%Q_s          <= 0.0_pREAL)) extmsg = trim(extmsg)//' Q_s'
      if (any(prm%tau_Peierls  <  0.0_pREAL)) extmsg = trim(extmsg)//' tau_Peierls'
      if (any(prm%B            <  0.0_pREAL)) extmsg = trim(extmsg)//' B'
      if (any(prm%d_caron      <  0.0_pREAL)) extmsg = trim(extmsg)//' d_caron(D_a,b_sl)'
      if (any(prm%f_at         <= 0.0_pREAL)) extmsg = trim(extmsg)//' f_at or b_sl'

    else slipActive
      rho_mob_0 = emptyRealArray
      rho_dip_0 = emptyRealArray
      allocate(prm%b_sl, &
               prm%d_caron, &
               prm%i_sl, &
               prm%f_at, &
               prm%tau_Peierls, &
               prm%Q_s, &
               prm%p, &
               prm%q, &
               prm%B, &
               prm%h, &
               prm%w, &
               prm%omega, &
               source = emptyRealArray)
      allocate(prm%forestProjection(0,0))
      allocate(prm%h_sl_sl         (0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_ID_phase == ph)
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%rho_mob = [startIndex,endIndex]
    stt%rho_mob => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_mob = spread(rho_mob_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%rho_dip = [startIndex,endIndex]
    stt%rho_dip => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_dip = spread(rho_dip_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    allocate(dst%Lambda_sl(prm%sum_N_sl,Nmembers), source=0.0_pREAL)
    allocate(dst%tau_pass (prm%sum_N_sl,Nmembers), source=0.0_pREAL)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_dislotungsten_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine dislotungsten_LpAndItsTangent(Lp,dLp_dMp, &
                                                     Mp,ph,en)
  real(pREAL), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress
  real(pREAL), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                     intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pREAL) :: &
    T                                                                                               !< temperature
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos,dot_gamma_neg, &
    ddot_gamma_dtau_pos,ddot_gamma_dtau_neg


  T = thermal_T(ph,en)
  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph))

    call kinetics(Mp,T,ph,en,dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg)
    do i = 1, prm%sum_N_sl
      Lp = Lp + (dot_gamma_pos(i)+dot_gamma_neg(i))*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_pos(i) * prm%P_sl(k,l,i) * prm%P_nS_pos(m,n,i) &
                         + ddot_gamma_dtau_neg(i) * prm%P_sl(k,l,i) * prm%P_nS_neg(m,n,i)
    end do

  end associate

end subroutine dislotungsten_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function dislotungsten_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos, dot_gamma_neg,&
    tau_pos,&
    tau_neg, &
    v_cl, &
    dot_rho_dip_formation, &
    dot_rho_dip_climb, &
    d_hat
  real(pREAL) :: &
    mu, T


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph), &
            dot_rho_mob => dotState(indexDotState(ph)%rho_mob(1):indexDotState(ph)%rho_mob(2)), &
            dot_rho_dip => dotState(indexDotState(ph)%rho_dip(1):indexDotState(ph)%rho_dip(2)), &
            dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)))

    mu = elastic_mu(ph,en,prm%isotropic_bound)
    T = thermal_T(ph,en)

    call kinetics(Mp,T,ph,en,&
                  dot_gamma_pos,dot_gamma_neg, &
                  tau_pos_out = tau_pos,tau_neg_out = tau_neg)

    dot_gamma_sl = abs(dot_gamma_pos+dot_gamma_neg)

    where(dEq0((tau_pos+tau_neg)*0.5_pREAL))
      dot_rho_dip_formation = 0.0_pREAL
      dot_rho_dip_climb     = 0.0_pREAL
    else where
      d_hat = math_clip(3.0_pREAL*mu*prm%b_sl/(16.0_pREAL*PI*abs(tau_pos+tau_neg)*0.5_pREAL), &
                        prm%d_caron, &                                                              ! lower limit
                        dst%Lambda_sl(:,en))                                                        ! upper limit
      dot_rho_dip_formation = merge(2.0_pREAL*(d_hat-prm%d_caron)*stt%rho_mob(:,en)*dot_gamma_sl/prm%b_sl, &
                                    0.0_pREAL, &
                                    prm%dipoleformation)
      v_cl = (3.0_pREAL*mu*prm%D_0*exp(-prm%Q_cl/(K_B*T))*prm%f_at/(TAU*K_B*T)) &
           * (1.0_pREAL/(d_hat+prm%d_caron))
      dot_rho_dip_climb = (4.0_pREAL*v_cl*stt%rho_dip(:,en))/(d_hat-prm%d_caron)                    ! ToDo: Discuss with Franz: Stress dependency?
    end where

    dot_rho_mob = dot_gamma_sl/(prm%b_sl*dst%Lambda_sl(:,en)) &                                     ! multiplication
                      - dot_rho_dip_formation &
                - (2.0_pREAL*prm%d_caron)/prm%b_sl*stt%rho_mob(:,en)*dot_gamma_sl                   ! Spontaneous annihilation of 2 edges
    dot_rho_dip = dot_rho_dip_formation &
                - (2.0_pREAL*prm%d_caron)/prm%b_sl*stt%rho_dip(:,en)*dot_gamma_sl &                 ! Spontaneous annihilation of an edge with a dipole
                      - dot_rho_dip_climb

  end associate

end function dislotungsten_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine dislotungsten_dependentState(ph,en)

  integer, intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    Lambda_sl_inv


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    dst%tau_pass(:,en) = elastic_mu(ph,en,prm%isotropic_bound)*prm%b_sl &
                       * sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,en)+stt%rho_dip(:,en)))

    Lambda_sl_inv = 1.0_pREAL/prm%D &
                  + sqrt(matmul(prm%forestProjection,stt%rho_mob(:,en)+stt%rho_dip(:,en)))/prm%i_sl
    dst%Lambda_sl(:,en) = Lambda_sl_inv**(-1.0_pREAL)

  end associate

end subroutine dislotungsten_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislotungsten_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_mob')
          call result_writeDataset(stt%rho_mob,group,trim(prm%output(ou)), &
                                   'mobile dislocation density','1/m²',prm%systems_sl)
        case('rho_dip')
          call result_writeDataset(stt%rho_dip,group,trim(prm%output(ou)), &
                                   'dislocation dipole density','1/m²',prm%systems_sl)
        case('gamma_sl')
          call result_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)
        case('Lambda_sl')
          call result_writeDataset(dst%Lambda_sl,group,trim(prm%output(ou)), &
                                   'mean free path for slip','m',prm%systems_sl)
        case('tau_pass')
          call result_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                   'threshold stress for slip','Pa',prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_dislotungsten_result


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,T,ph,en, &
                 dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg,tau_pos_out,tau_neg_out)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pREAL),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en

  real(pREAL),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos, &
    dot_gamma_neg
  real(pREAL),                  intent(out), optional, dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau_pos, &
    ddot_gamma_dtau_neg, &
    tau_pos_out, &
    tau_neg_out

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    StressRatio, &
    StressRatio_p,StressRatio_pminus1, &
    dvel, &
    tau_pos, tau_neg, tau_eff, &
    t_n, t_k, dtk,dtn
  integer :: j


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do j = 1, prm%sum_N_sl
      tau_pos(j) = math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,j))
      tau_neg(j) = math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,j))
    end do

    if (present(tau_pos_out)) tau_pos_out = tau_pos
    if (present(tau_neg_out)) tau_neg_out = tau_neg

    associate(BoltzmannRatio  => prm%Q_s/(K_B*T), &
              b_rho_half      => stt%rho_mob(:,en) * prm%b_sl * 0.5_pREAL, &
              effectiveLength => dst%Lambda_sl(:,en) - prm%w)

      tau_eff = abs(tau_pos)-dst%tau_pass(:,en)

      significantPositiveTau: where(tau_eff > tol_math_check)
        StressRatio = tau_eff/prm%tau_Peierls
        StressRatio_p       = StressRatio** prm%p
        StressRatio_pminus1 = StressRatio**(prm%p-1.0_pREAL)

        t_n = prm%b_sl*exp(BoltzmannRatio*(1.0_pREAL-StressRatio_p) ** prm%q) &
            / (prm%omega*effectiveLength)
        t_k = effectiveLength * prm%B /(2.0_pREAL*prm%b_sl*tau_eff)                                 ! corrected eq. (14)

        dot_gamma_pos = b_rho_half * sign(prm%h/(t_n + t_k),tau_pos)
      else where significantPositiveTau
        dot_gamma_pos = 0.0_pREAL
      end where significantPositiveTau

      if (present(ddot_gamma_dtau_pos)) then
        significantPositiveTau2: where(abs(tau_pos)-dst%tau_pass(:,en) > tol_math_check)
          dtn = -1.0_pREAL * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pREAL-StressRatio_p)**(prm%q - 1.0_pREAL) &
              * StressRatio_pminus1 / prm%tau_Peierls
          dtk = -1.0_pREAL * t_k / tau_pos

          dvel = -1.0_pREAL * prm%h * (dtk + dtn) / (t_n + t_k)**2

          ddot_gamma_dtau_pos = b_rho_half * dvel
        else where significantPositiveTau2
          ddot_gamma_dtau_pos = 0.0_pREAL
        end where significantPositiveTau2
      end if

      tau_eff = abs(tau_neg)-dst%tau_pass(:,en)

      significantNegativeTau: where(tau_eff > tol_math_check)
        StressRatio = tau_eff/prm%tau_Peierls
        StressRatio_p       = StressRatio** prm%p
        StressRatio_pminus1 = StressRatio**(prm%p-1.0_pREAL)

        t_n = prm%b_sl*exp(BoltzmannRatio*(1.0_pREAL-StressRatio_p) ** prm%q) &
            / (prm%omega*effectiveLength)
        t_k = effectiveLength * prm%B /(2.0_pREAL*prm%b_sl*tau_eff)                                 ! corrected eq. (14)

        dot_gamma_neg = b_rho_half * sign(prm%h/(t_n + t_k),tau_neg)
      else where significantNegativeTau
        dot_gamma_neg = 0.0_pREAL
      end where significantNegativeTau

      if (present(ddot_gamma_dtau_neg)) then
        significantNegativeTau2: where(abs(tau_neg)-dst%tau_pass(:,en) > tol_math_check)
          dtn = -1.0_pREAL * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pREAL-StressRatio_p)**(prm%q - 1.0_pREAL) &
              * StressRatio_pminus1 / prm%tau_Peierls
          dtk = -1.0_pREAL * t_k / tau_neg

          dvel = -1.0_pREAL * prm%h * (dtk + dtn) / (t_n + t_k)**2

          ddot_gamma_dtau_neg = b_rho_half * dvel
        else where significantNegativeTau2
          ddot_gamma_dtau_neg = 0.0_pREAL
        end where significantNegativeTau2
      end if

    end associate
  end associate

end subroutine kinetics

end submodule dislotungsten
