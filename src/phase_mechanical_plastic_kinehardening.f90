!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power law formulation for the shear rates
!! and a Voce-type kinematic hardening rule
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) kinehardening

  type :: tParameters
    real(pReal) :: &
      n           = 1.0_pReal, &                                                                    !< stress exponent for slip
      dot_gamma_0 = 1.0_pReal                                                                       !< reference shear strain rate for slip
    real(pReal),              allocatable, dimension(:) :: &
      h_0_f, &                                                                                      !< initial hardening rate of forward stress for each slip
      h_inf_f, &                                                                                    !< asymptotic hardening rate of forward stress for each slip
      h_0_b, &                                                                                      !< initial hardening rate of back stress for each slip
      h_inf_b, &                                                                                    !< asymptotic hardening rate of back stress for each slip
      xi_inf_f, &
      xi_inf_b
    real(pReal),              allocatable, dimension(:,:) :: &
      h_sl_sl                                                                                       !< slip resistance from slip activity
    real(pReal),              allocatable, dimension(:,:,:) :: &
      P, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl
    logical :: &
      nonSchmidActive = .false.
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl
  end type tParameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      xi, &
      chi, &
      gamma
  end type tIndexDotState

  type :: tKinehardeningState
    real(pReal), pointer, dimension(:,:) :: &
      xi, &                                                                                         !< resistance against plastic slip
      chi, &                                                                                        !< back stress
      chi_0, &                                                                                      !< back stress at last switch of stress sense
      gamma, &                                                                                      !< accumulated (absolute) shear
      gamma_0, &                                                                                    !< accumulated shear at last switch of stress sense
      sgn_gamma                                                                                     !< sense of acting shear stress (-1 or +1)
  end type tKinehardeningState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),         allocatable, dimension(:) :: param
  type(tIndexDotState),      allocatable, dimension(:) :: indexDotState
  type(tKinehardeningState), allocatable, dimension(:) :: state, deltaState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_kinehardening_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, o,  &
    Nmembers, &
    sizeState, sizeDeltaState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl
  real(pReal), dimension(:), allocatable :: &
    xi_0, &                                                                                         !< initial resistance against plastic flow
    a                                                                                               !< non-Schmid coefficients
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl

  myPlasticity = plastic_active('kinehardening')
  if(count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:kinehardening init  -+>>>'
  print'(/,a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)

  print'(/,1x,a)', 'J.A. Wollmershauser et al., International Journal of Fatigue 36:181–193, 2012'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijfatigue.2011.07.008'

  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  allocate(deltaState(phases%length))


  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get(ph)
    mech => phase%get('mechanical')
    pl => mech%get('plastic')

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = lattice_labels_slip(N_sl,phase_lattice(ph))
      prm%P = lattice_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal = emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%P_nS_pos = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%P_nS_neg = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%P_nS_pos = prm%P
        prm%P_nS_neg = prm%P
      end if
      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,pl%get_as1dFloat('h_sl-sl'), &
                                                   phase_lattice(ph))

      xi_0          = pl%get_as1dFloat('xi_0',       requiredSize=size(N_sl))
      prm%xi_inf_f  = pl%get_as1dFloat('xi_inf_f',   requiredSize=size(N_sl))
      prm%xi_inf_b  = pl%get_as1dFloat('xi_inf_b',   requiredSize=size(N_sl))
      prm%h_0_f     = pl%get_as1dFloat('h_0_f',      requiredSize=size(N_sl))
      prm%h_inf_f   = pl%get_as1dFloat('h_inf_f',    requiredSize=size(N_sl))
      prm%h_0_b     = pl%get_as1dFloat('h_0_b',      requiredSize=size(N_sl))
      prm%h_inf_b   = pl%get_as1dFloat('h_inf_b',    requiredSize=size(N_sl))

      prm%dot_gamma_0  = pl%get_asFloat('dot_gamma_0')
      prm%n            = pl%get_asFloat('n')

      ! expand: family => system
      xi_0          = math_expand(xi_0,            N_sl)
      prm%xi_inf_f  = math_expand(prm%xi_inf_f,    N_sl)
      prm%xi_inf_b  = math_expand(prm%xi_inf_b,    N_sl)
      prm%h_0_f     = math_expand(prm%h_0_f,       N_sl)
      prm%h_inf_f   = math_expand(prm%h_inf_f,     N_sl)
      prm%h_0_b     = math_expand(prm%h_0_b,       N_sl)
      prm%h_inf_b   = math_expand(prm%h_inf_b,     N_sl)

!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (    prm%dot_gamma_0  <= 0.0_pReal)   extmsg = trim(extmsg)//' dot_gamma_0'
      if (    prm%n            <= 0.0_pReal)   extmsg = trim(extmsg)//' n'
      if (any(xi_0             <= 0.0_pReal))  extmsg = trim(extmsg)//' xi_0'
      if (any(prm%xi_inf_f     <= 0.0_pReal))  extmsg = trim(extmsg)//' xi_inf_f'
      if (any(prm%xi_inf_b     <= 0.0_pReal))  extmsg = trim(extmsg)//' xi_inf_b'

    else slipActive
      xi_0 = emptyRealArray
      allocate(prm%xi_inf_f,prm%xi_inf_b,prm%h_0_f,prm%h_inf_f,prm%h_0_b,prm%h_inf_b,source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState   = size(['xi   ','chi  ', 'gamma']) * prm%sum_N_sl
    sizeDeltaState = size(['sgn_gamma',   'chi_0    ',    'gamma_0  ']) * prm%sum_N_sl
    sizeState = sizeDotState + sizeDeltaState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,sizeDeltaState)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%xi = [startIndex,endIndex]
    stt%xi => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi = spread(xi_0, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%chi = [startIndex,endIndex]
    stt%chi => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma = [startIndex,endIndex]
    stt%gamma => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    o = plasticState(ph)%offsetDeltaState
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%sgn_gamma => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%sgn_gamma => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%chi_0 => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%chi_0 => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%gamma_0 => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%gamma_0 => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do


end function plastic_kinehardening_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

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
    dot_gamma_pos,dot_gamma_neg, &
    ddot_gamma_dtau_pos,ddot_gamma_dtau_neg

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph))

  call kinetics(Mp,ph,en,dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg)

  do i = 1, prm%sum_N_sl
    Lp = Lp + (dot_gamma_pos(i)+dot_gamma_neg(i))*prm%P(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_pos(i) * prm%P(k,l,i) * prm%P_nS_pos(m,n,i) &
                       + ddot_gamma_dtau_neg(i) * prm%P(k,l,i) * prm%P_nS_neg(m,n,i)
  end do

  end associate

end subroutine kinehardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function plastic_kinehardening_dotState(Mp,ph,en) result(dotState)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pReal) :: &
    sumGamma
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos,dot_gamma_neg


  associate(prm => param(ph), stt => state(ph), &
            dot_xi => dotState(IndexDotState(ph)%xi(1):IndexDotState(ph)%xi(2)),&
            dot_chi => dotState(IndexDotState(ph)%chi(1):IndexDotState(ph)%chi(2)),&
            dot_gamma => dotState(IndexDotState(ph)%gamma(1):IndexDotState(ph)%gamma(2)))

    call kinetics(Mp,ph,en,dot_gamma_pos,dot_gamma_neg)
    dot_gamma = abs(dot_gamma_pos+dot_gamma_neg)
    sumGamma = sum(stt%gamma(:,en))


    dot_xi = matmul(prm%h_sl_sl,dot_gamma) &
                   * (  prm%h_inf_f &
                       + (prm%h_0_f - prm%h_inf_f + prm%h_0_f*prm%h_inf_f*sumGamma/prm%xi_inf_f) &
                       * exp(-sumGamma*prm%h_0_f/prm%xi_inf_f) &
                     )

    dot_chi = stt%sgn_gamma(:,en)*dot_gamma &
            * ( prm%h_inf_b &
               + (prm%h_0_b - prm%h_inf_b &
                 + prm%h_0_b*prm%h_inf_b/(prm%xi_inf_b+stt%chi_0(:,en))*(stt%gamma(:,en)-stt%gamma_0(:,en))&
               ) *exp(-(stt%gamma(:,en)-stt%gamma_0(:,en)) *prm%h_0_b/(prm%xi_inf_b+stt%chi_0(:,en))) &
             )

  end associate

end function plastic_kinehardening_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate (instantaneous) incremental change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_deltaState(Mp,ph,en)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos,dot_gamma_neg, &
    sgn_gamma


  associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph))

    call kinetics(Mp,ph,en,dot_gamma_pos,dot_gamma_neg)
    sgn_gamma = merge(state(ph)%sgn_gamma(:,en), &
                      sign(1.0_pReal,dot_gamma_pos+dot_gamma_neg), &
                      dEq0(dot_gamma_pos+dot_gamma_neg,1e-10_pReal))

    where(dNeq(sgn_gamma,stt%sgn_gamma(:,en),0.1_pReal)) ! ToDo sgn_gamma*stt%sgn_gamma(:,en)<0
      dlt%sgn_gamma (:,en) = sgn_gamma - stt%sgn_gamma(:,en)
      dlt%chi_0  (:,en) = abs(stt%chi(:,en)) - stt%chi_0(:,en)
      dlt%gamma_0(:,en) = stt%gamma(:,en) - stt%gamma_0(:,en)
    else where
      dlt%sgn_gamma (:,en) = 0.0_pReal
      dlt%chi_0  (:,en) = 0.0_pReal
      dlt%gamma_0(:,en) = 0.0_pReal
    end where

  end associate

end subroutine plastic_kinehardening_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case ('xi')
          call results_writeDataset(stt%xi,group,trim(prm%output(ou)), &
                                    'resistance against plastic slip','Pa',prm%systems_sl)
        case ('chi')
          call results_writeDataset(stt%chi,group,trim(prm%output(ou)), &
                                    'back stress','Pa',prm%systems_sl)
        case ('sgn(gamma)')
          call results_writeDataset(int(stt%sgn_gamma),group,trim(prm%output(ou)), &
                                    'sense of shear','1',prm%systems_sl)
        case ('chi_0')
          call results_writeDataset(stt%chi_0,group,trim(prm%output(ou)), &
                                    'back stress at last switch of stress sense','Pa',prm%systems_sl)
        case ('gamma_0')
          call results_writeDataset(stt%gamma_0,group,trim(prm%output(ou)), &
                                    'plastic shear at last switch of stress sense','1',prm%systems_sl)
        case ('gamma')
          call results_writeDataset(stt%gamma,group,trim(prm%output(ou)), &
                                    'plastic shear','1',prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_kinehardening_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details: Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,ph,en, &
                         dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_pos, &
    dot_gamma_neg
  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl), optional :: &
    ddot_gamma_dtau_pos, &
    ddot_gamma_dtau_neg

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau_pos, &
    tau_neg
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    do i = 1, prm%sum_N_sl
      tau_pos(i) =       math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,i)) - stt%chi(i,en)
      tau_neg(i) = merge(math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,i)) - stt%chi(i,en), &
                         0.0_pReal, prm%nonSchmidActive)
    end do

    where(dNeq0(tau_pos))
      dot_gamma_pos = prm%dot_gamma_0 * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &           ! 1/2 if non-Schmid active
                    * sign(abs(tau_pos/stt%xi(:,en))**prm%n,  tau_pos)
    else where
      dot_gamma_pos = 0.0_pReal
    end where

    where(dNeq0(tau_neg))
      dot_gamma_neg = prm%dot_gamma_0 * 0.5_pReal &                                                 ! only used if non-Schmid active, always 1/2
                    * sign(abs(tau_neg/stt%xi(:,en))**prm%n,  tau_neg)
    else where
      dot_gamma_neg = 0.0_pReal
    end where

    if (present(ddot_gamma_dtau_pos)) then
      where(dNeq0(dot_gamma_pos))
        ddot_gamma_dtau_pos = dot_gamma_pos*prm%n/tau_pos
      else where
        ddot_gamma_dtau_pos = 0.0_pReal
      end where
    end if
    if (present(ddot_gamma_dtau_neg)) then
      where(dNeq0(dot_gamma_neg))
        ddot_gamma_dtau_neg = dot_gamma_neg*prm%n/tau_neg
      else where
        ddot_gamma_dtau_neg = 0.0_pReal
      end where
    end if

  end associate

end subroutine kinetics

end submodule kinehardening
