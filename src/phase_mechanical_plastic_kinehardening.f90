!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power-law formulation for the shear rates
!! and a Voce-type kinematic hardening rule.
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) kinehardening

  type :: tParameters
    real(pREAL),              allocatable, dimension(:) :: &
      dot_gamma_0, &                                                                                !< reference shear strain rate for slip
      n, &                                                                                          !< stress exponent for slip
      h_0_xi, &                                                                                     !< initial hardening rate of forest stress per slip family
                                                                                                    !! θ_0,for
      h_0_chi, &                                                                                    !< initial hardening rate of back stress per slip family
                                                                                                    !! θ_0,bs
      h_inf_xi, &                                                                                   !< asymptotic hardening rate of forest stress per slip family
                                                                                                    !! θ_1,for
      h_inf_chi, &                                                                                  !< asymptotic hardening rate of back stress per slip family
                                                                                                    !! θ_1,bs
      xi_inf, &                                                                                     !< back-extrapolated forest stress from terminal linear hardening
      chi_inf                                                                                       !< back-extrapolated back stress from terminal linear hardening
    real(pREAL),              allocatable, dimension(:,:) :: &
      h_sl_sl                                                                                       !< slip resistance change per slip activity
    real(pREAL),              allocatable, dimension(:,:,:) :: &
      P, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
    character(len=:),       allocatable, dimension(:) :: &
      systems_sl
  end type tParameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      xi, &
      chi, &
      gamma
  end type tIndexDotState

  type :: tKinehardeningState
    real(pREAL), pointer, dimension(:,:) :: &
      xi, &                                                                                         !< forest stress
                                                                                                    !! τ_for
      chi, &                                                                                        !< back stress
                                                                                                    !! τ_bs
      chi_flip, &                                                                                   !< back stress at last reversal of stress sense
                                                                                                    !! χ_0
      gamma, &                                                                                      !< accumulated (absolute) shear
      gamma_flip, &                                                                                 !< accumulated shear at last reversal of stress sense
                                                                                                    !! γ_0
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
  real(pREAL), dimension(:), allocatable :: &
    xi_0                                                                                            !< initial forest stress
                                                                                                    !! τ_for,0
  real(pREAL), dimension(:,:), allocatable :: &
    a_nS                                                                                            !< non-Schmid coefficients
  character(len=:), allocatable :: &
    refs, &
    extmsg
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('kinehardening')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:kinehardening init  -+>>>'

  print'(/,1x,a)', 'J.A. Wollmershauser et al., International Journal of Fatigue 36:181–193, 2012'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijfatigue.2011.07.008'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  allocate(deltaState(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), &
              stt => state(ph), dlt => deltaState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech => phase%get_dict('mechanical')
    pl => mech%get_dict('plastic')

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(pl,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(pl)
#else
    prm%output = pl%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))
      prm%P = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        allocate(a_nS(3,size(pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray))),source=0.0_pREAL)
        a_nS(1,:) = pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray)
        prm%P_nS_pos = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=+1)
        prm%P_nS_neg = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=-1)
        deallocate(a_nS)
      else
        prm%P_nS_pos = +prm%P
        prm%P_nS_neg = -prm%P
      end if

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'), &
                                                   phase_lattice(ph))

      xi_0            = math_expand(pl%get_as1dReal('xi_0',        requiredSize=size(N_sl)),N_sl)
      prm%dot_gamma_0 = math_expand(pl%get_as1dReal('dot_gamma_0', requiredSize=size(N_sl)),N_sl)
      prm%n           = math_expand(pl%get_as1dReal('n',           requiredSize=size(N_sl)),N_sl)
      prm%xi_inf      = math_expand(pl%get_as1dReal('xi_inf',      requiredSize=size(N_sl)),N_sl)
      prm%chi_inf     = math_expand(pl%get_as1dReal('chi_inf',     requiredSize=size(N_sl)),N_sl)
      prm%h_0_xi      = math_expand(pl%get_as1dReal('h_0_xi',      requiredSize=size(N_sl)),N_sl)
      prm%h_0_chi     = math_expand(pl%get_as1dReal('h_0_chi',     requiredSize=size(N_sl)),N_sl)
      prm%h_inf_xi    = math_expand(pl%get_as1dReal('h_inf_xi',    requiredSize=size(N_sl)),N_sl)
      prm%h_inf_chi   = math_expand(pl%get_as1dReal('h_inf_chi',   requiredSize=size(N_sl)),N_sl)

!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (any(prm%dot_gamma_0  <= 0.0_pREAL))  extmsg = trim(extmsg)//' dot_gamma_0'
      if (any(prm%n            <= 0.0_pREAL))  extmsg = trim(extmsg)//' n'
      if (any(xi_0             <= 0.0_pREAL))  extmsg = trim(extmsg)//' xi_0'
      if (any(prm%xi_inf       <= 0.0_pREAL))  extmsg = trim(extmsg)//' xi_inf'
      if (any(prm%chi_inf      <= 0.0_pREAL))  extmsg = trim(extmsg)//' chi_inf'

    else slipActive
      xi_0 = emptyRealArray
      allocate(prm%dot_gamma_0, &
               prm%n, &
               prm%xi_inf, &
               prm%chi_inf, &
               prm%h_0_xi, &
               prm%h_0_chi, &
               prm%h_inf_xi, &
               prm%h_inf_chi, &
               source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_ID_phase == ph)
    sizeDotState   = prm%sum_N_sl * size(['xi   ',&
                                          'chi  ',&
                                          'gamma'])
    sizeDeltaState = prm%sum_N_sl * size(['sgn_gamma ',&
                                          'chi_flip  ',&
                                          'gamma_flip'])
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
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_xi',defaultVal=1.0_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%chi = [startIndex,endIndex]
    stt%chi => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_xi',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma = [startIndex,endIndex]
    stt%gamma => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    o = plasticState(ph)%offsetDeltaState
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%sgn_gamma => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%sgn_gamma => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%chi_flip => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%chi_flip => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%gamma_flip => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%gamma_flip => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do


end function plastic_kinehardening_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine kinehardening_LpAndItsTangent(Lp,dLp_dMp, Mp,ph,en)

  real(pREAL), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pREAL), dimension(3,3),     intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                         intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma, ddot_gamma_dtau
  real(pREAL), dimension(3,3,param(ph)%sum_N_sl) :: &
    P_nS


  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph))

    call kinetics(Mp,ph,en, dot_gamma,ddot_gamma_dtau)
    P_nS = merge(prm%P_nS_pos,prm%P_nS_neg, spread(spread(dot_gamma,1,3),2,3)>0.0_pREAL)            ! faster than 'merge' in loop
    do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma(i)*prm%P(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau(i) * prm%P(k,l,i) * P_nS(m,n,i)
    end do

  end associate

end subroutine kinehardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function plastic_kinehardening_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pREAL) :: &
    sumGamma


  associate(prm => param(ph), stt => state(ph), &
            dot_xi => dotState(IndexDotState(ph)%xi(1):IndexDotState(ph)%xi(2)),&
            dot_chi => dotState(IndexDotState(ph)%chi(1):IndexDotState(ph)%chi(2)),&
            dot_gamma => dotState(IndexDotState(ph)%gamma(1):IndexDotState(ph)%gamma(2)))

    call kinetics(Mp,ph,en, dot_gamma)
    sumGamma = sum(stt%gamma(:,en))
    dot_gamma = abs(dot_gamma)


    dot_xi = matmul(prm%h_sl_sl,dot_gamma) &
           * ( prm%h_inf_xi &
               + ( prm%h_0_xi  &
                 - prm%h_inf_xi * (1_pREAL -sumGamma*prm%h_0_xi/prm%xi_inf) ) &
               *                       exp(-sumGamma*prm%h_0_xi/prm%xi_inf) &
             )

    dot_chi = stt%sgn_gamma(:,en)*dot_gamma &
            * ( prm%h_inf_chi &
               + ( prm%h_0_chi &
                 - prm%h_inf_chi*(1_pREAL -(stt%gamma(:,en)-stt%gamma_flip(:,en))*prm%h_0_chi/(prm%chi_inf+stt%chi_flip(:,en))) ) &
               *                      exp(-(stt%gamma(:,en)-stt%gamma_flip(:,en))*prm%h_0_chi/(prm%chi_inf+stt%chi_flip(:,en))) &
              )

  end associate

end function plastic_kinehardening_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate (instantaneous) incremental change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_deltaState(Mp,ph,en)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma, &
    sgn_gamma


  associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph))

    call kinetics(Mp,ph,en, dot_gamma)
    sgn_gamma = merge(state(ph)%sgn_gamma(:,en), &
                      sign(1.0_pREAL,dot_gamma), &
                      dEq0(dot_gamma,1e-10_pREAL))

    where(dNeq(sgn_gamma,stt%sgn_gamma(:,en),0.1_pREAL)) ! ToDo sgn_gamma*stt%sgn_gamma(:,en)<0
      dlt%sgn_gamma (:,en) = sgn_gamma            - stt%sgn_gamma (:,en)
      dlt%chi_flip  (:,en) = abs(stt%chi  (:,en)) - stt%chi_flip  (:,en)
      dlt%gamma_flip(:,en) =     stt%gamma(:,en)  - stt%gamma_flip(:,en)
    else where
      dlt%sgn_gamma (:,en) = 0.0_pREAL
      dlt%chi_flip  (:,en) = 0.0_pREAL
      dlt%gamma_flip(:,en) = 0.0_pREAL
    end where

  end associate

end subroutine plastic_kinehardening_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case ('xi')
          call result_writeDataset(stt%xi,group,trim(prm%output(ou)), &
                                   'forest stress','Pa',prm%systems_sl)
        case ('chi')
          call result_writeDataset(stt%chi,group,trim(prm%output(ou)), &
                                   'back stress','Pa',prm%systems_sl)
        case ('sgn(gamma)')
          call result_writeDataset(int(stt%sgn_gamma),group,trim(prm%output(ou)), &
                                   'sense of shear','1',prm%systems_sl)
        case ('chi_flip')
          call result_writeDataset(stt%chi_flip,group,trim(prm%output(ou)), &
                                   'back stress at last reversal of stress sense','Pa',prm%systems_sl)
        case ('gamma_flip')
          call result_writeDataset(stt%gamma_flip,group,trim(prm%output(ou)), &
                                   'plastic shear at last reversal of stress sense','1',prm%systems_sl)
        case ('gamma')
          call result_writeDataset(stt%gamma,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_kinehardening_result


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details: Derivatives are calculated only optionally.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,ph,en, &
                         dot_gamma,ddot_gamma_dtau)

  real(pREAL), dimension(3,3),                           intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                                               intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl),           intent(out) :: &
    dot_gamma
  real(pREAL), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau_pos, &
    tau_neg
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    tau_pos = [(math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,i)) - stt%chi(i,en),i=1,prm%sum_N_sl)]
    tau_neg = [(math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,i)) + stt%chi(i,en),i=1,prm%sum_N_sl)]

    dot_gamma = merge(+1.0_pREAL,-1.0_pREAL, tau_pos>tau_neg) &
              * prm%dot_gamma_0  &
              * (max(tau_pos,tau_neg)/stt%xi(:,en))**prm%n

    if (present(ddot_gamma_dtau)) then
      where(dNeq0(dot_gamma))
        ddot_gamma_dtau = dot_gamma*prm%n/max(tau_pos,tau_neg)
      else where
        ddot_gamma_dtau = 0.0_pREAL
      end where
    end if

  end associate

end subroutine kinetics

end submodule kinehardening
