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
      interaction_slipslip                                                                          !< slip resistance from slip activity
    real(pReal),              allocatable, dimension(:,:,:) :: &
      P, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl
    logical :: &
      nonSchmidActive = .false.
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type :: tKinehardeningState
    real(pReal), pointer, dimension(:,:) :: &                                                       !< vectors along NipcMyInstance
      crss, &                                                                                       !< critical resolved stress
      crss_back, &                                                                                  !< critical resolved back stress
      sense, &                                                                                      !< sense of acting shear stress (-1 or +1)
      chi0, &                                                                                       !< backstress at last switch of stress sense
      gamma0, &                                                                                     !< accumulated shear at last switch of stress sense
      accshear                                                                                      !< accumulated (absolute) shear
  end type tKinehardeningState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),         allocatable, dimension(:) :: param
  type(tKinehardeningState), allocatable, dimension(:) :: &
    dotState, &
    deltaState, &
    state

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

  print'(/,a)', ' <<<+-  phase:mechanical:plastic:kinehardening init  -+>>>'
  print'(a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)


  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(dotState(phases%length))
  allocate(deltaState(phases%length))


  do ph = 1, phases%length
    if(.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), dot => dotState(ph), dlt => deltaState(ph), stt => state(ph))

    phase => phases%get(ph)
    mech  => phase%get('mechanical')
    pl  => mech%get('plastic')

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
      prm%P = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                        phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(trim(phase%get_asString('lattice')) == 'cI') then
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal = emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos  = prm%P
        prm%nonSchmid_neg  = prm%P
      endif
      prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(N_sl, &
                                                                pl%get_as1dFloat('h_sl_sl'), &
                                                                phase%get_asString('lattice'))

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

      !ToDo: Any sensible checks for theta?
    else slipActive
      xi_0 = emptyRealArray
      allocate(prm%xi_inf_f,prm%xi_inf_b,prm%h_0_f,prm%h_inf_f,prm%h_0_b,prm%h_inf_b,source=emptyRealArray)
      allocate(prm%interaction_SlipSlip(0,0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_phaseID == ph)
    sizeDotState   = size(['crss     ','crss_back', 'accshear ']) * prm%sum_N_sl !ToDo: adjust names like in material.yaml
    sizeDeltaState = size(['sense ',   'chi0  ',    'gamma0'   ]) * prm%sum_N_sl !ToDo: adjust names like in material.yaml
    sizeState = sizeDotState + sizeDeltaState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,sizeDeltaState)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%crss => plasticState(ph)%state   (startIndex:endIndex,:)
    stt%crss = spread(xi_0, 2, Nmembers)
    dot%crss => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%crss_back => plasticState(ph)%state   (startIndex:endIndex,:)
    dot%crss_back => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%accshear => plasticState(ph)%state   (startIndex:endIndex,:)
    dot%accshear => plasticState(ph)%dotState(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'
    ! global alias
    plasticState(ph)%slipRate => plasticState(ph)%dotState(startIndex:endIndex,:)

    o = plasticState(ph)%offsetDeltaState
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%sense => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%sense => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%chi0 => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%chi0 => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%gamma0 => plasticState(ph)%state     (startIndex  :endIndex  ,:)
    dlt%gamma0 => plasticState(ph)%deltaState(startIndex-o:endIndex-o,:)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(kinehardening)')

  enddo


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
    gdot_pos,gdot_neg, &
    dgdot_dtau_pos,dgdot_dtau_neg

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph))

  call kinetics(Mp,ph,en,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)

  do i = 1, prm%sum_N_sl
    Lp = Lp + (gdot_pos(i)+gdot_neg(i))*prm%P(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtau_pos(i) * prm%P(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + dgdot_dtau_neg(i) * prm%P(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo

  end associate

end subroutine kinehardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_dotState(Mp,ph,en)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pReal) :: &
    sumGamma
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    gdot_pos,gdot_neg


  associate(prm => param(ph), stt => state(ph),&
            dot => dotState(ph))

  call kinetics(Mp,ph,en,gdot_pos,gdot_neg)
  dot%accshear(:,en) = abs(gdot_pos+gdot_neg)
  sumGamma = sum(stt%accshear(:,en))


  dot%crss(:,en) = matmul(prm%interaction_SlipSlip,dot%accshear(:,en)) &
                 * (  prm%h_inf_f &
                     + (prm%h_0_f - prm%h_inf_f + prm%h_0_f*prm%h_inf_f*sumGamma/prm%xi_inf_f) &
                     * exp(-sumGamma*prm%h_0_f/prm%xi_inf_f) &
                   )

  dot%crss_back(:,en) = stt%sense(:,en)*dot%accshear(:,en) * &
           ( prm%h_inf_b + &
             (prm%h_0_b - prm%h_inf_b &
               + prm%h_0_b*prm%h_inf_b/(prm%xi_inf_b+stt%chi0(:,en))*(stt%accshear(:,en)-stt%gamma0(:,en))&
             ) *exp(-(stt%accshear(:,en)-stt%gamma0(:,en)) *prm%h_0_b/(prm%xi_inf_b+stt%chi0(:,en))) &
           )

  end associate

end subroutine plastic_kinehardening_dotState


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
    gdot_pos,gdot_neg, &
    sense

  associate(prm => param(ph), stt => state(ph), dlt => deltaState(ph))

  call kinetics(Mp,ph,en,gdot_pos,gdot_neg)
  sense = merge(state(ph)%sense(:,en), &                                                            ! keep existing...
                sign(1.0_pReal,gdot_pos+gdot_neg), &                                                ! ...or have a defined
                dEq0(gdot_pos+gdot_neg,1e-10_pReal))                                                ! current sense of shear direction


!--------------------------------------------------------------------------------------------------
! switch in sense of shear?
  where(dNeq(sense,stt%sense(:,en),0.1_pReal))
    dlt%sense (:,en) = sense - stt%sense(:,en)                                                      ! switch sense
    dlt%chi0  (:,en) = abs(stt%crss_back(:,en)) - stt%chi0(:,en)                                    ! remember current backstress magnitude
    dlt%gamma0(:,en) = stt%accshear(:,en) - stt%gamma0(:,en)                                        ! remember current accumulated shear
  else where
    dlt%sense (:,en) = 0.0_pReal
    dlt%chi0  (:,en) = 0.0_pReal
    dlt%gamma0(:,en) = 0.0_pReal
  end where

  end associate

end subroutine plastic_kinehardening_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(ph), stt => state(ph))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
     case ('xi')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%crss,group,trim(prm%output(o)), &
                                                    'resistance against plastic slip','Pa')
     case ('tau_b')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%crss_back,group,trim(prm%output(o)), &
                                                    'back stress against plastic slip','Pa')
     case ('sgn(gamma)')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%sense,group,trim(prm%output(o)), & ! ToDo: could be int
                                                    'sense of shear','1')
     case ('chi_0')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%chi0,group,trim(prm%output(o)), &
                                                    'tbd','Pa')
     case ('gamma_0')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%gamma0,group,trim(prm%output(o)), &
                                                    'tbd','1')
     case ('gamma')
       if(prm%sum_N_sl>0) call results_writeDataset(stt%accshear,group,trim(prm%output(o)), &
                                                    'plastic shear','1')
    end select
  enddo outputsLoop
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
                         gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                intent(in) :: &
    ph, &
    en

  real(pReal),                  intent(out), dimension(param(ph)%sum_N_sl) :: &
    gdot_pos, &
    gdot_neg
  real(pReal),                  intent(out), optional, dimension(param(ph)%sum_N_sl) :: &
    dgdot_dtau_pos, &
    dgdot_dtau_neg

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau_pos, &
    tau_neg
  integer :: i

  associate(prm => param(ph), stt => state(ph))

  do i = 1, prm%sum_N_sl
    tau_pos(i) =       math_tensordot(Mp,prm%nonSchmid_pos(1:3,1:3,i)) - stt%crss_back(i,en)
    tau_neg(i) = merge(math_tensordot(Mp,prm%nonSchmid_neg(1:3,1:3,i)) - stt%crss_back(i,en), &
                       0.0_pReal, prm%nonSchmidActive)
  enddo

  where(dNeq0(tau_pos))
    gdot_pos = prm%dot_gamma_0 * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &                  ! 1/2 if non-Schmid active
             * sign(abs(tau_pos/stt%crss(:,en))**prm%n,  tau_pos)
  else where
    gdot_pos = 0.0_pReal
  end where

  where(dNeq0(tau_neg))
    gdot_neg = prm%dot_gamma_0 * 0.5_pReal &                                                        ! only used if non-Schmid active, always 1/2
             * sign(abs(tau_neg/stt%crss(:,en))**prm%n,  tau_neg)
  else where
    gdot_neg = 0.0_pReal
  end where

  if (present(dgdot_dtau_pos)) then
    where(dNeq0(gdot_pos))
      dgdot_dtau_pos = gdot_pos*prm%n/tau_pos
    else where
      dgdot_dtau_pos = 0.0_pReal
    end where
  endif
  if (present(dgdot_dtau_neg)) then
    where(dNeq0(gdot_neg))
      dgdot_dtau_neg = gdot_neg*prm%n/tau_neg
    else where
      dgdot_dtau_neg = 0.0_pReal
    end where
  endif
  end associate

end subroutine kinetics

end submodule kinehardening
