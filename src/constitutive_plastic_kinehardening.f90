!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power law formulation for the shear rates
!! and a Voce-type kinematic hardening rule
!--------------------------------------------------------------------------------------------------
submodule(constitutive:constitutive_mech) plastic_kinehardening

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
      sum_N_sl, &                                                                                   !< total number of active slip system
      of_debug = 0
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
    Ninstances, &
    p, i, o,  &
    Nconstituents, &
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

  print'(/,a)', ' <<<+-  plastic_kinehardening init  -+>>>'

  myPlasticity = plastic_active('kinehardening')
  Ninstances = count(myPlasticity)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  allocate(param(Ninstances))
  allocate(state(Ninstances))
  allocate(dotState(Ninstances))
  allocate(deltaState(Ninstances))

  phases => config_material%get('phase')
  i = 0
  do p = 1, phases%length
    phase => phases%get(p)
    mech  => phase%get('mechanics')
    if(.not. myPlasticity(p)) cycle
    i = i + 1
    associate(prm => param(i), &
              dot => dotState(i), &
              dlt => deltaState(i), &
              stt => state(i))
    pl  => mech%get('plasticity')

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(pl)
#else
    prm%output = pl%get_asStrings('output',defaultVal=emptyStringArray)
#endif

#ifdef DEBUG
    if  (p==material_phaseAt(debugConstitutive%grain,debugConstitutive%element)) then
      prm%of_debug = material_phasememberAt(debugConstitutive%grain,debugConstitutive%ip,debugConstitutive%element)
    endif
#endif

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl         = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                        phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(trim(phase%get_asString('lattice')) == 'cI') then
        a = pl%get_asFloats('a_nonSchmid',defaultVal = emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(N_sl,a,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(N_sl,a,-1)
      else
        prm%nonSchmid_pos  = prm%P
        prm%nonSchmid_neg  = prm%P
      endif
      prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(N_sl, &
                                                                pl%get_asFloats('h_sl_sl'), &
                                                                phase%get_asString('lattice'))

      xi_0          = pl%get_asFloats('xi_0',       requiredSize=size(N_sl))
      prm%xi_inf_f  = pl%get_asFloats('xi_inf_f',   requiredSize=size(N_sl))
      prm%xi_inf_b  = pl%get_asFloats('xi_inf_b',   requiredSize=size(N_sl))
      prm%h_0_f     = pl%get_asFloats('h_0_f',      requiredSize=size(N_sl))
      prm%h_inf_f   = pl%get_asFloats('h_inf_f',    requiredSize=size(N_sl))
      prm%h_0_b     = pl%get_asFloats('h_0_b',      requiredSize=size(N_sl))
      prm%h_inf_b   = pl%get_asFloats('h_inf_b',    requiredSize=size(N_sl))

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
    Nconstituents = count(material_phaseAt == p) * discretization_nIPs
    sizeDotState   = size(['crss     ','crss_back', 'accshear ']) * prm%sum_N_sl!ToDo: adjust names, ask Philip
    sizeDeltaState = size(['sense ',   'chi0  ',    'gamma0'   ]) * prm%sum_N_sl !ToDo: adjust names
    sizeState = sizeDotState + sizeDeltaState

    call constitutive_allocateState(plasticState(p),Nconstituents,sizeState,sizeDotState,sizeDeltaState)

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%crss => plasticState(p)%state   (startIndex:endIndex,:)
    stt%crss = spread(xi_0, 2, Nconstituents)
    dot%crss => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%crss_back => plasticState(p)%state   (startIndex:endIndex,:)
    dot%crss_back => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_xi',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%accshear => plasticState(p)%state   (startIndex:endIndex,:)
    dot%accshear => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'
    ! global alias
    plasticState(p)%slipRate => plasticState(p)%dotState(startIndex:endIndex,:)

    o = plasticState(p)%offsetDeltaState
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%sense => plasticState(p)%state     (startIndex  :endIndex  ,:)
    dlt%sense => plasticState(p)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%chi0 => plasticState(p)%state     (startIndex  :endIndex  ,:)
    dlt%chi0 => plasticState(p)%deltaState(startIndex-o:endIndex-o,:)

    startIndex = endIndex + 1
    endIndex   = endIndex +  prm%sum_N_sl
    stt%gamma0 => plasticState(p)%state     (startIndex  :endIndex  ,:)
    dlt%gamma0 => plasticState(p)%deltaState(startIndex-o:endIndex-o,:)

    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(kinehardening)')

  enddo


end function plastic_kinehardening_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine plastic_kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)

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
    gdot_pos,gdot_neg, &
    dgdot_dtau_pos,dgdot_dtau_neg

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(instance))

  call kinetics(Mp,instance,of,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)

  do i = 1, prm%sum_N_sl
    Lp = Lp + (gdot_pos(i)+gdot_neg(i))*prm%P(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + dgdot_dtau_pos(i) * prm%P(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + dgdot_dtau_neg(i) * prm%P(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo

  end associate

end subroutine plastic_kinehardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_dotState(Mp,instance,of)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal) :: &
    sumGamma
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos,gdot_neg


  associate(prm => param(instance), stt => state(instance), dot => dotState(instance))

  call kinetics(Mp,instance,of,gdot_pos,gdot_neg)
  dot%accshear(:,of) = abs(gdot_pos+gdot_neg)
  sumGamma = sum(stt%accshear(:,of))


  dot%crss(:,of) = matmul(prm%interaction_SlipSlip,dot%accshear(:,of)) &
                 * (  prm%h_inf_f &
                     + (prm%h_0_f - prm%h_inf_f + prm%h_0_f*prm%h_inf_f*sumGamma/prm%xi_inf_f) &
                     * exp(-sumGamma*prm%h_0_f/prm%xi_inf_f) &
                   )

  dot%crss_back(:,of) = stt%sense(:,of)*dot%accshear(:,of) * &
           ( prm%h_inf_b + &
             (prm%h_0_b - prm%h_inf_b &
               + prm%h_0_b*prm%h_inf_b/(prm%xi_inf_b+stt%chi0(:,of))*(stt%accshear(:,of)-stt%gamma0(:,of))&
             ) *exp(-(stt%accshear(:,of)-stt%gamma0(:,of)) *prm%h_0_b/(prm%xi_inf_b+stt%chi0(:,of))) &
           )

  end associate

end subroutine plastic_kinehardening_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate (instantaneous) incremental change of microstructure.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_deltaState(Mp,instance,of)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    instance, &
    of

  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos,gdot_neg, &
    sense

  associate(prm => param(instance), stt => state(instance), dlt => deltaState(instance))

  call kinetics(Mp,instance,of,gdot_pos,gdot_neg)
  sense = merge(state(instance)%sense(:,of), &                                                      ! keep existing...
                sign(1.0_pReal,gdot_pos+gdot_neg), &                                                ! ...or have a defined
                dEq0(gdot_pos+gdot_neg,1e-10_pReal))                                                ! current sense of shear direction

#ifdef DEBUG
  if (debugConstitutive%extensive &
             .and. (of == prm%of_debug  .or. .not. debugConstitutive%selective)) then
    print*, '======= kinehardening delta state ======='
    print*, sense,state(instance)%sense(:,of)
  endif
#endif

!--------------------------------------------------------------------------------------------------
! switch in sense of shear?
  where(dNeq(sense,stt%sense(:,of),0.1_pReal))
    dlt%sense (:,of) = sense - stt%sense(:,of)                                                      ! switch sense
    dlt%chi0  (:,of) = abs(stt%crss_back(:,of)) - stt%chi0(:,of)                                    ! remember current backstress magnitude
    dlt%gamma0(:,of) = stt%accshear(:,of) - stt%gamma0(:,of)                                        ! remember current accumulated shear
  else where
    dlt%sense (:,of) = 0.0_pReal
    dlt%chi0  (:,of) = 0.0_pReal
    dlt%gamma0(:,of) = 0.0_pReal
  end where

  end associate

end subroutine plastic_kinehardening_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_kinehardening_results(instance,group)

  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
     case('xi')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%crss,trim(prm%output(o)), &
                                                    'resistance against plastic slip','Pa')
     case('tau_b')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%crss_back,trim(prm%output(o)), &
                                                    'back stress against plastic slip','Pa')
     case ('sgn(gamma)')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%sense,trim(prm%output(o)), & ! ToDo: could be int
                                                    'tbd','1')
     case ('chi_0')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%chi0,trim(prm%output(o)), &
                                                    'tbd','Pa')
     case ('gamma_0')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma0,trim(prm%output(o)), &
                                                    'tbd','1')
     case ('gamma')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%accshear,trim(prm%output(o)), &
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
pure subroutine kinetics(Mp,instance,of, &
                         gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                intent(in) :: &
    instance, &
    of

  real(pReal),                  intent(out), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos, &
    gdot_neg
  real(pReal),                  intent(out), optional, dimension(param(instance)%sum_N_sl) :: &
    dgdot_dtau_pos, &
    dgdot_dtau_neg

  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    tau_pos, &
    tau_neg
  integer :: i

  associate(prm => param(instance), stt => state(instance))

  do i = 1, prm%sum_N_sl
    tau_pos(i) =       math_tensordot(Mp,prm%nonSchmid_pos(1:3,1:3,i)) - stt%crss_back(i,of)
    tau_neg(i) = merge(math_tensordot(Mp,prm%nonSchmid_neg(1:3,1:3,i)) - stt%crss_back(i,of), &
                       0.0_pReal, prm%nonSchmidActive)
  enddo

  where(dNeq0(tau_pos))
    gdot_pos = prm%dot_gamma_0 * merge(0.5_pReal,1.0_pReal, prm%nonSchmidActive) &                  ! 1/2 if non-Schmid active
             * sign(abs(tau_pos/stt%crss(:,of))**prm%n,  tau_pos)
  else where
    gdot_pos = 0.0_pReal
  end where

  where(dNeq0(tau_neg))
    gdot_neg = prm%dot_gamma_0 * 0.5_pReal &                                                        ! only used if non-Schmid active, always 1/2
             * sign(abs(tau_neg/stt%crss(:,of))**prm%n,  tau_neg)
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

end submodule plastic_kinehardening
