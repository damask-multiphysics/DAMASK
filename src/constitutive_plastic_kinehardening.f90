!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power law formulation for the shear rates
!! and a Voce-type kinematic hardening rule
!--------------------------------------------------------------------------------------------------
submodule(constitutive) plastic_kinehardening

  type :: tParameters
    real(pReal) :: &
      gdot0, &                                                                                      !< reference shear strain rate for slip
      n                                                                                             !< stress exponent for slip
    real(pReal),              allocatable, dimension(:) :: &
      crss0, &                                                                                      !< initial critical shear stress for slip
      theta0, &                                                                                     !< initial hardening rate of forward stress for each slip
      theta1, &                                                                                     !< asymptotic hardening rate of forward stress for each slip
      theta0_b, &                                                                                   !< initial hardening rate of back stress for each slip
      theta1_b, &                                                                                   !< asymptotic hardening rate of back stress for each slip
      tau1, &
      tau1_b, &
      nonSchmidCoeff
    real(pReal),              allocatable, dimension(:,:) :: &
      interaction_slipslip                                                                          !< slip resistance from slip activity
    real(pReal),              allocatable, dimension(:,:,:) :: &
      P, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      of_debug = 0
    integer,                   allocatable, dimension(:) :: &
      Nslip                                                                                         !< number of active slip systems for each family
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
module subroutine plastic_kinehardening_init

  integer :: &
    Ninstance, &
    p, o, &
    NipcMyPhase, &
    sizeState, sizeDeltaState, sizeDotState, &
    startIndex, endIndex

  character(len=pStringLen) :: &
    extmsg = ''

  write(6,'(/,a)') ' <<<+-  plastic_'//PLASTICITY_KINEHARDENING_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_plasticity == PLASTICITY_KINEHARDENING_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))
  allocate(deltaState(Ninstance))

  do p = 1, size(phase_plasticityInstance)
    if (phase_plasticity(p) /= PLASTICITY_KINEHARDENING_ID) cycle
    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              dlt => deltaState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)),&
              config => config_phase(p))

    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

#ifdef DEBUG
    if  (p==material_phaseAt(debug_g,debug_e)) then
      prm%of_debug = material_phasememberAt(debug_g,debug_i,debug_e)
    endif
#endif

!--------------------------------------------------------------------------------------------------
! slip related parameters
    prm%Nslip    = config%getInts('nslip',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(prm%Nslip)
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P = lattice_SchmidMatrix_slip(prm%Nslip,config%getString('lattice_structure'),&
                                        config%getFloat('c/a',defaultVal=0.0_pReal))

      if(trim(config%getString('lattice_structure')) == 'bcc') then
        prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                               defaultVal = emptyRealArray)
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1)
      else
        prm%nonSchmid_pos  = prm%P
        prm%nonSchmid_neg  = prm%P
      endif
      prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(prm%Nslip, &
                                                                config%getFloats('interaction_slipslip'), &
                                                                config%getString('lattice_structure'))

      prm%crss0    = config%getFloats('crss0',    requiredSize=size(prm%Nslip))
      prm%tau1     = config%getFloats('tau1',     requiredSize=size(prm%Nslip))
      prm%tau1_b   = config%getFloats('tau1_b',   requiredSize=size(prm%Nslip))
      prm%theta0   = config%getFloats('theta0',   requiredSize=size(prm%Nslip))
      prm%theta1   = config%getFloats('theta1',   requiredSize=size(prm%Nslip))
      prm%theta0_b = config%getFloats('theta0_b', requiredSize=size(prm%Nslip))
      prm%theta1_b = config%getFloats('theta1_b', requiredSize=size(prm%Nslip))

      prm%gdot0    = config%getFloat('gdot0')
      prm%n        = config%getFloat('n_slip')

      ! expand: family => system
      prm%crss0    = math_expand(prm%crss0,   prm%Nslip)
      prm%tau1     = math_expand(prm%tau1,    prm%Nslip)
      prm%tau1_b   = math_expand(prm%tau1_b,  prm%Nslip)
      prm%theta0   = math_expand(prm%theta0,  prm%Nslip)
      prm%theta1   = math_expand(prm%theta1,  prm%Nslip)
      prm%theta0_b = math_expand(prm%theta0_b,prm%Nslip)
      prm%theta1_b = math_expand(prm%theta1_b,prm%Nslip)



!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (    prm%gdot0  <= 0.0_pReal)   extmsg = trim(extmsg)//' gdot0'
      if (    prm%n      <= 0.0_pReal)   extmsg = trim(extmsg)//' n_slip'
      if (any(prm%crss0  <= 0.0_pReal))  extmsg = trim(extmsg)//' crss0'
      if (any(prm%tau1   <= 0.0_pReal))  extmsg = trim(extmsg)//' tau1'
      if (any(prm%tau1_b <= 0.0_pReal))  extmsg = trim(extmsg)//' tau1_b'

      !ToDo: Any sensible checks for theta?

    endif slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
    sizeDotState   = size(['crss     ','crss_back', 'accshear ']) * prm%sum_N_sl
    sizeDeltaState = size(['sense ',   'chi0  ',    'gamma0'   ]) * prm%sum_N_sl
    sizeState = sizeDotState + sizeDeltaState

    call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,sizeDeltaState)

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and atol
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%crss => plasticState(p)%state   (startIndex:endIndex,:)
    stt%crss = spread(prm%crss0, 2, NipcMyPhase)
    dot%crss => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_resistance',defaultVal=1.0_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) <= 0.0_pReal)) &
      extmsg = trim(extmsg)//' atol_crss'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%crss_back => plasticState(p)%state   (startIndex:endIndex,:)
    dot%crss_back => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_resistance',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%accshear => plasticState(p)%state   (startIndex:endIndex,:)
    dot%accshear => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%atol(startIndex:endIndex) = config%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
    if(any(plasticState(p)%atol(startIndex:endIndex) <= 0.0_pReal)) &
      extmsg = trim(extmsg)//' atol_gamma'
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
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_KINEHARDENING_LABEL//')')

  enddo


end subroutine plastic_kinehardening_init


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
                 * (  prm%theta1 &
                     + (prm%theta0 - prm%theta1 + prm%theta0*prm%theta1*sumGamma/prm%tau1) &
                     * exp(-sumGamma*prm%theta0/prm%tau1) &
                   )

  dot%crss_back(:,of) = stt%sense(:,of)*dot%accshear(:,of) * &
           ( prm%theta1_b + &
             (prm%theta0_b - prm%theta1_b &
               + prm%theta0_b*prm%theta1_b/(prm%tau1_b+stt%chi0(:,of))*(stt%accshear(:,of)-stt%gamma0(:,of))&
             ) *exp(-(stt%accshear(:,of)-stt%gamma0(:,of)) *prm%theta0_b/(prm%tau1_b+stt%chi0(:,of))) &
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
  if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0 &
             .and. (of == prm%of_debug &
                    .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0)) then
    write(6,'(a)') '======= kinehardening delta state ======='
    write(6,*) sense,state(instance)%sense(:,of)
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
     case('resistance')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%crss,'xi_sl', &
                                                    'resistance against plastic slip','Pa')
     case('backstress')                                                                             ! ToDo: should be 'tau_back'
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%crss_back,'tau_back', &
                                                    'back stress against plastic slip','Pa')
     case ('sense')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%sense,'sense_of_shear', &
                                                    'tbd','1')
     case ('chi0')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%chi0,'chi0', &
                                                    'tbd','Pa')
     case ('gamma0')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma0,'gamma0', &
                                                    'tbd','1')
     case ('accumulatedshear')
       if(prm%sum_N_sl>0) call results_writeDataset(group,stt%accshear,'gamma_sl', &
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
  logical :: nonSchmidActive

  associate(prm => param(instance), stt => state(instance))

  nonSchmidActive = size(prm%nonSchmidCoeff) > 0

  do i = 1, prm%sum_N_sl
    tau_pos(i) =       math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i)) - stt%crss_back(i,of)
    tau_neg(i) = merge(math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,i)) - stt%crss_back(i,of), &
                       0.0_pReal, nonSchmidActive)
  enddo

  where(dNeq0(tau_pos))
    gdot_pos = prm%gdot0 * merge(0.5_pReal,1.0_pReal, nonSchmidActive) &                            ! 1/2 if non-Schmid active
             * sign(abs(tau_pos/stt%crss(:,of))**prm%n,  tau_pos)
  else where
    gdot_pos = 0.0_pReal
  end where

  where(dNeq0(tau_neg))
    gdot_neg = prm%gdot0 * 0.5_pReal &                                                              ! only used if non-Schmid active, always 1/2
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
