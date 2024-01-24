!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) phenopowerlaw

  type :: tParameters
    real(pREAL),               allocatable, dimension(:) :: &
      dot_gamma_0_sl, &                                                                             !< reference shear strain rate for slip
      dot_gamma_0_tw, &                                                                             !< reference shear strain rate for twin
      a_sl, &
      n_sl, &                                                                                       !< stress exponent for slip
      n_tw, &                                                                                       !< stress exponent for twin
      xi_inf_sl, &                                                                                  !< maximum critical shear stress for slip
      f_sat_sl_tw, &                                                                                !< push-up factor for slip saturation due to twinning
      c_1, &
      c_2, &
      c_3, &
      c_4, &
      h_0_sl_sl, &                                                                                  !< reference hardening slip - slip
      h_0_tw_sl, &                                                                                  !< reference hardening twin - slip
      h_0_tw_tw, &                                                                                  !< reference hardening twin - twin
      gamma_char                                                                                    !< characteristic shear for twins
    real(pREAL),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      h_sl_tw, &                                                                                    !< slip resistance from twin activity
      h_tw_sl, &                                                                                    !< twin resistance from slip activity
      h_tw_tw                                                                                       !< twin resistance from twin activity
    real(pREAL),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw                                                                                      !< total number of active twin systems
    character(len=pSTRLEN),    allocatable, dimension(:) :: &
      output
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl, &
      systems_tw
  end type tParameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      xi_sl, &
      xi_tw, &
      gamma_sl, &
      gamma_tw
  end type tIndexDotState

  type :: tPhenopowerlawState
    real(pREAL), pointer, dimension(:,:) :: &
      xi_sl, &
      xi_tw, &
      gamma_sl, &
      gamma_tw
  end type tPhenopowerlawState

!--------------------------------------------------------------------------------------------------
! containers for parameters, dot state index,  and state
  type(tParameters),         allocatable, dimension(:) :: param
  type(tIndexDotState),      allocatable, dimension(:) :: indexDotState
  type(tPhenopowerlawState), allocatable, dimension(:) :: state

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_phenopowerlaw_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl, &                                                                                         !< number of slip-systems for a given slip family
    N_tw                                                                                            !< number of twin-systems for a given twin family
  real(pREAL), dimension(:), allocatable :: &
    xi_0_sl, &                                                                                      !< initial critical shear stress for slip
    xi_0_tw                                                                                          !< initial critical shear stress for twin
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


  myPlasticity = plastic_active('phenopowerlaw')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:phenopowerlaw init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), stt => state(ph), &
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


    N_sl = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    N_tw = pl%get_as1dInt('N_tw',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    prm%sum_N_tw = sum(abs(N_tw))

!--------------------------------------------------------------------------------------------------
! slip related parameters
    slipActive: if (prm%sum_N_sl > 0) then
      prm%dot_gamma_0_sl = math_expand(pl%get_as1dReal('dot_gamma_0_sl',requiredSize=size(N_sl)), N_sl)
      prm%n_sl           = math_expand(pl%get_as1dReal('n_sl',          requiredSize=size(N_sl)), N_sl)
      prm%a_sl           = math_expand(pl%get_as1dReal('a_sl',          requiredSize=size(N_sl)), N_sl)
      prm%h_0_sl_sl      = math_expand(pl%get_as1dReal('h_0_sl-sl',     requiredSize=size(N_sl)), N_sl)
      xi_0_sl            = math_expand(pl%get_as1dReal('xi_0_sl',       requiredSize=size(N_sl)), N_sl)
      prm%xi_inf_sl      = math_expand(pl%get_as1dReal('xi_inf_sl',     requiredSize=size(N_sl)), N_sl)
      prm%c_1            = math_expand(pl%get_as1dReal('c_1',           requiredSize=size(N_sl), &
                                                                        defaultVal=misc_zeros(size(N_sl))), N_sl)
      prm%c_2            = math_expand(pl%get_as1dReal('c_2',           requiredSize=size(N_sl), &
                                                                        defaultVal=misc_ones(size(N_sl))), N_sl)
      prm%f_sat_sl_tw    = math_expand(pl%get_as1dReal('f_sat_sl-tw',   requiredSize=size(N_sl), &
                                                                        defaultVal=misc_zeros(size(N_sl))), N_sl)

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'),phase_lattice(ph))

      prm%P_sl = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        allocate(a_nS(3,size(pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray))),source=0.0_pREAL)
        a_nS(1,:) = pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray)
        prm%P_nS_pos = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=+1)
        prm%P_nS_neg = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=-1)
        deallocate(a_nS)
      else
        prm%P_nS_pos = +prm%P_sl
        prm%P_nS_neg = -prm%P_sl
      end if

      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))

      ! sanity checks
      if (any(prm%dot_gamma_0_sl <= 0.0_pREAL))   extmsg = trim(extmsg)//' dot_gamma_0_sl'
      if (any(prm%n_sl           <= 0.0_pREAL))   extmsg = trim(extmsg)//' n_sl'
      if (any(prm%a_sl           <= 0.0_pREAL))   extmsg = trim(extmsg)//' a_sl'
      if (any(xi_0_sl            <= 0.0_pREAL))   extmsg = trim(extmsg)//' xi_0_sl'
      if (any(prm%xi_inf_sl      <= 0.0_pREAL))   extmsg = trim(extmsg)//' xi_inf_sl'

    else slipActive
      xi_0_sl = emptyRealArray
      allocate(prm%dot_gamma_0_sl, &
               prm%a_sl, &
               prm%n_sl, &
               prm%xi_inf_sl, &
               prm%f_sat_sl_tw, &
               prm%c_1, &
               prm%c_2, &
               prm%h_0_sl_sl, &
               source=emptyRealArray)
      allocate(prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    twinActive: if (prm%sum_N_tw > 0) then
      prm%dot_gamma_0_tw = math_expand(pl%get_as1dReal('dot_gamma_0_tw', requiredSize=size(N_tw)), N_tw)
      prm%n_tw           = math_expand(pl%get_as1dReal('n_tw',           requiredSize=size(N_tw)), N_tw)
      prm%h_0_tw_tw      = math_expand(pl%get_as1dReal('h_0_tw-tw',      requiredSize=size(N_tw)), N_tw)
      xi_0_tw            = math_expand(pl%get_as1dReal('xi_0_tw',        requiredSize=size(N_tw)), N_tw)
      prm%c_3            = math_expand(pl%get_as1dReal('c_3',            requiredSize=size(N_tw), &
                                                                         defaultVal=misc_ones(size(N_tw))), N_tw)
      prm%c_4            = math_expand(pl%get_as1dReal('c_4',            requiredSize=size(N_tw), &
                                                                         defaultVal=misc_zeros(size(N_tw))), N_tw)

      prm%gamma_char = crystal_characteristicShear_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%h_tw_tw    = crystal_interaction_TwinByTwin(N_tw,pl%get_as1dReal('h_tw-tw'),phase_lattice(ph))

      prm%P_tw       = crystal_SchmidMatrix_twin(N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%systems_tw = crystal_labels_twin(N_tw,phase_lattice(ph))

      ! sanity checks
      if (any(prm%dot_gamma_0_tw <= 0.0_pREAL))   extmsg = trim(extmsg)//' dot_gamma_0_tw'
      if (any(prm%n_tw           <= 0.0_pREAL))   extmsg = trim(extmsg)//' n_tw'
      if (any(xi_0_tw            <= 0.0_pREAL))   extmsg = trim(extmsg)//' xi_0_tw'

    else twinActive
      xi_0_tw = emptyRealArray
      allocate(prm%dot_gamma_0_tw, &
               prm%n_tw, &
               prm%c_3, &
               prm%c_4, &
               prm%gamma_char, &
               prm%h_0_tw_sl, &
               prm%h_0_tw_tw, &
               source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
    end if twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
    slipAndTwinActive: if (prm%sum_N_sl > 0 .and. prm%sum_N_tw > 0) then
      prm%h_0_tw_sl = math_expand(pl%get_as1dReal('h_0_tw-sl',requiredSize=size(N_tw)), N_tw)
      prm%h_sl_tw    = crystal_interaction_SlipByTwin(N_sl,N_tw,pl%get_as1dReal('h_sl-tw'),phase_lattice(ph))
      prm%h_tw_sl    = crystal_interaction_TwinBySlip(N_tw,N_sl,pl%get_as1dReal('h_tw-sl'),phase_lattice(ph))
    else slipAndTwinActive
      allocate(prm%h_sl_tw(prm%sum_N_sl,prm%sum_N_tw))                                              ! at least one dimension is 0
      allocate(prm%h_tw_sl(prm%sum_N_tw,prm%sum_N_sl))                                              ! at least one dimension is 0
      prm%h_0_tw_sl = [(0.0_pREAL,i=1,size(N_tw))]
    end if slipAndTwinActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_ID_phase == ph)
    sizeDotState = size(['xi_sl   ','gamma_sl']) * prm%sum_N_sl &
                 + size(['xi_tw   ','gamma_tw']) * prm%sum_N_tw
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%xi_sl = [startIndex,endIndex]
    stt%xi_sl => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_sl = spread(xi_0_sl, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_xi',defaultVal=1.0_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_xi'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%xi_tw = [startIndex,endIndex]
    stt%xi_tw => plasticState(ph)%state(startIndex:endIndex,:)
    stt%xi_tw = spread(xi_0_tw, 2, Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_xi',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%gamma_tw = [startIndex,endIndex]
    stt%gamma_tw => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!> @details assumes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
pure module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pREAL), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress

  real(pREAL), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,               intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_sl
  real(pREAL), dimension(3,3,param(ph)%sum_N_sl) :: &
    P_nS
  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw,ddot_gamma_dtau_tw


  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph))

    call kinetics_sl(Mp,ph,en,dot_gamma_sl,ddot_gamma_dtau_sl)
    P_nS = merge(prm%P_nS_pos,prm%P_nS_neg, spread(spread(dot_gamma_sl,1,3),2,3)>0.0_pREAL)         ! faster than 'merge' in loop
    slipSystems: do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_sl(i) * prm%P_sl(k,l,i) * P_nS(m,n,i)
    end do slipSystems

    call kinetics_tw(Mp,ph,en,dot_gamma_tw,ddot_gamma_dtau_tw)
    twinSystems: do i = 1, prm%sum_N_tw
      Lp = Lp + dot_gamma_tw(i)*prm%P_tw(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_tw(i)*prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
    end do twinSystems

  end associate

end subroutine phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function phenopowerlaw_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  real(pREAL) :: &
    sumF
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    xi_sl_sat_offset, &
    left_SlipSlip


  associate(prm => param(ph), stt => state(ph), &
            dot_xi_sl => dotState(indexDotState(ph)%xi_sl(1):indexDotState(ph)%xi_sl(2)), &
            dot_xi_tw => dotState(indexDotState(ph)%xi_tw(1):indexDotState(ph)%xi_tw(2)), &
            dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)), &
            dot_gamma_tw => dotState(indexDotState(ph)%gamma_tw(1):indexDotState(ph)%gamma_tw(2)))

    call kinetics_sl(Mp,ph,en, dot_gamma_sl)
    call kinetics_tw(Mp,ph,en, dot_gamma_tw)
    dot_gamma_sl = abs(dot_gamma_sl)
    sumF = sum(stt%gamma_tw(:,en)/prm%gamma_char)

    xi_sl_sat_offset = prm%f_sat_sl_tw*sqrt(sumF)

    left_SlipSlip = sign(abs(1.0_pREAL - stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))**prm%a_sl, &
                             1.0_pREAL - stt%xi_sl(:,en) / (prm%xi_inf_sl+xi_sl_sat_offset))

    dot_xi_sl = prm%h_0_sl_sl * (1.0_pREAL + prm%c_1 * sumF**prm%c_2) &
              * left_SlipSlip &
              * matmul(prm%h_sl_sl,dot_gamma_sl) &
              + matmul(prm%h_sl_tw,dot_gamma_tw)

    dot_xi_tw = prm%h_0_tw_sl * sum(stt%gamma_sl(:,en))**prm%c_3 * matmul(prm%h_tw_sl,dot_gamma_sl) &
              + prm%h_0_tw_tw * sumF                   **prm%c_4 * matmul(prm%h_tw_tw,dot_gamma_tw)

  end associate

end function phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_phenopowerlaw_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('xi_sl')
          call result_writeDataset(stt%xi_sl,group,trim(prm%output(ou)), &
                                   'resistance against plastic slip','Pa',prm%systems_sl)
        case('gamma_sl')
          call result_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)

        case('xi_tw')
          call result_writeDataset(stt%xi_tw,group,trim(prm%output(ou)), &
                                   'resistance against twinning','Pa',prm%systems_tw)
        case('gamma_tw')
          call result_writeDataset(stt%gamma_tw,group,trim(prm%output(ou)), &
                                   'twinning shear','1',prm%systems_tw)

      end select

    end do

  end associate

end subroutine plastic_phenopowerlaw_result


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems and their derivatives with respect to resolved
!         stress.
!> @details Sign of dot_gamma_sl conveys sense of shear.
! Derivatives are calculated only optionally, hence, contrary to common convention,
! here the result (i.e. intent(out)) variables have to be put at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_sl(Mp,ph,en, &
                            dot_gamma_sl,ddot_gamma_dtau_sl)

  real(pREAL), dimension(3,3),                           intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                                               intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl),           intent(out) :: &
    dot_gamma_sl
  real(pREAL), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau_sl

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau_sl_pos, &
    tau_sl_neg
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    tau_sl_pos = [(math_tensordot(Mp,prm%P_nS_pos(1:3,1:3,i)),i=1,prm%sum_N_sl)]
    tau_sl_neg = [(math_tensordot(Mp,prm%P_nS_neg(1:3,1:3,i)),i=1,prm%sum_N_sl)]

    dot_gamma_sl = merge(+1.0_pREAL,-1.0_pREAL, tau_sl_pos>tau_sl_neg) &
                 * prm%dot_gamma_0_sl  &
                 * (max(tau_sl_pos,tau_sl_neg)/stt%xi_sl(:,en))**prm%n_sl

    if (present(ddot_gamma_dtau_sl)) then
      where(dNeq0(dot_gamma_sl))
        ddot_gamma_dtau_sl = dot_gamma_sl*prm%n_sl/max(tau_sl_pos,tau_sl_neg)
      else where
        ddot_gamma_dtau_sl = 0.0_pREAL
      end where
    end if

  end associate

end subroutine kinetics_sl


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved stress.
!         Twinning is assumed to take place only in an untwinned volume.
!> @details Derivatives are calculated and returned if corresponding output variables are present in the argument list.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_tw(Mp,ph,en,&
                            dot_gamma_tw,ddot_gamma_dtau_tw)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_tw), intent(out) :: &
    dot_gamma_tw
  real(pREAL), dimension(param(ph)%sum_N_tw), intent(out), optional :: &
    ddot_gamma_dtau_tw

  real(pREAL), dimension(param(ph)%sum_N_tw) :: &
    tau_tw
  integer :: i


  associate(prm => param(ph), stt => state(ph))

    tau_tw = [(math_tensordot(Mp,prm%P_tw(1:3,1:3,i)),i=1,prm%sum_N_tw)]

    where(tau_tw > 0.0_pREAL)
      dot_gamma_tw = (1.0_pREAL-sum(stt%gamma_tw(:,en)/prm%gamma_char)) &                           ! only twin in untwinned volume fraction
                   * prm%dot_gamma_0_tw*(tau_tw/stt%xi_tw(:,en))**prm%n_tw
    else where
      dot_gamma_tw = 0.0_pREAL
    end where

    if (present(ddot_gamma_dtau_tw)) then
      where(dNeq0(dot_gamma_tw))
        ddot_gamma_dtau_tw = dot_gamma_tw*prm%n_tw/tau_tw
      else where
        ddot_gamma_dtau_tw = 0.0_pREAL
      end where
    end if

  end associate

end subroutine kinetics_tw

end submodule phenopowerlaw
