!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Su Leen Wong, Max-Planck-Institut für Eisenforschung GmbH
!> @author Nan Jia, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) dislotwin

  type :: tParameters
    real(pReal) :: &
      Q_cl                = 1.0_pReal, &                                                            !< activation energy for dislocation climb
      omega               = 1.0_pReal, &                                                            !< frequency factor for dislocation climb
      D                   = 1.0_pReal, &                                                            !< grain size
      p_sb                = 1.0_pReal, &                                                            !< p-exponent in shear band velocity
      q_sb                = 1.0_pReal, &                                                            !< q-exponent in shear band velocity
      i_tw                = 1.0_pReal, &                                                            !< adjustment parameter to calculate MFP for twinning
      i_tr                = 1.0_pReal, &                                                            !< adjustment parameter to calculate MFP for transformation
      L_tw                = 1.0_pReal, &                                                            !< length of twin nuclei
      L_tr                = 1.0_pReal, &                                                            !< length of trans nuclei
      x_c                 = 1.0_pReal, &                                                            !< critical distance for formation of twin/trans nucleus
      V_cs                = 1.0_pReal, &                                                            !< cross slip volume
      xi_sb               = 1.0_pReal, &                                                            !< value for shearband resistance
      v_sb                = 1.0_pReal, &                                                            !< value for shearband velocity_0
      E_sb                = 1.0_pReal, &                                                            !< activation energy for shear bands
      h                   = 1.0_pReal, &                                                            !< stack height of hex nucleus
      gamma_char_tr       = sqrt(0.125_pReal), &                                                    !< Characteristic shear for transformation
      a_cF                = 1.0_pReal, &
      cOverA_hP           = 1.0_pReal, &
      V_mol               = 1.0_pReal, &
      rho                 = 1.0_pReal
    type(tPolynomial) :: &
      Gamma_sf, &                                                                                   !< stacking fault energy
      Delta_G                                                                                       !< free energy difference between austensite and martensite
    real(pReal),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< absolute length of Burgers vector [m] for each slip system
      b_tw, &                                                                                       !< absolute length of Burgers vector [m] for each twin system
      b_tr, &                                                                                       !< absolute length of Burgers vector [m] for each transformation system
      Q_sl,&                                                                                        !< activation energy for glide [J] for each slip system
      v_0, &                                                                                        !< dislocation velocity prefactor [m/s] for each slip system
      dot_N_0_tw, &                                                                                 !< twin nucleation rate [1/m³s] for each twin system
      t_tw, &                                                                                       !< twin thickness [m] for each twin system
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations for each slip system
      t_tr, &                                                                                       !< martensite lamellar thickness [m] for each trans system
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      r, &                                                                                          !< exponent in twin nucleation rate
      s, &                                                                                          !< exponent in trans nucleation rate
      tau_0, &                                                                                      !< strength due to elements in solid solution
      gamma_char_tw, &                                                                              !< characteristic shear for twins
      B, &                                                                                          !< drag coefficient
      d_caron                                                                                       !< distance of spontaneous annhihilation
    real(pReal),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< components of slip-slip interaction matrix
      h_sl_tw, &                                                                                    !< components of slip-twin interaction matrix
      h_tw_tw, &                                                                                    !< components of twin-twin interaction matrix
      h_sl_tr, &                                                                                    !< components of slip-trans interaction matrix
      h_tr_tr, &                                                                                    !< components of trans-trans interaction matrix
      n0_sl, &                                                                                      !< slip system normal
      forestProjection
    real(pReal),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_tw, &
      P_tr
    integer :: &
      sum_N_sl, &                                                                                   !< total number of active slip system
      sum_N_tw, &                                                                                   !< total number of active twin system
      sum_N_tr                                                                                      !< total number of active transformation system
    integer,                   allocatable, dimension(:)   :: &
      N_tw, &
      N_tr
    integer,                   allocatable, dimension(:,:) :: &
      fcc_twinNucleationSlipPair                                                                    ! ToDo: Better name? Is also use for trans
    character(len=:),          allocatable                 :: &
      lattice_tr
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
    logical :: &
      extendedDislocations, &                                                                       !< consider split into partials for climb calculation
      fccTwinTransNucleation, &                                                                     !< twinning and transformation models are for fcc
      omitDipoles                                                                                   !< flag controlling consideration of dipole formation
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl, &
      systems_tw
  end type tParameters                                                                              !< container type for internal constitutive parameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      rho_mob, &
      rho_dip, &
      gamma_sl, &
      f_tw, &
      f_tr
  end type tIndexDotState

  type :: tDislotwinState
    real(pReal),                  dimension(:,:),   pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl, &
      f_tw, &
      f_tr
  end type tDislotwinState

  type :: tDislotwinDependentState
    real(pReal),                  dimension(:,:),   allocatable :: &
      Lambda_sl, &                                                                                  !< mean free path between 2 obstacles seen by a moving dislocation
      Lambda_tw, &                                                                                  !< mean free path between 2 obstacles seen by a growing twin
      Lambda_tr, &                                                                                  !< mean free path between 2 obstacles seen by a growing martensite
      tau_pass                                                                                      !< threshold stress for slip
  end type tDislotwinDependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),              allocatable, dimension(:) :: param
  type(tIndexDotState),           allocatable, dimension(:) :: indexDotState
  type(tDislotwinState),          allocatable, dimension(:) :: state
  type(tDislotwinDependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislotwin_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,     dimension(:), allocatable :: &
    N_sl
  real(pReal) :: a_cF
  real(pReal), allocatable, dimension(:) :: &
    rho_mob_0, &                                                                                    !< initial unipolar dislocation density per slip system
    rho_dip_0                                                                                       !< initial dipole dislocation density per slip system
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislotwin')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:dislotwin init  -+>>>'
  print'(/,a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)

  print'(/,1x,a)', 'A. Ma and F. Roters, Acta Materialia 52(12):3603–3612, 2004'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2004.04.012'

  print'(/,1x,a)', 'F. Roters et al., Computational Materials Science 39:91–95, 2007'
  print'(  1x,a)', 'https://doi.org/10.1016/j.commatsci.2006.04.014'

  print'(/,1x,a)', 'S.L. Wong et al., Acta Materialia 118:140–151, 2016'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2016.07.032'


  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(indexDotState(phases%length))
  allocate(state(phases%length))
  allocate(dependentState(phases%length))

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), stt => state(ph), dst => dependentState(ph), &
              idx_dot => indexDotState(ph))

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
      prm%systems_sl = lattice_labels_slip(N_sl,phase_lattice(ph))
      prm%P_sl    = lattice_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%h_sl_sl = lattice_interaction_SlipBySlip(N_sl,pl%get_as1dFloat('h_sl-sl'),phase_lattice(ph))
      prm%forestProjection = lattice_forestProjection_edge(N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%forestProjection = transpose(prm%forestProjection)

      prm%n0_sl = lattice_slip_normal(N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%fccTwinTransNucleation = phase_lattice(ph) == 'cF' .and. (N_sl(1) == 12)
      if (prm%fccTwinTransNucleation) prm%fcc_twinNucleationSlipPair = lattice_CF_TWINNUCLEATIONSLIPPAIR

      rho_mob_0                = pl%get_as1dFloat('rho_mob_0',   requiredSize=size(N_sl))
      rho_dip_0                = pl%get_as1dFloat('rho_dip_0',   requiredSize=size(N_sl))
      prm%v_0                  = pl%get_as1dFloat('v_0',         requiredSize=size(N_sl))
      prm%b_sl                 = pl%get_as1dFloat('b_sl',        requiredSize=size(N_sl))
      prm%Q_sl                 = pl%get_as1dFloat('Q_sl',        requiredSize=size(N_sl))
      prm%i_sl                 = pl%get_as1dFloat('i_sl',        requiredSize=size(N_sl))
      prm%p                    = pl%get_as1dFloat('p_sl',        requiredSize=size(N_sl))
      prm%q                    = pl%get_as1dFloat('q_sl',        requiredSize=size(N_sl))
      prm%tau_0                = pl%get_as1dFloat('tau_0',       requiredSize=size(N_sl))
      prm%B                    = pl%get_as1dFloat('B',           requiredSize=size(N_sl), &
                                                  defaultVal=[(0.0_pReal, i=1,size(N_sl))])

      prm%Q_cl = pl%get_asFloat('Q_cl')

      prm%extendedDislocations = pl%get_asBool('extend_dislocations',defaultVal = .false.)
      prm%omitDipoles          = pl%get_asBool('omit_dipoles',defaultVal = .false.)

      ! multiplication factor according to crystal structure (nearest neighbors bcc vs fcc/hex)
      ! details: Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
      prm%omega = pl%get_asFloat('omega',  defaultVal = 1000.0_pReal) &
                * merge(12.0_pReal,8.0_pReal,any(phase_lattice(ph) == ['cF','hP']))

      ! expand: family => system
      rho_mob_0        = math_expand(rho_mob_0,       N_sl)
      rho_dip_0        = math_expand(rho_dip_0,       N_sl)
      prm%v_0          = math_expand(prm%v_0,         N_sl)
      prm%b_sl         = math_expand(prm%b_sl,        N_sl)
      prm%Q_sl         = math_expand(prm%Q_sl,        N_sl)
      prm%i_sl         = math_expand(prm%i_sl,        N_sl)
      prm%p            = math_expand(prm%p,           N_sl)
      prm%q            = math_expand(prm%q,           N_sl)
      prm%tau_0        = math_expand(prm%tau_0,       N_sl)
      prm%B            = math_expand(prm%B,           N_sl)
      prm%d_caron      = pl%get_asFloat('D_a') * prm%b_sl

      ! sanity checks
      if (    prm%Q_cl          <= 0.0_pReal)          extmsg = trim(extmsg)//' Q_cl'
      if (any(rho_mob_0         <  0.0_pReal))         extmsg = trim(extmsg)//' rho_mob_0'
      if (any(rho_dip_0         <  0.0_pReal))         extmsg = trim(extmsg)//' rho_dip_0'
      if (any(prm%v_0           <  0.0_pReal))         extmsg = trim(extmsg)//' v_0'
      if (any(prm%b_sl          <= 0.0_pReal))         extmsg = trim(extmsg)//' b_sl'
      if (any(prm%Q_sl           <= 0.0_pReal))        extmsg = trim(extmsg)//' Q_sl'
      if (any(prm%i_sl          <= 0.0_pReal))         extmsg = trim(extmsg)//' i_sl'
      if (any(prm%B             <  0.0_pReal))         extmsg = trim(extmsg)//' B'
      if (any(prm%d_caron       <  0.0_pReal))         extmsg = trim(extmsg)//' d_caron(D_a,b_sl)'
      if (any(prm%p<=0.0_pReal .or. prm%p>1.0_pReal))  extmsg = trim(extmsg)//' p_sl'
      if (any(prm%q< 1.0_pReal .or. prm%q>2.0_pReal))  extmsg = trim(extmsg)//' q_sl'
    else slipActive
      rho_mob_0 = emptyRealArray; rho_dip_0 = emptyRealArray
      allocate(prm%b_sl,prm%Q_sl,prm%v_0,prm%i_sl,prm%p,prm%q,prm%B,source=emptyRealArray)
      allocate(prm%forestProjection(0,0),prm%h_sl_sl(0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
    prm%N_tw = pl%get_as1dInt('N_tw', defaultVal=emptyIntArray)
    prm%sum_N_tw = sum(abs(prm%N_tw))
    twinActive: if (prm%sum_N_tw > 0) then
      prm%systems_tw = lattice_labels_twin(prm%N_tw,phase_lattice(ph))
      prm%P_tw  = lattice_SchmidMatrix_twin(prm%N_tw,phase_lattice(ph),phase_cOverA(ph))
      prm%h_tw_tw   = lattice_interaction_TwinByTwin(prm%N_tw,pl%get_as1dFloat('h_tw-tw'), &
                                                     phase_lattice(ph))

      prm%b_tw      = pl%get_as1dFloat('b_tw',     requiredSize=size(prm%N_tw))
      prm%t_tw      = pl%get_as1dFloat('t_tw',     requiredSize=size(prm%N_tw))
      prm%r         = pl%get_as1dFloat('p_tw',     requiredSize=size(prm%N_tw))

      prm%L_tw      = pl%get_asFloat('L_tw')
      prm%i_tw      = pl%get_asFloat('i_tw')

      prm%gamma_char_tw = lattice_characteristicShear_Twin(prm%N_tw,phase_lattice(ph),phase_cOverA(ph))

      ! expand: family => system
      prm%b_tw = math_expand(prm%b_tw,prm%N_tw)
      prm%t_tw = math_expand(prm%t_tw,prm%N_tw)
      prm%r    = math_expand(prm%r,prm%N_tw)

      ! sanity checks
      if (.not. prm%fccTwinTransNucleation)   extmsg = trim(extmsg)//' TWIP for non-fcc'
      if (    prm%L_tw          < 0.0_pReal)  extmsg = trim(extmsg)//' L_tw'
      if (    prm%i_tw          < 0.0_pReal)  extmsg = trim(extmsg)//' i_tw'
      if (any(prm%b_tw          < 0.0_pReal)) extmsg = trim(extmsg)//' b_tw'
      if (any(prm%t_tw          < 0.0_pReal)) extmsg = trim(extmsg)//' t_tw'
      if (any(prm%r             < 0.0_pReal)) extmsg = trim(extmsg)//' p_tw'
    else twinActive
      allocate(prm%gamma_char_tw,prm%b_tw,prm%t_tw,prm%r,source=emptyRealArray)
      allocate(prm%h_tw_tw(0,0))
    end if twinActive

!--------------------------------------------------------------------------------------------------
! transformation related parameters
    prm%N_tr = pl%get_as1dInt('N_tr', defaultVal=emptyIntArray)
    prm%sum_N_tr = sum(abs(prm%N_tr))
    transActive: if (prm%sum_N_tr > 0) then
      prm%b_tr = pl%get_as1dFloat('b_tr')
      prm%b_tr = math_expand(prm%b_tr,prm%N_tr)

      prm%i_tr       = pl%get_asFloat('i_tr')
      prm%Delta_G    = polynomial(pl%asDict(),'Delta_G','T')
      prm%L_tr       = pl%get_asFloat('L_tr')
      a_cF           = pl%get_asFloat('a_cF')
      prm%h          = 5.0_pReal * a_cF/sqrt(3.0_pReal)
      prm%cOverA_hP  = pl%get_asFloat('c/a_hP')
      prm%rho        = 4.0_pReal/(sqrt(3.0_pReal)*a_cF**2)/N_A
      prm%V_mol      = pl%get_asFloat('V_mol')
      prm%h_tr_tr = lattice_interaction_TransByTrans(prm%N_tr,pl%get_as1dFloat('h_tr-tr'),&
                                                     phase_lattice(ph))

      prm%P_tr    = lattice_SchmidMatrix_trans(prm%N_tr,'hP',prm%cOverA_hP)

      prm%t_tr = pl%get_as1dFloat('t_tr')
      prm%t_tr = math_expand(prm%t_tr,prm%N_tr)
      prm%s    = pl%get_as1dFloat('p_tr')
      prm%s    = math_expand(prm%s,prm%N_tr)

      ! sanity checks
      if (.not. prm%fccTwinTransNucleation)   extmsg = trim(extmsg)//' TRIP for non-fcc'
      if (    prm%L_tr          < 0.0_pReal)  extmsg = trim(extmsg)//' L_tr'
      if (    prm%V_mol         < 0.0_pReal)  extmsg = trim(extmsg)//' V_mol'
      if (    prm%i_tr          < 0.0_pReal)  extmsg = trim(extmsg)//' i_tr'
      if (any(prm%t_tr          < 0.0_pReal)) extmsg = trim(extmsg)//' t_tr'
      if (any(prm%s             < 0.0_pReal)) extmsg = trim(extmsg)//' p_tr'
    else transActive
      allocate(prm%s,prm%b_tr,prm%t_tr,source=emptyRealArray)
      allocate(prm%h_tr_tr(0,0))
    end if transActive

!--------------------------------------------------------------------------------------------------
! shearband related parameters
    prm%v_sb = pl%get_asFloat('v_sb',defaultVal=0.0_pReal)
    if (prm%v_sb > 0.0_pReal) then
      prm%xi_sb        = pl%get_asFloat('xi_sb')
      prm%E_sb         = pl%get_asFloat('Q_sb')
      prm%p_sb         = pl%get_asFloat('p_sb')
      prm%q_sb         = pl%get_asFloat('q_sb')

      ! sanity checks
      if (prm%xi_sb         <  0.0_pReal) extmsg = trim(extmsg)//' xi_sb'
      if (prm%E_sb          <  0.0_pReal) extmsg = trim(extmsg)//' Q_sb'
      if (prm%p_sb          <= 0.0_pReal) extmsg = trim(extmsg)//' p_sb'
      if (prm%q_sb          <= 0.0_pReal) extmsg = trim(extmsg)//' q_sb'
    end if

!--------------------------------------------------------------------------------------------------
! parameters required for several mechanisms and their interactions
    if (prm%sum_N_sl + prm%sum_N_tw + prm%sum_N_tw > 0) &
      prm%D = pl%get_asFloat('D')

    if (prm%sum_N_tw + prm%sum_N_tr > 0) then
      prm%x_c  = pl%get_asFloat('x_c')
      prm%V_cs = pl%get_asFloat('V_cs')
      if (prm%x_c  < 0.0_pReal)  extmsg = trim(extmsg)//' x_c'
      if (prm%V_cs < 0.0_pReal)  extmsg = trim(extmsg)//' V_cs'
    end if

    if (prm%sum_N_tw + prm%sum_N_tr > 0 .or. prm%extendedDislocations) &
      prm%Gamma_sf = polynomial(pl%asDict(),'Gamma_sf','T')

    slipAndTwinActive: if (prm%sum_N_sl * prm%sum_N_tw > 0) then
      prm%h_sl_tw = lattice_interaction_SlipByTwin(N_sl,prm%N_tw,pl%get_as1dFloat('h_sl-tw'), &
                                                   phase_lattice(ph))
      if (prm%fccTwinTransNucleation .and. size(prm%N_tw) /= 1) extmsg = trim(extmsg)//' N_tw: nucleation'
    endif slipAndTwinActive

    slipAndTransActive: if (prm%sum_N_sl * prm%sum_N_tr > 0) then
      prm%h_sl_tr = lattice_interaction_SlipByTrans(N_sl,prm%N_tr,pl%get_as1dFloat('h_sl-tr'), &
                                                    phase_lattice(ph))
      if (prm%fccTwinTransNucleation .and. size(prm%N_tr) /= 1) extmsg = trim(extmsg)//' N_tr: nucleation'
    endif slipAndTransActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers  = count(material_phaseID == ph)
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl &
                 + size(['f_tw'])                           * prm%sum_N_tw &
                 + size(['f_tr'])                           * prm%sum_N_tr
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%rho_mob = [startIndex,endIndex]
    stt%rho_mob=>plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_mob= spread(rho_mob_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_rho',defaultVal=1.0_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%rho_dip = [startIndex,endIndex]
    stt%rho_dip=>plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_dip= spread(rho_dip_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_rho',defaultVal=1.0_pReal)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl=>plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_gamma',defaultVal=1.0e-6_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_gamma'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tw
    idx_dot%f_tw = [startIndex,endIndex]
    stt%f_tw=>plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_f_tw',defaultVal=1.0e-6_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_f_tw'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_tr
    idx_dot%f_tr = [startIndex,endIndex]
    stt%f_tr=>plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_f_tr',defaultVal=1.0e-6_pReal)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_f_tr'

    allocate(dst%tau_pass (prm%sum_N_sl,Nmembers),source=0.0_pReal)
    allocate(dst%Lambda_sl(prm%sum_N_sl,Nmembers),source=0.0_pReal)
    allocate(dst%Lambda_tw(prm%sum_N_tw,Nmembers),source=0.0_pReal)
    allocate(dst%Lambda_tr(prm%sum_N_tr,Nmembers),source=0.0_pReal)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(dislotwin)')

  end do

end function plastic_dislotwin_init


!--------------------------------------------------------------------------------------------------
!> @brief Return the homogenized elasticity matrix.
!--------------------------------------------------------------------------------------------------
module function plastic_dislotwin_homogenizedC(ph,en) result(homogenizedC)

  integer,     intent(in) :: &
    ph, en
  real(pReal), dimension(6,6) :: &
    homogenizedC, &
    C
  real(pReal), dimension(:,:,:), allocatable :: &
    C66_tw, &
    C66_tr
  integer :: i
  real(pReal) :: f_matrix


  C = elastic_C66(ph,en)

  associate(prm => param(ph), stt => state(ph))

    f_matrix = 1.0_pReal &
             - sum(stt%f_tw(1:prm%sum_N_tw,en)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,en))

    homogenizedC = f_matrix * C

    twinActive: if (prm%sum_N_tw > 0) then
      C66_tw    = lattice_C66_twin(prm%N_tw,C,phase_lattice(ph),phase_cOverA(ph))
      do i=1,prm%sum_N_tw
        homogenizedC = homogenizedC &
                     + stt%f_tw(i,en)*C66_tw(1:6,1:6,i)
      end do
     end if twinActive

    transActive: if (prm%sum_N_tr > 0) then
      C66_tr    = lattice_C66_trans(prm%N_tr,C,'hP',prm%cOverA_hP)
      do i=1,prm%sum_N_tr
        homogenizedC = homogenizedC &
                     + stt%f_tr(i,en)*C66_tr(1:6,1:6,i)
      end do
    end if transActive

  end associate

end function plastic_dislotwin_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
module subroutine dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

  real(pReal), dimension(3,3),     intent(out) :: Lp
  real(pReal), dimension(3,3,3,3), intent(out) :: dLp_dMp
  real(pReal), dimension(3,3),     intent(in)  :: Mp
  integer,                         intent(in)  :: ph,en

  integer :: i,k,l,m,n
  real(pReal) :: &
    f_matrix,StressRatio_p,&
    E_kB_T, &
    ddot_gamma_dtau, &
    tau, &
    T
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_sl
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw,ddot_gamma_dtau_tw
  real(pReal), dimension(param(ph)%sum_N_tr) :: &
    dot_gamma_tr,ddot_gamma_dtau_tr
  real(pReal):: dot_gamma_sb
  real(pReal), dimension(3,3) :: eigVectors, P_sb
  real(pReal), dimension(3)   :: eigValues
  real(pReal), dimension(3,6), parameter :: &
    sb_sComposition = &
      reshape(real([&
         1, 0, 1, &
         1, 0,-1, &
         1, 1, 0, &
         1,-1, 0, &
         0, 1, 1, &
         0, 1,-1  &
         ],pReal),[ 3,6]), &
    sb_mComposition = &
      reshape(real([&
         1, 0,-1, &
         1, 0,+1, &
         1,-1, 0, &
         1, 1, 0, &
         0, 1,-1, &
         0, 1, 1  &
         ],pReal),[ 3,6])


  T = thermal_T(ph,en)
  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal

  associate(prm => param(ph), stt => state(ph))

    f_matrix = 1.0_pReal &
             - sum(stt%f_tw(1:prm%sum_N_tw,en)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,en))

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl,ddot_gamma_dtau_sl)
    slipContribution: do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_sl(i) * prm%P_sl(k,l,i) * prm%P_sl(m,n,i)
    end do slipContribution

    if (prm%sum_N_tw > 0) call kinetics_tw(Mp,T,dot_gamma_sl,ph,en,dot_gamma_tw,ddot_gamma_dtau_tw)
    twinContibution: do i = 1, prm%sum_N_tw
      Lp = Lp + dot_gamma_tw(i)*prm%P_tw(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_tw(i)* prm%P_tw(k,l,i)*prm%P_tw(m,n,i)
    end do twinContibution

    if (prm%sum_N_tr > 0) call kinetics_tr(Mp,T,dot_gamma_sl,ph,en,dot_gamma_tr,ddot_gamma_dtau_tr)
    transContibution: do i = 1, prm%sum_N_tr
      Lp = Lp + dot_gamma_tr(i)*prm%P_tr(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau_tr(i)* prm%P_tr(k,l,i)*prm%P_tr(m,n,i)
    end do transContibution

    Lp      = Lp      * f_matrix
    dLp_dMp = dLp_dMp * f_matrix

    shearBandingContribution: if (dNeq0(prm%v_sb)) then

      E_kB_T = prm%E_sb/(K_B*T)
      call math_eigh33(eigValues,eigVectors,Mp)                                                     ! is Mp symmetric by design?

      do i = 1,6
        P_sb = 0.5_pReal * math_outer(matmul(eigVectors,sb_sComposition(1:3,i)),&
                                      matmul(eigVectors,sb_mComposition(1:3,i)))
        tau = math_tensordot(Mp,P_sb)

        significantShearBandStress: if (abs(tau) > tol_math_check) then
          StressRatio_p = (abs(tau)/prm%xi_sb)**prm%p_sb
          dot_gamma_sb = sign(prm%v_sb*exp(-E_kB_T*(1-StressRatio_p)**prm%q_sb), tau)
          ddot_gamma_dtau = abs(dot_gamma_sb)*E_kB_T*prm%p_sb*prm%q_sb/prm%xi_sb &
                          * (abs(tau)/prm%xi_sb)**(prm%p_sb-1.0_pReal) &
                          * (1.0_pReal-StressRatio_p)**(prm%q_sb-1.0_pReal)

          Lp = Lp + dot_gamma_sb * P_sb
          forall (k=1:3,l=1:3,m=1:3,n=1:3) &
            dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                             + ddot_gamma_dtau * P_sb(k,l) * P_sb(m,n)
        end if significantShearBandStress
      end do

    end if shearBandingContribution

    end associate

end subroutine dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the rate of change of microstructure.
!--------------------------------------------------------------------------------------------------
module function dislotwin_dotState(Mp,ph,en) result(dotState)

  real(pReal), dimension(3,3),  intent(in):: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  integer :: i
  real(pReal) :: &
    f_matrix, &
    d_hat, &
    v_cl, &                                                                                         !< climb velocity
    tau, &
    sigma_cl, &                                                                                     !< climb stress
    b_d                                                                                             !< ratio of Burgers vector to stacking fault width
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_rho_dip_formation, &
    dot_rho_dip_climb, &
    dot_gamma_sl
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    dot_gamma_tw
  real(pReal), dimension(param(ph)%sum_N_tr) :: &
    dot_gamma_tr
  real(pReal) :: &
    mu, nu, &
    T


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph), &
            dot_rho_mob => dotState(indexDotState(ph)%rho_mob(1):indexDotState(ph)%rho_mob(2)), &
            dot_rho_dip => dotState(indexDotState(ph)%rho_dip(1):indexDotState(ph)%rho_dip(2)), &
            abs_dot_gamma_sl => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)), &
            dot_f_tw => dotState(indexDotState(ph)%f_tw(1):indexDotState(ph)%f_tw(2)), &
            dot_f_tr => dotState(indexDotState(ph)%f_tr(1):indexDotState(ph)%f_tr(2)))

    mu = elastic_mu(ph,en)
    nu = elastic_nu(ph,en)
    T = thermal_T(ph,en)

    f_matrix = 1.0_pReal &
             - sum(stt%f_tw(1:prm%sum_N_tw,en)) &
             - sum(stt%f_tr(1:prm%sum_N_tr,en))

    call kinetics_sl(Mp,T,ph,en,dot_gamma_sl)
    abs_dot_gamma_sl = abs(dot_gamma_sl)

    slipState: do i = 1, prm%sum_N_sl
      tau = math_tensordot(Mp,prm%P_sl(1:3,1:3,i))

      significantSlipStress: if (dEq0(tau) .or. prm%omitDipoles) then
        dot_rho_dip_formation(i) = 0.0_pReal
        dot_rho_dip_climb(i) = 0.0_pReal
      else significantSlipStress
        d_hat = 3.0_pReal*mu*prm%b_sl(i)/(16.0_pReal*PI*abs(tau))
        d_hat = math_clip(d_hat, right = dst%Lambda_sl(i,en))
        d_hat = math_clip(d_hat, left  = prm%d_caron(i))

        dot_rho_dip_formation(i) = 2.0_pReal*(d_hat-prm%d_caron(i))/prm%b_sl(i) &
                                 * stt%rho_mob(i,en)*abs_dot_gamma_sl(i)

        if (dEq(d_hat,prm%d_caron(i))) then
          dot_rho_dip_climb(i) = 0.0_pReal
        else
          ! Argon & Moffat, Acta Metallurgica, Vol. 29, pg 293 to 299, 1981
          sigma_cl = dot_product(prm%n0_sl(1:3,i),matmul(Mp,prm%n0_sl(1:3,i)))
          if (prm%extendedDislocations) then
            b_d = 24.0_pReal*PI*(1.0_pReal - nu)/(2.0_pReal + nu) * prm%Gamma_sf%at(T) / (mu*prm%b_sl(i))
          else
            b_d = 1.0_pReal
          end if
          v_cl = 2.0_pReal*prm%omega*b_d**2*exp(-prm%Q_cl/(K_B*T)) &
               * (exp(abs(sigma_cl)*prm%b_sl(i)**3/(K_B*T)) - 1.0_pReal)

          dot_rho_dip_climb(i) = 4.0_pReal*v_cl*stt%rho_dip(i,en) &
                               / (d_hat-prm%d_caron(i))
        end if
      end if significantSlipStress
    end do slipState

    dot_rho_mob = abs_dot_gamma_sl/(prm%b_sl*dst%Lambda_sl(:,en)) &
                - dot_rho_dip_formation &
                - 2.0_pReal*prm%d_caron/prm%b_sl * stt%rho_mob(:,en)*abs_dot_gamma_sl

    dot_rho_dip = dot_rho_dip_formation &
                - 2.0_pReal*prm%d_caron/prm%b_sl * stt%rho_dip(:,en)*abs_dot_gamma_sl &
                - dot_rho_dip_climb

    if (prm%sum_N_tw > 0) call kinetics_tw(Mp,T,abs_dot_gamma_sl,ph,en,dot_gamma_tw)
    dot_f_tw = f_matrix*dot_gamma_tw/prm%gamma_char_tw

    if (prm%sum_N_tr > 0) call kinetics_tr(Mp,T,abs_dot_gamma_sl,ph,en,dot_gamma_tr)
    dot_f_tr = f_matrix*dot_gamma_tr/prm%gamma_char_tr

  end associate

end function dislotwin_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine dislotwin_dependentState(ph,en)

  integer,       intent(in) :: &
    ph, &
    en

  real(pReal) :: &
    sumf_tw, sumf_tr
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    inv_lambda_sl
  real(pReal), dimension(param(ph)%sum_N_tw) :: &
    inv_lambda_tw_tw, &                                                                             !< 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
    f_over_t_tw
   real(pReal), dimension(param(ph)%sum_N_tr) :: &
    inv_lambda_tr_tr, &                                                                             !< 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite
    f_over_t_tr
  real(pReal) :: &
    mu


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en)
    sumf_tw = sum(stt%f_tw(1:prm%sum_N_tw,en))
    sumf_tr = sum(stt%f_tr(1:prm%sum_N_tr,en))

    !* rescaled volume fraction for topology
    f_over_t_tw = stt%f_tw(1:prm%sum_N_tw,en)/prm%t_tw                                              ! this is per system ...
    f_over_t_tr = sumf_tr/prm%t_tr                                                                  ! but this not
                                                                                                    ! ToDo ...Physically correct, but naming could be adjusted

    inv_lambda_sl = sqrt(matmul(prm%forestProjection,stt%rho_mob(:,en)+stt%rho_dip(:,en)))/prm%i_sl
    if (prm%sum_N_tw > 0 .and. prm%sum_N_sl > 0) &
      inv_lambda_sl = inv_lambda_sl + matmul(prm%h_sl_tw,f_over_t_tw)/(1.0_pReal-sumf_tw)
    if (prm%sum_N_tr > 0 .and. prm%sum_N_sl > 0) &
      inv_lambda_sl = inv_lambda_sl + matmul(prm%h_sl_tr,f_over_t_tr)/(1.0_pReal-sumf_tr)
    dst%Lambda_sl(:,en) = prm%D / (1.0_pReal+prm%D*inv_lambda_sl)

    inv_lambda_tw_tw = matmul(prm%h_tw_tw,f_over_t_tw)/(1.0_pReal-sumf_tw)
    dst%Lambda_tw(:,en) = prm%i_tw*prm%D/(1.0_pReal+prm%D*inv_lambda_tw_tw)

    inv_lambda_tr_tr = matmul(prm%h_tr_tr,f_over_t_tr)/(1.0_pReal-sumf_tr)
    dst%Lambda_tr(:,en) = prm%i_tr*prm%D/(1.0_pReal+prm%D*inv_lambda_tr_tr)

    !* threshold stress for dislocation motion
    dst%tau_pass(:,en) = mu*prm%b_sl* sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,en)+stt%rho_dip(:,en)))

  end associate

end subroutine dislotwin_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislotwin_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_mob')
          call results_writeDataset(stt%rho_mob,group,trim(prm%output(ou)), &
                                       'mobile dislocation density','1/m²',prm%systems_sl)
        case('rho_dip')
          call results_writeDataset(stt%rho_dip,group,trim(prm%output(ou)), &
                                       'dislocation dipole density','1/m²',prm%systems_sl)
        case('gamma_sl')
          call results_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                       'plastic shear','1',prm%systems_sl)
        case('Lambda_sl')
          call results_writeDataset(dst%Lambda_sl,group,trim(prm%output(ou)), &
                                       'mean free path for slip','m',prm%systems_sl)
        case('tau_pass')
          call results_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                       'passing stress for slip','Pa',prm%systems_sl)

        case('f_tw')
          call results_writeDataset(stt%f_tw,group,trim(prm%output(ou)), &
                                       'twinned volume fraction','m³/m³',prm%systems_tw)
        case('Lambda_tw')
          call results_writeDataset(dst%Lambda_tw,group,trim(prm%output(ou)), &
                                       'mean free path for twinning','m',prm%systems_tw)

        case('f_tr')
          if (prm%sum_N_tr>0) call results_writeDataset(stt%f_tr,group,trim(prm%output(ou)), &
                                                       'martensite volume fraction','m³/m³')

      end select

    end do

  end associate

end subroutine plastic_dislotwin_results


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_sl(Mp,T,ph,en, &
                            dot_gamma_sl,ddot_gamma_dtau_sl,tau_sl)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(param(ph)%sum_N_sl), intent(out) :: &
    dot_gamma_sl
  real(pReal), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau_sl, &
    tau_sl

  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    ddot_gamma_dtau
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau, &
    stressRatio, &
    StressRatio_p, &
    Q_kB_T, &
    v_wait_inverse, &                                                                               !< inverse of the effective velocity of a dislocation waiting at obstacles (unsigned)
    v_run_inverse, &                                                                                !< inverse of the velocity of a free moving dislocation (unsigned)
    dV_wait_inverse_dTau, &
    dV_run_inverse_dTau, &
    dV_dTau, &
    tau_eff                                                                                         !< effective resolved stress
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    tau = [(math_tensordot(Mp,prm%P_sl(1:3,1:3,i)),i = 1, prm%sum_N_sl)]

    tau_eff = abs(tau)-dst%tau_pass(:,en)

    significantStress: where(tau_eff > tol_math_check)
      stressRatio    = tau_eff/prm%tau_0
      StressRatio_p  = stressRatio** prm%p
      Q_kB_T = prm%Q_sl/(K_B*T)
      v_wait_inverse = exp(Q_kB_T*(1.0_pReal-StressRatio_p)** prm%q) &
                     / prm%v_0
      v_run_inverse  = prm%B/(tau_eff*prm%b_sl)

      dot_gamma_sl = sign(stt%rho_mob(:,en)*prm%b_sl/(v_wait_inverse+v_run_inverse),tau)

      dV_wait_inverse_dTau = -1.0_pReal * v_wait_inverse * prm%p * prm%q * Q_kB_T &
                           * (stressRatio**(prm%p-1.0_pReal)) &
                           * (1.0_pReal-StressRatio_p)**(prm%q-1.0_pReal) &
                           / prm%tau_0
      dV_run_inverse_dTau  = -1.0_pReal * v_run_inverse/tau_eff
      dV_dTau              = -1.0_pReal * (dV_wait_inverse_dTau+dV_run_inverse_dTau) &
                           / (v_wait_inverse+v_run_inverse)**2
      ddot_gamma_dtau = dV_dTau*stt%rho_mob(:,en)*prm%b_sl
    else where significantStress
      dot_gamma_sl    = 0.0_pReal
      ddot_gamma_dtau = 0.0_pReal
    end where significantStress

  end associate

  if (present(ddot_gamma_dtau_sl)) ddot_gamma_dtau_sl = ddot_gamma_dtau
  if (present(tau_sl))             tau_sl             = tau

end subroutine kinetics_sl


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on twin systems and their derivatives with respect to resolved
!         stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_tw(Mp,T,abs_dot_gamma_sl,ph,en,&
                            dot_gamma_tw,ddot_gamma_dtau_tw)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(param(ph)%sum_N_sl), intent(in) :: &
    abs_dot_gamma_sl
  real(pReal), dimension(param(ph)%sum_N_tw), intent(out) :: &
    dot_gamma_tw
  real(pReal), dimension(param(ph)%sum_N_tw), optional, intent(out) :: &
    ddot_gamma_dtau_tw

  real :: &
    tau, tau_r, tau_hat, &
    dot_N_0, &
    x0, V, &
    Gamma_sf, &
    mu, nu, &
    P_ncs, dP_ncs_dtau, &
    P, dP_dtau
  integer, dimension(2) :: &
    s
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en)
    nu = elastic_nu(ph,en)
    Gamma_sf = prm%Gamma_sf%at(T)

    tau_hat = 3.0_pReal*prm%b_tw(1)*mu/prm%L_tw &
            + Gamma_sf/(3.0_pReal*prm%b_tw(1))
    x0 = mu*prm%b_sl(1)**2*(2.0_pReal+nu)/(Gamma_sf*8.0_pReal*PI*(1.0_pReal-nu))
    tau_r = mu*prm%b_sl(1)/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%x_c)+cos(PI/3.0_pReal)/x0)

    do i = 1, prm%sum_N_tw
      tau = math_tensordot(Mp,prm%P_tw(1:3,1:3,i))

      if (tau > tol_math_check .and. tau < tau_r) then
        P = exp(-(tau_hat/tau)**prm%r(i))
        dP_dTau = prm%r(i) * (tau_hat/tau)**prm%r(i)/tau * P

        s = prm%fcc_twinNucleationSlipPair(1:2,i)
        dot_N_0 = sum(abs_dot_gamma_sl(s(2:1:-1))*(stt%rho_mob(s,en)+stt%rho_dip(s,en)))/(prm%L_tw*3.0_pReal)

        P_ncs = 1.0_pReal-exp(-prm%V_cs/(K_B*T)*(tau_r-tau))
        dP_ncs_dtau = prm%V_cs / (K_B * T) * (P_ncs - 1.0_pReal)

        V = PI/4.0_pReal*dst%Lambda_tw(i,en)**2*prm%t_tw(i)
        dot_gamma_tw(i) = V*dot_N_0*P_ncs*P*prm%gamma_char_tw(i)
        if (present(ddot_gamma_dtau_tw)) &
          ddot_gamma_dtau_tw(i) = V*dot_N_0*(P*dP_ncs_dtau + P_ncs*dP_dtau)*prm%gamma_char_tw(i)
      else
        dot_gamma_tw(i) = 0.0_pReal
        if (present(ddot_gamma_dtau_tw)) ddot_gamma_dtau_tw(i) = 0.0_pReal
      end if
    end do

  end associate

end subroutine kinetics_tw


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on transformation systems and their derivatives with respect to
!         resolved stress.
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_tr(Mp,T,abs_dot_gamma_sl,ph,en,&
                            dot_gamma_tr,ddot_gamma_dtau_tr)

  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    ph, &
    en
  real(pReal), dimension(param(ph)%sum_N_sl), intent(in) :: &
    abs_dot_gamma_sl
  real(pReal), dimension(param(ph)%sum_N_tr), intent(out) :: &
    dot_gamma_tr
  real(pReal), dimension(param(ph)%sum_N_tr), optional, intent(out) :: &
    ddot_gamma_dtau_tr

  real :: &
    tau, tau_r, tau_hat, &
    dot_N_0, &
    x0, V, &
    Gamma_sf, &
    mu, nu, &
    P_ncs, dP_ncs_dtau, &
    P, dP_dtau
  integer, dimension(2) :: &
    s
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    mu = elastic_mu(ph,en)
    nu = elastic_nu(ph,en)
    Gamma_sf = prm%Gamma_sf%at(T)

    tau_hat = 3.0_pReal*prm%b_tr(1)*mu/prm%L_tr &
            + (Gamma_sf + (prm%h/prm%V_mol - 2.0_pReal*prm%rho)*prm%Delta_G%at(T))/(3.0_pReal*prm%b_tr(1))
    x0 = mu*prm%b_sl(1)**2*(2.0_pReal+nu)/(Gamma_sf*8.0_pReal*PI*(1.0_pReal-nu))
    tau_r = mu*prm%b_sl(1)/(2.0_pReal*PI)*(1.0_pReal/(x0+prm%x_c)+cos(PI/3.0_pReal)/x0)

    do i = 1, prm%sum_N_tr
      tau = math_tensordot(Mp,prm%P_tr(1:3,1:3,i))

      if (tau > tol_math_check .and. tau < tau_r) then
        P = exp(-(tau_hat/tau)**prm%s(i))
        dP_dTau = prm%s(i) * (tau_hat/tau)**prm%s(i)/tau * P

        s = prm%fcc_twinNucleationSlipPair(1:2,i)
        dot_N_0 = sum(abs_dot_gamma_sl(s(2:1:-1))*(stt%rho_mob(s,en)+stt%rho_dip(s,en)))/(prm%L_tr*3.0_pReal)

        P_ncs = 1.0_pReal-exp(-prm%V_cs/(K_B*T)*(tau_r-tau))
        dP_ncs_dtau = prm%V_cs / (K_B * T) * (P_ncs - 1.0_pReal)

        V = PI/4.0_pReal*dst%Lambda_tr(i,en)**2*prm%t_tr(i)
        dot_gamma_tr(i) = V*dot_N_0*P_ncs*P*prm%gamma_char_tr
        if (present(ddot_gamma_dtau_tr)) &
          ddot_gamma_dtau_tr(i) = V*dot_N_0*(P*dP_ncs_dtau + P_ncs*dP_dtau)*prm%gamma_char_tr
      else
        dot_gamma_tr(i) = 0.0_pReal
        if (present(ddot_gamma_dtau_tr)) ddot_gamma_dtau_tr(i) = 0.0_pReal
      end if
    end do

  end associate

end subroutine kinetics_tr

end submodule dislotwin
