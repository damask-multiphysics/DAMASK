!--------------------------------------------------------------------------------------------------
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for plasticity including dislocation flux
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) nonlocal
  use geometry_plastic_nonlocal, only: &
    nCellNeighbors  => geometry_plastic_nonlocal_nIPneighbors, &
    IPneighborhood  => geometry_plastic_nonlocal_IPneighborhood, &
    IPvolume0       => geometry_plastic_nonlocal_IPvolume0, &
    IParea0         => geometry_plastic_nonlocal_IParea0, &
    IPareaNormal0   => geometry_plastic_nonlocal_IPareaNormal0, &
    geometry_plastic_nonlocal_disable

  type :: tGeometry
    real(pREAL), dimension(:), allocatable :: v_0
    real(pREAL), dimension(:,:), allocatable :: a_0, x_0
    real(pREAL), dimension(:,:,:), allocatable :: n_0
    integer, dimension(:,:,:), allocatable :: IPneighborhood
  end type tGeometry

  type(tGeometry), dimension(:), allocatable :: geom

  ! storage order of dislocation types
  integer, dimension(*), parameter :: &
    sgl = [1,2,3,4,5,6,7,8]                                                                         !< signed (single)
  integer, dimension(*), parameter :: &
    edg = [1,2,5,6,9], &                                                                            !< edge
    scr = [3,4,7,8,10]                                                                              !< screw
  integer, dimension(*), parameter :: &
    mob = [1,2,3,4], &                                                                              !< mobile
    imm = [5,6,7,8]                                                                                 !< immobile (blocked)
  integer, dimension(*), parameter :: &
    dip = [9,10], &                                                                                 !< dipole
    imm_edg = imm(1:2), &                                                                           !< immobile edge
    imm_scr = imm(3:4)                                                                              !< immobile screw
  integer, parameter :: &
    mob_edg_pos = 1, &                                                                              !< mobile edge positive
    mob_edg_neg = 2, &                                                                              !< mobile edge negative
    mob_scr_pos = 3, &                                                                              !< mobile screw positive
    mob_scr_neg = 4                                                                                 !< mobile screw positive

  real(pREAL), dimension(:,:,:,:,:,:), allocatable :: &
    compatibility                                                                                   !< slip system compatibility between en and my neighbors

  type :: tInitialParameters                                                                        !< container type for internal constitutive parameters
    real(pREAL) :: &
      sigma_rho_u, &                                                                                !< standard deviation of scatter in initial dislocation density
      random_rho_u, &
      random_rho_u_binning
    real(pREAL), dimension(:), allocatable :: &
      rho_u_ed_pos_0, &                                                                             !< initial edge_pos dislocation density
      rho_u_ed_neg_0, &                                                                             !< initial edge_neg dislocation density
      rho_u_sc_pos_0, &                                                                             !< initial screw_pos dislocation density
      rho_u_sc_neg_0, &                                                                             !< initial screw_neg dislocation density
      rho_d_ed_0, &                                                                                 !< initial edge dipole dislocation density
      rho_d_sc_0                                                                                    !< initial screw dipole dislocation density
    integer,     dimension(:), allocatable :: &
      N_sl
  end type tInitialParameters

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pREAL) :: &
      V_at, &                                                                                       !< atomic volume
      D_0, &                                                                                        !< prefactor for self-diffusion coefficient
      Q_cl, &                                                                                       !< activation enthalpy for diffusion
      atol_rho, &                                                                                   !< absolute tolerance for dislocation density in state integration
      rho_significant, &                                                                            !< density considered significant
      rho_min, &                                                                                    !< number of dislocations considered significant
      w, &                                                                                          !< width of a doubkle kink in multiples of the Burgers vector length b
      Q_sol, &                                                                                      !< activation energy for solid solution in J
      f_sol, &                                                                                      !< solid solution obstacle size in multiples of the Burgers vector length
      c_sol, &                                                                                      !< concentration of solid solution in atomic parts
      p, &                                                                                          !< parameter for kinetic law (Kocks,Argon,Ashby)
      q, &                                                                                          !< parameter for kinetic law (Kocks,Argon,Ashby)
      B, &                                                                                          !< drag coefficient in Pa s
      nu_a, &                                                                                       !< attack frequency in Hz
      chi_surface, &                                                                                !< transmissivity at free surface
      chi_GB, &                                                                                     !< transmissivity at grain boundary (identified by different texture)
      C_CFL, &                                                                                      !< safety factor for CFL flux condition
      f_ed_mult, &                                                                                  !< factor that determines how much edge dislocations contribute to multiplication (0...1)
      f_F, &
      f_ed, &
      mu, &
      nu
    real(pREAL), dimension(:),      allocatable :: &
      i_sl, &                                                                                       !< mean free path prefactor for each
      b_sl                                                                                          !< absolute length of Burgers vector [m]
    real(pREAL), dimension(:,:),   allocatable :: &
      slip_normal, &
      slip_direction, &
      slip_transverse, &
      minDipoleHeight, &                                                                            ! minimum stable dipole height edge and screw
      peierlsstress, &                                                                              ! edge and screw
      h_sl_sl ,&                                                                                    !< coefficients for slip-slip interaction
      forestProjection_Edge, &                                                                      !< matrix of forest projections of edge dislocations
      forestProjection_Screw                                                                        !< matrix of forest projections of screw dislocations
    real(pREAL), dimension(:,:,:), allocatable :: &
      P_sl, &                                                                                       !< Schmid contribution
      P_nS_pos, &
      P_nS_neg                                                                                      !< combined projection of Schmid and non-Schmid contributions to the resolved shear stress (only for screws)
    integer :: &
      sum_N_sl = 0
    integer,     dimension(:),     allocatable :: &
      colinearSystem                                                                                !< colinear system to the active slip system (only valid for fcc!)
    character(len=:),              allocatable :: &
      isotropic_bound
    character(len=pSTRLEN), dimension(:), allocatable :: &
      output
    logical :: &
      shortRangeStressCorrection                                                                    !< use of short range stress correction by excess density gradient term
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl
  end type tParameters

  type :: tNonlocalDependentState
    real(pREAL), allocatable, dimension(:,:) :: &
      tau_pass, &
      tau_back, &
      rho_forest, &
      max_dipole_height
    real(pREAL), allocatable, dimension(:,:,:,:,:) :: &
      compatibility
  end type tNonlocalDependentState

  type :: tNonlocalState
    real(pREAL), pointer, dimension(:,:) :: &
      rho, &                                                                                        ! < all dislocations
        rho_sgl, &
          rho_sgl_mob, &
            rho_sgl_mob_edg_pos, &
            rho_sgl_mob_edg_neg, &
            rho_sgl_mob_scr_pos, &
            rho_sgl_mob_scr_neg, &
          rho_sgl_imm, &
            rho_sgl_imm_edg_pos, &
            rho_sgl_imm_edg_neg, &
            rho_sgl_imm_scr_pos, &
            rho_sgl_imm_scr_neg, &
        rho_dip, &
          rho_dip_edg, &
          rho_dip_scr, &
      gamma, &
      v, &
        v_edg_pos, &
        v_edg_neg, &
        v_scr_pos, &
        v_scr_neg
  end type tNonlocalState

 type(tNonlocalState), allocatable, dimension(:) :: &
    deltaState, &
    dotState, &
    state, &
    state0

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters
  type(tNonlocalDependentState), dimension(:), allocatable :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_nonlocal_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, &
    Nmembers, &
    sizeState, sizeDotState, sizeDeltaState, &
    s1, s2
  real(pREAL), dimension(:,:), allocatable :: &
    a_nS                                                                                            !< non-Schmid coefficients
  character(len=:), allocatable :: &
    refs, &
    extmsg
  type(tInitialParameters) :: &
    ini
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl

  myPlasticity = plastic_active('nonlocal')
  if (count(myPlasticity) == 0) then
    call geometry_plastic_nonlocal_disable()
    return
  end if

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:nonlocal init  -+>>>'

  print'(/,1x,a)', 'C. Reuber et al., Acta Materialia 71:333–348, 2014'
  print'(  1x,a)', 'https://doi.org/10.1016/j.actamat.2014.03.012'

  print'(/,1x,a)', 'C. Kords, Dissertation RWTH Aachen, 2014'
  print'(  1x,a)', 'http://publications.rwth-aachen.de/record/229993'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(geom(phases%length))
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(state0(phases%length))
  allocate(dotState(phases%length))
  allocate(deltaState(phases%length))
  allocate(dependentState(phases%length))
  extmsg = ''

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph),  dot => dotState(ph),   stt => state(ph), &
              st0 => state0(ph), del => deltaState(ph), dst => dependentState(ph))

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

    plasticState(ph)%nonlocal = pl%get_asBool('flux',defaultVal=.True.)
    prm%isotropic_bound = pl%get_asStr('isotropic_bound',defaultVal='isostrain')
    prm%atol_rho = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)

    ini%N_sl     = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(ini%N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%systems_sl = crystal_labels_slip(ini%N_sl,phase_lattice(ph))
      prm%P_sl = crystal_SchmidMatrix_slip(ini%N_sl,phase_lattice(ph), phase_cOverA(ph))

      if (phase_lattice(ph) == 'cI') then
        allocate(a_nS(3,size(pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray))),source=0.0_pREAL)          ! anticipating parameters for all three families
        a_nS(1,:) = pl%get_as1dReal('a_nonSchmid_110',defaultVal=emptyRealArray)
        prm%P_nS_pos = crystal_SchmidMatrix_slip(ini%N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=+1)
        prm%P_nS_neg = crystal_SchmidMatrix_slip(ini%N_sl,phase_lattice(ph),phase_cOverA(ph),nonSchmidCoefficients=a_nS,sense=-1)
        deallocate(a_nS)
      else
        prm%P_nS_pos = +prm%P_sl
        prm%P_nS_neg = -prm%P_sl
      end if

      prm%h_sl_sl = crystal_interaction_SlipBySlip(ini%N_sl,pl%get_as1dReal('h_sl-sl'), &
                                                   phase_lattice(ph))

      prm%forestProjection_edge  = crystal_forestProjection_edge (ini%N_sl,phase_lattice(ph),&
                                                                  phase_cOverA(ph))
      prm%forestProjection_screw = crystal_forestProjection_screw(ini%N_sl,phase_lattice(ph),&
                                                                  phase_cOverA(ph))

      prm%slip_direction  = crystal_slip_direction (ini%N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%slip_transverse = crystal_slip_transverse(ini%N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%slip_normal     = crystal_slip_normal    (ini%N_sl,phase_lattice(ph),phase_cOverA(ph))

      ! collinear systems (only for octahedral slip systems in fcc)
      allocate(prm%colinearSystem(prm%sum_N_sl), source = -1)
      do s1 = 1, prm%sum_N_sl
        do s2 = 1, prm%sum_N_sl
           if (all(dEq0 (math_cross(prm%slip_direction(1:3,s1),prm%slip_direction(1:3,s2)))) .and. &
               any(dNeq0(math_cross(prm%slip_normal   (1:3,s1),prm%slip_normal   (1:3,s2))))) &
             prm%colinearSystem(s1) = s2
        end do
      end do

      ini%rho_u_ed_pos_0  = pl%get_as1dReal('rho_u_ed_pos_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_ed_neg_0  = pl%get_as1dReal('rho_u_ed_neg_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_sc_pos_0  = pl%get_as1dReal('rho_u_sc_pos_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_sc_neg_0  = pl%get_as1dReal('rho_u_sc_neg_0',   requiredSize=size(ini%N_sl))
      ini%rho_d_ed_0      = pl%get_as1dReal('rho_d_ed_0',       requiredSize=size(ini%N_sl))
      ini%rho_d_sc_0      = pl%get_as1dReal('rho_d_sc_0',       requiredSize=size(ini%N_sl))

      prm%i_sl            = math_expand(pl%get_as1dReal('i_sl', requiredSize=size(ini%N_sl)),ini%N_sl)
      prm%b_sl            = math_expand(pl%get_as1dReal('b_sl', requiredSize=size(ini%N_sl)),ini%N_sl)

      allocate(prm%minDipoleHeight(prm%sum_N_sl,2))
      prm%minDipoleHeight(:,1) = math_expand(pl%get_as1dReal('d_ed', requiredSize=size(ini%N_sl)),ini%N_sl)
      prm%minDipoleHeight(:,2) = math_expand(pl%get_as1dReal('d_sc', requiredSize=size(ini%N_sl)),ini%N_sl)

      allocate(prm%peierlsstress(prm%sum_N_sl,2))
      prm%peierlsstress(:,1)   = math_expand(pl%get_as1dReal('tau_Peierls_ed', requiredSize=size(ini%N_sl)),ini%N_sl)
      prm%peierlsstress(:,2)   = math_expand(pl%get_as1dReal('tau_Peierls_sc', requiredSize=size(ini%N_sl)),ini%N_sl)

      prm%rho_significant = pl%get_asReal('rho_significant')
      prm%rho_min         = pl%get_asReal('rho_min', 0.0_pREAL)
      prm%C_CFL           = pl%get_asReal('C_CFL',defaultVal=2.0_pREAL)

      prm%V_at            = pl%get_asReal('V_at')
      prm%D_0             = pl%get_asReal('D_0')
      prm%Q_cl            = pl%get_asReal('Q_cl')
      prm%f_F             = pl%get_asReal('f_F')
      prm%f_ed            = pl%get_asReal('f_ed')
      prm%w               = pl%get_asReal('w')
      prm%Q_sol           = pl%get_asReal('Q_sol')
      prm%f_sol           = pl%get_asReal('f_sol')
      prm%c_sol           = pl%get_asReal('c_sol')

      prm%p               = pl%get_asReal('p_sl')
      prm%q               = pl%get_asReal('q_sl')
      prm%B               = pl%get_asReal('B')
      prm%nu_a            = pl%get_asReal('nu_a')

      ! ToDo: discuss logic
      ini%sigma_rho_u     = pl%get_asReal('sigma_rho_u')
      ini%random_rho_u    = pl%get_asReal('random_rho_u',defaultVal= 0.0_pREAL)
      if (pl%contains('random_rho_u')) &
        ini%random_rho_u_binning = pl%get_asReal('random_rho_u_binning',defaultVal=0.0_pREAL) !ToDo: useful default?
     ! if (rhoSglRandom(instance) < 0.0_pREAL) &
     ! if (rhoSglRandomBinning(instance) <= 0.0_pREAL) &

      prm%chi_surface                 = pl%get_asReal('chi_surface',defaultVal=1.0_pREAL)
      prm%chi_GB                      = pl%get_asReal('chi_GB',     defaultVal=-1.0_pREAL)
      prm%f_ed_mult                   = pl%get_asReal('f_ed_mult')
      prm%shortRangeStressCorrection  = pl%get_asBool('short_range_stress_correction', defaultVal = .false.)


!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (any(prm%b_sl             <  0.0_pREAL)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%i_sl             <= 0.0_pREAL)) extmsg = trim(extmsg)//' i_sl'

      if (any(ini%rho_u_ed_pos_0   <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_u_ed_pos_0'
      if (any(ini%rho_u_ed_neg_0   <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_u_ed_neg_0'
      if (any(ini%rho_u_sc_pos_0   <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_u_sc_pos_0'
      if (any(ini%rho_u_sc_neg_0   <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_u_sc_neg_0'
      if (any(ini%rho_d_ed_0       <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_d_ed_0'
      if (any(ini%rho_d_sc_0       <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_d_sc_0'

      if (any(prm%peierlsstress    <  0.0_pREAL)) extmsg = trim(extmsg)//' tau_peierls'
      if (any(prm%minDipoleHeight  <  0.0_pREAL)) extmsg = trim(extmsg)//' d_ed or d_sc'

      if (prm%B                    <  0.0_pREAL)  extmsg = trim(extmsg)//' B'
      if (prm%Q_cl                 <  0.0_pREAL)  extmsg = trim(extmsg)//' Q_cl'
      if (prm%nu_a                <=  0.0_pREAL)  extmsg = trim(extmsg)//' nu_a'
      if (prm%w                   <=  0.0_pREAL)  extmsg = trim(extmsg)//' w'
      if (prm%D_0                 <   0.0_pREAL)  extmsg = trim(extmsg)//' D_0'
      if (prm%V_at                <=  0.0_pREAL)  extmsg = trim(extmsg)//' V_at'                   ! ToDo: in dislotungsten, the atomic volume is given as a factor

      if (prm%rho_min             <   0.0_pREAL)  extmsg = trim(extmsg)//' rho_min'
      if (prm%rho_significant     <   0.0_pREAL)  extmsg = trim(extmsg)//' rho_significant'
      if (prm%atol_rho            <   0.0_pREAL)  extmsg = trim(extmsg)//' atol_rho'
      if (prm%C_CFL               <   0.0_pREAL)  extmsg = trim(extmsg)//' C_CFL'

      if (prm%p    <= 0.0_pREAL .or. prm%p > 1.0_pREAL) extmsg = trim(extmsg)//' p_sl'
      if (prm%q    <  1.0_pREAL .or. prm%q > 2.0_pREAL) extmsg = trim(extmsg)//' q_sl'

      if (prm%f_F  < 0.0_pREAL  .or. prm%f_F  > 1.0_pREAL) &
                                                  extmsg = trim(extmsg)//' f_F'
      if (prm%f_ed <  0.0_pREAL .or. prm%f_ed > 1.0_pREAL) &
                                                  extmsg = trim(extmsg)//' f_ed'

      if (prm%Q_sol               <=  0.0_pREAL)  extmsg = trim(extmsg)//' Q_sol'
      if (prm%f_sol               <=  0.0_pREAL)  extmsg = trim(extmsg)//' f_sol'
      if (prm%c_sol               <=  0.0_pREAL)  extmsg = trim(extmsg)//' c_sol'

      if (prm%chi_GB              >   1.0_pREAL)  extmsg = trim(extmsg)//' chi_GB'
      if (prm%chi_surface  < 0.0_pREAL .or. prm%chi_surface  > 1.0_pREAL) &
                                                  extmsg = trim(extmsg)//' chi_surface'

      if (prm%f_ed_mult    < 0.0_pREAL .or. prm%f_ed_mult    > 1.0_pREAL) &
                                                  extmsg = trim(extmsg)//' f_ed_mult'

    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers  = count(material_ID_phase == ph)
    sizeDotState = size(['rho_sgl_mob_edg_pos', 'rho_sgl_mob_edg_neg', &
                         'rho_sgl_mob_scr_pos', 'rho_sgl_mob_scr_neg', &
                         'rho_sgl_imm_edg_pos', 'rho_sgl_imm_edg_neg', &
                         'rho_sgl_imm_scr_pos', 'rho_sgl_imm_scr_neg', &
                         'rho_dip_edg        ', 'rho_dip_scr        ', &
                         'gamma              ' ]) * prm%sum_N_sl                                    !< "basic" microstructural state variables that are independent from other state variables
    sizeState          = sizeDotState &
                       + size([ 'velocityEdgePos     ','velocityEdgeNeg     ', &
                                'velocityScrewPos    ','velocityScrewNeg    ']) * prm%sum_N_sl      !< other dependent state variables that are not updated by microstructure
    sizeDeltaState            = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,sizeDeltaState,0)     ! ToDo: state structure does not follow convention

    allocate(geom(ph)%v_0(Nmembers))
    allocate(geom(ph)%a_0(nCellNeighbors,Nmembers))
    allocate(geom(ph)%x_0(3,Nmembers))
    allocate(geom(ph)%n_0(3,nCellNeighbors,Nmembers))
    allocate(geom(ph)%IPneighborhood(3,nCellNeighbors,Nmembers))
    call storeGeometry(ph)

    if (plasticState(ph)%nonlocal .and. .not. allocated(IPneighborhood)) &
      call IO_error(212,ext_msg='IPneighborhood does not exist')

    st0%rho => plasticState(ph)%state0                             (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    stt%rho => plasticState(ph)%state                              (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    dot%rho => plasticState(ph)%dotState                           (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    del%rho => plasticState(ph)%deltaState                         (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    plasticState(ph)%atol(1:10*prm%sum_N_sl) = prm%atol_rho

      stt%rho_sgl => plasticState(ph)%state                        (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
      dot%rho_sgl => plasticState(ph)%dotState                     (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
      del%rho_sgl => plasticState(ph)%deltaState                   (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)

        st0%rho_sgl_mob => plasticState(ph)%state0                 (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
        stt%rho_sgl_mob => plasticState(ph)%state                  (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
        dot%rho_sgl_mob => plasticState(ph)%dotState               (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
        del%rho_sgl_mob => plasticState(ph)%deltaState             (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)

            stt%rho_sgl_mob_edg_pos => plasticState(ph)%state      (0*prm%sum_N_sl+1: 1*prm%sum_N_sl,:)
            dot%rho_sgl_mob_edg_pos => plasticState(ph)%dotState   (0*prm%sum_N_sl+1: 1*prm%sum_N_sl,:)
            del%rho_sgl_mob_edg_pos => plasticState(ph)%deltaState (0*prm%sum_N_sl+1: 1*prm%sum_N_sl,:)

            stt%rho_sgl_mob_edg_neg => plasticState(ph)%state      (1*prm%sum_N_sl+1: 2*prm%sum_N_sl,:)
            dot%rho_sgl_mob_edg_neg => plasticState(ph)%dotState   (1*prm%sum_N_sl+1: 2*prm%sum_N_sl,:)
            del%rho_sgl_mob_edg_neg => plasticState(ph)%deltaState (1*prm%sum_N_sl+1: 2*prm%sum_N_sl,:)

            stt%rho_sgl_mob_scr_pos => plasticState(ph)%state      (2*prm%sum_N_sl+1: 3*prm%sum_N_sl,:)
            dot%rho_sgl_mob_scr_pos => plasticState(ph)%dotState   (2*prm%sum_N_sl+1: 3*prm%sum_N_sl,:)
            del%rho_sgl_mob_scr_pos => plasticState(ph)%deltaState (2*prm%sum_N_sl+1: 3*prm%sum_N_sl,:)

            stt%rho_sgl_mob_scr_neg => plasticState(ph)%state      (3*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
            dot%rho_sgl_mob_scr_neg => plasticState(ph)%dotState   (3*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
            del%rho_sgl_mob_scr_neg => plasticState(ph)%deltaState (3*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)

        stt%rho_sgl_imm => plasticState(ph)%state                  (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
        dot%rho_sgl_imm => plasticState(ph)%dotState               (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
        del%rho_sgl_imm => plasticState(ph)%deltaState             (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)

            stt%rho_sgl_imm_edg_pos => plasticState(ph)%state      (4*prm%sum_N_sl+1: 5*prm%sum_N_sl,:)
            dot%rho_sgl_imm_edg_pos => plasticState(ph)%dotState   (4*prm%sum_N_sl+1: 5*prm%sum_N_sl,:)
            del%rho_sgl_imm_edg_pos => plasticState(ph)%deltaState (4*prm%sum_N_sl+1: 5*prm%sum_N_sl,:)

            stt%rho_sgl_imm_edg_neg => plasticState(ph)%state      (5*prm%sum_N_sl+1: 6*prm%sum_N_sl,:)
            dot%rho_sgl_imm_edg_neg => plasticState(ph)%dotState   (5*prm%sum_N_sl+1: 6*prm%sum_N_sl,:)
            del%rho_sgl_imm_edg_neg => plasticState(ph)%deltaState (5*prm%sum_N_sl+1: 6*prm%sum_N_sl,:)

            stt%rho_sgl_imm_scr_pos => plasticState(ph)%state      (6*prm%sum_N_sl+1: 7*prm%sum_N_sl,:)
            dot%rho_sgl_imm_scr_pos => plasticState(ph)%dotState   (6*prm%sum_N_sl+1: 7*prm%sum_N_sl,:)
            del%rho_sgl_imm_scr_pos => plasticState(ph)%deltaState (6*prm%sum_N_sl+1: 7*prm%sum_N_sl,:)

            stt%rho_sgl_imm_scr_neg => plasticState(ph)%state      (7*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
            dot%rho_sgl_imm_scr_neg => plasticState(ph)%dotState   (7*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
            del%rho_sgl_imm_scr_neg => plasticState(ph)%deltaState (7*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)

      stt%rho_dip => plasticState(ph)%state                        (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
      dot%rho_dip => plasticState(ph)%dotState                     (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
      del%rho_dip => plasticState(ph)%deltaState                   (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)

        stt%rho_dip_edg => plasticState(ph)%state                  (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)
        dot%rho_dip_edg => plasticState(ph)%dotState               (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)
        del%rho_dip_edg => plasticState(ph)%deltaState             (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)

        stt%rho_dip_scr => plasticState(ph)%state                  (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
        dot%rho_dip_scr => plasticState(ph)%dotState               (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
        del%rho_dip_scr => plasticState(ph)%deltaState             (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)

    stt%gamma => plasticState(ph)%state                      (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    dot%gamma => plasticState(ph)%dotState                   (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    del%gamma => plasticState(ph)%deltaState                 (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    plasticState(ph)%atol(10*prm%sum_N_sl+1:11*prm%sum_N_sl )  = pl%get_asReal('atol_gamma', defaultVal = 1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(10*prm%sum_N_sl+1:11*prm%sum_N_sl) < 0.0_pREAL)) &
      extmsg = trim(extmsg)//' atol_gamma'

    stt%v          => plasticState(ph)%state                 (11*prm%sum_N_sl + 1:15*prm%sum_N_sl,1:Nmembers)
    st0%v          => plasticState(ph)%state0                (11*prm%sum_N_sl + 1:15*prm%sum_N_sl,1:Nmembers)
        stt%v_edg_pos  => plasticState(ph)%state             (11*prm%sum_N_sl + 1:12*prm%sum_N_sl,1:Nmembers)
        stt%v_edg_neg  => plasticState(ph)%state             (12*prm%sum_N_sl + 1:13*prm%sum_N_sl,1:Nmembers)
        stt%v_scr_pos  => plasticState(ph)%state             (13*prm%sum_N_sl + 1:14*prm%sum_N_sl,1:Nmembers)
        stt%v_scr_neg  => plasticState(ph)%state             (14*prm%sum_N_sl + 1:15*prm%sum_N_sl,1:Nmembers)

    allocate(dst%tau_pass(prm%sum_N_sl,Nmembers),source=0.0_pREAL)
    allocate(dst%tau_back(prm%sum_N_sl,Nmembers),source=0.0_pREAL)
    allocate(dst%rho_forest(prm%sum_N_sl,Nmembers),source=0.0_pREAL)
    allocate(dst%max_dipole_height(2*prm%sum_N_sl,Nmembers),source=0.0_pREAL)                       ! edge and screw
    allocate(dst%compatibility(2,maxval(param%sum_N_sl),maxval(param%sum_N_sl),nCellNeighbors,Nmembers),source=0.0_pREAL)
    end associate

    if (Nmembers > 0) call stateInit(ini,ph,Nmembers)

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

  allocate(compatibility(2,maxval(param%sum_N_sl),maxval(param%sum_N_sl),nCellNeighbors,&
                         discretization_nIPs,discretization_Nelems), source=0.0_pREAL)

end function plastic_nonlocal_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates quantities characterizing the microstructure
!--------------------------------------------------------------------------------------------------
module subroutine nonlocal_dependentState(ph, en)

  integer, intent(in) :: &
    ph, &
    en

  integer :: &
    en_nbr, &                                                                                       ! neighbor phase entry
    el_nbr, &                                                                                       ! element number of neighboring material point
    ip_nbr, &                                                                                       ! integration point of neighboring material point
    c, &                                                                                            ! index of dilsocation character (edge, screw)
    s, &                                                                                            ! slip system index
    dir, &
    n
  real(pREAL) :: &
    FVsize, &
    nRealNeighbors, &                                                                               ! number of really existing neighbors
    mu, &
    nu
  integer, dimension(2) :: &
    neighbors
  real(pREAL), dimension(2) :: &
    rhoExcessGradient, &
    rhoExcessGradient_over_rho, &
    rho_sum
  real(pREAL), dimension(3) :: &
    rhoExcessDifferences, &
    normal_latticeConf
  real(pREAL), dimension(3,3) :: &
    invFe, &                                                                                        !< inverse of elastic deformation gradient
    invFp, &                                                                                        !< inverse of plastic deformation gradient
    connections, &
    invConnections
  real(pREAL), dimension(3,nCellNeighbors) :: &
    connection_latticeConf
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    rho_edg_delta, &
    rho_scr_delta
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rho_0, &
    rho_0_nbr
  real(pREAL), dimension(param(ph)%sum_N_sl,param(ph)%sum_N_sl) :: &
    myInteractionMatrix                                                                             ! corrected slip interaction matrix
  real(pREAL), dimension(param(ph)%sum_N_sl,nCellNeighbors) :: &
    rho_edg_delta_nbr, &
    rho_scr_delta_nbr
  real(pREAL), dimension(2,maxval(param%sum_N_sl),nCellNeighbors) :: &
    rho_delta_nbr, &                                                                                ! excess density at neighboring material point
    rho_sum_nbr                                                                                     ! total density at neighboring material point
  real(pREAL), dimension(3,param(ph)%sum_N_sl,2) :: &
    m                                                                                               ! direction of dislocation motion

  associate(prm => param(ph),dst => dependentState(ph), stt => state(ph))

  mu = elastic_mu(ph,en,prm%isotropic_bound)
  nu = elastic_nu(ph,en,prm%isotropic_bound)
  rho = getRho(ph,en)

  dst%rho_forest(:,en) = matmul(prm%forestProjection_Edge, sum(abs(rho(:,edg)),2)) &
                       + matmul(prm%forestProjection_Screw,sum(abs(rho(:,scr)),2))


  ! coefficients are corrected for the line tension effect
  ! (see Kubin,Devincre,Hoc; 2008; Modeling dislocation storage rates and mean free paths in face-centered cubic crystals)
  if (any(phase_lattice(ph) == ['cI','cF'])) then
    myInteractionMatrix = prm%h_sl_sl &
                        * spread((  1.0_pREAL - prm%f_F &
                                   + prm%f_F &
                                   * log(0.35_pREAL * prm%b_sl * sqrt(max(dst%rho_forest(:,en),prm%rho_significant))) &
                                   / log(0.35_pREAL * prm%b_sl * 1e6_pREAL))**2,2,prm%sum_N_sl)
  else
    myInteractionMatrix = prm%h_sl_sl
  end if

  dst%tau_pass(:,en) = mu * prm%b_sl &
                     * sqrt(matmul(myInteractionMatrix,sum(abs(rho),2)))

!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

  !#################################################################################################
  ! ToDo: MD: this is most likely only correct for F_i = I
  !#################################################################################################

  rho_0 = getRho0(ph,en)
  if (plasticState(ph)%nonlocal .and. prm%shortRangeStressCorrection) then
    invFp = math_inv33(phase_mechanical_Fp(ph)%data(1:3,1:3,en))
    invFe = math_inv33(phase_mechanical_Fe(ph)%data(1:3,1:3,en))

    rho_edg_delta = rho_0(:,mob_edg_pos) - rho_0(:,mob_edg_neg)
    rho_scr_delta = rho_0(:,mob_scr_pos) - rho_0(:,mob_scr_neg)

    FVsize = geom(ph)%v_0(en)**(1.0_pREAL/3.0_pREAL)

    !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities

    nRealNeighbors = 0.0_pREAL
    rho_sum_nbr = 0.0_pREAL
    do n = 1,nCellNeighbors
      el_nbr = geom(ph)%IPneighborhood(1,n,en)
      ip_nbr = geom(ph)%IPneighborhood(2,n,en)

      if (el_nbr > 0 .and. ip_nbr > 0) then
        if (material_ID_phase(1,(el_nbr-1)*discretization_nIPs + ip_nbr) == ph) then
            en_nbr = material_entry_phase(1,(el_nbr-1)*discretization_nIPs + ip_nbr)
            nRealNeighbors = nRealNeighbors + 1.0_pREAL
            rho_0_nbr = getRho0(ph,en_nbr)

            rho_edg_delta_nbr(:,n) = rho_0_nbr(:,mob_edg_pos) - rho_0_nbr(:,mob_edg_neg)
            rho_scr_delta_nbr(:,n) = rho_0_nbr(:,mob_scr_pos) - rho_0_nbr(:,mob_scr_neg)

            rho_sum_nbr(1,:,n) = sum(abs(rho_0_nbr(:,edg)),2)
            rho_sum_nbr(2,:,n) = sum(abs(rho_0_nbr(:,scr)),2)

            connection_latticeConf(1:3,n) = matmul(invFe, geom(ph)%x_0(1:3,en_nbr) &
                                          - geom(ph)%x_0(1:3,en))
            normal_latticeConf = matmul(transpose(invFp), geom(ph)%n_0(1:3,n,en))
            if (math_inner(normal_latticeConf,connection_latticeConf(1:3,n)) < 0.0_pREAL) &         ! neighboring connection points in opposite direction to face normal: must be periodic image
              connection_latticeConf(1:3,n) = normal_latticeConf * geom(ph)%v_0(en)/geom(ph)%a_0(n,en)  ! instead take the surface normal scaled with the diameter of the cell
        else
          ! local neighbor or different lattice structure or different constitution instance -> use central values instead
          connection_latticeConf(1:3,n) = 0.0_pREAL
          rho_edg_delta_nbr(:,n) = rho_edg_delta
          rho_scr_delta_nbr(:,n) = rho_scr_delta
        end if
      else
        ! free surface -> use central values instead
        connection_latticeConf(1:3,n) = 0.0_pREAL
        rho_edg_delta_nbr(:,n) = rho_edg_delta
        rho_scr_delta_nbr(:,n) = rho_scr_delta
      end if
    end do

    rho_delta_nbr(1,:,:) = rho_edg_delta_nbr
    rho_delta_nbr(2,:,:) = rho_scr_delta_nbr

    !* loop through the slip systems and calculate the dislocation gradient by
    !* 1. interpolation of the excess density in the neighorhood
    !* 2. interpolation of the dead dislocation density in the central volume
    m(1:3,:,1) =  prm%slip_direction
    m(1:3,:,2) = -prm%slip_transverse

    do s = 1,prm%sum_N_sl

      ! gradient from interpolation of neighboring excess density ...
      do c = 1,2
        do dir = 1,3
          neighbors(1) = 2 * dir - 1
          neighbors(2) = 2 * dir
          connections(dir,1:3) = connection_latticeConf(1:3,neighbors(1)) &
                               - connection_latticeConf(1:3,neighbors(2))
          rhoExcessDifferences(dir) = rho_delta_nbr(c,s,neighbors(1)) &
                                    - rho_delta_nbr(c,s,neighbors(2))
        end do
        invConnections = math_inv33(connections)
        if (all(dEq0(invConnections))) call IO_error(-1,ext_msg='back stress calculation: inversion error')

        rhoExcessGradient(c) = math_inner(m(1:3,s,c), matmul(invConnections,rhoExcessDifferences))
      end do

        ! ... plus gradient from deads ...
      rhoExcessGradient(1) = rhoExcessGradient(1) + sum(rho(s,imm_edg)) / FVsize
      rhoExcessGradient(2) = rhoExcessGradient(2) + sum(rho(s,imm_scr)) / FVsize

        ! ... normalized with the total density ...
      rho_sum(1) = (sum(abs(rho(s,edg))) + sum(rho_sum_nbr(1,s,:))) / (1.0_pREAL + nRealNeighbors)
      rho_sum(2) = (sum(abs(rho(s,scr))) + sum(rho_sum_nbr(2,s,:))) / (1.0_pREAL + nRealNeighbors)

      rhoExcessGradient_over_rho = 0.0_pREAL
      where(rho_sum > 0.0_pREAL) rhoExcessGradient_over_rho = rhoExcessGradient / rho_sum

        ! ... gives the local stress correction when multiplied with a factor
      dst%tau_back(s,en) = - mu * prm%b_sl(s) / (2.0_pREAL * PI) &
                         * (  rhoExcessGradient_over_rho(1) / (1.0_pREAL - nu) &
                            + rhoExcessGradient_over_rho(2))
    end do
  end if

 end associate

end subroutine nonlocal_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
module subroutine nonlocal_LpAndItsTangent(Lp,dLp_dMp, &
                                                   Mp,ph,en)
  real(pREAL), dimension(3,3), intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(3,3), intent(in) :: &
    Mp
                                                                                                    !< derivative of Lp with respect to Mp
  integer :: &
    i, j, k, l, &
    t, &                                                                                            !< dislocation type
    s                                                                                               !< index of my current slip system
  real(pREAL), dimension(param(ph)%sum_N_sl,8) :: &
    rho_sgl                                                                                          !< single dislocation densities (including blocked)
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    rho
  real(pREAL), dimension(param(ph)%sum_N_sl,4) :: &
    v, &                                                                                            !< velocity
    tauNS, &                                                                                        !< resolved shear stress including non Schmid and backstress terms
    dv_dtau, &                                                                                      !< velocity derivative with respect to the shear stress
    dv_dtauNS                                                                                       !< velocity derivative with respect to the shear stress
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau, &                                                                                          !< resolved shear stress including backstress terms
    dot_gamma                                                                                       !< shear rate
  real(pREAL) :: &
    Temperature                                                                                     !< temperature


  Temperature = thermal_T(ph,en)
  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph),dst=>dependentState(ph),stt=>state(ph))

    !*** shortcut to state variables
    rho = getRho(ph,en)
    rho_sgl = rho(:,sgl)

    do s = 1,prm%sum_N_sl
      tau(s) = math_tensordot(Mp, prm%P_sl(1:3,1:3,s))
      tauNS(s,1) = tau(s)
      tauNS(s,2) = tau(s)
      if (tau(s) > 0.0_pREAL) then
        tauNS(s,3) = math_tensordot(Mp, +prm%P_nS_pos(1:3,1:3,s))
        tauNS(s,4) = math_tensordot(Mp, -prm%P_nS_neg(1:3,1:3,s))
      else
        tauNS(s,3) = math_tensordot(Mp, +prm%P_nS_neg(1:3,1:3,s))
        tauNS(s,4) = math_tensordot(Mp, -prm%P_nS_pos(1:3,1:3,s))
      end if
    end do
    tauNS = tauNS + spread(dst%tau_back(:,en),2,4)
    tau   = tau   + dst%tau_back(:,en)

    ! edges
    call kinetics(v(:,1), dv_dtau(:,1), dv_dtauNS(:,1), &
                  tau, tauNS(:,1), dst%tau_pass(:,en),1,Temperature, ph)
    v(:,2)         = v(:,1)
    dv_dtau(:,2)   = dv_dtau(:,1)
    dv_dtauNS(:,2) = dv_dtauNS(:,1)

    !screws
    do t = 3,4
      call kinetics(v(:,t), dv_dtau(:,t), dv_dtauNS(:,t), &
                    tau, tauNS(:,t), dst%tau_pass(:,en),2,Temperature, ph)
    end do

    stt%v(:,en) = pack(v,.true.)

    !*** Bauschinger effect
    forall (s = 1:prm%sum_N_sl, t = 5:8, rho_sgl(s,t) * v(s,t-4) < 0.0_pREAL) &
      rho_sgl(s,t-4) = rho_sgl(s,t-4) + abs(rho_sgl(s,t))

    dot_gamma = sum(rho_sgl(:,1:4) * v, 2) * prm%b_sl

    do s = 1,prm%sum_N_sl
      Lp = Lp + dot_gamma(s) * prm%P_sl(1:3,1:3,s)
      forall (i=1:3,j=1:3,k=1:3,l=1:3) &
        dLp_dMp(i,j,k,l) = dLp_dMp(i,j,k,l) &
            + prm%P_sl(i,j,s) * prm%P_sl(k,l,s) &
            * sum(rho_sgl(s,1:4) * dv_dtau(s,1:4)) * prm%b_sl(s) &
            + prm%P_sl(i,j,s) &
            * (+ prm%P_nS_pos(k,l,s) * rho_sgl(s,3) * dv_dtauNS(s,3) &
               - prm%P_nS_neg(k,l,s) * rho_sgl(s,4) * dv_dtauNS(s,4))  * prm%b_sl(s)
    end do

  end associate

end subroutine nonlocal_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_deltaState(Mp,ph,en)

  real(pREAL), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< MandelStress
  integer, intent(in) :: &
    ph, &
    en

  integer :: &
    c, &                                                                                            ! character of dislocation
    t, &                                                                                            ! type of dislocation
    s                                                                                               ! index of my current slip system
  real(pREAL) :: &
    mu, &
    nu
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    deltaRhoRemobilization, &                                                                       ! density increment by remobilization
    deltaRhoDipole2SingleStress                                                                     ! density increment by dipole dissociation (by stress change)
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    rho                                                                                             ! current  dislocation densities
  real(pREAL), dimension(param(ph)%sum_N_sl,4) :: &
    v                                                                                               ! dislocation glide velocity
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau                                                                                             ! current resolved shear stress
  real(pREAL), dimension(param(ph)%sum_N_sl,2) :: &
    rho_dip, &                                                                                      ! current dipole dislocation densities (screw and edge dipoles)
    dUpper, &                                                                                       ! current maximum stable dipole distance for edges and screws
    dUpperOld, &                                                                                    ! old maximum stable dipole distance for edges and screws
    deltaDUpper                                                                                     ! change in maximum stable dipole distance for edges and screws


  associate(prm => param(ph),dst => dependentState(ph),del => deltaState(ph), stt=>state(ph))

  mu = elastic_mu(ph,en,prm%isotropic_bound)
  nu = elastic_nu(ph,en,prm%isotropic_bound)

  !*** shortcut to state variables
  v = reshape(stt%v(:,en),[prm%sum_N_sl,4])
  dUpperOld = reshape(dst%max_dipole_height(:,en),[prm%sum_N_sl,2])

  rho =  getRho(ph,en)
  rho_dip = rho(:,dip)

  !****************************************************************************
  !*** dislocation remobilization (bauschinger effect)
  where(rho(:,imm) * v < 0.0_pREAL)
    deltaRhoRemobilization(:,mob) = abs(rho(:,imm))
    deltaRhoRemobilization(:,imm) =   - rho(:,imm)
    rho(:,mob) = rho(:,mob) + abs(rho(:,imm))
    rho(:,imm) = 0.0_pREAL
  elsewhere
    deltaRhoRemobilization(:,mob) = 0.0_pREAL
    deltaRhoRemobilization(:,imm) = 0.0_pREAL
  endwhere
  deltaRhoRemobilization(:,dip) = 0.0_pREAL

  !****************************************************************************
  !*** calculate dipole formation and dissociation by stress change

  !*** calculate limits for stable dipole height
  do s = 1,prm%sum_N_sl
    tau(s) = math_tensordot(Mp, prm%P_sl(1:3,1:3,s)) +dst%tau_back(s,en)
    if (abs(tau(s)) < 1.0e-15_pREAL) tau(s) = 1.0e-15_pREAL
  end do

  dUpper(:,1) = mu * prm%b_sl/(8.0_pREAL * PI * (1.0_pREAL - nu) * abs(tau))
  dUpper(:,2) = mu * prm%b_sl/(4.0_pREAL * PI * abs(tau))

  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(:,1) = min(1.0_pREAL/sqrt(sum(abs(rho(:,edg)),2)),dUpper(:,1))
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(:,2) = min(1.0_pREAL/sqrt(sum(abs(rho(:,scr)),2)),dUpper(:,2))

  dUpper = max(dUpper,prm%minDipoleHeight)
  deltaDUpper = dUpper - dUpperOld


  !*** dissociation by stress increase
  deltaRhoDipole2SingleStress = 0.0_pREAL
  forall (c=1:2, s=1:prm%sum_N_sl, deltaDUpper(s,c) < 0.0_pREAL .and. &
                                          dNeq0(dUpperOld(s,c) - prm%minDipoleHeight(s,c))) &
    deltaRhoDipole2SingleStress(s,8+c) = rho_dip(s,c) * deltaDUpper(s,c) &
                                       / (dUpperOld(s,c) - prm%minDipoleHeight(s,c))

  forall (t=1:4) deltaRhoDipole2SingleStress(:,t) = -0.5_pREAL * deltaRhoDipole2SingleStress(:,(t-1)/2+9)
  dst%max_dipole_height(:,en) = pack(dUpper,.true.)

  plasticState(ph)%deltaState(:,en) = 0.0_pREAL
  del%rho(:,en) = reshape(deltaRhoRemobilization + deltaRhoDipole2SingleStress, [10*prm%sum_N_sl])

  end associate

end subroutine plastic_nonlocal_deltaState


!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
module subroutine nonlocal_dotState(Mp,timestep, &
                                    ph,en)

  real(pREAL), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< MandelStress
  real(pREAL), intent(in) :: &
    timestep                                                                                        !< substepped crystallite time increment
  integer, intent(in) :: &
    ph, &
    en

  integer ::  &
    c, &                                                                                            !< character of dislocation
    s                                                                                               !< index of my current slip system
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rhoDot, &                                                                                       !< density evolution
    rhoDotMultiplication, &                                                                         !< density evolution by multiplication
    rhoDotSingle2DipoleGlide, &                                                                     !< density evolution by dipole formation (by glide)
    rhoDotAthermalAnnihilation, &                                                                   !< density evolution by athermal annihilation
    rhoDotThermalAnnihilation                                                                       !< density evolution by thermal annihilation
  real(pREAL), dimension(param(ph)%sum_N_sl,8) :: &
    rho_sgl                                                                                         !< current single dislocation densities (positive/negative screw and edge without dipoles)
  real(pREAL), dimension(param(ph)%sum_N_sl,4) :: &
    v, &                                                                                            !< current dislocation glide velocity
    dot_gamma                                                                                       !< shear rates
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    tau, &                                                                                          !< current resolved shear stress
    v_climb                                                                                         !< climb velocity of edge dipoles
  real(pREAL), dimension(param(ph)%sum_N_sl,2) :: &
    rho_dip, &                                                                                      !< current dipole dislocation densities (screw and edge dipoles)
    dLower, &                                                                                       !< minimum stable dipole distance for edges and screws
    dUpper                                                                                          !< current maximum stable dipole distance for edges and screws
  real(pREAL) :: &
    D_SD, &
    mu, &
    nu, Temperature

  if (timestep <= 0.0_pREAL) then
    plasticState(ph)%dotState = 0.0_pREAL
    return
  end if

  associate(prm => param(ph), dst => dependentState(ph), dot => dotState(ph), &
            stt => state(ph), st0 => state0(ph))

  mu = elastic_mu(ph,en,prm%isotropic_bound)
  nu = elastic_nu(ph,en,prm%isotropic_bound)
  Temperature = thermal_T(ph,en)

  tau = 0.0_pREAL
  dot_gamma = 0.0_pREAL

  rho = getRho(ph,en)
  rho_sgl = rho(:,sgl)
  rho_dip = rho(:,dip)

  v = reshape(stt%v(:,en),[prm%sum_N_sl,4])
  dot_gamma = rho_sgl(:,1:4) * v * spread(prm%b_sl,2,4)


  ! limits for stable dipole height
  do s = 1,prm%sum_N_sl
    tau(s) = math_tensordot(Mp, prm%P_sl(1:3,1:3,s)) + dst%tau_back(s,en)
    if (abs(tau(s)) < 1.0e-15_pREAL) tau(s) = 1.0e-15_pREAL
  end do

  dLower = prm%minDipoleHeight
  dUpper(:,1) = mu * prm%b_sl/(8.0_pREAL * PI * (1.0_pREAL - nu) * abs(tau))
  dUpper(:,2) = mu * prm%b_sl/(4.0_pREAL * PI * abs(tau))

  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(:,1) = min(1.0_pREAL/sqrt(sum(abs(rho(:,edg)),2)),dUpper(:,1))
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(:,2) = min(1.0_pREAL/sqrt(sum(abs(rho(:,scr)),2)),dUpper(:,2))

  dUpper = max(dUpper,dLower)


  ! dislocation multiplication
  rhoDotMultiplication = 0.0_pREAL
  isBCC: if (phase_lattice(ph) == 'cI') then
    forall (s = 1:prm%sum_N_sl, sum(abs(v(s,1:4))) > 0.0_pREAL)
      rhoDotMultiplication(s,1:2) = sum(abs(dot_gamma(s,3:4))) / prm%b_sl(s) &                      ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(dst%rho_forest(s,en)) / prm%i_sl(s) ! &                    ! mean free path
                                  ! * 2.0_pREAL * sum(abs(v(s,3:4))) / sum(abs(v(s,1:4)))           ! ratio of screw to overall velocity determines edge generation
      rhoDotMultiplication(s,3:4) = sum(abs(dot_gamma(s,3:4))) /prm%b_sl(s) &                       ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(dst%rho_forest(s,en)) / prm%i_sl(s) ! &                    ! mean free path
                                  ! * 2.0_pREAL * sum(abs(v(s,1:2))) / sum(abs(v(s,1:4)))           ! ratio of edge to overall velocity determines screw generation
    endforall

  else isBCC
    rhoDotMultiplication(:,1:4) = spread( &
          (sum(abs(dot_gamma(:,1:2)),2) * prm%f_ed_mult + sum(abs(dot_gamma(:,3:4)),2)) &
        * sqrt(dst%rho_forest(:,en)) / prm%i_sl / prm%b_sl, 2, 4)                                   ! eq. 3.26
  end if isBCC


  !****************************************************************************
  ! dipole formation and annihilation

  ! formation by glide
  do c = 1,2
    rhoDotSingle2DipoleGlide(:,2*c-1) = -2.0_pREAL * dUpper(:,c) / prm%b_sl &
                                                   * (      rho_sgl(:,2*c-1)  * abs(dot_gamma(:,2*c)) &   ! negative mobile --> positive mobile
                                                      +     rho_sgl(:,2*c)    * abs(dot_gamma(:,2*c-1)) & ! positive mobile --> negative mobile
                                                      + abs(rho_sgl(:,2*c+4)) * abs(dot_gamma(:,2*c-1)))  ! positive mobile --> negative immobile

    rhoDotSingle2DipoleGlide(:,2*c) = -2.0_pREAL * dUpper(:,c) / prm%b_sl &
                                                 * (      rho_sgl(:,2*c-1)  * abs(dot_gamma(:,2*c)) &     ! negative mobile --> positive mobile
                                                    +     rho_sgl(:,2*c)    * abs(dot_gamma(:,2*c-1)) &   ! positive mobile --> negative mobile
                                                    + abs(rho_sgl(:,2*c+3)) * abs(dot_gamma(:,2*c)))      ! negative mobile --> positive immobile

    rhoDotSingle2DipoleGlide(:,2*c+3) = -2.0_pREAL * dUpper(:,c) / prm%b_sl &
                                                   * rho_sgl(:,2*c+3) * abs(dot_gamma(:,2*c))             ! negative mobile --> positive immobile

    rhoDotSingle2DipoleGlide(:,2*c+4) = -2.0_pREAL * dUpper(:,c) / prm%b_sl &
                                                   * rho_sgl(:,2*c+4) * abs(dot_gamma(:,2*c-1))           ! positive mobile --> negative immobile

    rhoDotSingle2DipoleGlide(:,c+8) = abs(rhoDotSingle2DipoleGlide(:,2*c+3)) &
                                    + abs(rhoDotSingle2DipoleGlide(:,2*c+4)) &
                                    - rhoDotSingle2DipoleGlide(:,2*c-1) &
                                    - rhoDotSingle2DipoleGlide(:,2*c)
  end do


  ! athermal annihilation
  rhoDotAthermalAnnihilation = 0.0_pREAL
  forall (c=1:2) &
    rhoDotAthermalAnnihilation(:,c+8) = -2.0_pREAL * dLower(:,c) / prm%b_sl &
       * (  2.0_pREAL * (rho_sgl(:,2*c-1) * abs(dot_gamma(:,2*c)) + rho_sgl(:,2*c) * abs(dot_gamma(:,2*c-1))) & ! was single hitting single
          + 2.0_pREAL * (abs(rho_sgl(:,2*c+3)) * abs(dot_gamma(:,2*c)) + abs(rho_sgl(:,2*c+4)) * abs(dot_gamma(:,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
          + rho_dip(:,c) * (abs(dot_gamma(:,2*c-1)) + abs(dot_gamma(:,2*c))))                                  ! single knocks dipole constituent

  ! annihilated screw dipoles leave edge jogs behind on the colinear system
  if (phase_lattice(ph) == 'cF') &
    forall (s = 1:prm%sum_N_sl, prm%colinearSystem(s) > 0) &
      rhoDotAthermalAnnihilation(prm%colinearSystem(s),1:2) = - rhoDotAthermalAnnihilation(s,10) &
        * 0.25_pREAL * sqrt(dst%rho_forest(s,en)) * (dUpper(s,2) + dLower(s,2)) * prm%f_ed


  ! thermally activated annihilation of edge dipoles by climb
  rhoDotThermalAnnihilation = 0.0_pREAL
  D_SD = prm%D_0 * exp(-prm%Q_cl / (K_B * Temperature))                                              ! eq. 3.53
  v_climb = D_SD * mu * prm%V_at &
          / (PI * (1.0_pREAL-nu) * (dUpper(:,1) + dLower(:,1)) * K_B * Temperature)              ! eq. 3.54
  forall (s = 1:prm%sum_N_sl, dUpper(s,1) > dLower(s,1)) &
    rhoDotThermalAnnihilation(s,9) = max(- 4.0_pREAL * rho_dip(s,1) * v_climb(s) / (dUpper(s,1) - dLower(s,1)), &
                                         - rho_dip(s,1) / timestep - rhoDotAthermalAnnihilation(s,9) &
                                                                  - rhoDotSingle2DipoleGlide(s,9))  ! make sure that we do not annihilate more dipoles than we have

  rhoDot = rhoDotFlux(timestep, ph,en) &
         + rhoDotMultiplication &
         + rhoDotSingle2DipoleGlide &
         + rhoDotAthermalAnnihilation &
         + rhoDotThermalAnnihilation


  if (    any(rho(:,mob) + rhoDot(:,1:4)  * timestep < -prm%atol_rho) &
     .or. any(rho(:,dip) + rhoDot(:,9:10) * timestep < -prm%atol_rho)) then
    plasticState(ph)%dotState = IEEE_value(1.0_pREAL,IEEE_quiet_NaN)
  else
    dot%rho(:,en) = pack(rhoDot,.true.)
    dot%gamma(:,en) = sum(dot_gamma,2)
  end if

  end associate

end subroutine nonlocal_dotState


!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
#if __INTEL_COMPILER >= 2020
non_recursive function rhoDotFlux(timestep,ph,en)
#else
function rhoDotFlux(timestep,ph,en)
#endif
  real(pREAL), intent(in) :: &
    timestep                                                                                        !< substepped crystallite time increment
  integer, intent(in) :: &
    ph, &
    en

  integer ::  &
    ns, &                                                                                           !< short notation for the total number of active slip systems
    c, &                                                                                            !< character of dislocation
    n, &                                                                                            !< index of my current neighbor
    el_nbr, &                                                                                       !< element number of my neighbor
    ip_nbr, &                                                                                       !< integration point of my neighbor
    n_nbr, &                                                                                        !< neighbor index pointing to en when looking from my neighbor
    opposite_neighbor, &                                                                            !< index of my opposite neighbor
    opposite_ip, &                                                                                  !< ip of my opposite neighbor
    opposite_el, &                                                                                  !< element index of my opposite neighbor
    opposite_n, &                                                                                   !< neighbor index pointing to en when looking from my opposite neighbor
    t, &                                                                                            !< type of dislocation
    en_nbr,&                                                                                        !< neighbor phase entry
    ph_nbr,&                                                                                        !< neighbor phase ID
    topp, &                                                                                         !< type of dislocation with opposite sign to t
    s                                                                                               !< index of my current slip system
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rho_0, &                                                                                        !< dislocation density at beginning of time step
    rhoDotFlux                                                                                      !< density evolution by flux
  real(pREAL), dimension(param(ph)%sum_N_sl,4) :: &
    rho_0_sgl_mob, &                                                                                !< mobile dislocation densities of neighboring ip (positive/negative screw and edge)
    rho_0_sgl_mob_nbr, &                                                                            !< mobile dislocation densities of neighboring ip (positive/negative screw and edge)
    v, &                                                                                            !< dislocation glide velocity
    v_0, &
    v_0_nbr, &                                                                                      !< dislocation glide velocity of enighboring ip
    dot_gamma                                                                                       !< shear rates
  real(pREAL), dimension(3,param(ph)%sum_N_sl,4) :: &
    m                                                                                               !< direction of dislocation motion
  real(pREAL), dimension(3,3) :: &
    my_F, &                                                                                         !< my total deformation gradient
    F_nbr, &                                                                                        !< total deformation gradient of my neighbor
    my_Fe, &                                                                                        !< my elastic deformation gradient
    F_e_nbr, &                                                                                      !< elastic deformation gradient of my neighbor
    Favg                                                                                            !< average total deformation gradient of en and my neighbor
  real(pREAL), dimension(3) :: &
    normal_neighbor2me, &                                                                           !< interface normal pointing from my neighbor to en in neighbor's lattice configuration
    normal_neighbor2me_defConf, &                                                                   !< interface normal pointing from my neighbor to en in shared deformed configuration
    normal_me2neighbor, &                                                                           !< interface normal pointing from en to my neighbor in my lattice configuration
    normal_me2neighbor_defConf                                                                      !< interface normal pointing from en to my neighbor in shared deformed configuration
  real(pREAL) :: &
    a, &                                                                                            !< area of the current interface
    transmissivity, &                                                                               !< overall transmissivity of dislocation flux to neighboring material point
    lineLength                                                                                      !< dislocation line length leaving the current interface


  associate(prm => param(ph), &
            dst => dependentState(ph), &
            stt => state(ph), &
            st0 => state0(ph))

  ns = prm%sum_N_sl

  dot_gamma = 0.0_pREAL

  rho = getRho(ph,en)
  rho_0 = getRho0(ph,en)
  rho_0_sgl_mob = rho_0(:,mob)

  v = reshape(stt%v(:,en),[prm%sum_N_sl,4])                                                         !ToDo: MD: I think we should use state0 here
  dot_gamma = rho(:,mob) * v * spread(prm%b_sl,2,4)

  v_0 = reshape(st0%v(:,en),[prm%sum_N_sl,4])

  !****************************************************************************
  !*** calculate dislocation fluxes (only for nonlocal plasticity)
  rhoDotFlux = 0.0_pREAL
  if (plasticState(ph)%nonlocal) then

    !*** check CFL (Courant-Friedrichs-Lewy) condition for flux
    if (any( abs(dot_gamma) > 0.0_pREAL &                                                           ! any active slip system ...
            .and. prm%C_CFL * abs(v_0) * timestep &
                > geom(ph)%v_0(en)/ maxval(geom(ph)%a_0(:,en)))) then                               ! ...with velocity above critical value (we use the reference volume and area for simplicity here)
      rhoDotFlux = IEEE_value(1.0_pREAL,IEEE_quiet_NaN)                                             ! enforce cutback
      return
    end if


    !*** be aware of the definition of slip_transverse = slip_direction x slip_normal !!!
    !*** opposite sign to our t vector in the (s,t,n) triplet !!!

    m(1:3,:,1) =  prm%slip_direction
    m(1:3,:,2) = -prm%slip_direction
    m(1:3,:,3) = -prm%slip_transverse
    m(1:3,:,4) =  prm%slip_transverse

    my_F = phase_mechanical_F(ph)%data(1:3,1:3,en)
    my_Fe = matmul(my_F, math_inv33(phase_mechanical_Fp(ph)%data(1:3,1:3,en)))

    neighbors: do n = 1,nCellNeighbors

      el_nbr = geom(ph)%IPneighborhood(1,n,en)
      ip_nbr = geom(ph)%IPneighborhood(2,n,en)
      n_nbr  = geom(ph)%IPneighborhood(3,n,en)
      ph_nbr = material_ID_phase(1,(el_nbr-1)*discretization_nIPs + ip_nbr)
      en_nbr = material_entry_phase(1,(el_nbr-1)*discretization_nIPs + ip_nbr)

      opposite_neighbor = n + mod(n,2) - mod(n+1,2)
      opposite_el = geom(ph)%IPneighborhood(1,opposite_neighbor,en)
      opposite_ip = geom(ph)%IPneighborhood(2,opposite_neighbor,en)
      opposite_n  = geom(ph)%IPneighborhood(3,opposite_neighbor,en)

      if (n_nbr > 0) then                                                                           ! if neighbor exists, average deformation gradient
        F_nbr = phase_mechanical_F(ph_nbr)%data(1:3,1:3,en_nbr)
        F_e_nbr = matmul(F_nbr, math_inv33(phase_mechanical_Fp(ph_nbr)%data(1:3,1:3,en_nbr)))
        Favg = 0.5_pREAL * (my_F + F_nbr)
      else                                                                                          ! if no neighbor, take my value as average
        Favg = my_F
      end if

      v_0_nbr = 0.0_pREAL        ! needed for check of sign change in flux density below

      !* FLUX FROM MY NEIGHBOR TO ME
      !* This is only considered, if I have a neighbor of nonlocal plasticity
      !* (also nonlocal constitutive law with local properties) that is at least a little bit
      !* compatible.
      !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of
      !* my neighbor's interface.
      !* The entering flux from my neighbor will be distributed on my slip systems according to the
      !* compatibility
      if (n_nbr > 0) then
      if (mechanical_plasticity_type(ph_nbr) == MECHANICAL_PLASTICITY_NONLOCAL .and. &
          any(dependentState(ph)%compatibility(:,:,:,n,en) > 0.0_pREAL)) then

        ! ToDo MD: Not sure if ns is correct here, but I think that compatibility is 0 if different phase
        v_0_nbr = reshape(state0(ph_nbr)%v(:,en_nbr),[ns,4])
        rho_0_sgl_mob_nbr = max(reshape(state0(ph_nbr)%rho_sgl_mob(:,en_nbr),[ns,4]),0.0_pREAL)

        where (rho_0_sgl_mob_nbr * geom(ph_nbr)%v_0(en_nbr) ** 0.667_pREAL < prm%rho_min &
          .or. rho_0_sgl_mob_nbr < prm%rho_significant) &
          rho_0_sgl_mob_nbr = 0.0_pREAL
        normal_neighbor2me_defConf = math_det33(Favg) * matmul(math_inv33(transpose(Favg)), &
                                     geom(ph_nbr)%n_0(1:3,n_nbr,en_nbr))                            ! normal of the interface in (average) deformed configuration (pointing neighbor => en)
        normal_neighbor2me = matmul(transpose(F_e_nbr), normal_neighbor2me_defConf) &
                           / math_det33(F_e_nbr)                                                    ! interface normal in the lattice configuration of my neighbor
        a = geom(ph_nbr)%a_0(n_nbr,en_nbr) * norm2(normal_neighbor2me)
        normal_neighbor2me = normal_neighbor2me / norm2(normal_neighbor2me)                         ! normalize the surface normal to unit length
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            topp = t + mod(t,2) - mod(t+1,2)
            if (v_0_nbr(s,t) * math_inner(m(1:3,s,t), normal_neighbor2me) > 0.0_pREAL &             ! flux from my neighbor to en == entering flux for en
                .and. v_0(s,t) * v_0_nbr(s,t) >= 0.0_pREAL ) then                                   ! ... only if no sign change in flux density
              lineLength = rho_0_sgl_mob_nbr(s,t) * v_0_nbr(s,t) &
                         * math_inner(m(1:3,s,t), normal_neighbor2me) * a                           ! positive line length that wants to enter through this interface
              where (dependentState(ph)%compatibility(c,:,s,n,en) > 0.0_pREAL) &
                rhoDotFlux(:,t) = rhoDotFlux(1:ns,t) &
                                + lineLength/geom(ph)%v_0(en)*dependentState(ph)%compatibility(c,:,s,n,en)**2        ! transferring to equally signed mobile dislocation type
              where (dependentState(ph)%compatibility(c,:,s,n,en) < 0.0_pREAL) &
                rhoDotFlux(:,topp) = rhoDotFlux(:,topp) &
                                   + lineLength/geom(ph)%v_0(en)*dependentState(ph)%compatibility(c,:,s,n,en)**2     ! transferring to opposite signed mobile dislocation type

            end if
          end do
        end do
      end if; end if


      !* FLUX FROM ME TO MY NEIGHBOR
      !* This is not considered, if my opposite neighbor has a different constitutive law than nonlocal (still considered for nonlocal law with local properties).
      !* Then, we assume, that the opposite(!) neighbor sends an equal amount of dislocations to en.
      !* So the net flux in the direction of my neighbor is equal to zero:
      !*    leaving flux to neighbor == entering flux from opposite neighbor
      !* In case of reduced transmissivity, part of the leaving flux is stored as dead dislocation density.
      !* That means for an interface of zero transmissivity the leaving flux is fully converted to dead dislocations.
      if (opposite_n > 0) then
      if (mechanical_plasticity_type(ph_nbr) == MECHANICAL_PLASTICITY_NONLOCAL) then

        normal_me2neighbor_defConf = math_det33(Favg) &
                                   * matmul(math_inv33(transpose(Favg)),geom(ph)%n_0(1:3,n,en))  ! normal of the interface in (average) deformed configuration (pointing en => neighbor)
        normal_me2neighbor = matmul(transpose(my_Fe), normal_me2neighbor_defConf) &
                           / math_det33(my_Fe)                                                      ! interface normal in my lattice configuration
        a = geom(ph)%a_0(n,en) * norm2(normal_me2neighbor)
        normal_me2neighbor = normal_me2neighbor / norm2(normal_me2neighbor)                         ! normalize the surface normal to unit length
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            if (v_0(s,t) * math_inner(m(1:3,s,t), normal_me2neighbor) > 0.0_pREAL ) then             ! flux from en to my neighbor == leaving flux for en (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
              if (v_0(s,t) * v_0_nbr(s,t) >= 0.0_pREAL) then                                         ! no sign change in flux density
                transmissivity = sum(dependentState(ph)%compatibility(c,:,s,n,en)**2)               ! overall transmissivity from this slip system to my neighbor
              else                                                                                  ! sign change in flux density means sign change in stress which does not allow for dislocations to arive at the neighbor
                transmissivity = 0.0_pREAL
              end if
              lineLength = rho_0_sgl_mob(s,t) * v_0(s,t) &
                         * math_inner(m(1:3,s,t), normal_me2neighbor) * a                           ! positive line length of mobiles that wants to leave through this interface
              rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / geom(ph)%v_0(en)                     ! subtract dislocation flux from current type
              rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) &
                                + lineLength / geom(ph)%v_0(en) * (1.0_pREAL - transmissivity) &
                                * sign(1.0_pREAL, v_0(s,t))                                          ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
            end if
          end do
        end do
      end if; end if

    end do neighbors
  end if

  end associate

end function rhoDotFlux


!--------------------------------------------------------------------------------------------------
!> @brief Compatibility update
!> @detail Compatibility is defined as normalized product of signed cosine of the angle between the slip
! plane normals and signed cosine of the angle between the slip directions. Only the largest values
! that sum up to a total of 1 are considered, all others are set to zero.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_updateCompatibility(orientation,ph,en)

  type(tRotationContainer), dimension(:), intent(in) :: &
    orientation                                                                                     ! crystal orientation
  integer, intent(in) :: &
    ph, en

  integer :: &
    n, &                                                                                            ! neighbor index
    el_nbr, &                                                                                   ! element index of my neighbor
    ip_nbr, &                                                                                   ! integration point index of my neighbor
    ce_nbr, &
    en_nbr, &
    ph_nbr, &
    ns, &                                                                                           ! number of active slip systems
    s1, &                                                                                           ! slip system index (en)
    s2                                                                                              ! slip system index (my neighbor)
  real(pREAL), dimension(2,param(ph)%sum_N_sl,param(ph)%sum_N_sl,nCellNeighbors) :: &
    my_compatibility                                                                                ! my_compatibility for current element and ip
  real(pREAL) :: &
    my_compatibilitySum, &
    thresholdValue, &
    nThresholdValues
  logical, dimension(param(ph)%sum_N_sl) :: &
    belowThreshold
  type(tRotation) :: mis


  associate(prm => param(ph))
    ns = prm%sum_N_sl

    !*** start out fully compatible
    my_compatibility = 0.0_pREAL
    forall(s1 = 1:ns) my_compatibility(:,s1,s1,:) = 1.0_pREAL

    neighbors: do n = 1,nCellNeighbors
      el_nbr = geom(ph)%IPneighborhood(1,n,en)
      ip_nbr = geom(ph)%IPneighborhood(2,n,en)
      ce_nbr = (el_nbr-1)*discretization_nIPs + ip_nbr
      en_nbr = material_entry_phase(1,ce_nbr)
      ph_nbr = material_ID_phase(1,ce_nbr)

      if (ce_nbr <= 0) then                                                                         ! free surface
        forall(s1 = 1:ns) my_compatibility(:,s1,s1,n) = sqrt(prm%chi_surface)
      elseif (ph_nbr /= ph) then                                                                    ! phase boundary
        if (plasticState(ph_nbr)%nonlocal .and. plasticState(ph)%nonlocal) &
          forall(s1 = 1:ns) my_compatibility(:,s1,s1,n) = 0.0_pREAL
      elseif (prm%chi_GB >= 0.0_pREAL) then                                                         ! grain boundary
        if (any(dNeq(phase_O_0(ph)%data(en)%asQuaternion(), &
                     phase_O_0(ph_nbr)%data(en_nbr)%asQuaternion())) .and. &
            plasticState(ph_nbr)%nonlocal) &
          forall(s1 = 1:ns) my_compatibility(:,s1,s1,n) = sqrt(prm%chi_GB)
      else
        !* GRAIN BOUNDARY ?
        !* Compatibility defined by relative orientation of slip systems:
        !* The my_compatibility value is defined as the product of the slip normal projection and the slip direction projection.
        !* Its sign is always positive for screws, for edges it has the same sign as the slip normal projection.
        !* Since the sum for each slip system can easily exceed one (which would result in a transmissivity larger than one),
        !* only values above or equal to a certain threshold value are considered. This threshold value is chosen, such that
        !* the number of compatible slip systems is minimized with the sum of the original compatibility values exceeding one.
        !* Finally the smallest compatibility value is decreased until the sum is exactly equal to one.
        !* All values below the threshold are set to zero.
        mis = orientation(ph)%data(en)%misorientation(orientation(ph_nbr)%data(en_nbr))
        mySlipSystems: do s1 = 1,ns
          neighborSlipSystems: do s2 = 1,ns
            my_compatibility(1,s2,s1,n) =     math_inner(prm%slip_normal(1:3,s1), &
                                                         mis%rotate(prm%slip_normal(1:3,s2))) &
                                        * abs(math_inner(prm%slip_direction(1:3,s1), &
                                                         mis%rotate(prm%slip_direction(1:3,s2))))
            my_compatibility(2,s2,s1,n) = abs(math_inner(prm%slip_normal(1:3,s1), &
                                                         mis%rotate(prm%slip_normal(1:3,s2)))) &
                                        * abs(math_inner(prm%slip_direction(1:3,s1), &
                                                         mis%rotate(prm%slip_direction(1:3,s2))))
          end do neighborSlipSystems

          my_compatibilitySum = 0.0_pREAL
          belowThreshold = .true.
          do while (my_compatibilitySum < 1.0_pREAL .and. any(belowThreshold))
            thresholdValue = maxval(my_compatibility(2,:,s1,n), belowThreshold)                     ! screws always positive
            nThresholdValues = real(count(my_compatibility(2,:,s1,n) >= thresholdValue),pREAL)
            where (my_compatibility(2,:,s1,n) >= thresholdValue) belowThreshold = .false.
            if (my_compatibilitySum + thresholdValue * nThresholdValues > 1.0_pREAL) &
              where (abs(my_compatibility(:,:,s1,n)) >= thresholdValue) &
                my_compatibility(:,:,s1,n) = sign((1.0_pREAL - my_compatibilitySum)/nThresholdValues,&
                                                   my_compatibility(:,:,s1,n))
            my_compatibilitySum = my_compatibilitySum + nThresholdValues * thresholdValue
          end do

          where(belowThreshold) my_compatibility(1,:,s1,n) = 0.0_pREAL
          where(belowThreshold) my_compatibility(2,:,s1,n) = 0.0_pREAL

        end do mySlipSystems
      end if

    end do neighbors

    dependentState(ph)%compatibility(:,:,:,:,en) = my_compatibility

  end associate

end subroutine plastic_nonlocal_updateCompatibility


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_result(ph,group)

  integer,         intent(in) :: ph
  character(len=*),intent(in) :: group

  integer :: ou

  associate(prm => param(ph),dst => dependentState(ph),stt=>state(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_u_ed_pos')
          call result_writeDataset(stt%rho_sgl_mob_edg_pos,group,trim(prm%output(ou)), &
                                   'positive mobile edge density','1/m²', prm%systems_sl)
        case('rho_b_ed_pos')
          call result_writeDataset(stt%rho_sgl_imm_edg_pos,group,trim(prm%output(ou)), &
                                   'positive immobile edge density','1/m²', prm%systems_sl)
        case('rho_u_ed_neg')
          call result_writeDataset(stt%rho_sgl_mob_edg_neg,group,trim(prm%output(ou)), &
                                   'negative mobile edge density','1/m²', prm%systems_sl)
        case('rho_b_ed_neg')
          call result_writeDataset(stt%rho_sgl_imm_edg_neg,group,trim(prm%output(ou)), &
                                   'negative immobile edge density','1/m²', prm%systems_sl)
        case('rho_d_ed')
          call result_writeDataset(stt%rho_dip_edg,group,trim(prm%output(ou)), &
                                   'edge dipole density','1/m²', prm%systems_sl)
        case('rho_u_sc_pos')
          call result_writeDataset(stt%rho_sgl_mob_scr_pos,group,trim(prm%output(ou)), &
                                   'positive mobile screw density','1/m²', prm%systems_sl)
        case('rho_b_sc_pos')
          call result_writeDataset(stt%rho_sgl_imm_scr_pos,group,trim(prm%output(ou)), &
                                   'positive immobile screw density','1/m²', prm%systems_sl)
        case('rho_u_sc_neg')
          call result_writeDataset(stt%rho_sgl_mob_scr_neg,group,trim(prm%output(ou)), &
                                   'negative mobile screw density','1/m²', prm%systems_sl)
        case('rho_b_sc_neg')
          call result_writeDataset(stt%rho_sgl_imm_scr_neg,group,trim(prm%output(ou)), &
                                   'negative immobile screw density','1/m²', prm%systems_sl)
        case('rho_d_sc')
          call result_writeDataset(stt%rho_dip_scr,group,trim(prm%output(ou)), &
                                   'screw dipole density','1/m²', prm%systems_sl)
        case('rho_f')
          call result_writeDataset(dst%rho_forest,group,trim(prm%output(ou)), &
                                   'forest density','1/m²', prm%systems_sl)
        case('v_ed_pos')
          call result_writeDataset(stt%v_edg_pos,group,trim(prm%output(ou)), &
                                   'positive edge velocity','m/s', prm%systems_sl)
        case('v_ed_neg')
          call result_writeDataset(stt%v_edg_neg,group,trim(prm%output(ou)), &
                                   'negative edge velocity','m/s', prm%systems_sl)
        case('v_sc_pos')
          call result_writeDataset(stt%v_scr_pos,group,trim(prm%output(ou)), &
                                   'positive srew velocity','m/s', prm%systems_sl)
        case('v_sc_neg')
          call result_writeDataset(stt%v_scr_neg,group,trim(prm%output(ou)), &
                                   'negative screw velocity','m/s', prm%systems_sl)
        case('gamma')
          call result_writeDataset(stt%gamma,group,trim(prm%output(ou)), &
                                   'plastic shear','1', prm%systems_sl)
        case('tau_pass')
          call result_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                   'passing stress for slip','Pa', prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_nonlocal_result


!--------------------------------------------------------------------------------------------------
!> @brief populates the initial dislocation density
!--------------------------------------------------------------------------------------------------
subroutine stateInit(ini,phase,Nentries)

  type(tInitialParameters) :: &
    ini
  integer,intent(in) :: &
    phase, &
    Nentries

  integer :: &
    e, &
    f, &
    from, &
    upto, &
    s
  real(pREAL), dimension(2) :: &
    rnd
  real(pREAL) :: &
    meanDensity, &
    totalVolume, &
    densityBinning, &
    minimumIpVolume


  associate(stt => state(phase))

    if (ini%random_rho_u > 0.0_pREAL) then ! randomly distribute dislocation segments on random slip system and of random type in the volume
      totalVolume     = sum(geom(phase)%v_0)
      minimumIPVolume = minval(geom(phase)%v_0)
      densityBinning  = ini%random_rho_u_binning / minimumIpVolume ** (2.0_pREAL / 3.0_pREAL)

      ! fill random material points with dislocation segments until the desired overall density is reached
      meanDensity = 0.0_pREAL
      do while(meanDensity < ini%random_rho_u)
        call random_number(rnd)
        e = nint(rnd(1)*real(Nentries,pREAL) + 0.5_pREAL)
        s = nint(rnd(2)*real(sum(ini%N_sl),pREAL)*4.0_pREAL + 0.5_pREAL)
        meanDensity = meanDensity + densityBinning * geom(phase)%v_0(e) / totalVolume
        stt%rho_sgl_mob(s,e) = densityBinning
      end do
    else                                ! homogeneous distribution with noise
      do f = 1,size(ini%N_sl,1)
        from = 1 + sum(ini%N_sl(1:f-1))
        upto = sum(ini%N_sl(1:f))
        call math_normal(stt%rho_sgl_mob_edg_pos(from:upto,:),ini%rho_u_ed_pos_0(f),ini%sigma_rho_u)
        call math_normal(stt%rho_sgl_mob_edg_neg(from:upto,:),ini%rho_u_ed_neg_0(f),ini%sigma_rho_u)
        call math_normal(stt%rho_sgl_mob_scr_pos(from:upto,:),ini%rho_u_sc_pos_0(f),ini%sigma_rho_u)
        call math_normal(stt%rho_sgl_mob_scr_neg(from:upto,:),ini%rho_u_sc_neg_0(f),ini%sigma_rho_u)
        stt%rho_dip_edg(from:upto,:) = ini%rho_d_ed_0(f)
        stt%rho_dip_scr(from:upto,:) = ini%rho_d_sc_0(f)
      end do
    end if

  end associate

end subroutine stateInit


!--------------------------------------------------------------------------------------------------
!> @brief calculates kinetics
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(v, dv_dtau, dv_dtauNS, tau, tauNS, tauThreshold, c, T, ph)

  integer, intent(in) :: &
    c, &                                                                                            !< dislocation character (1:edge, 2:screw)
    ph
  real(pREAL), dimension(param(ph)%sum_N_sl), intent(in) :: &
    tau, &                                                                                          !< resolved external shear stress (without non Schmid effects)
    tauNS, &                                                                                        !< resolved external shear stress (including non Schmid effects)
    tauThreshold                                                                                    !< threshold shear stress
  real(pREAL), intent(in) :: &
    T                                                                                               !< T
  real(pREAL), dimension(param(ph)%sum_N_sl), intent(out) ::  &
    v, &                                                                                            !< velocity
    dv_dtau, &                                                                                      !< velocity derivative with respect to resolved shear stress (without non Schmid contributions)
    dv_dtauNS                                                                                       !< velocity derivative with respect to resolved shear stress (including non Schmid contributions)

  integer :: &
    s                                                                                               !< index of my current slip system
  real(pREAL) :: &
    tauRel_P, &
    tauRel_S, &
    tauEff, &                                                                                       !< effective shear stress
    tPeierls, &                                                                                     !< waiting time in front of a peierls barriers
    tSolidSolution, &                                                                               !< waiting time in front of a solid solution obstacle
    dtPeierls_dtau, &                                                                               !< derivative with respect to resolved shear stress
    dtSolidSolution_dtau, &                                                                         !< derivative with respect to resolved shear stress
    lambda_S, &                                                                                     !< mean free distance between two solid solution obstacles
    lambda_P, &                                                                                     !< mean free distance between two Peierls barriers
    activationVolume_P, &                                                                           !< volume that needs to be activated to overcome barrier
    activationVolume_S, &                                                                           !< volume that needs to be activated to overcome barrier
    activationEnergy_P, &                                                                           !< energy that is needed to overcome barrier
    criticalStress_P, &                                                                             !< maximum obstacle strength
    criticalStress_S                                                                                !< maximum obstacle strength


  v = 0.0_pREAL
  dv_dtau = 0.0_pREAL
  dv_dtauNS = 0.0_pREAL

  associate(prm => param(ph))

    do s = 1,prm%sum_N_sl
      if (abs(tau(s)) > tauThreshold(s)) then

        !* Peierls contribution
        tauEff = max(0.0_pREAL, abs(tauNS(s)) - tauThreshold(s))
        lambda_P = prm%b_sl(s)
        activationVolume_P = prm%w * prm%b_sl(s)**3
        criticalStress_P = prm%peierlsStress(s,c)
        activationEnergy_P = criticalStress_P * activationVolume_P
        tauRel_P = min(1.0_pREAL, tauEff / criticalStress_P)
        tPeierls = 1.0_pREAL / prm%nu_a &
                 * exp(activationEnergy_P / (K_B * T) &
                       * (1.0_pREAL - tauRel_P**prm%p)**prm%q)
        dtPeierls_dtau = merge(tPeierls * prm%p * prm%q * activationVolume_P / (K_B * T) &
                               * (1.0_pREAL - tauRel_P**prm%p)**(prm%q-1.0_pREAL) * tauRel_P**(prm%p-1.0_pREAL), &
                               0.0_pREAL, &
                               tauEff < criticalStress_P)

        ! Contribution from solid solution strengthening
        tauEff = abs(tau(s)) - tauThreshold(s)
        lambda_S = prm%b_sl(s) / sqrt(prm%c_sol)
        activationVolume_S = prm%f_sol * prm%b_sl(s)**3 / sqrt(prm%c_sol)
        criticalStress_S = prm%Q_sol / activationVolume_S
        tauRel_S = min(1.0_pREAL, tauEff / criticalStress_S)
        tSolidSolution = 1.0_pREAL /  prm%nu_a &
                       * exp(prm%Q_sol / (K_B * T)* (1.0_pREAL - tauRel_S**prm%p)**prm%q)
        dtSolidSolution_dtau = merge(tSolidSolution * prm%p * prm%q * activationVolume_S / (K_B * T) &
                                     * (1.0_pREAL - tauRel_S**prm%p)**(prm%q-1.0_pREAL)* tauRel_S**(prm%p-1.0_pREAL), &
                                     0.0_pREAL, &
                                     tauEff < criticalStress_S)

        !* viscous glide velocity
        tauEff = abs(tau(s)) - tauThreshold(s)


        v(s) = sign(1.0_pREAL,tau(s)) &
             / (tPeierls / lambda_P + tSolidSolution / lambda_S + prm%B /(prm%b_sl(s) * tauEff))
        dv_dtau(s)   = v(s)**2 * (dtSolidSolution_dtau / lambda_S + prm%B / (prm%b_sl(s) * tauEff**2))
        dv_dtauNS(s) = v(s)**2 * dtPeierls_dtau / lambda_P

      end if
    end do

  end associate

end subroutine kinetics


!--------------------------------------------------------------------------------------------------
!> @brief returns copy of current dislocation densities from state
!> @details raw values is rectified
!--------------------------------------------------------------------------------------------------
pure function getRho(ph,en) result(rho)

  integer, intent(in) :: ph, en
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: rho


  associate(prm => param(ph))

    rho = reshape(state(ph)%rho(:,en),[prm%sum_N_sl,10])

    ! ensure positive densities (not for imm, they have a sign)
    rho(:,mob) = max(rho(:,mob),0.0_pREAL)
    rho(:,dip) = max(rho(:,dip),0.0_pREAL)

    where(abs(rho) < max(prm%rho_min/geom(ph)%v_0(en)**(2.0_pREAL/3.0_pREAL),prm%rho_significant)) &
      rho = 0.0_pREAL

  end associate

end function getRho


!--------------------------------------------------------------------------------------------------
!> @brief returns copy of current dislocation densities from state
!> @details raw values is rectified
!--------------------------------------------------------------------------------------------------
pure function getRho0(ph,en) result(rho_0)

  integer, intent(in) :: ph, en
  real(pREAL), dimension(param(ph)%sum_N_sl,10) :: rho_0


  associate(prm => param(ph))

    rho_0 = reshape(state0(ph)%rho(:,en),[prm%sum_N_sl,10])

    ! ensure positive densities (not for imm, they have a sign)
    rho_0(:,mob) = max(rho_0(:,mob),0.0_pREAL)
    rho_0(:,dip) = max(rho_0(:,dip),0.0_pREAL)

    where (abs(rho_0) < max(prm%rho_min/geom(ph)%v_0(en)**(2.0_pREAL/3.0_pREAL),prm%rho_significant)) &
      rho_0 = 0.0_pREAL

  end associate

end function getRho0


!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine storeGeometry(ph)

  integer, intent(in) :: ph

  integer :: ce, nCell
  real(pREAL), dimension(:), allocatable :: v_0
  real(pREAL), dimension(:,:), allocatable :: a_0, x_0
  real(pREAL), dimension(:,:,:), allocatable :: n_0
  integer, dimension(:,:,:), allocatable :: neighborhood


  nCell = product(shape(IPVolume0))

  v_0 = reshape(IPVolume0,[nCell])
  a_0 = reshape(IPArea0,[nCellNeighbors,nCell])
  x_0 = reshape(discretization_IPcoords,[3,nCell])
  n_0 = reshape(IPAreaNormal0,[3,nCellNeighbors,nCell])
  neighborhood = reshape(IPneighborhood,[3,nCellNeighbors,nCell])

  do ce = 1, size(material_entry_homogenization,1)
    if (material_ID_phase(1,ce) == ph) then
      geom(ph)%v_0(material_entry_phase(1,ce)) = v_0(ce)
      geom(ph)%a_0(:,material_entry_phase(1,ce)) = a_0(:,ce)
      geom(ph)%x_0(:,material_entry_phase(1,ce)) = x_0(:,ce)
      geom(ph)%n_0(:,:,material_entry_phase(1,ce)) = n_0(:,:,ce)
      geom(ph)%IPneighborhood(:,:,material_entry_phase(1,ce)) = neighborhood(:,:,ce)
    end if
  end do

end subroutine storeGeometry

end submodule nonlocal
