!--------------------------------------------------------------------------------------------------
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for plasticity including dislocation flux
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) nonlocal
  use geometry_plastic_nonlocal, only: &
    nIPneighbors    => geometry_plastic_nonlocal_nIPneighbors, &
    IPneighborhood  => geometry_plastic_nonlocal_IPneighborhood, &
    IPvolume        => geometry_plastic_nonlocal_IPvolume0, &
    IParea          => geometry_plastic_nonlocal_IParea0, &
    IPareaNormal    => geometry_plastic_nonlocal_IPareaNormal0, &
    geometry_plastic_nonlocal_disable

  type :: tGeometry
    real(pReal), dimension(:), allocatable :: V_0
  end type tGeometry

  type(tGeometry), dimension(:), allocatable :: geom

  real(pReal), parameter :: &
    kB = 1.38e-23_pReal                                                                             !< Boltzmann constant in J/Kelvin

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

  ! BEGIN DEPRECATED
  integer, dimension(:,:,:), allocatable :: &
    iRhoU, &                                                                                        !< state indices for unblocked density
    iV, &                                                                                           !< state indices for dislocation velocities
    iD                                                                                              !< state indices for stable dipole height
  !END DEPRECATED

  real(pReal), dimension(:,:,:,:,:,:), allocatable :: &
    compatibility                                                                                   !< slip system compatibility between en and my neighbors

  type :: tInitialParameters                                                                        !< container type for internal constitutive parameters
    real(pReal) :: &
      sigma_rho_u, &                                                                                !< standard deviation of scatter in initial dislocation density
      random_rho_u, &
      random_rho_u_binning
    real(pReal), dimension(:), allocatable :: &
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
    real(pReal) :: &
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
      eta, &                                                                                        !< viscosity for dislocation glide in Pa s
      nu_a, &                                                                                       !< attack frequency in Hz
      chi_surface, &                                                                                !< transmissivity at free surface
      chi_GB, &                                                                                     !< transmissivity at grain boundary (identified by different texture)
      f_c, &                                                                                        !< safety factor for CFL flux condition
      f_ed_mult, &                                                                                  !< factor that determines how much edge dislocations contribute to multiplication (0...1)
      f_F, &
      f_ed, &
      mu, &
      nu
    real(pReal), dimension(:),      allocatable :: &
      d_ed, &                                                                                       !< minimum stable edge dipole height
      d_sc, &                                                                                       !< minimum stable screw dipole height
      tau_Peierls_ed, &
      tau_Peierls_sc, &
      i_sl, &                                                                                       !< mean free path prefactor for each
      b_sl                                                                                          !< absolute length of Burgers vector [m]
    real(pReal), dimension(:,:),   allocatable :: &
      slip_normal, &
      slip_direction, &
      slip_transverse, &
      minDipoleHeight, &                                                                            ! edge and screw
      peierlsstress, &                                                                              ! edge and screw
      h_sl_sl ,&                                                                                    !< coefficients for slip-slip interaction
      forestProjection_Edge, &                                                                      !< matrix of forest projections of edge dislocations
      forestProjection_Screw                                                                        !< matrix of forest projections of screw dislocations
    real(pReal), dimension(:,:,:), allocatable :: &
      Schmid, &                                                                                     !< Schmid contribution
      nonSchmid_pos, &
      nonSchmid_neg                                                                                 !< combined projection of Schmid and non-Schmid contributions to the resolved shear stress (only for screws)
    integer :: &
      sum_N_sl
    integer,     dimension(:),     allocatable :: &
      colinearSystem                                                                                !< colinear system to the active slip system (only valid for fcc!)
    character(len=pStringLen), dimension(:), allocatable :: &
      output
    logical :: &
      shortRangeStressCorrection, &                                                                 !< use of short range stress correction by excess density gradient term
      nonSchmidActive = .false.
  end type tParameters

  type :: tNonlocalMicrostructure
    real(pReal), allocatable, dimension(:,:) :: &
     tau_pass, &
     tau_Back
  end type tNonlocalMicrostructure

  type :: tNonlocalState
    real(pReal), pointer, dimension(:,:) :: &
      rho, &                                                                                        ! < all dislocations
        rhoSgl, &
          rhoSglMobile, &                       ! iRhoU
            rho_sgl_mob_edg_pos, &
            rho_sgl_mob_edg_neg, &
            rho_sgl_mob_scr_pos, &
            rho_sgl_mob_scr_neg, &
          rhoSglImmobile, &
            rho_sgl_imm_edg_pos, &
            rho_sgl_imm_edg_neg, &
            rho_sgl_imm_scr_pos, &
            rho_sgl_imm_scr_neg, &
        rhoDip, &
          rho_dip_edg, &
          rho_dip_scr, &
        rho_forest, &
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

  type(tNonlocalMicrostructure), dimension(:), allocatable :: microstructure

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_nonlocal_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    Ninstances, &
    ph, &
    Nmembers, &
    sizeState, sizeDotState, sizeDependentState, sizeDeltaState, &
    s1, s2, &
    s, t, l
  real(pReal), dimension(:), allocatable :: &
    a
  character(len=pStringLen) :: &
    extmsg  = ''
  type(tInitialParameters) :: &
    ini
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl

  myPlasticity = plastic_active('nonlocal')
  Ninstances = count(myPlasticity)
  if(Ninstances == 0) then
    call geometry_plastic_nonlocal_disable
    return
  endif

  print'(/,a)', ' <<<+-  phase:mechanical:plastic:nonlocal init  -+>>>'
  print'(a,i0)', ' # phases: ',Ninstances; flush(IO_STDOUT)

  print*, 'C. Reuber et al., Acta Materialia 71:333–348, 2014'
  print*, 'https://doi.org/10.1016/j.actamat.2014.03.012'//IO_EOL

  print*, 'C. Kords, Dissertation RWTH Aachen, 2014'
  print*, 'http://publications.rwth-aachen.de/record/229993'


  phases => config_material%get('phase')
  allocate(geom(phases%length))

  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(state0(phases%length))
  allocate(dotState(phases%length))
  allocate(deltaState(phases%length))
  allocate(microstructure(phases%length))

  do ph = 1, phases%length
    if(.not. myPlasticity(ph)) cycle

    associate(prm => param(ph),  dot => dotState(ph),   stt => state(ph), &
              st0 => state0(ph), del => deltaState(ph), dst => microstructure(ph))

    phase => phases%get(ph)
    mech  => phase%get('mechanical')
    pl  => mech%get('plastic')

    phase_localPlasticity(ph) = .not. pl%contains('nonlocal')

#if defined (__GFORTRAN__)
    prm%output = output_as1dString(pl)
#else
    prm%output = pl%get_as1dString('output',defaultVal=emptyStringArray)
#endif

    prm%atol_rho   = pl%get_asFloat('atol_rho',defaultVal=1.0e4_pReal)

    ! This data is read in already in lattice
    prm%mu = lattice_mu(ph)
    prm%nu = lattice_nu(ph)

    ini%N_sl     = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(ini%N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%Schmid = lattice_SchmidMatrix_slip(ini%N_sl,phase%get_asString('lattice'),&
                                             phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      if(trim(phase%get_asString('lattice')) == 'cI') then
        a = pl%get_as1dFloat('a_nonSchmid',defaultVal = emptyRealArray)
        if(size(a) > 0) prm%nonSchmidActive = .true.
        prm%nonSchmid_pos = lattice_nonSchmidMatrix(ini%N_sl,a,+1)
        prm%nonSchmid_neg = lattice_nonSchmidMatrix(ini%N_sl,a,-1)
      else
        prm%nonSchmid_pos = prm%Schmid
        prm%nonSchmid_neg = prm%Schmid
      endif

      prm%h_sl_sl = lattice_interaction_SlipBySlip(ini%N_sl, &
                                                   pl%get_as1dFloat('h_sl_sl'), &
                                                   phase%get_asString('lattice'))

      prm%forestProjection_edge  = lattice_forestProjection_edge (ini%N_sl,phase%get_asString('lattice'),&
                                                                  phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%forestProjection_screw = lattice_forestProjection_screw(ini%N_sl,phase%get_asString('lattice'),&
                                                                  phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      prm%slip_direction  = lattice_slip_direction (ini%N_sl,phase%get_asString('lattice'),&
                                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_transverse = lattice_slip_transverse(ini%N_sl,phase%get_asString('lattice'),&
                                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_normal     = lattice_slip_normal    (ini%N_sl,phase%get_asString('lattice'),&
                                                    phase%get_asFloat('c/a',defaultVal=0.0_pReal))

      ! collinear systems (only for octahedral slip systems in fcc)
      allocate(prm%colinearSystem(prm%sum_N_sl), source = -1)
      do s1 = 1, prm%sum_N_sl
        do s2 = 1, prm%sum_N_sl
           if (all(dEq0 (math_cross(prm%slip_direction(1:3,s1),prm%slip_direction(1:3,s2)))) .and. &
               any(dNeq0(math_cross(prm%slip_normal   (1:3,s1),prm%slip_normal   (1:3,s2))))) &
             prm%colinearSystem(s1) = s2
        enddo
      enddo

      ini%rho_u_ed_pos_0  = pl%get_as1dFloat('rho_u_ed_pos_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_ed_neg_0  = pl%get_as1dFloat('rho_u_ed_neg_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_sc_pos_0  = pl%get_as1dFloat('rho_u_sc_pos_0',   requiredSize=size(ini%N_sl))
      ini%rho_u_sc_neg_0  = pl%get_as1dFloat('rho_u_sc_neg_0',   requiredSize=size(ini%N_sl))
      ini%rho_d_ed_0      = pl%get_as1dFloat('rho_d_ed_0',       requiredSize=size(ini%N_sl))
      ini%rho_d_sc_0      = pl%get_as1dFloat('rho_d_sc_0',       requiredSize=size(ini%N_sl))

      prm%i_sl            = pl%get_as1dFloat('i_sl',             requiredSize=size(ini%N_sl))
      prm%b_sl            = pl%get_as1dFloat('b_sl',             requiredSize=size(ini%N_sl))

      prm%i_sl = math_expand(prm%i_sl,ini%N_sl)
      prm%b_sl = math_expand(prm%b_sl,ini%N_sl)

      prm%d_ed            = pl%get_as1dFloat('d_ed',             requiredSize=size(ini%N_sl))
      prm%d_sc            = pl%get_as1dFloat('d_sc',             requiredSize=size(ini%N_sl))
      prm%d_ed            = math_expand(prm%d_ed,ini%N_sl)
      prm%d_sc            = math_expand(prm%d_sc,ini%N_sl)
      allocate(prm%minDipoleHeight(prm%sum_N_sl,2))
      prm%minDipoleHeight(:,1)  = prm%d_ed
      prm%minDipoleHeight(:,2)  = prm%d_sc

      prm%tau_Peierls_ed    = pl%get_as1dFloat('tau_Peierls_ed', requiredSize=size(ini%N_sl))
      prm%tau_Peierls_sc    = pl%get_as1dFloat('tau_Peierls_sc', requiredSize=size(ini%N_sl))
      prm%tau_Peierls_ed    = math_expand(prm%tau_Peierls_ed,ini%N_sl)
      prm%tau_Peierls_sc    = math_expand(prm%tau_Peierls_sc,ini%N_sl)
      allocate(prm%peierlsstress(prm%sum_N_sl,2))
      prm%peierlsstress(:,1)    = prm%tau_Peierls_ed
      prm%peierlsstress(:,2)    = prm%tau_Peierls_sc

      prm%rho_significant       = pl%get_asFloat('rho_significant')
      prm%rho_min               = pl%get_asFloat('rho_min', 0.0_pReal)
      prm%f_c                   = pl%get_asFloat('f_c',defaultVal=2.0_pReal)

      prm%V_at                  = pl%get_asFloat('V_at')
      prm%D_0                   = pl%get_asFloat('D_0')
      prm%Q_cl                  = pl%get_asFloat('Q_cl')
      prm%f_F                   = pl%get_asFloat('f_F')
      prm%f_ed                  = pl%get_asFloat('f_ed') !,'edgejogs'
      prm%w                     = pl%get_asFloat('w')
      prm%Q_sol                 = pl%get_asFloat('Q_sol')
      prm%f_sol                 = pl%get_asFloat('f_sol')
      prm%c_sol                 = pl%get_asFloat('c_sol')

      prm%p                     = pl%get_asFloat('p_sl')
      prm%q                     = pl%get_asFloat('q_sl')
      prm%eta                   = pl%get_asFloat('eta')
      prm%nu_a                  = pl%get_asFloat('nu_a')

      ! ToDo: discuss logic
      ini%sigma_rho_u           = pl%get_asFloat('sigma_rho_u')
      ini%random_rho_u          = pl%get_asFloat('random_rho_u',defaultVal= 0.0_pReal)
      if (pl%contains('random_rho_u')) &
        ini%random_rho_u_binning      = pl%get_asFloat('random_rho_u_binning',defaultVal=0.0_pReal) !ToDo: useful default?
     ! if (rhoSglRandom(instance) < 0.0_pReal) &
     ! if (rhoSglRandomBinning(instance) <= 0.0_pReal) &

      prm%chi_surface                 = pl%get_asFloat('chi_surface',defaultVal=1.0_pReal)
      prm%chi_GB                      = pl%get_asFloat('chi_GB',     defaultVal=-1.0_pReal)
      prm%f_ed_mult                   = pl%get_asFloat('f_ed_mult')
      prm%shortRangeStressCorrection  = pl%get_asBool('short_range_stress_correction', defaultVal = .false.)


!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (any(prm%b_sl             <  0.0_pReal)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%i_sl             <= 0.0_pReal)) extmsg = trim(extmsg)//' i_sl'

      if (any(ini%rho_u_ed_pos_0   <  0.0_pReal)) extmsg = trim(extmsg)//' rho_u_ed_pos_0'
      if (any(ini%rho_u_ed_neg_0   <  0.0_pReal)) extmsg = trim(extmsg)//' rho_u_ed_neg_0'
      if (any(ini%rho_u_sc_pos_0   <  0.0_pReal)) extmsg = trim(extmsg)//' rho_u_sc_pos_0'
      if (any(ini%rho_u_sc_neg_0   <  0.0_pReal)) extmsg = trim(extmsg)//' rho_u_sc_neg_0'
      if (any(ini%rho_d_ed_0       <  0.0_pReal)) extmsg = trim(extmsg)//' rho_d_ed_0'
      if (any(ini%rho_d_sc_0       <  0.0_pReal)) extmsg = trim(extmsg)//' rho_d_sc_0'

      if (any(prm%peierlsstress    <  0.0_pReal)) extmsg = trim(extmsg)//' tau_peierls'
      if (any(prm%minDipoleHeight  <  0.0_pReal)) extmsg = trim(extmsg)//' d_ed or d_sc'

      if (prm%eta                 <=  0.0_pReal)  extmsg = trim(extmsg)//' eta'
      if (prm%Q_cl                <=  0.0_pReal)  extmsg = trim(extmsg)//' Q_cl'
      if (prm%nu_a                <=  0.0_pReal)  extmsg = trim(extmsg)//' nu_a'
      if (prm%w                   <=  0.0_pReal)  extmsg = trim(extmsg)//' w'
      if (prm%D_0                 <   0.0_pReal)  extmsg = trim(extmsg)//' D_0'
      if (prm%V_at                <=  0.0_pReal)  extmsg = trim(extmsg)//' V_at'                   ! ToDo: in dislotungsten, the atomic volume is given as a factor

      if (prm%rho_min             <   0.0_pReal)  extmsg = trim(extmsg)//' rho_min'
      if (prm%rho_significant     <   0.0_pReal)  extmsg = trim(extmsg)//' rho_significant'
      if (prm%atol_rho            <   0.0_pReal)  extmsg = trim(extmsg)//' atol_rho'
      if (prm%f_c                 <   0.0_pReal)  extmsg = trim(extmsg)//' f_c'

      if (prm%p    <= 0.0_pReal .or. prm%p > 1.0_pReal) extmsg = trim(extmsg)//' p_sl'
      if (prm%q    <  1.0_pReal .or. prm%q > 2.0_pReal) extmsg = trim(extmsg)//' q_sl'

      if (prm%f_F  < 0.0_pReal  .or. prm%f_F  > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' f_F'
      if (prm%f_ed <  0.0_pReal .or. prm%f_ed > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' f_ed'

      if (prm%Q_sol               <=  0.0_pReal)  extmsg = trim(extmsg)//' Q_sol'
      if (prm%f_sol               <=  0.0_pReal)  extmsg = trim(extmsg)//' f_sol'
      if (prm%c_sol               <=  0.0_pReal)  extmsg = trim(extmsg)//' c_sol'

      if (prm%chi_GB              >   1.0_pReal)  extmsg = trim(extmsg)//' chi_GB'
      if (prm%chi_surface  < 0.0_pReal .or. prm%chi_surface  > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' chi_surface'

      if (prm%f_ed_mult    < 0.0_pReal .or. prm%f_ed_mult    > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' f_ed_mult'

    endif slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers  = count(material_phaseID == ph)
    sizeDotState = size([   'rhoSglEdgePosMobile   ','rhoSglEdgeNegMobile   ', &
                            'rhoSglScrewPosMobile  ','rhoSglScrewNegMobile  ', &
                            'rhoSglEdgePosImmobile ','rhoSglEdgeNegImmobile ', &
                            'rhoSglScrewPosImmobile','rhoSglScrewNegImmobile', &
                            'rhoDipEdge            ','rhoDipScrew           ', &
                            'gamma                 ' ]) * prm%sum_N_sl                              !< "basic" microstructural state variables that are independent from other state variables
    sizeDependentState = size([ 'rhoForest   ']) * prm%sum_N_sl                                     !< microstructural state variables that depend on other state variables
    sizeState          = sizeDotState + sizeDependentState &
                       + size([ 'velocityEdgePos     ','velocityEdgeNeg     ', &
                                'velocityScrewPos    ','velocityScrewNeg    ', &
                                'maxDipoleHeightEdge ','maxDipoleHeightScrew' ]) * prm%sum_N_sl     !< other dependent state variables that are not updated by microstructure
    sizeDeltaState            = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,sizeDeltaState)

    allocate(geom(ph)%V_0(Nmembers))
    call storeGeometry(ph)

    plasticState(ph)%nonlocal         = pl%get_asBool('nonlocal')
    if(plasticState(ph)%nonlocal .and. .not. allocated(IPneighborhood)) &
      call IO_error(212,ext_msg='IPneighborhood does not exist')


    plasticState(ph)%offsetDeltaState = 0                                                            ! ToDo: state structure does not follow convention

    st0%rho => plasticState(ph)%state0                             (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    stt%rho => plasticState(ph)%state                              (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    dot%rho => plasticState(ph)%dotState                           (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    del%rho => plasticState(ph)%deltaState                         (0*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
    plasticState(ph)%atol(1:10*prm%sum_N_sl) = prm%atol_rho

      stt%rhoSgl => plasticState(ph)%state                         (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
      dot%rhoSgl => plasticState(ph)%dotState                      (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
      del%rhoSgl => plasticState(ph)%deltaState                    (0*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)

        stt%rhoSglMobile => plasticState(ph)%state                 (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
        dot%rhoSglMobile => plasticState(ph)%dotState              (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)
        del%rhoSglMobile => plasticState(ph)%deltaState            (0*prm%sum_N_sl+1: 4*prm%sum_N_sl,:)

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

        stt%rhoSglImmobile => plasticState(ph)%state               (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
        dot%rhoSglImmobile => plasticState(ph)%dotState            (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)
        del%rhoSglImmobile => plasticState(ph)%deltaState          (4*prm%sum_N_sl+1: 8*prm%sum_N_sl,:)

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

      stt%rhoDip => plasticState(ph)%state                         (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
      dot%rhoDip => plasticState(ph)%dotState                      (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
      del%rhoDip => plasticState(ph)%deltaState                    (8*prm%sum_N_sl+1:10*prm%sum_N_sl,:)

        stt%rho_dip_edg => plasticState(ph)%state                  (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)
        dot%rho_dip_edg => plasticState(ph)%dotState               (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)
        del%rho_dip_edg => plasticState(ph)%deltaState             (8*prm%sum_N_sl+1: 9*prm%sum_N_sl,:)

        stt%rho_dip_scr => plasticState(ph)%state                  (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
        dot%rho_dip_scr => plasticState(ph)%dotState               (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)
        del%rho_dip_scr => plasticState(ph)%deltaState             (9*prm%sum_N_sl+1:10*prm%sum_N_sl,:)

    stt%gamma => plasticState(ph)%state                      (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    dot%gamma => plasticState(ph)%dotState                   (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    del%gamma => plasticState(ph)%deltaState                 (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)
    plasticState(ph)%atol(10*prm%sum_N_sl+1:11*prm%sum_N_sl )  = pl%get_asFloat('atol_gamma', defaultVal = 1.0e-2_pReal)
    if(any(plasticState(ph)%atol(10*prm%sum_N_sl+1:11*prm%sum_N_sl) < 0.0_pReal)) &
      extmsg = trim(extmsg)//' atol_gamma'
    plasticState(ph)%slipRate => plasticState(ph)%dotState    (10*prm%sum_N_sl + 1:11*prm%sum_N_sl,1:Nmembers)

    stt%rho_forest => plasticState(ph)%state                 (11*prm%sum_N_sl + 1:12*prm%sum_N_sl,1:Nmembers)
    stt%v          => plasticState(ph)%state                 (12*prm%sum_N_sl + 1:16*prm%sum_N_sl,1:Nmembers)
        stt%v_edg_pos  => plasticState(ph)%state             (12*prm%sum_N_sl + 1:13*prm%sum_N_sl,1:Nmembers)
        stt%v_edg_neg  => plasticState(ph)%state             (13*prm%sum_N_sl + 1:14*prm%sum_N_sl,1:Nmembers)
        stt%v_scr_pos  => plasticState(ph)%state             (14*prm%sum_N_sl + 1:15*prm%sum_N_sl,1:Nmembers)
        stt%v_scr_neg  => plasticState(ph)%state             (15*prm%sum_N_sl + 1:16*prm%sum_N_sl,1:Nmembers)

    allocate(dst%tau_pass(prm%sum_N_sl,Nmembers),source=0.0_pReal)
    allocate(dst%tau_back(prm%sum_N_sl,Nmembers),source=0.0_pReal)
    end associate

    if (Nmembers > 0) call stateInit(ini,ph,Nmembers)

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(nonlocal)')

  enddo

  allocate(compatibility(2,maxval(param%sum_N_sl),maxval(param%sum_N_sl),nIPneighbors,&
                         discretization_nIPs,discretization_Nelems), source=0.0_pReal)

! BEGIN DEPRECATED----------------------------------------------------------------------------------
  allocate(iRhoU(maxval(param%sum_N_sl),4,phases%length), source=0)
  allocate(iV(maxval(param%sum_N_sl),4,phases%length),    source=0)
  allocate(iD(maxval(param%sum_N_sl),2,phases%length),    source=0)

  do ph = 1, phases%length

    if(.not. myPlasticity(ph)) cycle

    phase => phases%get(ph)
    Nmembers = count(material_phaseID == ph)
    l = 0
    do t = 1,4
      do s = 1,param(ph)%sum_N_sl
        l = l + 1
        iRhoU(s,t,ph) = l
      enddo
    enddo
    l = l + (4+2+1+1)*param(ph)%sum_N_sl ! immobile(4), dipole(2), shear, forest
    do t = 1,4
      do s = 1,param(ph)%sum_N_sl
        l = l + 1
        iV(s,t,ph) = l
      enddo
    enddo
    do t = 1,2
      do s = 1,param(ph)%sum_N_sl
        l = l + 1
        iD(s,t,ph) = l
      enddo
    enddo
    if (iD(param(ph)%sum_N_sl,2,ph) /= plasticState(ph)%sizeState) &
      error stop 'state indices not properly set (nonlocal)'
  enddo

end function plastic_nonlocal_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates quantities characterizing the microstructure
!--------------------------------------------------------------------------------------------------
module subroutine nonlocal_dependentState(ph, en, ip, el)

  integer, intent(in) :: &
    ph, &
    en, &
    ip, &
    el

  integer :: &
    no, &                                                                                           !< neighbor offset
    neighbor_el, &                                                                                  ! element number of neighboring material point
    neighbor_ip, &                                                                                  ! integration point of neighboring material point
    c, &                                                                                            ! index of dilsocation character (edge, screw)
    s, &                                                                                            ! slip system index
    dir, &
    n
  real(pReal) :: &
    FVsize, &
    nRealNeighbors                                                                                  ! number of really existing neighbors
  integer, dimension(2) :: &
    neighbors
  real(pReal), dimension(2) :: &
    rhoExcessGradient, &
    rhoExcessGradient_over_rho, &
    rhoTotal
  real(pReal), dimension(3) :: &
    rhoExcessDifferences, &
    normal_latticeConf
  real(pReal), dimension(3,3) :: &
    invFe, &                                                                                        !< inverse of elastic deformation gradient
    invFp, &                                                                                        !< inverse of plastic deformation gradient
    connections, &
    invConnections
  real(pReal), dimension(3,nIPneighbors) :: &
    connection_latticeConf
  real(pReal), dimension(2,param(ph)%sum_N_sl) :: &
    rhoExcess
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    rho_edg_delta, &
    rho_scr_delta
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rho0, &
    rho_neighbor0
  real(pReal), dimension(param(ph)%sum_N_sl,param(ph)%sum_N_sl) :: &
    myInteractionMatrix                                                                             ! corrected slip interaction matrix
  real(pReal), dimension(param(ph)%sum_N_sl,nIPneighbors) :: &
    rho_edg_delta_neighbor, &
    rho_scr_delta_neighbor
  real(pReal), dimension(2,maxval(param%sum_N_sl),nIPneighbors) :: &
    neighbor_rhoExcess, &                                                                           ! excess density at neighboring material point
    neighbor_rhoTotal                                                                               ! total density at neighboring material point
  real(pReal), dimension(3,param(ph)%sum_N_sl,2) :: &
    m                                                                                               ! direction of dislocation motion

  associate(prm => param(ph),dst => microstructure(ph), stt => state(ph))

  rho = getRho(ph,en)

  stt%rho_forest(:,en) = matmul(prm%forestProjection_Edge, sum(abs(rho(:,edg)),2)) &
                       + matmul(prm%forestProjection_Screw,sum(abs(rho(:,scr)),2))


  ! coefficients are corrected for the line tension effect
  ! (see Kubin,Devincre,Hoc; 2008; Modeling dislocation storage rates and mean free paths in face-centered cubic crystals)
  if (any(phase_lattice(ph) == ['cI','cF'])) then
    myInteractionMatrix = prm%h_sl_sl &
                        * spread((  1.0_pReal - prm%f_F &
                                   + prm%f_F &
                                   * log(0.35_pReal * prm%b_sl * sqrt(max(stt%rho_forest(:,en),prm%rho_significant))) &
                                   / log(0.35_pReal * prm%b_sl * 1e6_pReal))** 2.0_pReal,2,prm%sum_N_sl)
  else
    myInteractionMatrix = prm%h_sl_sl
  endif

  dst%tau_pass(:,en) = prm%mu * prm%b_sl &
                     * sqrt(matmul(myInteractionMatrix,sum(abs(rho),2)))

!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

  !#################################################################################################
  ! ToDo: MD: this is most likely only correct for F_i = I
  !#################################################################################################

  rho0 = getRho0(ph,en)
  if (.not. phase_localPlasticity(ph) .and. prm%shortRangeStressCorrection) then
    invFp = math_inv33(phase_mechanical_Fp(ph)%data(1:3,1:3,en))
    invFe = math_inv33(phase_mechanical_Fe(ph)%data(1:3,1:3,en))

    rho_edg_delta = rho0(:,mob_edg_pos) - rho0(:,mob_edg_neg)
    rho_scr_delta = rho0(:,mob_scr_pos) - rho0(:,mob_scr_neg)

    rhoExcess(1,:) = rho_edg_delta
    rhoExcess(2,:) = rho_scr_delta

    FVsize = geom(ph)%V_0(en) ** (1.0_pReal/3.0_pReal)

    !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities

    nRealNeighbors = 0.0_pReal
    neighbor_rhoTotal = 0.0_pReal
    do n = 1,nIPneighbors
      neighbor_el = IPneighborhood(1,n,ip,el)
      neighbor_ip = IPneighborhood(2,n,ip,el)
      no = material_phasememberAt(1,neighbor_ip,neighbor_el)
      if (neighbor_el > 0 .and. neighbor_ip > 0) then
        if (material_phaseAt(1,neighbor_el) == ph) then

            nRealNeighbors = nRealNeighbors + 1.0_pReal
            rho_neighbor0 = getRho0(ph,no)

            rho_edg_delta_neighbor(:,n) = rho_neighbor0(:,mob_edg_pos) - rho_neighbor0(:,mob_edg_neg)
            rho_scr_delta_neighbor(:,n) = rho_neighbor0(:,mob_scr_pos) - rho_neighbor0(:,mob_scr_neg)

            neighbor_rhoTotal(1,:,n) = sum(abs(rho_neighbor0(:,edg)),2)
            neighbor_rhoTotal(2,:,n) = sum(abs(rho_neighbor0(:,scr)),2)

            connection_latticeConf(1:3,n) = matmul(invFe, discretization_IPcoords(1:3,neighbor_el+neighbor_ip-1) &
                                          - discretization_IPcoords(1:3,el+neighbor_ip-1))
            normal_latticeConf = matmul(transpose(invFp), IPareaNormal(1:3,n,ip,el))
            if (math_inner(normal_latticeConf,connection_latticeConf(1:3,n)) < 0.0_pReal) &         ! neighboring connection points in opposite direction to face normal: must be periodic image
              connection_latticeConf(1:3,n) = normal_latticeConf * IPvolume(ip,el)/IParea(n,ip,el)  ! instead take the surface normal scaled with the diameter of the cell
        else
          ! local neighbor or different lattice structure or different constitution instance -> use central values instead
          connection_latticeConf(1:3,n) = 0.0_pReal
          rho_edg_delta_neighbor(:,n) = rho_edg_delta
          rho_scr_delta_neighbor(:,n) = rho_scr_delta
        endif
      else
        ! free surface -> use central values instead
        connection_latticeConf(1:3,n) = 0.0_pReal
        rho_edg_delta_neighbor(:,n) = rho_edg_delta
        rho_scr_delta_neighbor(:,n) = rho_scr_delta
      endif
    enddo

    neighbor_rhoExcess(1,:,:) = rho_edg_delta_neighbor
    neighbor_rhoExcess(2,:,:) = rho_scr_delta_neighbor

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
          rhoExcessDifferences(dir) = neighbor_rhoExcess(c,s,neighbors(1)) &
                                    - neighbor_rhoExcess(c,s,neighbors(2))
        enddo
        invConnections = math_inv33(connections)
        if (all(dEq0(invConnections))) call IO_error(-1,ext_msg='back stress calculation: inversion error')

        rhoExcessGradient(c) = math_inner(m(1:3,s,c), matmul(invConnections,rhoExcessDifferences))
      enddo

        ! ... plus gradient from deads ...
      rhoExcessGradient(1) = rhoExcessGradient(1) + sum(rho(s,imm_edg)) / FVsize
      rhoExcessGradient(2) = rhoExcessGradient(2) + sum(rho(s,imm_scr)) / FVsize

        ! ... normalized with the total density ...
      rhoTotal(1) = (sum(abs(rho(s,edg))) + sum(neighbor_rhoTotal(1,s,:))) / (1.0_pReal + nRealNeighbors)
      rhoTotal(2) = (sum(abs(rho(s,scr))) + sum(neighbor_rhoTotal(2,s,:))) / (1.0_pReal + nRealNeighbors)

      rhoExcessGradient_over_rho = 0.0_pReal
      where(rhoTotal > 0.0_pReal) rhoExcessGradient_over_rho = rhoExcessGradient / rhoTotal

        ! ... gives the local stress correction when multiplied with a factor
      dst%tau_back(s,en) = - prm%mu * prm%b_sl(s) / (2.0_pReal * PI) &
                         * (  rhoExcessGradient_over_rho(1) / (1.0_pReal - prm%nu) &
                            + rhoExcessGradient_over_rho(2))
    enddo
  endif

#ifdef DEBUG
  if (debugConstitutive%extensive &
      .and. ((debugConstitutive%element == el .and. debugConstitutive%ip == ip)&
             .or. .not. debugConstitutive%selective)) then
    print'(/,a,i8,1x,i2,1x,i1,/)', '<< CONST >> nonlocal_microstructure at el ip ',el,ip
    print'(a,/,12x,12(e10.3,1x))', '<< CONST >> rhoForest', stt%rho_forest(:,en)
    print'(a,/,12x,12(f10.5,1x))', '<< CONST >> tauThreshold / MPa', dst%tau_pass(:,en)*1e-6
    print'(a,/,12x,12(f10.5,1x),/)', '<< CONST >> tauBack / MPa', dst%tau_back(:,en)*1e-6
  endif
#endif

 end associate

end subroutine nonlocal_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
module subroutine nonlocal_LpAndItsTangent(Lp,dLp_dMp, &
                                                   Mp,Temperature,ph,en)
  real(pReal), dimension(3,3), intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp
  integer, intent(in) :: &
    ph, &
    en
  real(pReal), intent(in) :: &
    Temperature                                                                                     !< temperature

  real(pReal), dimension(3,3), intent(in) :: &
    Mp
                                                                                                    !< derivative of Lp with respect to Mp
  integer :: &
    ns, &                                                                                           !< short notation for the total number of active slip systems
    i, &
    j, &
    k, &
    l, &
    t, &                                                                                            !< dislocation type
    s                                                                                               !< index of my current slip system
  real(pReal), dimension(param(ph)%sum_N_sl,8) :: &
    rhoSgl                                                                                          !< single dislocation densities (including blocked)
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    rho
  real(pReal), dimension(param(ph)%sum_N_sl,4) :: &
    v, &                                                                                            !< velocity
    tauNS, &                                                                                        !< resolved shear stress including non Schmid and backstress terms
    dv_dtau, &                                                                                      !< velocity derivative with respect to the shear stress
    dv_dtauNS                                                                                       !< velocity derivative with respect to the shear stress
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau, &                                                                                          !< resolved shear stress including backstress terms
    gdotTotal                                                                                       !< shear rate

  associate(prm => param(ph),dst=>microstructure(ph),&
  stt=>state(ph))
  ns = prm%sum_N_sl

  !*** shortcut to state variables
  rho = getRho(ph,en)
  rhoSgl = rho(:,sgl)

  do s = 1,ns
    tau(s) = math_tensordot(Mp, prm%Schmid(1:3,1:3,s))
    tauNS(s,1) = tau(s)
    tauNS(s,2) = tau(s)
    if (tau(s) > 0.0_pReal) then
      tauNS(s,3) = math_tensordot(Mp, +prm%nonSchmid_pos(1:3,1:3,s))
      tauNS(s,4) = math_tensordot(Mp, -prm%nonSchmid_neg(1:3,1:3,s))
    else
      tauNS(s,3) = math_tensordot(Mp, +prm%nonSchmid_neg(1:3,1:3,s))
      tauNS(s,4) = math_tensordot(Mp, -prm%nonSchmid_pos(1:3,1:3,s))
    endif
  enddo
  tauNS = tauNS + spread(dst%tau_back(:,en),2,4)
  tau   = tau   + dst%tau_back(:,en)

  ! edges
  call kinetics(v(:,1), dv_dtau(:,1), dv_dtauNS(:,1), &
                tau, tauNS(:,1), dst%tau_pass(:,en),1,Temperature, ph)
  v(:,2)         = v(:,1)
  dv_dtau(:,2)   = dv_dtau(:,1)
  dv_dtauNS(:,2) = dv_dtauNS(:,1)

  !screws
  if (prm%nonSchmidActive) then
    do t = 3,4
      call kinetics(v(:,t), dv_dtau(:,t), dv_dtauNS(:,t), &
                    tau, tauNS(:,t), dst%tau_pass(:,en),2,Temperature, ph)
    enddo
  else
    v(:,3:4)         = spread(v(:,1),2,2)
    dv_dtau(:,3:4)   = spread(dv_dtau(:,1),2,2)
    dv_dtauNS(:,3:4) = spread(dv_dtauNS(:,1),2,2)
  endif

  stt%v(:,en) = pack(v,.true.)

  !*** Bauschinger effect
  forall (s = 1:ns, t = 5:8, rhoSgl(s,t) * v(s,t-4) < 0.0_pReal) &
    rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(rhoSgl(s,t))

  gdotTotal = sum(rhoSgl(:,1:4) * v, 2) * prm%b_sl

  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal
  do s = 1,ns
    Lp = Lp + gdotTotal(s) * prm%Schmid(1:3,1:3,s)
    forall (i=1:3,j=1:3,k=1:3,l=1:3) &
      dLp_dMp(i,j,k,l) = dLp_dMp(i,j,k,l) &
          + prm%Schmid(i,j,s) * prm%Schmid(k,l,s) &
          * sum(rhoSgl(s,1:4) * dv_dtau(s,1:4)) * prm%b_sl(s) &
          + prm%Schmid(i,j,s) &
          * (+ prm%nonSchmid_pos(k,l,s) * rhoSgl(s,3) * dv_dtauNS(s,3) &
             - prm%nonSchmid_neg(k,l,s) * rhoSgl(s,4) * dv_dtauNS(s,4))  * prm%b_sl(s)
  enddo

  end associate

end subroutine nonlocal_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_deltaState(Mp,ph,en)

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< MandelStress
  integer, intent(in) :: &
    ph, &
    en

  integer :: &
    ns, &                                                                                           ! short notation for the total number of active slip systems
    c, &                                                                                            ! character of dislocation
    t, &                                                                                            ! type of dislocation
    s                                                                                               ! index of my current slip system
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    deltaRhoRemobilization, &                                                                       ! density increment by remobilization
    deltaRhoDipole2SingleStress                                                                     ! density increment by dipole dissociation (by stress change)
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    rho                                                                                             ! current  dislocation densities
  real(pReal), dimension(param(ph)%sum_N_sl,4) :: &
    v                                                                                               ! dislocation glide velocity
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau                                                                                             ! current resolved shear stress
  real(pReal), dimension(param(ph)%sum_N_sl,2) :: &
    rhoDip, &                                                                                       ! current dipole dislocation densities (screw and edge dipoles)
    dUpper, &                                                                                       ! current maximum stable dipole distance for edges and screws
    dUpperOld, &                                                                                    ! old maximum stable dipole distance for edges and screws
    deltaDUpper                                                                                     ! change in maximum stable dipole distance for edges and screws

  associate(prm => param(ph),dst => microstructure(ph),del => deltaState(ph))
  ns = prm%sum_N_sl

  !*** shortcut to state variables
  forall (s = 1:ns, t = 1:4) v(s,t) = plasticState(ph)%state(iV(s,t,ph),en)
  forall (s = 1:ns, c = 1:2) dUpperOld(s,c) = plasticState(ph)%state(iD(s,c,ph),en)

  rho =  getRho(ph,en)
  rhoDip = rho(:,dip)

  !****************************************************************************
  !*** dislocation remobilization (bauschinger effect)
  where(rho(:,imm) * v < 0.0_pReal)
    deltaRhoRemobilization(:,mob) = abs(rho(:,imm))
    deltaRhoRemobilization(:,imm) =   - rho(:,imm)
    rho(:,mob) = rho(:,mob) + abs(rho(:,imm))
    rho(:,imm) = 0.0_pReal
  elsewhere
    deltaRhoRemobilization(:,mob) = 0.0_pReal
    deltaRhoRemobilization(:,imm) = 0.0_pReal
  endwhere
  deltaRhoRemobilization(:,dip) = 0.0_pReal

  !****************************************************************************
  !*** calculate dipole formation and dissociation by stress change

  !*** calculate limits for stable dipole height
  do s = 1,prm%sum_N_sl
    tau(s) = math_tensordot(Mp, prm%Schmid(1:3,1:3,s)) +dst%tau_back(s,en)
    if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
  enddo

  dUpper(:,1) = prm%mu * prm%b_sl/(8.0_pReal * PI * (1.0_pReal - prm%nu) * abs(tau))
  dUpper(:,2) = prm%mu * prm%b_sl/(4.0_pReal * PI * abs(tau))

  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(:,1) = min(1.0_pReal/sqrt(sum(abs(rho(:,edg)),2)),dUpper(:,1))
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(:,2) = min(1.0_pReal/sqrt(sum(abs(rho(:,scr)),2)),dUpper(:,2))

  dUpper = max(dUpper,prm%minDipoleHeight)
  deltaDUpper = dUpper - dUpperOld


  !*** dissociation by stress increase
  deltaRhoDipole2SingleStress = 0.0_pReal
  forall (c=1:2, s=1:ns, deltaDUpper(s,c) < 0.0_pReal .and. &
                                          dNeq0(dUpperOld(s,c) - prm%minDipoleHeight(s,c))) &
    deltaRhoDipole2SingleStress(s,8+c) = rhoDip(s,c) * deltaDUpper(s,c) &
                                       / (dUpperOld(s,c) - prm%minDipoleHeight(s,c))

  forall (t=1:4) deltaRhoDipole2SingleStress(:,t) = -0.5_pReal * deltaRhoDipole2SingleStress(:,(t-1)/2+9)
  forall (s = 1:ns, c = 1:2) plasticState(ph)%state(iD(s,c,ph),en) = dUpper(s,c)

  plasticState(ph)%deltaState(:,en) = 0.0_pReal
  del%rho(:,en) = reshape(deltaRhoRemobilization + deltaRhoDipole2SingleStress, [10*ns])

  end associate

end subroutine plastic_nonlocal_deltaState


!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
module subroutine nonlocal_dotState(Mp, Temperature,timestep, &
                                            ph,en,ip,el)

  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< MandelStress
  real(pReal), intent(in) :: &
    Temperature, &                                                                                  !< temperature
    timestep                                                                                        !< substepped crystallite time increment
  integer, intent(in) :: &
    ph, &
    en, &
    ip, &                                                                                           !< current integration point
    el                                                                                              !< current element number

  integer ::  &
    ns, &                                                                                           !< short notation for the total number of active slip systems
    c, &                                                                                            !< character of dislocation
    t, &                                                                                            !< type of dislocation
    s                                                                                               !< index of my current slip system
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rho0, &                                                                                         !< dislocation density at beginning of time step
    rhoDot, &                                                                                       !< density evolution
    rhoDotMultiplication, &                                                                         !< density evolution by multiplication
    rhoDotSingle2DipoleGlide, &                                                                     !< density evolution by dipole formation (by glide)
    rhoDotAthermalAnnihilation, &                                                                   !< density evolution by athermal annihilation
    rhoDotThermalAnnihilation                                                                       !< density evolution by thermal annihilation
  real(pReal), dimension(param(ph)%sum_N_sl,8) :: &
    rhoSgl, &                                                                                       !< current single dislocation densities (positive/negative screw and edge without dipoles)
    my_rhoSgl0                                                                                      !< single dislocation densities of central ip (positive/negative screw and edge without dipoles)
  real(pReal), dimension(param(ph)%sum_N_sl,4) :: &
    v, &                                                                                            !< current dislocation glide velocity
    v0, &
    gdot                                                                                            !< shear rates
  real(pReal), dimension(param(ph)%sum_N_sl) :: &
    tau, &                                                                                          !< current resolved shear stress
    vClimb                                                                                          !< climb velocity of edge dipoles
  real(pReal), dimension(param(ph)%sum_N_sl,2) :: &
    rhoDip, &                                                                                       !< current dipole dislocation densities (screw and edge dipoles)
    dLower, &                                                                                       !< minimum stable dipole distance for edges and screws
    dUpper                                                                                          !< current maximum stable dipole distance for edges and screws
  real(pReal) :: &
    selfDiffusion                                                                                   !< self diffusion

  if (timestep <= 0.0_pReal) then
    plasticState(ph)%dotState = 0.0_pReal
    return
  endif

  associate(prm => param(ph), &
            dst => microstructure(ph), &
            dot => dotState(ph), &
            stt => state(ph))
  ns = prm%sum_N_sl

  tau = 0.0_pReal
  gdot = 0.0_pReal

  rho  = getRho(ph,en)
  rhoSgl = rho(:,sgl)
  rhoDip = rho(:,dip)
  rho0 = getRho0(ph,en)
  my_rhoSgl0 = rho0(:,sgl)

  forall (s = 1:ns, t = 1:4) v(s,t) = plasticState(ph)%state(iV(s,t,ph),en)
  gdot = rhoSgl(:,1:4) * v * spread(prm%b_sl,2,4)

#ifdef DEBUG
    if (debugConstitutive%basic &
        .and. ((debugConstitutive%element == el .and. debugConstitutive%ip == ip) &
               .or. .not. debugConstitutive%selective)) then
      print'(a,/,10(12x,12(e12.5,1x),/))', '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
      print'(a,/,4(12x,12(e12.5,1x),/))', '<< CONST >> gdot / 1/s',gdot
    endif
#endif

  !****************************************************************************
  !*** limits for stable dipole height
  do s = 1,ns
    tau(s) = math_tensordot(Mp, prm%Schmid(1:3,1:3,s)) + dst%tau_back(s,en)
    if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
  enddo

  dLower = prm%minDipoleHeight
  dUpper(:,1) = prm%mu * prm%b_sl/(8.0_pReal * PI * (1.0_pReal - prm%nu) * abs(tau))
  dUpper(:,2) = prm%mu * prm%b_sl/(4.0_pReal * PI * abs(tau))

  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(:,1) = min(1.0_pReal/sqrt(sum(abs(rho(:,edg)),2)),dUpper(:,1))
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(:,2) = min(1.0_pReal/sqrt(sum(abs(rho(:,scr)),2)),dUpper(:,2))

  dUpper = max(dUpper,dLower)

  !****************************************************************************
  !*** dislocation multiplication
  rhoDotMultiplication = 0.0_pReal
  isBCC: if (phase_lattice(ph) == 'cI') then
    forall (s = 1:ns, sum(abs(v(s,1:4))) > 0.0_pReal)
      rhoDotMultiplication(s,1:2) = sum(abs(gdot(s,3:4))) / prm%b_sl(s) &                           ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(stt%rho_forest(s,en)) / prm%i_sl(s) ! &                    ! mean free path
                                  ! * 2.0_pReal * sum(abs(v(s,3:4))) / sum(abs(v(s,1:4)))           ! ratio of screw to overall velocity determines edge generation
      rhoDotMultiplication(s,3:4) = sum(abs(gdot(s,3:4))) /prm%b_sl(s) &                            ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(stt%rho_forest(s,en)) / prm%i_sl(s) ! &                    ! mean free path
                                  ! * 2.0_pReal * sum(abs(v(s,1:2))) / sum(abs(v(s,1:4)))           ! ratio of edge to overall velocity determines screw generation
    endforall

  else isBCC
    rhoDotMultiplication(:,1:4) = spread( &
          (sum(abs(gdot(:,1:2)),2) * prm%f_ed_mult + sum(abs(gdot(:,3:4)),2)) &
        * sqrt(stt%rho_forest(:,en)) / prm%i_sl / prm%b_sl, 2, 4)
  endif isBCC

  forall (s = 1:ns, t = 1:4) v0(s,t) = plasticState(ph)%state0(iV(s,t,ph),en)


  !****************************************************************************
  !*** calculate dipole formation and annihilation

  !*** formation by glide
  do c = 1,2
    rhoDotSingle2DipoleGlide(:,2*c-1) = -2.0_pReal * dUpper(:,c) / prm%b_sl &
                                                   * (      rhoSgl(:,2*c-1)  * abs(gdot(:,2*c)) &   ! negative mobile --> positive mobile
                                                      +     rhoSgl(:,2*c)    * abs(gdot(:,2*c-1)) & ! positive mobile --> negative mobile
                                                      + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1)))  ! positive mobile --> negative immobile

    rhoDotSingle2DipoleGlide(:,2*c) = -2.0_pReal * dUpper(:,c) / prm%b_sl &
                                                 * (      rhoSgl(:,2*c-1)  * abs(gdot(:,2*c)) &     ! negative mobile --> positive mobile
                                                    +     rhoSgl(:,2*c)    * abs(gdot(:,2*c-1)) &   ! positive mobile --> negative mobile
                                                    + abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)))      ! negative mobile --> positive immobile

    rhoDotSingle2DipoleGlide(:,2*c+3) = -2.0_pReal * dUpper(:,c) / prm%b_sl &
                                                   * rhoSgl(:,2*c+3) * abs(gdot(:,2*c))             ! negative mobile --> positive immobile

    rhoDotSingle2DipoleGlide(:,2*c+4) = -2.0_pReal * dUpper(:,c) / prm%b_sl &
                                                   * rhoSgl(:,2*c+4) * abs(gdot(:,2*c-1))           ! positive mobile --> negative immobile

    rhoDotSingle2DipoleGlide(:,c+8) = abs(rhoDotSingle2DipoleGlide(:,2*c+3)) &
                                    + abs(rhoDotSingle2DipoleGlide(:,2*c+4)) &
                                    - rhoDotSingle2DipoleGlide(:,2*c-1) &
                                    - rhoDotSingle2DipoleGlide(:,2*c)
  enddo


  !*** athermal annihilation
  rhoDotAthermalAnnihilation = 0.0_pReal
  forall (c=1:2) &
    rhoDotAthermalAnnihilation(:,c+8) = -2.0_pReal * dLower(:,c) / prm%b_sl &
       * (  2.0_pReal * (rhoSgl(:,2*c-1) * abs(gdot(:,2*c)) + rhoSgl(:,2*c) * abs(gdot(:,2*c-1))) & ! was single hitting single
          + 2.0_pReal * (abs(rhoSgl(:,2*c+3)) * abs(gdot(:,2*c)) + abs(rhoSgl(:,2*c+4)) * abs(gdot(:,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
          + rhoDip(:,c) * (abs(gdot(:,2*c-1)) + abs(gdot(:,2*c))))                                  ! single knocks dipole constituent

  ! annihilated screw dipoles leave edge jogs behind on the colinear system
  if (phase_lattice(ph) == 'cF') &
    forall (s = 1:ns, prm%colinearSystem(s) > 0) &
      rhoDotAthermalAnnihilation(prm%colinearSystem(s),1:2) = - rhoDotAthermalAnnihilation(s,10) &
        * 0.25_pReal * sqrt(stt%rho_forest(s,en)) * (dUpper(s,2) + dLower(s,2)) * prm%f_ed


  !*** thermally activated annihilation of edge dipoles by climb
  rhoDotThermalAnnihilation = 0.0_pReal
  selfDiffusion = prm%D_0 * exp(-prm%Q_cl / (kB * Temperature))
  vClimb = prm%V_at * selfDiffusion * prm%mu &
         / ( kB * Temperature  * PI * (1.0_pReal-prm%nu) * (dUpper(:,1) + dLower(:,1)))
  forall (s = 1:ns, dUpper(s,1) > dLower(s,1)) &
    rhoDotThermalAnnihilation(s,9) = max(- 4.0_pReal * rhoDip(s,1) * vClimb(s) / (dUpper(s,1) - dLower(s,1)), &
                                         - rhoDip(s,1) / timestep - rhoDotAthermalAnnihilation(s,9) &
                                                                  - rhoDotSingle2DipoleGlide(s,9))  ! make sure that we do not annihilate more dipoles than we have

  rhoDot = rhoDotFlux(timestep, ph,en,ip,el) &
         + rhoDotMultiplication &
         + rhoDotSingle2DipoleGlide &
         + rhoDotAthermalAnnihilation &
         + rhoDotThermalAnnihilation


  if (    any(rho(:,mob) + rhoDot(:,1:4)  * timestep < -prm%atol_rho) &
     .or. any(rho(:,dip) + rhoDot(:,9:10) * timestep < -prm%atol_rho)) then
#ifdef DEBUG
    if (debugConstitutive%extensive) then
      print'(a,i5,a,i2)', '<< CONST >> evolution rate leads to negative density at el ',el,' ip ',ip
      print'(a)', '<< CONST >> enforcing cutback !!!'
    endif
#endif
    plasticState(ph)%dotState = IEEE_value(1.0_pReal,IEEE_quiet_NaN)
  else
    dot%rho(:,en) = pack(rhoDot,.true.)
    dot%gamma(:,en) = sum(gdot,2)
  endif

  end associate

end subroutine nonlocal_dotState


!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
function rhoDotFlux(timestep,ph,en,ip,el)

  real(pReal), intent(in) :: &
    timestep                                                                                        !< substepped crystallite time increment
  integer, intent(in) :: &
    ph, &
    en, &
    ip, &                                                                                           !< current integration point
    el                                                                                              !< current element number

  integer ::  &
    neighbor_ph, &                                                                                  !< phase of my neighbor's plasticity
    ns, &                                                                                           !< short notation for the total number of active slip systems
    c, &                                                                                            !< character of dislocation
    n, &                                                                                            !< index of my current neighbor
    neighbor_el, &                                                                                  !< element number of my neighbor
    neighbor_ip, &                                                                                  !< integration point of my neighbor
    neighbor_n, &                                                                                   !< neighbor index pointing to en when looking from my neighbor
    opposite_neighbor, &                                                                            !< index of my opposite neighbor
    opposite_ip, &                                                                                  !< ip of my opposite neighbor
    opposite_el, &                                                                                  !< element index of my opposite neighbor
    opposite_n, &                                                                                   !< neighbor index pointing to en when looking from my opposite neighbor
    t, &                                                                                            !< type of dislocation
    no,&                                                                                            !< neighbor offset shortcut
    np,&                                                                                            !< neighbor phase shortcut
    topp, &                                                                                         !< type of dislocation with opposite sign to t
    s                                                                                               !< index of my current slip system
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: &
    rho, &
    rho0, &                                                                                         !< dislocation density at beginning of time step
    rhoDotFlux                                                                                      !< density evolution by flux
  real(pReal), dimension(param(ph)%sum_N_sl,8) :: &
    rhoSgl, &                                                                                       !< current single dislocation densities (positive/negative screw and edge without dipoles)
    neighbor_rhoSgl0, &                                                                             !< current single dislocation densities of neighboring ip (positive/negative screw and edge without dipoles)
    my_rhoSgl0                                                                                      !< single dislocation densities of central ip (positive/negative screw and edge without dipoles)
  real(pReal), dimension(param(ph)%sum_N_sl,4) :: &
    v, &                                                                                            !< current dislocation glide velocity
    v0, &
    neighbor_v0, &                                                                                  !< dislocation glide velocity of enighboring ip
    gdot                                                                                            !< shear rates
  real(pReal), dimension(3,param(ph)%sum_N_sl,4) :: &
    m                                                                                               !< direction of dislocation motion
  real(pReal), dimension(3,3) :: &
    my_F, &                                                                                         !< my total deformation gradient
    neighbor_F, &                                                                                   !< total deformation gradient of my neighbor
    my_Fe, &                                                                                        !< my elastic deformation gradient
    neighbor_Fe, &                                                                                  !< elastic deformation gradient of my neighbor
    Favg                                                                                            !< average total deformation gradient of en and my neighbor
  real(pReal), dimension(3) :: &
    normal_neighbor2me, &                                                                           !< interface normal pointing from my neighbor to en in neighbor's lattice configuration
    normal_neighbor2me_defConf, &                                                                   !< interface normal pointing from my neighbor to en in shared deformed configuration
    normal_me2neighbor, &                                                                           !< interface normal pointing from en to my neighbor in my lattice configuration
    normal_me2neighbor_defConf                                                                      !< interface normal pointing from en to my neighbor in shared deformed configuration
  real(pReal) :: &
    area, &                                                                                         !< area of the current interface
    transmissivity, &                                                                               !< overall transmissivity of dislocation flux to neighboring material point
    lineLength                                                                                      !< dislocation line length leaving the current interface


  associate(prm => param(ph), &
            dst => microstructure(ph), &
            dot => dotState(ph), &
            stt => state(ph))
  ns = prm%sum_N_sl

  gdot = 0.0_pReal

  rho  = getRho(ph,en)
  rhoSgl = rho(:,sgl)
  rho0 = getRho0(ph,en)
  my_rhoSgl0 = rho0(:,sgl)

  forall (s = 1:ns, t = 1:4) v(s,t) = plasticState(ph)%state(iV(s,t,ph),en)                   !ToDo: MD: I think we should use state0 here
  gdot = rhoSgl(:,1:4) * v * spread(prm%b_sl,2,4)


  forall (s = 1:ns, t = 1:4) v0(s,t) = plasticState(ph)%state0(iV(s,t,ph),en)

  !****************************************************************************
  !*** calculate dislocation fluxes (only for nonlocal plasticity)
  rhoDotFlux = 0.0_pReal
  if (.not. phase_localPlasticity(ph)) then

    !*** check CFL (Courant-Friedrichs-Lewy) condition for flux
    if (any( abs(gdot) > 0.0_pReal &                                                                ! any active slip system ...
            .and. prm%f_c * abs(v0) * timestep &
                > IPvolume(ip,el) / maxval(IParea(:,ip,el)))) then                                  ! ...with velocity above critical value (we use the reference volume and area for simplicity here)
#ifdef DEBUG
    if (debugConstitutive%extensive) then
      print'(a,i5,a,i2)', '<< CONST >> CFL condition not fullfilled at el ',el,' ip ',ip
      print'(a,e10.3,a,e10.3)', '<< CONST >> velocity is at  ', &
        maxval(abs(v0), abs(gdot) > 0.0_pReal &
                       .and.  prm%f_c * abs(v0) * timestep &
                             > IPvolume(ip,el) / maxval(IParea(:,ip,el))), &
        ' at a timestep of ',timestep
      print*, '<< CONST >> enforcing cutback !!!'
    endif
#endif
      rhoDotFlux = IEEE_value(1.0_pReal,IEEE_quiet_NaN)                                             ! enforce cutback
      return
    endif


    !*** be aware of the definition of slip_transverse = slip_direction x slip_normal !!!
    !*** opposite sign to our t vector in the (s,t,n) triplet !!!

    m(1:3,:,1) =  prm%slip_direction
    m(1:3,:,2) = -prm%slip_direction
    m(1:3,:,3) = -prm%slip_transverse
    m(1:3,:,4) =  prm%slip_transverse

    my_F = phase_mechanical_F(ph)%data(1:3,1:3,en)
    my_Fe = matmul(my_F, math_inv33(phase_mechanical_Fp(ph)%data(1:3,1:3,en)))

    neighbors: do n = 1,nIPneighbors

      neighbor_el = IPneighborhood(1,n,ip,el)
      neighbor_ip = IPneighborhood(2,n,ip,el)
      neighbor_n  = IPneighborhood(3,n,ip,el)
      np = material_phaseAt(1,neighbor_el)
      no = material_phasememberAt(1,neighbor_ip,neighbor_el)

      opposite_neighbor = n + mod(n,2) - mod(n+1,2)
      opposite_el = IPneighborhood(1,opposite_neighbor,ip,el)
      opposite_ip = IPneighborhood(2,opposite_neighbor,ip,el)
      opposite_n  = IPneighborhood(3,opposite_neighbor,ip,el)

      if (neighbor_n > 0) then                                                                      ! if neighbor exists, average deformation gradient
        neighbor_ph = material_phaseAt(1,neighbor_el)
        neighbor_F = phase_mechanical_F(np)%data(1:3,1:3,no)
        neighbor_Fe = matmul(neighbor_F, math_inv33(phase_mechanical_Fp(np)%data(1:3,1:3,no)))
        Favg = 0.5_pReal * (my_F + neighbor_F)
      else                                                                                          ! if no neighbor, take my value as average
        Favg = my_F
      endif

      neighbor_v0 = 0.0_pReal        ! needed for check of sign change in flux density below

      !* FLUX FROM MY NEIGHBOR TO ME
      !* This is only considered, if I have a neighbor of nonlocal plasticity
      !* (also nonlocal constitutive law with local properties) that is at least a little bit
      !* compatible.
      !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of
      !* my neighbor's interface.
      !* The entering flux from my neighbor will be distributed on my slip systems according to the
      !* compatibility
      if (neighbor_n > 0) then
      if (phase_plasticity(material_phaseAt(1,neighbor_el)) == PLASTICITY_NONLOCAL_ID .and. &
          any(compatibility(:,:,:,n,ip,el) > 0.0_pReal)) then

        forall (s = 1:ns, t = 1:4)
          neighbor_v0(s,t) =          plasticState(np)%state0(iV   (s,t,neighbor_ph),no)
          neighbor_rhoSgl0(s,t) = max(plasticState(np)%state0(iRhoU(s,t,neighbor_ph),no),0.0_pReal)
        endforall

        where (neighbor_rhoSgl0 * IPvolume(neighbor_ip,neighbor_el) ** 0.667_pReal < prm%rho_min &
          .or. neighbor_rhoSgl0 < prm%rho_significant) &
          neighbor_rhoSgl0 = 0.0_pReal
        normal_neighbor2me_defConf = math_det33(Favg) * matmul(math_inv33(transpose(Favg)), &
                                     IPareaNormal(1:3,neighbor_n,neighbor_ip,neighbor_el))          ! normal of the interface in (average) deformed configuration (pointing neighbor => en)
        normal_neighbor2me = matmul(transpose(neighbor_Fe), normal_neighbor2me_defConf) &
                           / math_det33(neighbor_Fe)                                                ! interface normal in the lattice configuration of my neighbor
        area = IParea(neighbor_n,neighbor_ip,neighbor_el) * norm2(normal_neighbor2me)
        normal_neighbor2me = normal_neighbor2me / norm2(normal_neighbor2me)                         ! normalize the surface normal to unit length
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            topp = t + mod(t,2) - mod(t+1,2)
            if (neighbor_v0(s,t) * math_inner(m(1:3,s,t), normal_neighbor2me) > 0.0_pReal &         ! flux from my neighbor to en == entering flux for en
                .and. v0(s,t) * neighbor_v0(s,t) >= 0.0_pReal ) then                                ! ... only if no sign change in flux density
              lineLength = neighbor_rhoSgl0(s,t) * neighbor_v0(s,t) &
                         * math_inner(m(1:3,s,t), normal_neighbor2me) * area                        ! positive line length that wants to enter through this interface
              where (compatibility(c,:,s,n,ip,el) > 0.0_pReal) &
                rhoDotFlux(:,t) = rhoDotFlux(1:ns,t) &
                                + lineLength/IPvolume(ip,el)*compatibility(c,:,s,n,ip,el)**2.0_pReal    ! transferring to equally signed mobile dislocation type
              where (compatibility(c,:,s,n,ip,el) < 0.0_pReal) &
                rhoDotFlux(:,topp) = rhoDotFlux(:,topp) &
                                   + lineLength/IPvolume(ip,el)*compatibility(c,:,s,n,ip,el)**2.0_pReal ! transferring to opposite signed mobile dislocation type

            endif
          enddo
        enddo
      endif; endif


      !* FLUX FROM ME TO MY NEIGHBOR
      !* This is not considered, if my opposite neighbor has a different constitutive law than nonlocal (still considered for nonlocal law with local properties).
      !* Then, we assume, that the opposite(!) neighbor sends an equal amount of dislocations to en.
      !* So the net flux in the direction of my neighbor is equal to zero:
      !*    leaving flux to neighbor == entering flux from opposite neighbor
      !* In case of reduced transmissivity, part of the leaving flux is stored as dead dislocation density.
      !* That means for an interface of zero transmissivity the leaving flux is fully converted to dead dislocations.
      if (opposite_n > 0) then
      if (phase_plasticity(material_phaseAt(1,opposite_el)) == PLASTICITY_NONLOCAL_ID) then

        normal_me2neighbor_defConf = math_det33(Favg) &
                                   * matmul(math_inv33(transpose(Favg)),IPareaNormal(1:3,n,ip,el))  ! normal of the interface in (average) deformed configuration (pointing en => neighbor)
        normal_me2neighbor = matmul(transpose(my_Fe), normal_me2neighbor_defConf) &
                           / math_det33(my_Fe)                                                      ! interface normal in my lattice configuration
        area = IParea(n,ip,el) * norm2(normal_me2neighbor)
        normal_me2neighbor = normal_me2neighbor / norm2(normal_me2neighbor)                         ! normalize the surface normal to unit length
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            if (v0(s,t) * math_inner(m(1:3,s,t), normal_me2neighbor) > 0.0_pReal ) then             ! flux from en to my neighbor == leaving flux for en (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
              if (v0(s,t) * neighbor_v0(s,t) >= 0.0_pReal) then                                     ! no sign change in flux density
                transmissivity = sum(compatibility(c,:,s,n,ip,el)**2.0_pReal)                       ! overall transmissivity from this slip system to my neighbor
              else                                                                                  ! sign change in flux density means sign change in stress which does not allow for dislocations to arive at the neighbor
                transmissivity = 0.0_pReal
              endif
              lineLength = my_rhoSgl0(s,t) * v0(s,t) &
                         * math_inner(m(1:3,s,t), normal_me2neighbor) * area                        ! positive line length of mobiles that wants to leave through this interface
              rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / IPvolume(ip,el)                      ! subtract dislocation flux from current type
              rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) &
                                + lineLength / IPvolume(ip,el) * (1.0_pReal - transmissivity) &
                                * sign(1.0_pReal, v0(s,t))                                          ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
            endif
          enddo
        enddo
      endif; endif

    enddo neighbors
  endif

  end associate

end function rhoDotFlux


!--------------------------------------------------------------------------------------------------
!> @brief Compatibility update
!> @detail Compatibility is defined as normalized product of signed cosine of the angle between the slip
! plane normals and signed cosine of the angle between the slip directions. Only the largest values
! that sum up to a total of 1 are considered, all others are set to zero.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_updateCompatibility(orientation,ph,i,e)

  type(tRotationContainer), dimension(:), intent(in) :: &
    orientation                                                                                     ! crystal orientation
  integer, intent(in) :: &
    ph, &
    i, &
    e

  integer :: &
    n, &                                                                                            ! neighbor index
    en, &
    neighbor_e, &                                                                                   ! element index of my neighbor
    neighbor_i, &                                                                                   ! integration point index of my neighbor
    neighbor_me, &
    neighbor_phase, &
    ns, &                                                                                           ! number of active slip systems
    s1, &                                                                                           ! slip system index (en)
    s2                                                                                              ! slip system index (my neighbor)
  real(pReal), dimension(2,param(ph)%sum_N_sl,param(ph)%sum_N_sl,nIPneighbors) :: &
    my_compatibility                                                                                ! my_compatibility for current element and ip
  real(pReal) :: &
    my_compatibilitySum, &
    thresholdValue, &
    nThresholdValues
  logical, dimension(param(ph)%sum_N_sl) :: &
    belowThreshold
  type(rotation) :: mis

  associate(prm => param(ph))
  ns = prm%sum_N_sl

  en = material_phaseMemberAt(1,i,e)
  !*** start out fully compatible
  my_compatibility = 0.0_pReal
  forall(s1 = 1:ns) my_compatibility(:,s1,s1,:) = 1.0_pReal

  neighbors: do n = 1,nIPneighbors
    neighbor_e = IPneighborhood(1,n,i,e)
    neighbor_i = IPneighborhood(2,n,i,e)
    neighbor_me = material_phaseMemberAt(1,neighbor_i,neighbor_e)
    neighbor_phase = material_phaseAt(1,neighbor_e)

    if (neighbor_e <= 0 .or. neighbor_i <= 0) then
      !* FREE SURFACE
      !* Set surface transmissivity to the value specified in the material.config
      forall(s1 = 1:ns) my_compatibility(:,s1,s1,n) = sqrt(prm%chi_surface)
    elseif (neighbor_phase /= ph) then
      !* PHASE BOUNDARY
      !* If we encounter a different nonlocal phase at the neighbor,
      !* we consider this to be a real "physical" phase boundary, so completely incompatible.
      !* If one of the two phases has a local plasticity law,
      !* we do not consider this to be a phase boundary, so completely compatible.
      if (.not. phase_localPlasticity(neighbor_phase) .and. .not. phase_localPlasticity(ph)) &
        forall(s1 = 1:ns) my_compatibility(:,s1,s1,n) = 0.0_pReal
    elseif (prm%chi_GB >= 0.0_pReal) then
      !* GRAIN BOUNDARY !
      !* fixed transmissivity for adjacent ips with different texture (only if explicitly given in material.config)
      if (any(dNeq(phase_orientation0(ph)%data(en)%asQuaternion(), &
                   phase_orientation0(neighbor_phase)%data(neighbor_me)%asQuaternion())) .and. &
          (.not. phase_localPlasticity(neighbor_phase))) &
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
      mis = orientation(ph)%data(en)%misorientation(orientation(neighbor_phase)%data(neighbor_me))
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
        enddo neighborSlipSystems

        my_compatibilitySum = 0.0_pReal
        belowThreshold = .true.
        do while (my_compatibilitySum < 1.0_pReal .and. any(belowThreshold))
          thresholdValue = maxval(my_compatibility(2,:,s1,n), belowThreshold)                       ! screws always positive
          nThresholdValues = real(count(my_compatibility(2,:,s1,n) >= thresholdValue),pReal)
          where (my_compatibility(2,:,s1,n) >= thresholdValue) belowThreshold = .false.
          if (my_compatibilitySum + thresholdValue * nThresholdValues > 1.0_pReal) &
            where (abs(my_compatibility(:,:,s1,n)) >= thresholdValue) &
              my_compatibility(:,:,s1,n) = sign((1.0_pReal - my_compatibilitySum)/nThresholdValues,&
                                                 my_compatibility(:,:,s1,n))
          my_compatibilitySum = my_compatibilitySum + nThresholdValues * thresholdValue
        enddo

        where(belowThreshold) my_compatibility(1,:,s1,n) = 0.0_pReal
        where(belowThreshold) my_compatibility(2,:,s1,n) = 0.0_pReal

      enddo mySlipSystems
    endif

  enddo neighbors

  compatibility(:,:,:,:,i,e) = my_compatibility

  end associate

end subroutine plastic_nonlocal_updateCompatibility


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine plastic_nonlocal_results(ph,group)

  integer,         intent(in) :: ph
  character(len=*),intent(in) :: group

  integer :: o

  associate(prm => param(ph),dst => microstructure(ph),stt=>state(ph))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))
      case('rho_u_ed_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_mob_edg_pos,group,trim(prm%output(o)), &
                                                     'positive mobile edge density','1/m²')
      case('rho_b_ed_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_imm_edg_pos,group,trim(prm%output(o)), &
                                                     'positive immobile edge density','1/m²')
      case('rho_u_ed_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_mob_edg_neg,group,trim(prm%output(o)), &
                                                     'negative mobile edge density','1/m²')
      case('rho_b_ed_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_imm_edg_neg,group,trim(prm%output(o)), &
                                                     'negative immobile edge density','1/m²')
      case('rho_d_ed')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_dip_edg,group,trim(prm%output(o)), &
                                                     'edge dipole density','1/m²')
      case('rho_u_sc_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_mob_scr_pos,group,trim(prm%output(o)), &
                                                     'positive mobile screw density','1/m²')
      case('rho_b_sc_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_imm_scr_pos,group,trim(prm%output(o)), &
                                                     'positive immobile screw density','1/m²')
      case('rho_u_sc_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_mob_scr_neg,group,trim(prm%output(o)), &
                                                     'negative mobile screw density','1/m²')
      case('rho_b_sc_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_sgl_imm_scr_neg,group,trim(prm%output(o)), &
                                                     'negative immobile screw density','1/m²')
      case('rho_d_sc')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_dip_scr,group,trim(prm%output(o)), &
                                                     'screw dipole density','1/m²')
      case('rho_f')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%rho_forest,group,trim(prm%output(o)), &
                                                     'forest density','1/m²')
      case('v_ed_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%v_edg_pos,group,trim(prm%output(o)), &
                                                     'positive edge velocity','m/s')
      case('v_ed_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%v_edg_neg,group,trim(prm%output(o)), &
                                                     'negative edge velocity','m/s')
      case('v_sc_pos')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%v_scr_pos,group,trim(prm%output(o)), &
                                                     'positive srew velocity','m/s')
      case('v_sc_neg')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%v_scr_neg,group,trim(prm%output(o)), &
                                                     'negative screw velocity','m/s')
      case('gamma')
        if(prm%sum_N_sl>0) call results_writeDataset(stt%gamma,group,trim(prm%output(o)), &
                                                     'plastic shear','1')
      case('tau_pass')
        if(prm%sum_N_sl>0) call results_writeDataset(dst%tau_pass,group,trim(prm%output(o)), &
                                                     'passing stress for slip','Pa')
    end select
  enddo outputsLoop
  end associate

end subroutine plastic_nonlocal_results


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
  real(pReal), dimension(2) :: &
    noise, &
    rnd
  real(pReal) :: &
    meanDensity, &
    totalVolume, &
    densityBinning, &
    minimumIpVolume


  associate(stt => state(phase))

    if (ini%random_rho_u > 0.0_pReal) then ! randomly distribute dislocation segments on random slip system and of random type in the volume
      totalVolume     = sum(geom(phase)%V_0)
      minimumIPVolume = minval(geom(phase)%V_0)
      densityBinning  = ini%random_rho_u_binning / minimumIpVolume ** (2.0_pReal / 3.0_pReal)

      ! fill random material points with dislocation segments until the desired overall density is reached
      meanDensity = 0.0_pReal
      do while(meanDensity < ini%random_rho_u)
        call random_number(rnd)
        e = nint(rnd(1)*real(Nentries,pReal) + 0.5_pReal)
        s = nint(rnd(2)*real(sum(ini%N_sl),pReal)*4.0_pReal + 0.5_pReal)
        meanDensity = meanDensity + densityBinning * geom(phase)%V_0(e) / totalVolume
        stt%rhoSglMobile(s,e) = densityBinning
      enddo
    else                                ! homogeneous distribution with noise
      do e = 1, Nentries
        do f = 1,size(ini%N_sl,1)
          from = 1 + sum(ini%N_sl(1:f-1))
          upto = sum(ini%N_sl(1:f))
          do s = from,upto
            noise = [math_sampleGaussVar(0.0_pReal, ini%sigma_rho_u), &
                     math_sampleGaussVar(0.0_pReal, ini%sigma_rho_u)]
            stt%rho_sgl_mob_edg_pos(s,e) = ini%rho_u_ed_pos_0(f) + noise(1)
            stt%rho_sgl_mob_edg_neg(s,e) = ini%rho_u_ed_neg_0(f) + noise(1)
            stt%rho_sgl_mob_scr_pos(s,e) = ini%rho_u_sc_pos_0(f) + noise(2)
            stt%rho_sgl_mob_scr_neg(s,e) = ini%rho_u_sc_neg_0(f) + noise(2)
          enddo
          stt%rho_dip_edg(from:upto,e)   = ini%rho_d_ed_0(f)
          stt%rho_dip_scr(from:upto,e)   = ini%rho_d_sc_0(f)
        enddo
      enddo
    endif

  end associate

end subroutine stateInit


!--------------------------------------------------------------------------------------------------
!> @brief calculates kinetics
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(v, dv_dtau, dv_dtauNS, tau, tauNS, tauThreshold, c, Temperature, ph)

  integer, intent(in) :: &
    c, &                                                                                            !< dislocation character (1:edge, 2:screw)
    ph
  real(pReal), intent(in) :: &
    Temperature                                                                                     !< temperature
  real(pReal), dimension(param(ph)%sum_N_sl), intent(in) :: &
    tau, &                                                                                          !< resolved external shear stress (without non Schmid effects)
    tauNS, &                                                                                        !< resolved external shear stress (including non Schmid effects)
    tauThreshold                                                                                    !< threshold shear stress
  real(pReal), dimension(param(ph)%sum_N_sl), intent(out) ::  &
    v, &                                                                                            !< velocity
    dv_dtau, &                                                                                      !< velocity derivative with respect to resolved shear stress (without non Schmid contributions)
    dv_dtauNS                                                                                       !< velocity derivative with respect to resolved shear stress (including non Schmid contributions)

  integer :: &
    s                                                                                               !< index of my current slip system
  real(pReal) :: &
    tauRel_P, &
    tauRel_S, &
    tauEff, &                                                                                       !< effective shear stress
    tPeierls, &                                                                                     !< waiting time in front of a peierls barriers
    tSolidSolution, &                                                                               !< waiting time in front of a solid solution obstacle
    vViscous, &                                                                                     !< viscous glide velocity
    dtPeierls_dtau, &                                                                               !< derivative with respect to resolved shear stress
    dtSolidSolution_dtau, &                                                                         !< derivative with respect to resolved shear stress
    meanfreepath_S, &                                                                               !< mean free travel distance for dislocations between two solid solution obstacles
    meanfreepath_P, &                                                                               !< mean free travel distance for dislocations between two Peierls barriers
    jumpWidth_P, &                                                                                  !< depth of activated area
    jumpWidth_S, &                                                                                  !< depth of activated area
    activationLength_P, &                                                                           !< length of activated dislocation line
    activationLength_S, &                                                                           !< length of activated dislocation line
    activationVolume_P, &                                                                           !< volume that needs to be activated to overcome barrier
    activationVolume_S, &                                                                           !< volume that needs to be activated to overcome barrier
    activationEnergy_P, &                                                                           !< energy that is needed to overcome barrier
    activationEnergy_S, &                                                                           !< energy that is needed to overcome barrier
    criticalStress_P, &                                                                             !< maximum obstacle strength
    criticalStress_S, &                                                                             !< maximum obstacle strength
    mobility                                                                                        !< dislocation mobility

  associate(prm => param(ph))
  v = 0.0_pReal
  dv_dtau = 0.0_pReal
  dv_dtauNS = 0.0_pReal

  do s = 1,prm%sum_N_sl
    if (abs(tau(s)) > tauThreshold(s)) then

      !* Peierls contribution
      !* Effective stress includes non Schmid constributions
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
      tauEff = max(0.0_pReal, abs(tauNS(s)) - tauThreshold(s))                                      ! ensure that the effective stress is positive
      meanfreepath_P = prm%b_sl(s)
      jumpWidth_P = prm%b_sl(s)
      activationLength_P = prm%w *prm%b_sl(s)
      activationVolume_P = activationLength_P * jumpWidth_P * prm%b_sl(s)
      criticalStress_P = prm%peierlsStress(s,c)
      activationEnergy_P = criticalStress_P * activationVolume_P
      tauRel_P = min(1.0_pReal, tauEff / criticalStress_P)                                          ! ensure that the activation probability cannot become greater than one
      tPeierls = 1.0_pReal / prm%nu_a &
               * exp(activationEnergy_P / (kB * Temperature) &
                     * (1.0_pReal - tauRel_P**prm%p)**prm%q)
      if (tauEff < criticalStress_P) then
        dtPeierls_dtau = tPeierls * prm%p * prm%q * activationVolume_P / (kB * Temperature) &
                       * (1.0_pReal - tauRel_P**prm%p)**(prm%q-1.0_pReal) * tauRel_P**(prm%p-1.0_pReal)
      else
        dtPeierls_dtau = 0.0_pReal
      endif

      !* Contribution from solid solution strengthening
      !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
      tauEff = abs(tau(s)) - tauThreshold(s)
      meanfreepath_S = prm%b_sl(s) / sqrt(prm%c_sol)
      jumpWidth_S = prm%f_sol * prm%b_sl(s)
      activationLength_S = prm%b_sl(s) / sqrt(prm%c_sol)
      activationVolume_S = activationLength_S * jumpWidth_S * prm%b_sl(s)
      activationEnergy_S = prm%Q_sol
      criticalStress_S = activationEnergy_S / activationVolume_S
      tauRel_S = min(1.0_pReal, tauEff / criticalStress_S)                                          ! ensure that the activation probability cannot become greater than one
      tSolidSolution = 1.0_pReal /  prm%nu_a &
                     * exp(activationEnergy_S / (kB * Temperature)* (1.0_pReal - tauRel_S**prm%p)**prm%q)
      if (tauEff < criticalStress_S) then
        dtSolidSolution_dtau = tSolidSolution * prm%p * prm%q * activationVolume_S / (kB * Temperature) &
                             * (1.0_pReal - tauRel_S**prm%p)**(prm%q-1.0_pReal)* tauRel_S**(prm%p-1.0_pReal)
      else
        dtSolidSolution_dtau = 0.0_pReal
      endif

      !* viscous glide velocity
      tauEff = abs(tau(s)) - tauThreshold(s)
      mobility = prm%b_sl(s) / prm%eta
      vViscous = mobility * tauEff

      !* Mean velocity results from waiting time at peierls barriers and solid solution obstacles with respective meanfreepath of
      !* free flight at glide velocity in between.
      !* adopt sign from resolved stress
      v(s) = sign(1.0_pReal,tau(s)) &
           / (tPeierls / meanfreepath_P + tSolidSolution / meanfreepath_S + 1.0_pReal / vViscous)
      dv_dtau(s)   = v(s)**2.0_pReal * (dtSolidSolution_dtau / meanfreepath_S + mobility /vViscous**2.0_pReal)
      dv_dtauNS(s) = v(s)**2.0_pReal * dtPeierls_dtau / meanfreepath_P
    endif
  enddo

  end associate

end subroutine kinetics


!--------------------------------------------------------------------------------------------------
!> @brief returns copy of current dislocation densities from state
!> @details raw values is rectified
!--------------------------------------------------------------------------------------------------
pure function getRho(ph,en) result(rho)

  integer, intent(in) :: ph, en
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: rho


  associate(prm => param(ph))

    rho = reshape(state(ph)%rho(:,en),[prm%sum_N_sl,10])

    ! ensure positive densities (not for imm, they have a sign)
    rho(:,mob) = max(rho(:,mob),0.0_pReal)
    rho(:,dip) = max(rho(:,dip),0.0_pReal)

    where(abs(rho) < max(prm%rho_min/geom(ph)%V_0(en)**(2.0_pReal/3.0_pReal),prm%rho_significant)) &
      rho = 0.0_pReal

  end associate

end function getRho


!--------------------------------------------------------------------------------------------------
!> @brief returns copy of current dislocation densities from state
!> @details raw values is rectified
!--------------------------------------------------------------------------------------------------
pure function getRho0(ph,en) result(rho0)

  integer, intent(in) :: ph, en
  real(pReal), dimension(param(ph)%sum_N_sl,10) :: rho0


  associate(prm => param(ph))

    rho0 = reshape(state0(ph)%rho(:,en),[prm%sum_N_sl,10])

    ! ensure positive densities (not for imm, they have a sign)
    rho0(:,mob) = max(rho0(:,mob),0.0_pReal)
    rho0(:,dip) = max(rho0(:,dip),0.0_pReal)

    where (abs(rho0) < max(prm%rho_min/geom(ph)%V_0(en)**(2.0_pReal/3.0_pReal),prm%rho_significant)) &
      rho0 = 0.0_pReal

  end associate

end function getRho0


subroutine storeGeometry(ph)

  integer, intent(in) :: ph

  integer :: ce, co
  real(pReal), dimension(:), allocatable :: V


  V = reshape(IPvolume,[product(shape(IPvolume))])
  do ce = 1, size(material_homogenizationEntry,1)
    do co = 1, homogenization_maxNconstituents
      if (material_phaseID(co,ce) == ph) geom(ph)%V_0(material_phaseEntry(co,ce)) = V(ce)
    enddo
  enddo

end subroutine

end submodule nonlocal
