!--------------------------------------------------------------------------------------------------
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for plasticity including dislocation flux
!--------------------------------------------------------------------------------------------------
module plastic_nonlocal
  use prec
  use IO
  use math
  use debug
  use material
  use lattice
  use rotations
  use config
  use lattice
  use discretization
  use geometry_plastic_nonlocal, only: &
    nIPneighbors    => geometry_plastic_nonlocal_nIPneighbors, &
    IPneighborhood  => geometry_plastic_nonlocal_IPneighborhood, &
    IPvolume        => geometry_plastic_nonlocal_IPvolume0, &
    IParea          => geometry_plastic_nonlocal_IParea0, &
    IPareaNormal    => geometry_plastic_nonlocal_IPareaNormal0
    
  implicit none
  private
  real(pReal), parameter, private :: &
    KB = 1.38e-23_pReal                                                                             !< Physical parameter, Boltzmann constant in J/Kelvin
 
  integer, dimension(:,:), allocatable, target, public :: &
    plastic_nonlocal_sizePostResult                                                                 !< size of each post result output
  
  character(len=64), dimension(:,:), allocatable, target, public :: &
    plastic_nonlocal_output                                                                         !< name of each post result output
 

  ! storage order of dislocation types
  integer, dimension(8), parameter :: &
    sgl = [1,2,3,4,5,6,7,8]                                                                         !< signed (single)
  integer, dimension(5), parameter :: &
    edg = [1,2,5,6,9], &                                                                            !< edge
    scr = [3,4,7,8,10]                                                                              !< screw
  integer, dimension(4), parameter :: &
    mob = [1,2,3,4], &                                                                              !< mobile
    imm = [5,6,7,8]                                                                                 !< immobile (blocked)
  integer, dimension(2), parameter :: &
    dip = [9,10], &                                                                                 !< dipole
    imm_edg = imm(1:2), &                                                                           !< immobile edge
    imm_scr = imm(3:4)                                                                              !< immobile screw
  integer, parameter :: &
    mob_edg_pos = 1, &                                                                              !< mobile edge positive
    mob_edg_neg = 2, &                                                                              !< mobile edge negative
    mob_scr_pos = 3, &                                                                              !< mobile screw positive
    mob_scr_neg = 4                                                                                 !< mobile screw positive

  ! BEGIN DEPRECATES
  integer, dimension(:,:,:), allocatable, private :: &
    iRhoU, &                                                                                        !< state indices for unblocked density
    iRhoB, &                                                                                        !< state indices for blocked density
    iRhoD, &                                                                                        !< state indices for dipole density
    iV, &                                                                                           !< state indices for dislcation velocities
    iD                                                                                              !< state indices for stable dipole height
  integer, dimension(:), allocatable, private, protected :: &
    totalNslip                                                                                      !< total number of active slip systems for each instance
  !END DEPRECATED
  
  real(pReal), dimension(:,:,:,:,:,:), allocatable, private :: &
    compatibility                                                                                    !< slip system compatibility between me and my neighbors
  
  enum, bind(c) 
    enumerator :: &
      undefined_ID, &
      rho_sgl_mob_edg_pos_ID, &
      rho_sgl_mob_edg_neg_ID, &
      rho_sgl_mob_scr_pos_ID, &
      rho_sgl_mob_scr_neg_ID, &
      rho_sgl_imm_edg_pos_ID, &
      rho_sgl_imm_edg_neg_ID, &
      rho_sgl_imm_scr_pos_ID, &
      rho_sgl_imm_scr_neg_ID, &
      rho_dip_edg_ID, &
      rho_dip_scr_ID, &
      rho_forest_ID, &
      resolvedstress_back_ID, &
      resistance_ID, &
      rho_dot_sgl_ID, &
      rho_dot_sgl_mobile_ID, &
      rho_dot_dip_ID, &
      velocity_edge_pos_ID, &
      velocity_edge_neg_ID, &
      velocity_screw_pos_ID, &
      velocity_screw_neg_ID, &
      accumulatedshear_ID
  end enum
  
  type, private :: tParameters                                                                      !< container type for internal constitutive parameters
    real(pReal) :: &
      atomicVolume, &                                                                               !< atomic volume
      Dsd0, &                                                                                       !< prefactor for self-diffusion coefficient
      selfDiffusionEnergy, &                                                                        !< activation enthalpy for diffusion
      aTolRho, &                                                                                    !< absolute tolerance for dislocation density in state integration
      aTolShear, &                                                                                  !< absolute tolerance for accumulated shear in state integration
      significantRho, &                                                                             !< density considered significant
      significantN, &                                                                               !< number of dislocations considered significant
      doublekinkwidth, &                                                                            !< width of a doubkle kink in multiples of the burgers vector length b
      solidSolutionEnergy, &                                                                        !< activation energy for solid solution in J
      solidSolutionSize, &                                                                          !< solid solution obstacle size in multiples of the burgers vector length
      solidSolutionConcentration, &                                                                 !< concentration of solid solution in atomic parts
      p, &                                                                                          !< parameter for kinetic law (Kocks,Argon,Ashby)
      q, &                                                                                          !< parameter for kinetic law (Kocks,Argon,Ashby)
      viscosity, &                                                                                  !< viscosity for dislocation glide in Pa s
      fattack, &                                                                                    !< attack frequency in Hz
      rhoSglScatter, &                                                                              !< standard deviation of scatter in initial dislocation density
      surfaceTransmissivity, &                                                                      !< transmissivity at free surface
      grainboundaryTransmissivity, &                                                                !< transmissivity at grain boundary (identified by different texture)
      CFLfactor, &                                                                                  !< safety factor for CFL flux condition
      fEdgeMultiplication, &                                                                        !< factor that determines how much edge dislocations contribute to multiplication (0...1)
      rhoSglRandom, &
      rhoSglRandomBinning, &
      linetensionEffect, &
      edgeJogFactor, &
      mu, &
      nu
    real(pReal), dimension(:), allocatable     :: &
      minDipoleHeight_edge, &                                                                       !< minimum stable edge dipole height
      minDipoleHeight_screw, &                                                                      !< minimum stable screw dipole height
      peierlsstress_edge, &
      peierlsstress_screw, &
      rhoSglEdgePos0, &                                                                             !< initial edge_pos dislocation density
      rhoSglEdgeNeg0, &                                                                             !< initial edge_neg dislocation density
      rhoSglScrewPos0, &                                                                            !< initial screw_pos dislocation density
      rhoSglScrewNeg0, &                                                                            !< initial screw_neg dislocation density
      rhoDipEdge0, &                                                                                !< initial edge dipole dislocation density
      rhoDipScrew0,&                                                                                !< initial screw dipole dislocation density
      lambda0, &                                                                                    !< mean free path prefactor for each
      burgers                                                                                       !< absolute length of burgers vector [m]
    real(pReal), dimension(:,:), allocatable     :: &
      slip_normal, &
      slip_direction, &
      slip_transverse, &
      minDipoleHeight, &                                                                            ! edge and screw
      peierlsstress, &                                                                              ! edge and screw
      interactionSlipSlip ,&                                                                        !< coefficients for slip-slip interaction
      forestProjection_Edge, &                                                                      !< matrix of forest projections of edge dislocations
      forestProjection_Screw                                                                        !< matrix of forest projections of screw dislocations
    real(pReal), dimension(:), allocatable, private :: &
      nonSchmidCoeff
    real(pReal), dimension(:,:,:), allocatable, private :: &
      Schmid, &                                                                                     !< Schmid contribution
      nonSchmid_pos, &
      nonSchmid_neg                                                                                 !< combined projection of Schmid and non-Schmid contributions to the resolved shear stress (only for screws)
    integer :: &
      totalNslip
    integer, dimension(:) ,allocatable , public:: &
      Nslip,&
      colinearSystem                                                                                !< colinear system to the active slip system (only valid for fcc!)
 
    logical, private :: &
      shortRangeStressCorrection, &                                                                 !< flag indicating the use of the short range stress correction by a excess density gradient term
      probabilisticMultiplication
    
    integer(kind(undefined_ID)),  dimension(:), allocatable :: &
      outputID                                                                                      !< ID of each post result output
      
  end type tParameters
   
  type, private :: tNonlocalMicrostructure
    real(pReal), allocatable, dimension(:,:) :: &
     tau_Threshold, & 
     tau_Back 
  end type tNonlocalMicrostructure
  
  type, private :: tNonlocalState
    real(pReal), pointer, dimension(:,:) :: &
      rho, &                                                                                        ! < all dislocations
        rhoSgl, &
          rhoSglMobile, &                       ! iRhoU                               
            rho_sgl_mob_edg_pos, &
            rho_sgl_mob_edg_neg, &
            rho_sgl_mob_scr_pos, &
            rho_sgl_mob_scr_neg, &
          rhoSglImmobile, &                     ! iRhoB
            rho_sgl_imm_edg_pos, &
            rho_sgl_imm_edg_neg, &
            rho_sgl_imm_scr_pos, &
            rho_sgl_imm_scr_neg, &
        rhoDip, &                               ! iRhoD
          rho_dip_edg, &
          rho_dip_scr, &
        rho_forest, &
      accumulatedshear, &
      v
  end type tNonlocalState
  
  type(tNonlocalState), allocatable, dimension(:), private :: &
    deltaState, &
    dotState, &
    state
 
  type(tParameters), dimension(:), allocatable, private :: param                                    !< containers of constitutive parameters (len Ninstance)
 
  type(tNonlocalMicrostructure), dimension(:), allocatable, private :: microstructure

  integer(kind(undefined_ID)), dimension(:,:),   allocatable, private :: & 
    plastic_nonlocal_outputID                                                                       !< ID of each post result output
  
  public :: &
    plastic_nonlocal_init, &
    plastic_nonlocal_dependentState, &
    plastic_nonlocal_LpAndItsTangent, &
    plastic_nonlocal_dotState, &
    plastic_nonlocal_deltaState, &
    plastic_nonlocal_updateCompatibility, &
    plastic_nonlocal_postResults, &
    plastic_nonlocal_results
  
  private :: &
    plastic_nonlocal_kinetics

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_init

  character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  integer,                dimension(0), parameter :: emptyIntArray    = [integer::]
  real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]

  integer :: &
    sizeState, sizeDotState,sizeDependentState, sizeDeltaState, &
    maxNinstances, &
    p, i, &
    l, &
    s1, s2, &
    s, &
    t, &
    c

  integer(kind(undefined_ID)) :: &
    outputID
  character(len=512) :: &
    extmsg    = '', &
    structure
  character(len=65536), dimension(:), allocatable :: outputs
  integer :: NofMyPhase 
 
  write(6,'(/,a)') ' <<<+-  constitutive_'//PLASTICITY_NONLOCAL_label//' init  -+>>>'

  write(6,'(/,a)') ' Reuber et al., Acta Materialia 71:333–348, 2014'
  write(6,'(a)')   ' https://doi.org/10.1016/j.actamat.2014.03.012'

  write(6,'(/,a)') ' Kords, Dissertation RWTH Aachen, 2014'
  write(6,'(a)')   ' http://publications.rwth-aachen.de/record/229993'

  maxNinstances = count(phase_plasticity == PLASTICITY_NONLOCAL_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',maxNinstances

  allocate(param(maxNinstances))
  allocate(state(maxNinstances))
  allocate(dotState(maxNinstances))
  allocate(deltaState(maxNinstances))
  allocate(microstructure(maxNinstances))

  allocate(plastic_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstances), source=0)
  allocate(plastic_nonlocal_output(maxval(phase_Noutput), maxNinstances))
           plastic_nonlocal_output = ''
  allocate(plastic_nonlocal_outputID(maxval(phase_Noutput), maxNinstances), source=undefined_ID)
  allocate(totalNslip(maxNinstances), source=0)


  do p=1, size(config_phase)
    if (phase_plasticity(p) /= PLASTICITY_NONLOCAL_ID) cycle

    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)), &
              del => deltaState(phase_plasticityInstance(p)), &
              dst => microstructure(phase_plasticityInstance(p)), &
              config => config_phase(p))

    prm%aTolRho    = config%getFloat('atol_rho',   defaultVal=0.0_pReal)
    prm%aTolShear  = config%getFloat('atol_shear', defaultVal=0.0_pReal)
   
    structure      = config%getString('lattice_structure')

    ! This data is read in already in lattice
    prm%mu = lattice_mu(p)
    prm%nu = lattice_nu(p)


    prm%Nslip      = config%getInts('nslip',defaultVal=emptyIntArray)
    prm%totalNslip = sum(prm%Nslip)
    slipActive: if (prm%totalNslip > 0) then
      prm%Schmid = lattice_SchmidMatrix_slip(prm%Nslip,config%getString('lattice_structure'),&
                                             config%getFloat('c/a',defaultVal=0.0_pReal))

      if(trim(config%getString('lattice_structure')) == 'bcc') then
        prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                               defaultVal = emptyRealArray)
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1)
      else
        prm%nonSchmid_pos  = prm%Schmid
        prm%nonSchmid_neg  = prm%Schmid
      endif

      prm%interactionSlipSlip = lattice_interaction_SlipBySlip(prm%Nslip, &
                                                               config%getFloats('interaction_slipslip'), &
                                                               config%getString('lattice_structure'))

      prm%forestProjection_edge  = lattice_forestProjection_edge (prm%Nslip,config%getString('lattice_structure'),&
                                                                  config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%forestProjection_screw = lattice_forestProjection_screw(prm%Nslip,config%getString('lattice_structure'),&
                                                                  config%getFloat('c/a',defaultVal=0.0_pReal))

      prm%slip_direction  = lattice_slip_direction (prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_transverse = lattice_slip_transverse(prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_normal     = lattice_slip_normal    (prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))

      ! collinear systems (only for octahedral slip systems in fcc)
      allocate(prm%colinearSystem(prm%totalNslip), source = -1)
      do s1 = 1, prm%totalNslip
        do s2 = 1, prm%totalNslip
           if (all(dEq0 (math_cross(prm%slip_direction(1:3,s1),prm%slip_direction(1:3,s2)))) .and. &
               any(dNeq0(math_cross(prm%slip_normal   (1:3,s1),prm%slip_normal   (1:3,s2))))) &
             prm%colinearSystem(s1) = s2
        enddo
      enddo
    
      prm%rhoSglEdgePos0  = config%getFloats('rhosgledgepos0',   requiredSize=size(prm%Nslip))
      prm%rhoSglEdgeNeg0  = config%getFloats('rhosgledgeneg0',   requiredSize=size(prm%Nslip))
      prm%rhoSglScrewPos0 = config%getFloats('rhosglscrewpos0',  requiredSize=size(prm%Nslip))
      prm%rhoSglScrewNeg0 = config%getFloats('rhosglscrewneg0',  requiredSize=size(prm%Nslip))
      prm%rhoDipEdge0     = config%getFloats('rhodipedge0',      requiredSize=size(prm%Nslip))
      prm%rhoDipScrew0    = config%getFloats('rhodipscrew0',     requiredSize=size(prm%Nslip))
   
      prm%lambda0         = config%getFloats('lambda0',          requiredSize=size(prm%Nslip))
      prm%burgers         = config%getFloats('burgers',          requiredSize=size(prm%Nslip))
      
      prm%lambda0 = math_expand(prm%lambda0,prm%Nslip)
      prm%burgers = math_expand(prm%burgers,prm%Nslip)
      
      prm%minDipoleHeight_edge  = config%getFloats('minimumdipoleheightedge',  requiredSize=size(prm%Nslip))
      prm%minDipoleHeight_screw = config%getFloats('minimumdipoleheightscrew', requiredSize=size(prm%Nslip))
      prm%minDipoleHeight_edge  = math_expand(prm%minDipoleHeight_edge,prm%Nslip)
      prm%minDipoleHeight_screw =  math_expand(prm%minDipoleHeight_screw,prm%Nslip)
      allocate(prm%minDipoleHeight(prm%totalNslip,2))
      prm%minDipoleHeight(:,1)  = prm%minDipoleHeight_edge
      prm%minDipoleHeight(:,2)  = prm%minDipoleHeight_screw
      
      prm%peierlsstress_edge    = config%getFloats('peierlsstressedge',        requiredSize=size(prm%Nslip))
      prm%peierlsstress_screw   = config%getFloats('peierlsstressscrew',       requiredSize=size(prm%Nslip))
      prm%peierlsstress_edge    = math_expand(prm%peierlsstress_edge,prm%Nslip)
      prm%peierlsstress_screw   = math_expand(prm%peierlsstress_screw,prm%Nslip)
      allocate(prm%peierlsstress(prm%totalNslip,2))
      prm%peierlsstress(:,1)    = prm%peierlsstress_edge
      prm%peierlsstress(:,2)    = prm%peierlsstress_screw

      prm%significantRho              = config%getFloat('significantrho')
      prm%significantN                = config%getFloat('significantn', 0.0_pReal)
      prm%CFLfactor                   = config%getFloat('cflfactor',defaultVal=2.0_pReal)
   
      prm%atomicVolume                = config%getFloat('atomicvolume')
      prm%Dsd0                        = config%getFloat('selfdiffusionprefactor') !,'dsd0'
      prm%selfDiffusionEnergy         = config%getFloat('selfdiffusionenergy') !,'qsd'
      prm%linetensionEffect           = config%getFloat('linetension')
      prm%edgeJogFactor               = config%getFloat('edgejog')!,'edgejogs'
      prm%doublekinkwidth             = config%getFloat('doublekinkwidth')
      prm%solidSolutionEnergy         = config%getFloat('solidsolutionenergy')
      prm%solidSolutionSize           = config%getFloat('solidsolutionsize')
      prm%solidSolutionConcentration  = config%getFloat('solidsolutionconcentration')
   
      prm%p                           = config%getFloat('p')
      prm%q                           = config%getFloat('q')
      prm%viscosity                   = config%getFloat('viscosity')
      prm%fattack                     = config%getFloat('attackfrequency')

      ! ToDo: discuss logic
      prm%rhoSglScatter               = config%getFloat('rhosglscatter')
      prm%rhoSglRandom                = config%getFloat('rhosglrandom',0.0_pReal)
      if (config%keyExists('/rhosglrandom/')) & 
        prm%rhoSglRandomBinning       = config%getFloat('rhosglrandombinning',0.0_pReal) !ToDo: useful default?
     ! if (rhoSglRandom(instance) < 0.0_pReal) &
     ! if (rhoSglRandomBinning(instance) <= 0.0_pReal) &
   
      prm%surfaceTransmissivity       = config%getFloat('surfacetransmissivity',defaultVal=1.0_pReal)
      prm%grainboundaryTransmissivity = config%getFloat('grainboundarytransmissivity',defaultVal=-1.0_pReal)
      prm%fEdgeMultiplication         = config%getFloat('edgemultiplication')
      prm%shortRangeStressCorrection  = config%keyExists('/shortrangestresscorrection/')
   
!--------------------------------------------------------------------------------------------------
!  sanity checks
      if (any(prm%burgers          <  0.0_pReal)) extmsg = trim(extmsg)//' burgers'
      if (any(prm%lambda0          <= 0.0_pReal)) extmsg = trim(extmsg)//' lambda0'
      
      if (any(prm%rhoSglEdgePos0   <  0.0_pReal)) extmsg = trim(extmsg)//' rhoSglEdgePos0'
      if (any(prm%rhoSglEdgeNeg0   <  0.0_pReal)) extmsg = trim(extmsg)//' rhoSglEdgeNeg0'
      if (any(prm%rhoSglScrewPos0  <  0.0_pReal)) extmsg = trim(extmsg)//' rhoSglScrewPos0'
      if (any(prm%rhoSglScrewNeg0  <  0.0_pReal)) extmsg = trim(extmsg)//' rhoSglScrewNeg0'
      if (any(prm%rhoDipEdge0      <  0.0_pReal)) extmsg = trim(extmsg)//' rhoDipEdge0'
      if (any(prm%rhoDipScrew0     <  0.0_pReal)) extmsg = trim(extmsg)//' rhoDipScrew0'
      
      if (any(prm%peierlsstress    <  0.0_pReal)) extmsg = trim(extmsg)//' peierlsstress'
      if (any(prm%minDipoleHeight  <  0.0_pReal)) extmsg = trim(extmsg)//' minDipoleHeight'
          
      if (prm%viscosity           <= 0.0_pReal)   extmsg = trim(extmsg)//' viscosity'
      if (prm%selfDiffusionEnergy <= 0.0_pReal)   extmsg = trim(extmsg)//' selfDiffusionEnergy'
      if (prm%fattack             <= 0.0_pReal)   extmsg = trim(extmsg)//' fattack'
      if (prm%doublekinkwidth     <= 0.0_pReal)   extmsg = trim(extmsg)//' doublekinkwidth'
      if (prm%Dsd0                < 0.0_pReal)    extmsg = trim(extmsg)//' Dsd0'
      if (prm%atomicVolume        <= 0.0_pReal)   extmsg = trim(extmsg)//' atomicVolume'            ! ToDo: in disloUCLA/dislotwin, the atomic volume is given as a factor
      
      if (prm%significantN         < 0.0_pReal)   extmsg = trim(extmsg)//' significantN'
      if (prm%significantrho       < 0.0_pReal)   extmsg = trim(extmsg)//' significantrho'
      if (prm%atolshear           <= 0.0_pReal)   extmsg = trim(extmsg)//' atolshear'
      if (prm%atolrho             <= 0.0_pReal)   extmsg = trim(extmsg)//' atolrho'
      if (prm%CFLfactor            < 0.0_pReal)   extmsg = trim(extmsg)//' CFLfactor'
      
      if (prm%p <= 0.0_pReal .or. prm%p > 1.0_pReal) extmsg = trim(extmsg)//' p' 
      if (prm%q <  1.0_pReal .or. prm%q > 2.0_pReal) extmsg = trim(extmsg)//' q' 
      
      if (prm%linetensionEffect <  0.0_pReal .or. prm%linetensionEffect  > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' linetensionEffect'
      if (prm%edgeJogFactor     <  0.0_pReal .or. prm%edgeJogFactor      > 1.0_pReal) &
                                                  extmsg = trim(extmsg)//' edgeJogFactor'
                                                  
      if (prm%solidSolutionEnergy        <= 0.0_pReal) extmsg = trim(extmsg)//' solidSolutionEnergy'
      if (prm%solidSolutionSize          <= 0.0_pReal) extmsg = trim(extmsg)//' solidSolutionSize'
      if (prm%solidSolutionConcentration <= 0.0_pReal) extmsg = trim(extmsg)//' solidSolutionConcentration'

      if (prm%grainboundaryTransmissivity  > 1.0_pReal) extmsg = trim(extmsg)//' grainboundaryTransmissivity' 
      if (prm%surfaceTransmissivity  <  0.0_pReal .or. prm%surfaceTransmissivity  > 1.0_pReal) &
                                                        extmsg = trim(extmsg)//' surfaceTransmissivity' 
 
      if (prm%fEdgeMultiplication  <  0.0_pReal .or. prm%fEdgeMultiplication  > 1.0_pReal) &
      extmsg = trim(extmsg)//' fEdgeMultiplication' 

    endif slipActive

    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))
    do i=1, size(outputs)
      outputID = undefined_ID
      select case(trim(outputs(i)))
        case ('rho_sgl_edge_pos_mobile')
          outputID = merge(rho_sgl_mob_edg_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_edge_neg_mobile')
          outputID = merge(rho_sgl_mob_edg_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_screw_pos_mobile')
          outputID = merge(rho_sgl_mob_scr_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_screw_neg_mobile')
          outputID = merge(rho_sgl_mob_scr_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_edge_pos_immobile')
          outputID = merge(rho_sgl_imm_edg_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_edge_neg_immobile')
          outputID = merge(rho_sgl_imm_edg_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_screw_pos_immobile')
          outputID = merge(rho_sgl_imm_scr_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_sgl_screw_neg_immobile')
          outputID = merge(rho_sgl_imm_scr_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_dip_edge')
          outputID = merge(rho_dip_edg_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_dip_screw')
          outputID = merge(rho_dip_scr_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_forest')
          outputID = merge(rho_forest_ID,undefined_ID,prm%totalNslip>0)
        case ('resolvedstress_back')
          outputID = merge(resolvedstress_back_ID,undefined_ID,prm%totalNslip>0)
        case ('resistance')
          outputID = merge(resistance_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_dot_sgl')
          outputID = merge(rho_dot_sgl_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_dot_sgl_mobile')
          outputID = merge(rho_dot_sgl_mobile_ID,undefined_ID,prm%totalNslip>0)
        case ('rho_dot_dip')
          outputID = merge(rho_dot_dip_ID,undefined_ID,prm%totalNslip>0)
        case ('velocity_edge_pos')
          outputID = merge(velocity_edge_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('velocity_edge_neg')
          outputID = merge(velocity_edge_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('velocity_screw_pos')
          outputID = merge(velocity_screw_pos_ID,undefined_ID,prm%totalNslip>0)
        case ('velocity_screw_neg')
          outputID = merge(velocity_screw_neg_ID,undefined_ID,prm%totalNslip>0)
        case ('accumulatedshear','accumulated_shear')
          outputID = merge(accumulatedshear_ID,undefined_ID,prm%totalNslip>0)
      end select

      if (outputID /= undefined_ID) then
        plastic_nonlocal_output(i,phase_plasticityInstance(p)) = outputs(i)
        plastic_nonlocal_sizePostResult(i,phase_plasticityInstance(p)) = prm%totalNslip
        prm%outputID = [prm%outputID , outputID]
      endif

    enddo

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NofMyPhase   = count(material_phaseAt==p) * discretization_nIP
    sizeDotState =     size([   'rhoSglEdgePosMobile   ','rhoSglEdgeNegMobile   ', &
                                'rhoSglScrewPosMobile  ','rhoSglScrewNegMobile  ', &
                                'rhoSglEdgePosImmobile ','rhoSglEdgeNegImmobile ', &
                                'rhoSglScrewPosImmobile','rhoSglScrewNegImmobile', &
                                'rhoDipEdge            ','rhoDipScrew           ', &
                                'accumulatedshear      ' ]) * prm%totalNslip                        !< "basic" microstructural state variables that are independent from other state variables
    sizeDependentState = size([ 'rhoForest   ']) * prm%totalNslip                                   !< microstructural state variables that depend on other state variables
    sizeState          = sizeDotState + sizeDependentState &
                       + size([ 'velocityEdgePos     ','velocityEdgeNeg     ', & 
                                'velocityScrewPos    ','velocityScrewNeg    ', &
                                'maxDipoleHeightEdge ','maxDipoleHeightScrew' ]) * prm%totalNslip   !< other dependent state variables that are not updated by microstructure
    sizeDeltaState            = sizeDotState
    
    call material_allocatePlasticState(p,NofMyPhase,sizeState,sizeDotState,sizeDeltaState, &
                                       prm%totalNslip,0,0)    
    plasticState(p)%nonlocal = .true.
    plasticState(p)%offsetDeltaState = 0                                                            ! ToDo: state structure does not follow convention
    plasticState(p)%sizePostResults = sum(plastic_nonlocal_sizePostResult(:,phase_plasticityInstance(p)))
    
    totalNslip(phase_plasticityInstance(p)) =  prm%totalNslip
    
    stt%rho => plasticState(p)%state                              (0*prm%totalNslip+1:10*prm%totalNslip,:)
    dot%rho => plasticState(p)%dotState                           (0*prm%totalNslip+1:10*prm%totalNslip,:)
    del%rho => plasticState(p)%deltaState                         (0*prm%totalNslip+1:10*prm%totalNslip,:)
    plasticState(p)%aTolState(1:10*prm%totalNslip) = prm%aTolRho
    
      stt%rhoSgl => plasticState(p)%state                         (0*prm%totalNslip+1: 8*prm%totalNslip,:)
      dot%rhoSgl => plasticState(p)%dotState                      (0*prm%totalNslip+1: 8*prm%totalNslip,:)
      del%rhoSgl => plasticState(p)%deltaState                    (0*prm%totalNslip+1: 8*prm%totalNslip,:)
      
        stt%rhoSglMobile => plasticState(p)%state                 (0*prm%totalNslip+1: 4*prm%totalNslip,:)
        dot%rhoSglMobile => plasticState(p)%dotState              (0*prm%totalNslip+1: 4*prm%totalNslip,:)
        del%rhoSglMobile => plasticState(p)%deltaState            (0*prm%totalNslip+1: 4*prm%totalNslip,:)
          
            stt%rho_sgl_mob_edg_pos => plasticState(p)%state      (0*prm%totalNslip+1: 1*prm%totalNslip,:)
            dot%rho_sgl_mob_edg_pos => plasticState(p)%dotState   (0*prm%totalNslip+1: 1*prm%totalNslip,:)
            del%rho_sgl_mob_edg_pos => plasticState(p)%deltaState (0*prm%totalNslip+1: 1*prm%totalNslip,:)
            
            stt%rho_sgl_mob_edg_neg => plasticState(p)%state      (1*prm%totalNslip+1: 2*prm%totalNslip,:)
            dot%rho_sgl_mob_edg_neg => plasticState(p)%dotState   (1*prm%totalNslip+1: 2*prm%totalNslip,:)
            del%rho_sgl_mob_edg_neg => plasticState(p)%deltaState (1*prm%totalNslip+1: 2*prm%totalNslip,:)
          
            stt%rho_sgl_mob_scr_pos => plasticState(p)%state     (2*prm%totalNslip+1: 3*prm%totalNslip,:)
            dot%rho_sgl_mob_scr_pos => plasticState(p)%dotState  (2*prm%totalNslip+1: 3*prm%totalNslip,:)
            del%rho_sgl_mob_scr_pos => plasticState(p)%deltaState  (2*prm%totalNslip+1: 3*prm%totalNslip,:)

            stt%rho_sgl_mob_scr_neg => plasticState(p)%state     (3*prm%totalNslip+1: 4*prm%totalNslip,:)
            dot%rho_sgl_mob_scr_neg => plasticState(p)%dotState  (3*prm%totalNslip+1: 4*prm%totalNslip,:)
            del%rho_sgl_mob_scr_neg => plasticState(p)%deltaState  (3*prm%totalNslip+1: 4*prm%totalNslip,:)
            
        stt%rhoSglImmobile => plasticState(p)%state               (4*prm%totalNslip+1: 8*prm%totalNslip,:)
        dot%rhoSglImmobile => plasticState(p)%dotState            (4*prm%totalNslip+1: 8*prm%totalNslip,:)
        del%rhoSglImmobile => plasticState(p)%deltaState            (4*prm%totalNslip+1: 8*prm%totalNslip,:)
          
            stt%rho_sgl_imm_edg_pos => plasticState(p)%state    (4*prm%totalNslip+1: 5*prm%totalNslip,:)
            dot%rho_sgl_imm_edg_pos => plasticState(p)%dotState (4*prm%totalNslip+1: 5*prm%totalNslip,:)
            del%rho_sgl_imm_edg_pos => plasticState(p)%deltaState (4*prm%totalNslip+1: 5*prm%totalNslip,:)
            
            stt%rho_sgl_imm_edg_neg => plasticState(p)%state    (5*prm%totalNslip+1: 6*prm%totalNslip,:)
            dot%rho_sgl_imm_edg_neg => plasticState(p)%dotState (5*prm%totalNslip+1: 6*prm%totalNslip,:)
            del%rho_sgl_imm_edg_neg => plasticState(p)%deltaState (5*prm%totalNslip+1: 6*prm%totalNslip,:)
            
            stt%rho_sgl_imm_scr_pos => plasticState(p)%state   (6*prm%totalNslip+1: 7*prm%totalNslip,:)
            dot%rho_sgl_imm_scr_pos => plasticState(p)%dotState(6*prm%totalNslip+1: 7*prm%totalNslip,:)
            del%rho_sgl_imm_scr_pos => plasticState(p)%deltaState(6*prm%totalNslip+1: 7*prm%totalNslip,:)
            
            stt%rho_sgl_imm_scr_neg => plasticState(p)%state   (7*prm%totalNslip+1: 8*prm%totalNslip,:)
            dot%rho_sgl_imm_scr_neg => plasticState(p)%dotState(7*prm%totalNslip+1: 8*prm%totalNslip,:)
            del%rho_sgl_imm_scr_neg => plasticState(p)%deltaState(7*prm%totalNslip+1: 8*prm%totalNslip,:)
    
      stt%rhoDip => plasticState(p)%state                         (8*prm%totalNslip+1:10*prm%totalNslip,:)
      dot%rhoDip => plasticState(p)%dotState                      (8*prm%totalNslip+1:10*prm%totalNslip,:)
      del%rhoDip => plasticState(p)%deltaState                      (8*prm%totalNslip+1:10*prm%totalNslip,:)
      
        stt%rho_dip_edg => plasticState(p)%state                   (8*prm%totalNslip+1: 9*prm%totalNslip,:)
        dot%rho_dip_edg => plasticState(p)%dotState                (8*prm%totalNslip+1: 9*prm%totalNslip,:)
        del%rho_dip_edg => plasticState(p)%deltaState              (8*prm%totalNslip+1: 9*prm%totalNslip,:)
        
        stt%rho_dip_scr => plasticState(p)%state                  (9*prm%totalNslip+1:10*prm%totalNslip,:)
        dot%rho_dip_scr => plasticState(p)%dotState               (9*prm%totalNslip+1:10*prm%totalNslip,:)
        del%rho_dip_scr => plasticState(p)%deltaState               (9*prm%totalNslip+1:10*prm%totalNslip,:)
     
    stt%accumulatedshear => plasticState(p)%state       (10*prm%totalNslip + 1:11*prm%totalNslip ,1:NofMyPhase)
    dot%accumulatedshear => plasticState(p)%dotState   (10*prm%totalNslip + 1:11*prm%totalNslip ,1:NofMyPhase)
    del%accumulatedshear => plasticState(p)%deltaState  (10*prm%totalNslip + 1:11*prm%totalNslip ,1:NofMyPhase)
    plasticState(p)%aTolState(10*prm%totalNslip + 1:11*prm%totalNslip )  = prm%aTolShear
    plasticState(p)%slipRate => plasticState(p)%dotState(10*prm%totalNslip + 1:11*prm%totalNslip ,1:NofMyPhase)
    plasticState(p)%accumulatedSlip => plasticState(p)%state (10*prm%totalNslip + 1:11*prm%totalNslip ,1:NofMyPhase)
     
    stt%rho_forest => plasticState(p)%state       (11*prm%totalNslip + 1:12*prm%totalNslip ,1:NofMyPhase)
    stt%v          => plasticState(p)%state       (12*prm%totalNslip + 1:16*prm%totalNslip ,1:NofMyPhase)

    allocate(dst%tau_Threshold(prm%totalNslip,NofMyPhase),source=0.0_pReal)
    allocate(dst%tau_Back(prm%totalNslip,NofMyPhase),source=0.0_pReal)
    end associate
  

    if (NofMyPhase > 0) call stateInit(p,NofMyPhase)
    plasticState(p)%state0 = plasticState(p)%state

  enddo
 
  allocate(compatibility(2,maxval(totalNslip),maxval(totalNslip),nIPneighbors,&
                         discretization_nIP,discretization_nElem), source=0.0_pReal)
                                            
! BEGIN DEPRECATED----------------------------------------------------------------------------------
  allocate(iRhoU(maxval(totalNslip),4,maxNinstances), source=0)
  allocate(iRhoB(maxval(totalNslip),4,maxNinstances), source=0)
  allocate(iRhoD(maxval(totalNslip),2,maxNinstances), source=0)
  allocate(iV(maxval(totalNslip),4,maxNinstances),    source=0)
  allocate(iD(maxval(totalNslip),2,maxNinstances),    source=0)

  initializeInstances: do p = 1, size(phase_plasticity)
    NofMyPhase = count(material_phaseAt==p) * discretization_nIP
    myPhase2: if (phase_plasticity(p) == PLASTICITY_NONLOCAL_ID) then

      !*** determine indices to state array

      l = 0
      do t = 1,4
        do s = 1,param(phase_plasticityInstance(p))%totalNslip
          l = l + 1
          iRhoU(s,t,phase_plasticityInstance(p)) = l
        enddo
      enddo
      do t = 1,4
        do s = 1,param(phase_plasticityInstance(p))%totalNslip
          l = l + 1
          iRhoB(s,t,phase_plasticityInstance(p)) = l
        enddo
      enddo
      do c = 1,2
        do s = 1,param(phase_plasticityInstance(p))%totalNslip
          l = l + 1
          iRhoD(s,c,phase_plasticityInstance(p)) = l
        enddo
      enddo
      
      l = l + param(phase_plasticityInstance(p))%totalNslip ! shear(rates)
      l = l + param(phase_plasticityInstance(p))%totalNslip ! rho_forest

      do t = 1,4
        do s = 1,param(phase_plasticityInstance(p))%totalNslip
          l = l + 1
          iV(s,t,phase_plasticityInstance(p)) = l
        enddo
      enddo
      do c = 1,2
        do s = 1,param(phase_plasticityInstance(p))%totalNslip
          l = l + 1
          iD(s,c,phase_plasticityInstance(p)) = l
        enddo
      enddo
      if (iD(param(phase_plasticityInstance(p))%totalNslip,2,phase_plasticityInstance(p)) /= plasticState(p)%sizeState) &
        call IO_error(0, ext_msg = 'state indices not properly set ('//PLASTICITY_NONLOCAL_label//')')
      
      
    endif myPhase2
   
  enddo initializeInstances
! END DEPRECATED------------------------------------------------------------------------------------

   
  contains
  !--------------------------------------------------------------------------------------------------
  !> @brief populates the initial dislocation density
  !--------------------------------------------------------------------------------------------------
  subroutine stateInit(phase,NofMyPhase)
      
   integer,intent(in) ::&
     phase, &
     NofMyPhase
   integer :: &
     e, &
     i, &
     f, &
     from, &
     upto, &
     s, &
     instance, &
     phasemember
   real(pReal), dimension(2) :: &
     noise, &
     rnd   
   real(pReal) :: &
     meanDensity, &
     totalVolume, &
     densityBinning, &
     minimumIpVolume
   real(pReal), dimension(NofMyPhase) :: &
     volume
  
  
   instance = phase_plasticityInstance(phase)
   associate(prm => param(instance), stt => state(instance))
  
   ! randomly distribute dislocation segments on random slip system and of random type in the volume 
   if (prm%rhoSglRandom > 0.0_pReal) then
  
    ! get the total volume of the instance
    do e = 1,discretization_nElem
      do i = 1,discretization_nIP
        if (material_phaseAt(1,e) == phase) volume(material_phasememberAt(1,i,e)) = IPvolume(i,e)
      enddo
    enddo
    totalVolume = sum(volume)
    minimumIPVolume = minval(volume)
    densityBinning = prm%rhoSglRandomBinning / minimumIpVolume ** (2.0_pReal / 3.0_pReal)
  
    ! subsequently fill random ips with dislocation segments until we reach the desired overall density
    meanDensity = 0.0_pReal
    do while(meanDensity < prm%rhoSglRandom)
      call random_number(rnd)
      phasemember = nint(rnd(1)*real(NofMyPhase,pReal) + 0.5_pReal)
      s           = nint(rnd(2)*real(prm%totalNslip,pReal)*4.0_pReal + 0.5_pReal)
      meanDensity = meanDensity + densityBinning * volume(phasemember) / totalVolume
      stt%rhoSglMobile(s,phasemember) = densityBinning
    enddo
   ! homogeneous distribution of density with some noise
   else
     do e = 1, NofMyPhase
       do f = 1,size(prm%Nslip,1)
         from = 1 + sum(prm%Nslip(1:f-1))
         upto = sum(prm%Nslip(1:f))
         do s = from,upto
           noise = [math_sampleGaussVar(0.0_pReal, prm%rhoSglScatter), &
                    math_sampleGaussVar(0.0_pReal, prm%rhoSglScatter)]
           stt%rho_sgl_mob_edg_pos(s,e) = prm%rhoSglEdgePos0(f)  + noise(1)
           stt%rho_sgl_mob_edg_neg(s,e) = prm%rhoSglEdgeNeg0(f)  + noise(1)
           stt%rho_sgl_mob_scr_pos(s,e) = prm%rhoSglScrewPos0(f) + noise(2)
           stt%rho_sgl_mob_scr_neg(s,e) = prm%rhoSglScrewNeg0(f) + noise(2)
         enddo
         stt%rho_dip_edg(from:upto,e)   = prm%rhoDipEdge0(f)
         stt%rho_dip_scr(from:upto,e)   = prm%rhoDipScrew0(f)
       enddo
     enddo
    endif
    
    end associate
  
  end subroutine stateInit

end subroutine plastic_nonlocal_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates quantities characterizing the microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_dependentState(Fe, Fp, ip, el)
  
  integer, intent(in) :: &
    ip, &
    el
  real(pReal), dimension(3,3), intent(in) :: &
    Fe, &
    Fp
  
  integer :: &
    ph, &                                                                                           !< phase
    of, &                                                                                           !< offset
    no, &                                                                                           !< neighbor offset
    ns, &
    neighbor_el, &                                                                                  ! element number of neighboring material point
    neighbor_ip, &                                                                                  ! integration point of neighboring material point
    instance, &                                                                                     ! my instance of this plasticity
    neighbor_instance, &                                                                            ! instance of this plasticity of neighboring material point
    c, &                                                                                            ! index of dilsocation character (edge, screw)
    s, &                                                                                            ! slip system index
    dir, &
    n
  real(pReal) :: &
    FVsize, &
    correction, &
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
  real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    rhoExcess
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    rho_edg_delta, &
    rho_scr_delta
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),10) :: &
    rho, &
    rho_neighbor
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))), &
                         totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    myInteractionMatrix                                                                             ! corrected slip interaction matrix
                                  
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),nIPneighbors) :: &
    rho_edg_delta_neighbor, &
    rho_scr_delta_neighbor
  real(pReal), dimension(2,maxval(totalNslip),nIPneighbors) :: &
    neighbor_rhoExcess, &                                                                           ! excess density at neighboring material point
    neighbor_rhoTotal                                                                               ! total density at neighboring material point
  real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),2) :: &
    m                                                                                               ! direction of dislocation motion
  
  ph = material_phaseAt(1,el)
  of = material_phasememberAt(1,ip,el)
  instance = phase_plasticityInstance(ph)

  associate(prm => param(instance),dst => microstructure(instance), stt => state(instance))
  
  ns = prm%totalNslip
  
  rho = getRho(instance,of,ip,el)
  
  stt%rho_forest(:,of) = matmul(prm%forestProjection_Edge, sum(abs(rho(:,edg)),2)) &
                       + matmul(prm%forestProjection_Screw,sum(abs(rho(:,scr)),2))        
  
  
  ! coefficients are corrected for the line tension effect 
  ! (see Kubin,Devincre,Hoc; 2008; Modeling dislocation storage rates and mean free paths in face-centered cubic crystals)
  if (lattice_structure(ph) ==  LATTICE_bcc_ID .or. lattice_structure(ph) == LATTICE_fcc_ID) then     ! only fcc and bcc
    do s = 1,ns
      correction = (  1.0_pReal - prm%linetensionEffect &
                    + prm%linetensionEffect &
                    * log(0.35_pReal * prm%burgers(s) * sqrt(max(stt%rho_forest(s,of),prm%significantRho))) &
                    / log(0.35_pReal * prm%burgers(s) * 1e6_pReal)) ** 2.0_pReal
      myInteractionMatrix(1:ns,s) = correction * prm%interactionSlipSlip(1:ns,s) 
    enddo
  else
    myInteractionMatrix = prm%interactionSlipSlip
  endif

  forall (s = 1:ns) &
    dst%tau_threshold(s,of) = prm%mu * prm%burgers(s) &
                  * sqrt(dot_product(sum(abs(rho),2), myInteractionMatrix(1:ns,s)))


!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local plasticity

  !#################################################################################################
  ! ToDo: MD: this is most likely only correct for F_i = I
  !#################################################################################################


  if (.not. phase_localPlasticity(ph) .and. prm%shortRangeStressCorrection) then
    invFe = math_inv33(Fe)
    invFp = math_inv33(Fp)
  
    rho_edg_delta = rho(:,mob_edg_pos) - rho(:,mob_edg_neg)
    rho_scr_delta = rho(:,mob_scr_pos) - rho(:,mob_scr_neg)
  
    rhoExcess(1,1:ns) = rho_edg_delta
    rhoExcess(2,1:ns) = rho_scr_delta
  
    FVsize = IPvolume(ip,el) ** (1.0_pReal/3.0_pReal)
    
    !* loop through my neighborhood and get the connection vectors (in lattice frame) and the excess densities
    
    nRealNeighbors = 0.0_pReal
    neighbor_rhoTotal = 0.0_pReal
    do n = 1,nIPneighbors
      neighbor_el = IPneighborhood(1,n,ip,el)
      neighbor_ip = IPneighborhood(2,n,ip,el)
      no = material_phasememberAt(1,neighbor_ip,neighbor_el)
      if (neighbor_el > 0 .and. neighbor_ip > 0) then
        neighbor_instance = phase_plasticityInstance(material_phaseAt(1,neighbor_el))
        if (neighbor_instance == instance) then
  
            nRealNeighbors = nRealNeighbors + 1.0_pReal
            rho_neighbor = getRho(instance,no,neighbor_ip,neighbor_el)
  
            rho_edg_delta_neighbor(:,n) = rho_neighbor(:,mob_edg_pos) - rho_neighbor(:,mob_edg_neg)
            rho_scr_delta_neighbor(:,n) = rho_neighbor(:,mob_scr_pos) - rho_neighbor(:,mob_scr_neg)
             
            neighbor_rhoTotal(1,:,n) = sum(abs(rho_neighbor(:,edg)),2)
            neighbor_rhoTotal(2,:,n) = sum(abs(rho_neighbor(:,scr)),2)
  
            connection_latticeConf(1:3,n) = matmul(invFe, discretization_IPcoords(1:3,neighbor_el+neighbor_ip-1) &
                                          - discretization_IPcoords(1:3,el+neighbor_ip-1))
            normal_latticeConf = matmul(transpose(invFp), IPareaNormal(1:3,n,ip,el))
            if (math_inner(normal_latticeConf,connection_latticeConf(1:3,n)) < 0.0_pReal) &        ! neighboring connection points in opposite direction to face normal: must be periodic image
              connection_latticeConf(1:3,n) = normal_latticeConf * IPvolume(ip,el)/IParea(n,ip,el) ! instead take the surface normal scaled with the diameter of the cell
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
    m(1:3,1:ns,1) =  prm%slip_direction
    m(1:3,1:ns,2) = -prm%slip_transverse
  
    do s = 1,ns
      
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
      where(rhoTotal > 0.0_pReal) &
        rhoExcessGradient_over_rho = rhoExcessGradient / rhoTotal
      
        ! ... gives the local stress correction when multiplied with a factor
      dst%tau_back(s,of) = - prm%mu * prm%burgers(s) / (2.0_pReal * pi) &
                         * (rhoExcessGradient_over_rho(1) / (1.0_pReal - prm%nu) &
                         + rhoExcessGradient_over_rho(2))
  
    enddo
  endif

#ifdef DEBUG
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0 &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0)) then
    write(6,'(/,a,i8,1x,i2,1x,i1,/)') '<< CONST >> nonlocal_microstructure at el ip ',el,ip
    write(6,'(a,/,12x,12(e10.3,1x))') '<< CONST >> rhoForest', stt%rho_forest(:,of)
    write(6,'(a,/,12x,12(f10.5,1x))') '<< CONST >> tauThreshold / MPa', dst%tau_threshold(:,of)*1e-6
    write(6,'(a,/,12x,12(f10.5,1x),/)') '<< CONST >> tauBack / MPa', dst%tau_back(:,of)*1e-6
  endif
#endif

 end associate

end subroutine plastic_nonlocal_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief calculates kinetics 
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_kinetics(v, dv_dtau, dv_dtauNS, tau, tauNS, &
                                     tauThreshold, c, Temperature, instance, of)
  integer, intent(in) :: &
    c, &                                                                                            !< dislocation character (1:edge, 2:screw)
    instance, of
  real(pReal), intent(in) :: &
    Temperature                                                                                     !< temperature
  real(pReal), dimension(param(instance)%totalNslip), intent(in) :: &
    tau, &                                                                                          !< resolved external shear stress (without non Schmid effects)
    tauNS, &                                                                                        !< resolved external shear stress (including non Schmid effects)
    tauThreshold                                                                                    !< threshold shear stress
  
  real(pReal), dimension(param(instance)%totalNslip), intent(out) ::  &
    v, &                                                                                            !< velocity
    dv_dtau, &                                                                                      !< velocity derivative with respect to resolved shear stress (without non Schmid contributions)
    dv_dtauNS                                                                                       !< velocity derivative with respect to resolved shear stress (including non Schmid contributions)
  
  integer :: &
    ns, &                                                                                           !< short notation for the total number of active slip systems
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
  
  associate(prm => param(instance))
  ns = prm%totalNslip
  v = 0.0_pReal
  dv_dtau = 0.0_pReal
  dv_dtauNS = 0.0_pReal
  
  
  if (Temperature > 0.0_pReal) then
    do s = 1,ns
      if (abs(tau(s)) > tauThreshold(s)) then
  
        !* Peierls contribution
        !* Effective stress includes non Schmid constributions
        !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
        
        tauEff = max(0.0_pReal, abs(tauNS(s)) - tauThreshold(s))                                    ! ensure that the effective stress is positive
        meanfreepath_P = prm%burgers(s)
        jumpWidth_P = prm%burgers(s)
        activationLength_P = prm%doublekinkwidth *prm%burgers(s)
        activationVolume_P = activationLength_P * jumpWidth_P * prm%burgers(s)
        criticalStress_P = prm%peierlsStress(s,c)
        activationEnergy_P = criticalStress_P * activationVolume_P
        tauRel_P = min(1.0_pReal, tauEff / criticalStress_P)                                        ! ensure that the activation probability cannot become greater than one
        tPeierls = 1.0_pReal / prm%fattack &
                 * exp(activationEnergy_P / (KB * Temperature) &
                       * (1.0_pReal - tauRel_P**prm%p)**prm%q)
        if (tauEff < criticalStress_P) then
          dtPeierls_dtau = tPeierls * prm%p * prm%q * activationVolume_P / (KB * Temperature) &
                         * (1.0_pReal - tauRel_P**prm%p)**(prm%q-1.0_pReal) &
                                      * tauRel_P**(prm%p-1.0_pReal) 
        else
          dtPeierls_dtau = 0.0_pReal
        endif
  
  
        !* Contribution from solid solution strengthening
        !* The derivative only gives absolute values; the correct sign is taken care of in the formula for the derivative of the velocity
  
        tauEff = abs(tau(s)) - tauThreshold(s)
        meanfreepath_S = prm%burgers(s) / sqrt(prm%solidSolutionConcentration)
        jumpWidth_S = prm%solidSolutionSize * prm%burgers(s)
        activationLength_S = prm%burgers(s) / sqrt(prm%solidSolutionConcentration)
        activationVolume_S = activationLength_S * jumpWidth_S * prm%burgers(s)
        activationEnergy_S = prm%solidSolutionEnergy
        criticalStress_S = activationEnergy_S / activationVolume_S
        tauRel_S = min(1.0_pReal, tauEff / criticalStress_S)                                        ! ensure that the activation probability cannot become greater than one
        tSolidSolution = 1.0_pReal /  prm%fattack &
                       * exp(activationEnergy_S / (KB * Temperature) &
                             * (1.0_pReal - tauRel_S**prm%p)**prm%q)
        if (tauEff < criticalStress_S) then
          dtSolidSolution_dtau = tSolidSolution * prm%p * prm%q &
                               * activationVolume_S / (KB * Temperature) &
                               * (1.0_pReal - tauRel_S**prm%p)**(prm%q-1.0_pReal) &
                                              * tauRel_S**(prm%p-1.0_pReal) 
        else
          dtSolidSolution_dtau = 0.0_pReal
        endif
  
  
        !* viscous glide velocity
        
        tauEff = abs(tau(s)) - tauThreshold(s)
        mobility = prm%burgers(s) / prm%viscosity
        vViscous = mobility * tauEff
  
  
        !* Mean velocity results from waiting time at peierls barriers and solid solution obstacles with respective meanfreepath of 
        !* free flight at glide velocity in between.       
        !* adopt sign from resolved stress
  
        v(s) = sign(1.0_pReal,tau(s)) &
             / (tPeierls / meanfreepath_P + tSolidSolution / meanfreepath_S + 1.0_pReal / vViscous)
        dv_dtau(s) = v(s) * v(s) * (dtSolidSolution_dtau / meanfreepath_S &
                                   + mobility / (vViscous * vViscous))
        dv_dtauNS(s) = v(s) * v(s) * dtPeierls_dtau / meanfreepath_P 
      endif
    enddo
  endif
    

#ifdef DEBUGTODO
  write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tauThreshold / MPa', tauThreshold * 1e-6_pReal
  write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tau / MPa', tau * 1e-6_pReal
  write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> tauNS / MPa', tauNS * 1e-6_pReal
  write(6,'(a,/,12x,12(f12.5,1x))') '<< CONST >> v / mm/s', v * 1e3
  write(6,'(a,/,12x,12(e12.5,1x))') '<< CONST >> dv_dtau', dv_dtau
  write(6,'(a,/,12x,12(e12.5,1x))') '<< CONST >> dv_dtauNS', dv_dtauNS
#endif

  end associate

end subroutine plastic_nonlocal_kinetics


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_LpAndItsTangent(Lp, dLp_dMp, &
                                            Mp, Temperature, volume, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< current integration point
    el                                                                                              !< current element number
  real(pReal), intent(in) :: &
    Temperature, &                                                                                  !< temperature
    volume                                                                                          !< volume of the materialpoint
  real(pReal), dimension(3,3), intent(in) :: &
    Mp
  real(pReal), dimension(3,3), intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to Tstar (9x9 matrix)
  
  
  integer :: &
    instance, &                                                                                     !< current instance of this plasticity
    ns, &                                                                                           !< short notation for the total number of active slip systems
    i, &
    j, &
    k, &
    l, &
    ph, &                                                                                           !phase number
    of, &                                                                                           !offset
    t, &                                                                                            !< dislocation type
    s                                                                                               !< index of my current slip system
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),8) :: &
    rhoSgl                                                                                          !< single dislocation densities (including blocked)
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),10) :: &
    rho
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),4) :: &
    v, &                                                                                            !< velocity
    tauNS, &                                                                                        !< resolved shear stress including non Schmid and backstress terms
    dv_dtau, &                                                                                      !< velocity derivative with respect to the shear stress
    dv_dtauNS                                                                                       !< velocity derivative with respect to the shear stress
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    tau, &                                                                                          !< resolved shear stress including backstress terms
    gdotTotal                                                                                       !< shear rate
  
  !*** shortcut for mapping
  ph = material_phaseAt(1,el)
  of = material_phasememberAt(1,ip,el)
  
  instance = phase_plasticityInstance(ph)
  associate(prm => param(instance),dst=>microstructure(instance))
  ns = prm%totalNslip
  
  !*** shortcut to state variables 
  rho = getRho(instance,of,ip,el)
  rhoSgl = rho(:,sgl)
  
  
  !*** get resolved shear stress
  !*** for screws possible non-schmid contributions are also taken into account
  do s = 1,ns
    tau(s) = math_mul33xx33(Mp, prm%Schmid(1:3,1:3,s))
    tauNS(s,1) = tau(s)
    tauNS(s,2) = tau(s)
    if (tau(s) > 0.0_pReal) then
      tauNS(s,3) = math_mul33xx33(Mp, +prm%nonSchmid_pos(1:3,1:3,s))
      tauNS(s,4) = math_mul33xx33(Mp, -prm%nonSchmid_neg(1:3,1:3,s))
    else
      tauNS(s,3) = math_mul33xx33(Mp, +prm%nonSchmid_neg(1:3,1:3,s))
      tauNS(s,4) = math_mul33xx33(Mp, -prm%nonSchmid_pos(1:3,1:3,s))
    endif
  enddo
  forall (t = 1:4) &
    tauNS(1:ns,t) = tauNS(1:ns,t) + dst%tau_back(:,of)
  tau = tau + dst%tau_back(:,of)
  
  
  !*** get dislocation velocity and its tangent and store the velocity in the state array
  
  ! edges 
  call plastic_nonlocal_kinetics(v(1:ns,1), dv_dtau(1:ns,1), dv_dtauNS(1:ns,1), &
                                      tau(1:ns), tauNS(1:ns,1), dst%tau_Threshold(1:ns,of), &
                                        1, Temperature, instance, of)
  v(1:ns,2) = v(1:ns,1)
  dv_dtau(1:ns,2) = dv_dtau(1:ns,1)
  dv_dtauNS(1:ns,2) = dv_dtauNS(1:ns,1)
  
  !screws
  if (size(prm%nonSchmidCoeff) == 0) then
    forall(t = 3:4)
      v(1:ns,t) = v(1:ns,1)
      dv_dtau(1:ns,t) = dv_dtau(1:ns,1)
      dv_dtauNS(1:ns,t) = dv_dtauNS(1:ns,1)
    endforall
  else
    do t = 3,4
      call plastic_nonlocal_kinetics(v(1:ns,t), dv_dtau(1:ns,t), dv_dtauNS(1:ns,t), &
                                          tau(1:ns), tauNS(1:ns,t), dst%tau_Threshold(1:ns,of), &
                                            2 , Temperature, instance, of)
    enddo
  endif
  
  
  !*** store velocity in state
  forall (t = 1:4) &
    plasticState(ph)%state(iV(1:ns,t,instance),of) = v(1:ns,t)
  
  !*** Bauschinger effect
  forall (s = 1:ns, t = 5:8, rhoSgl(s,t) * v(s,t-4) < 0.0_pReal) &
    rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(rhoSgl(s,t))
  
  
  gdotTotal = sum(rhoSgl(1:ns,1:4) * v, 2) * prm%burgers(1:ns)
  
  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal
  
  do s = 1,ns
    Lp = Lp + gdotTotal(s) * prm%Schmid(1:3,1:3,s)
    forall (i=1:3,j=1:3,k=1:3,l=1:3) &
      dLp_dMp(i,j,k,l) = dLp_dMp(i,j,k,l) &
          + prm%Schmid(i,j,s) * prm%Schmid(k,l,s) &
          * sum(rhoSgl(s,1:4) * dv_dtau(s,1:4)) * prm%burgers(s) &
          + prm%Schmid(i,j,s) &
          * ( prm%nonSchmid_pos(k,l,s) * rhoSgl(s,3) * dv_dtauNS(s,3) & 
            - prm%nonSchmid_neg(k,l,s) * rhoSgl(s,4) * dv_dtauNS(s,4))  * prm%burgers(s)
  enddo
  
  
  end associate
 
end subroutine plastic_nonlocal_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_deltaState(Mp,ip,el)
  
  integer, intent(in) :: &
    ip, &
    el
  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< MandelStress

  integer :: &
    ph, &                                                                                           !< phase
    of, &                                                                                           !< offset
    instance, &                                                                                     ! current instance of this plasticity
    ns, &                                                                                           ! short notation for the total number of active slip systems
    c, &                                                                                            ! character of dislocation
    t, &                                                                                            ! type of dislocation
    s                                                                                               ! index of my current slip system
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),10) :: &
    deltaRhoRemobilization, &                                                                       ! density increment by remobilization
    deltaRhoDipole2SingleStress                                                                     ! density increment by dipole dissociation (by stress change)
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),10) :: &
    rho                                                                                             ! current  dislocation densities
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),4) :: &
    v                                                                                               ! dislocation glide velocity
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    tau                                                                                             ! current resolved shear stress
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),2) :: &
    rhoDip, &                                                                                       ! current dipole dislocation densities (screw and edge dipoles)
    dUpper, &                                                                                       ! current maximum stable dipole distance for edges and screws
    dUpperOld, &                                                                                    ! old maximum stable dipole distance for edges and screws
    deltaDUpper                                                                                     ! change in maximum stable dipole distance for edges and screws
  
   ph = material_phaseAt(1,el)
   of = material_phasememberAt(1,ip,el)
   instance = phase_plasticityInstance(ph)

   associate(prm => param(instance),dst => microstructure(instance),del => deltaState(instance))
   ns = totalNslip(instance)
   
  !*** shortcut to state variables 
  forall (s = 1:ns, t = 1:4) &
    v(s,t) = plasticState(ph)%state(iV(s,t,instance),of)
  forall (s = 1:ns, c = 1:2) &
    dUpperOld(s,c) = plasticState(ph)%state(iD(s,c,instance),of)
  
  rho =  getRho(instance,of,ip,el)
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
  do s = 1,prm%totalNslip
    tau(s) = math_mul33xx33(Mp, prm%Schmid(1:3,1:3,s)) +dst%tau_back(s,of)
    if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
  enddo
  
  dUpper(1:ns,1) = prm%mu * prm%burgers/(8.0_pReal * PI * (1.0_pReal - prm%nu) * abs(tau))
  dUpper(1:ns,2) = prm%mu * prm%burgers/(4.0_pReal * PI * abs(tau))
  
  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(1:ns,1) = min(1.0_pReal/sqrt(sum(abs(rho(:,edg)),2)),dUpper(1:ns,1))
  
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(1:ns,2) = min(1.0_pReal/sqrt(sum(abs(rho(:,scr)),2)),dUpper(1:ns,2))
  
  dUpper = max(dUpper,prm%minDipoleHeight)
  deltaDUpper = dUpper - dUpperOld
  
  
  !*** dissociation by stress increase
  deltaRhoDipole2SingleStress = 0.0_pReal
  forall (c=1:2, s=1:ns, deltaDUpper(s,c) < 0.0_pReal .and. &
                                          dNeq0(dUpperOld(s,c) - prm%minDipoleHeight(s,c))) &
    deltaRhoDipole2SingleStress(s,8+c) = rhoDip(s,c) * deltaDUpper(s,c) &
                                             / (dUpperOld(s,c) - prm%minDipoleHeight(s,c))
  
  forall (t=1:4) &
    deltaRhoDipole2SingleStress(1:ns,t) = -0.5_pReal &
                                    * deltaRhoDipole2SingleStress(1:ns,(t-1)/2+9)
  
  forall (s = 1:ns, c = 1:2) &
    plasticState(ph)%state(iD(s,c,instance),of) = dUpper(s,c) 
  
  plasticState(ph)%deltaState(:,of) = 0.0_pReal
  del%rho(:,of) = reshape(deltaRhoRemobilization + deltaRhoDipole2SingleStress, [10*ns])
  
#ifdef DEBUG
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0 &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0 )) then
    write(6,'(a,/,8(12x,12(e12.5,1x),/))') '<< CONST >> dislocation remobilization', deltaRhoRemobilization(1:ns,1:8)
    write(6,'(a,/,10(12x,12(e12.5,1x),/),/)') '<< CONST >> dipole dissociation by stress increase', deltaRhoDipole2SingleStress
  endif
#endif

  end associate

end subroutine plastic_nonlocal_deltaState


!---------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!---------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_dotState(Mp, Fe, Fp, Temperature, &
                                     timestep,ip,el)
  
  integer, intent(in) :: &
    ip, &                                                                                           !< current integration point
    el                                                                                              !< current element number
  real(pReal), intent(in) :: &
    Temperature, &                                                                                  !< temperature
    timestep                                                                                        !< substepped crystallite time increment
  real(pReal), dimension(3,3), intent(in) ::&
    Mp                                                                                              !< MandelStress
  real(pReal), dimension(3,3,homogenization_maxNgrains,discretization_nIP,discretization_nElem), intent(in) :: &
    Fe, &                                                                                           !< elastic deformation gradient
    Fp                                                                                              !< plastic deformation gradient
  
  integer ::  &
    ph, &                              
    instance, &                                                                                     !< current instance of this plasticity
    neighbor_instance, &                                                                            !< instance of my neighbor's plasticity
    ns, &                                                                                           !< short notation for the total number of active slip systems
    c, &                                                                                            !< character of dislocation
    n, &                                                                                            !< index of my current neighbor
    neighbor_el, &                                                                                  !< element number of my neighbor
    neighbor_ip, &                                                                                  !< integration point of my neighbor
    neighbor_n, &                                                                                   !< neighbor index pointing to me when looking from my neighbor
    opposite_neighbor, &                                                                            !< index of my opposite neighbor
    opposite_ip, &                                                                                  !< ip of my opposite neighbor
    opposite_el, &                                                                                  !< element index of my opposite neighbor
    opposite_n, &                                                                                   !< neighbor index pointing to me when looking from my opposite neighbor
    t, &                                                                                            !< type of dislocation
    o,&                                                                                             !< offset shortcut
    no,&                                                                                            !< neighbor offset shortcut
    p,&                                                                                             !< phase shortcut
    np,&                                                                                            !< neighbor phase shortcut
    topp, &                                                                                         !< type of dislocation with opposite sign to t
    s                                                                                               !< index of my current slip system
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),10) :: &
    rho, &
    rhoDot, &                                                                                       !< density evolution
    rhoDotMultiplication, &                                                                         !< density evolution by multiplication
    rhoDotFlux, &                                                                                   !< density evolution by flux
    rhoDotSingle2DipoleGlide, &                                                                     !< density evolution by dipole formation (by glide)
    rhoDotAthermalAnnihilation, &                                                                   !< density evolution by athermal annihilation
    rhoDotThermalAnnihilation                                                                       !< density evolution by thermal annihilation
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),8) :: &
    rhoSgl, &                                                                                       !< current single dislocation densities (positive/negative screw and edge without dipoles)
    neighbor_rhoSgl, &                                                                              !< current single dislocation densities of neighboring ip (positive/negative screw and edge without dipoles)
    my_rhoSgl                                                                                       !< single dislocation densities of central ip (positive/negative screw and edge without dipoles)
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),4) :: &
    v, &                                                                                            !< current dislocation glide velocity
    my_v, &                                                                                         !< dislocation glide velocity of central ip
    neighbor_v, &                                                                                   !< dislocation glide velocity of enighboring ip
    gdot                                                                                            !< shear rates
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el)))) :: &
    tau, &                                                                                          !< current resolved shear stress
    vClimb                                                                                          !< climb velocity of edge dipoles
                                              
  real(pReal), dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),2) :: &
    rhoDip, &                                                                                       !< current dipole dislocation densities (screw and edge dipoles)
    dLower, &                                                                                       !< minimum stable dipole distance for edges and screws
    dUpper                                                                                          !< current maximum stable dipole distance for edges and screws
  real(pReal), dimension(3,totalNslip(phase_plasticityInstance(material_phaseAt(1,el))),4) :: &
    m                                                                                               !< direction of dislocation motion
  real(pReal), dimension(3,3) :: &
    my_F, &                                                                                         !< my total deformation gradient
    neighbor_F, &                                                                                   !< total deformation gradient of my neighbor
    my_Fe, &                                                                                        !< my elastic deformation gradient
    neighbor_Fe, &                                                                                  !< elastic deformation gradient of my neighbor
    Favg                                                                                            !< average total deformation gradient of me and my neighbor
  real(pReal), dimension(3) :: &
    normal_neighbor2me, &                                                                           !< interface normal pointing from my neighbor to me in neighbor's lattice configuration
    normal_neighbor2me_defConf, &                                                                   !< interface normal pointing from my neighbor to me in shared deformed configuration
    normal_me2neighbor, &                                                                           !< interface normal pointing from me to my neighbor in my lattice configuration
    normal_me2neighbor_defConf                                                                      !< interface normal pointing from me to my neighbor in shared deformed configuration
  real(pReal) :: &
    area, &                                                                                         !< area of the current interface
    transmissivity, &                                                                               !< overall transmissivity of dislocation flux to neighboring material point
    lineLength, &                                                                                   !< dislocation line length leaving the current interface
    selfDiffusion                                                                                   !< self diffusion
  
  logical :: &
    considerEnteringFlux, &
    considerLeavingFlux
                                              
  p = material_phaseAt(1,el)
  o = material_phasememberAt(1,ip,el)
  
  if (timestep <= 0.0_pReal) then
    plasticState(p)%dotState = 0.0_pReal
    return
  endif
  
  ph = material_phaseAt(1,el)
  instance = phase_plasticityInstance(ph)
  associate(prm => param(instance), &
            dst => microstructure(instance), &
            dot => dotState(instance), &
            stt => state(instance))
  ns = totalNslip(instance)
  
  tau = 0.0_pReal
  gdot = 0.0_pReal
  
  rho = getRho(instance,o,ip,el)
  rhoSgl = rho(:,sgl)
  rhoDip = rho(:,dip)
  
  forall (s = 1:ns, t = 1:4)
    v(s,t) = plasticState(p)%state(iV   (s,t,instance),o)
  endforall
  
  
  !****************************************************************************
  !*** Calculate shear rate
  
  forall (t = 1:4) &
    gdot(1:ns,t) = rhoSgl(1:ns,t) *  prm%burgers(1:ns) * v(1:ns,t)
  
#ifdef DEBUG
    if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0 &
        .and. ((debug_e == el .and. debug_i == ip)&
               .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0 )) then
      write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
      write(6,'(a,/,4(12x,12(e12.5,1x),/))') '<< CONST >> gdot / 1/s',gdot
    endif
#endif
  
  
  
  !****************************************************************************
  !*** calculate limits for stable dipole height
  
  do s = 1,ns   ! loop over slip systems
    tau(s) = math_mul33xx33(Mp, prm%Schmid(1:3,1:3,s)) + dst%tau_back(s,o)
    if (abs(tau(s)) < 1.0e-15_pReal) tau(s) = 1.0e-15_pReal
  enddo
  
  dLower = prm%minDipoleHeight
  dUpper(1:ns,1) = prm%mu * prm%burgers/(8.0_pReal * PI * (1.0_pReal - prm%nu) * abs(tau))
  dUpper(1:ns,2) = prm%mu * prm%burgers/(4.0_pReal * PI * abs(tau))
  
  where(dNeq0(sqrt(sum(abs(rho(:,edg)),2)))) &
    dUpper(1:ns,1) = min(1.0_pReal/sqrt(sum(abs(rho(:,edg)),2)),dUpper(1:ns,1))
 
  where(dNeq0(sqrt(sum(abs(rho(:,scr)),2)))) &
    dUpper(1:ns,2) = min(1.0_pReal/sqrt(sum(abs(rho(:,scr)),2)),dUpper(1:ns,2))

  dUpper = max(dUpper,dLower)
  
  !****************************************************************************
  !*** calculate dislocation multiplication
  rhoDotMultiplication = 0.0_pReal
  isBCC: if (lattice_structure(ph) == LATTICE_bcc_ID) then
    forall (s = 1:ns, sum(abs(v(s,1:4))) > 0.0_pReal)
      rhoDotMultiplication(s,1:2) = sum(abs(gdot(s,3:4))) / prm%burgers(s) &                        ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(stt%rho_forest(s,o)) / prm%lambda0(s) ! &                  ! mean free path
                                  ! * 2.0_pReal * sum(abs(v(s,3:4))) / sum(abs(v(s,1:4)))           ! ratio of screw to overall velocity determines edge generation
      rhoDotMultiplication(s,3:4) = sum(abs(gdot(s,3:4))) /prm%burgers(s) &                         ! assuming double-cross-slip of screws to be decisive for multiplication
                                  * sqrt(stt%rho_forest(s,o)) / prm%lambda0(s) ! &                  ! mean free path
                                  ! * 2.0_pReal * sum(abs(v(s,1:2))) / sum(abs(v(s,1:4)))           ! ratio of edge to overall velocity determines screw generation
    endforall
  
  else isBCC
    rhoDotMultiplication(1:ns,1:4) = spread( &
          (sum(abs(gdot(1:ns,1:2)),2) * prm%fEdgeMultiplication + sum(abs(gdot(1:ns,3:4)),2)) &
        * sqrt(stt%rho_forest(:,o)) / prm%lambda0 / prm%burgers(1:ns), 2, 4)
  endif isBCC
  
  
  !****************************************************************************
  !*** calculate dislocation fluxes (only for nonlocal plasticity)
  rhoDotFlux = 0.0_pReal
  if (.not. phase_localPlasticity(material_phaseAt(1,el))) then
  
    !*** check CFL (Courant-Friedrichs-Lewy) condition for flux
    if (any( abs(gdot) > 0.0_pReal &                                                                ! any active slip system ...
            .and. prm%CFLfactor * abs(v) * timestep &
                > IPvolume(ip,el) / maxval(IParea(:,ip,el)))) then                                  ! ...with velocity above critical value (we use the reference volume and area for simplicity here)
#ifdef DEBUG
    if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0) then 
      write(6,'(a,i5,a,i2)') '<< CONST >> CFL condition not fullfilled at el ',el,' ip ',ip
      write(6,'(a,e10.3,a,e10.3)') '<< CONST >> velocity is at  ', &
        maxval(abs(v), abs(gdot) > 0.0_pReal &
                       .and.  prm%CFLfactor * abs(v) * timestep &
                             > IPvolume(ip,el) / maxval(IParea(:,ip,el))), &
        ' at a timestep of ',timestep
      write(6,'(a)') '<< CONST >> enforcing cutback !!!'
    endif
#endif
      plasticState(p)%dotState = IEEE_value(1.0_pReal,IEEE_quiet_NaN)                               ! -> return NaN and, hence, enforce cutback
      return
    endif
  
  
    !*** be aware of the definition of slip_transverse = slip_direction x slip_normal !!!
    !*** opposite sign to our p vector in the (s,p,n) triplet !!!
    
    m(1:3,1:ns,1) =  prm%slip_direction
    m(1:3,1:ns,2) = -prm%slip_direction
    m(1:3,1:ns,3) = -prm%slip_transverse
    m(1:3,1:ns,4) =  prm%slip_transverse
    
    my_Fe = Fe(1:3,1:3,1,ip,el)
    my_F = matmul(my_Fe, Fp(1:3,1:3,1,ip,el))
    
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
        neighbor_instance = phase_plasticityInstance(material_phaseAt(1,neighbor_el))
        neighbor_Fe = Fe(1:3,1:3,1,neighbor_ip,neighbor_el)
        neighbor_F = matmul(neighbor_Fe, Fp(1:3,1:3,1,neighbor_ip,neighbor_el))
        Favg = 0.5_pReal * (my_F + neighbor_F)
      else                                                                                          ! if no neighbor, take my value as average
        Favg = my_F
      endif
      
  
      !* FLUX FROM MY NEIGHBOR TO ME
      !* This is only considered, if I have a neighbor of nonlocal plasticity 
      !* (also nonlocal constitutive law with local properties) that is at least a little bit 
      !* compatible.
      !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of
      !* my neighbor's interface.
      !* The entering flux from my neighbor will be distributed on my slip systems according to the 
      !* compatibility
      
      considerEnteringFlux = .false.
      neighbor_v = 0.0_pReal        ! needed for check of sign change in flux density below 
      neighbor_rhoSgl = 0.0_pReal
      if (neighbor_n > 0) then
        if (phase_plasticity(material_phaseAt(1,neighbor_el)) == PLASTICITY_NONLOCAL_ID &
            .and. any(compatibility(:,:,:,n,ip,el) > 0.0_pReal)) &
          considerEnteringFlux = .true.
      endif
      
      enteringFlux: if (considerEnteringFlux) then
        forall (s = 1:ns, t = 1:4)
          neighbor_v(s,t) =          plasticState(np)%state(iV   (s,t,neighbor_instance),no)
          neighbor_rhoSgl(s,t) = max(plasticState(np)%state(iRhoU(s,t,neighbor_instance),no), &
                                                                                            0.0_pReal)
        endforall
  
        where (neighbor_rhoSgl * IPvolume(neighbor_ip,neighbor_el) ** 0.667_pReal < prm%significantN &
          .or. neighbor_rhoSgl < prm%significantRho) &
          neighbor_rhoSgl = 0.0_pReal
        normal_neighbor2me_defConf = math_det33(Favg) * matmul(math_inv33(transpose(Favg)), &
                                     IPareaNormal(1:3,neighbor_n,neighbor_ip,neighbor_el))          ! calculate the normal of the interface in (average) deformed configuration (now pointing from my neighbor to me!!!)
        normal_neighbor2me = matmul(transpose(neighbor_Fe), normal_neighbor2me_defConf) &
                           / math_det33(neighbor_Fe)                                                ! interface normal in the lattice configuration of my neighbor
        area = IParea(neighbor_n,neighbor_ip,neighbor_el) * norm2(normal_neighbor2me)
        normal_neighbor2me = normal_neighbor2me / norm2(normal_neighbor2me)                         ! normalize the surface normal to unit length
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            topp = t + mod(t,2) - mod(t+1,2)
            if (neighbor_v(s,t) * math_inner(m(1:3,s,t), normal_neighbor2me) > 0.0_pReal &          ! flux from my neighbor to me == entering flux for me
                .and. v(s,t) * neighbor_v(s,t) >= 0.0_pReal ) then                                  ! ... only if no sign change in flux density  
              lineLength = neighbor_rhoSgl(s,t) * neighbor_v(s,t) &
                         * math_inner(m(1:3,s,t), normal_neighbor2me) * area                        ! positive line length that wants to enter through this interface
              where (compatibility(c,1:ns,s,n,ip,el) > 0.0_pReal) &                                 ! positive compatibility...
                rhoDotFlux(1:ns,t) = rhoDotFlux(1:ns,t) &
                                   + lineLength / IPvolume(ip,el) &                                 ! ... transferring to equally signed mobile dislocation type
                                   * compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
              where (compatibility(c,1:ns,s,n,ip,el) < 0.0_pReal) &                                 ! ..negative compatibility...
                rhoDotFlux(1:ns,topp) = rhoDotFlux(1:ns,topp) &
                                      + lineLength / IPvolume(ip,el) &                              ! ... transferring to opposite signed mobile dislocation type
                                      * compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
            endif
          enddo
        enddo
      endif enteringFlux
     
   
      !* FLUX FROM ME TO MY NEIGHBOR
      !* This is not considered, if my opposite neighbor has a different constitutive law than nonlocal (still considered for nonlocal law with local properties). 
      !* Then, we assume, that the opposite(!) neighbor sends an equal amount of dislocations to me.
      !* So the net flux in the direction of my neighbor is equal to zero:
      !*    leaving flux to neighbor == entering flux from opposite neighbor
      !* In case of reduced transmissivity, part of the leaving flux is stored as dead dislocation density.
      !* That means for an interface of zero transmissivity the leaving flux is fully converted to dead dislocations.
      
      considerLeavingFlux = .true.
      if (opposite_n > 0) then
        if (phase_plasticity(material_phaseAt(1,opposite_el)) /= PLASTICITY_NONLOCAL_ID) &
          considerLeavingFlux = .false.
      endif
  
      leavingFlux: if (considerLeavingFlux) then
        my_rhoSgl = rhoSgl
        my_v = v
        
        normal_me2neighbor_defConf = math_det33(Favg) &
                                   * matmul(math_inv33(transpose(Favg)), & 
                                                             IPareaNormal(1:3,n,ip,el))             ! calculate the normal of the interface in (average) deformed configuration (pointing from me to my neighbor!!!)
        normal_me2neighbor = matmul(transpose(my_Fe), normal_me2neighbor_defConf) &
                           / math_det33(my_Fe)                                                      ! interface normal in my lattice configuration
        area = IParea(n,ip,el) * norm2(normal_me2neighbor)
        normal_me2neighbor = normal_me2neighbor / norm2(normal_me2neighbor)                         ! normalize the surface normal to unit length    
        do s = 1,ns
          do t = 1,4
            c = (t + 1) / 2
            if (my_v(s,t) * math_inner(m(1:3,s,t), normal_me2neighbor) > 0.0_pReal ) then           ! flux from me to my neighbor == leaving flux for me (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
              if (my_v(s,t) * neighbor_v(s,t) >= 0.0_pReal) then                                    ! no sign change in flux density
                transmissivity = sum(compatibility(c,1:ns,s,n,ip,el)**2.0_pReal)                    ! overall transmissivity from this slip system to my neighbor
              else                                                                                  ! sign change in flux density means sign change in stress which does not allow for dislocations to arive at the neighbor
                transmissivity = 0.0_pReal
              endif
              lineLength = my_rhoSgl(s,t) * my_v(s,t) &
                         * math_inner(m(1:3,s,t), normal_me2neighbor) * area                       ! positive line length of mobiles that wants to leave through this interface
              rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / IPvolume(ip,el)                     ! subtract dislocation flux from current type
              rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) &
                                     + lineLength / IPvolume(ip,el) * (1.0_pReal - transmissivity) &
                                     * sign(1.0_pReal, my_v(s,t))                                   ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
            endif
          enddo
        enddo
      endif leavingFlux
      
    enddo neighbors
  endif
  
  
  
  !****************************************************************************
  !*** calculate dipole formation and annihilation
  
  !*** formation by glide
  
  do c = 1,2
    rhoDotSingle2DipoleGlide(1:ns,2*c-1) = -2.0_pReal * dUpper(1:ns,c) / prm%burgers(1:ns) &
                                                      * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) & ! negative mobile --> positive mobile
                                                         + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &   ! positive mobile --> negative mobile
                                                         + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) ! positive mobile --> negative immobile
  
    rhoDotSingle2DipoleGlide(1:ns,2*c) = -2.0_pReal * dUpper(1:ns,c) / prm%burgers(1:ns) &
                                                    * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &   ! negative mobile --> positive mobile
                                                       + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) & ! positive mobile --> negative mobile
                                                       + abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c))) ! negative mobile --> positive immobile
  
    rhoDotSingle2DipoleGlide(1:ns,2*c+3) = -2.0_pReal * dUpper(1:ns,c) / prm%burgers(1:ns) &
                                                      * rhoSgl(1:ns,2*c+3) * abs(gdot(1:ns,2*c))    ! negative mobile --> positive immobile
  
    rhoDotSingle2DipoleGlide(1:ns,2*c+4) = -2.0_pReal * dUpper(1:ns,c) / prm%burgers(1:ns)&
                                                      * rhoSgl(1:ns,2*c+4) * abs(gdot(1:ns,2*c-1))  ! positive mobile --> negative immobile
  
    rhoDotSingle2DipoleGlide(1:ns,c+8) = - rhoDotSingle2DipoleGlide(1:ns,2*c-1) &
                                         - rhoDotSingle2DipoleGlide(1:ns,2*c) &
                                         + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+3)) &
                                         + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+4))
  enddo
  
  
  !*** athermal annihilation
  
  rhoDotAthermalAnnihilation = 0.0_pReal
  
  forall (c=1:2) &  
    rhoDotAthermalAnnihilation(1:ns,c+8) = -2.0_pReal * dLower(1:ns,c) / prm%burgers(1:ns) &
       * (  2.0_pReal * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &        ! was single hitting single
     + 2.0_pReal * (abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
     + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                                 ! single knocks dipole constituent
  ! annihilated screw dipoles leave edge jogs behind on the colinear system
  
  if (lattice_structure(ph) == LATTICE_fcc_ID) &
    forall (s = 1:ns, prm%colinearSystem(s) > 0) &
      rhoDotAthermalAnnihilation(prm%colinearSystem(s),1:2) = - rhoDotAthermalAnnihilation(s,10) &
        * 0.25_pReal * sqrt(stt%rho_forest(s,o)) * (dUpper(s,2) + dLower(s,2)) * prm%edgeJogFactor
  
    
    
  !*** thermally activated annihilation of edge dipoles by climb
  
  rhoDotThermalAnnihilation = 0.0_pReal
  selfDiffusion = prm%Dsd0 * exp(-prm%selfDiffusionEnergy / (KB * Temperature))
  vClimb =  prm%atomicVolume * selfDiffusion / ( KB * Temperature ) &
            * prm%mu / ( 2.0_pReal * PI * (1.0_pReal-prm%nu) ) &
            * 2.0_pReal / ( dUpper(1:ns,1) + dLower(1:ns,1) )
  forall (s = 1:ns, dUpper(s,1) > dLower(s,1)) &
    rhoDotThermalAnnihilation(s,9) = max(- 4.0_pReal * rhoDip(s,1) * vClimb(s) / (dUpper(s,1) - dLower(s,1)), &
                                         - rhoDip(s,1) / timestep - rhoDotAthermalAnnihilation(s,9) &
                                                                  - rhoDotSingle2DipoleGlide(s,9))    ! make sure that we do not annihilate more dipoles than we have

  rhoDot = 0.0_pReal
  rhoDot = rhoDotFlux &
         + rhoDotMultiplication &
         + rhoDotSingle2DipoleGlide &
         + rhoDotAthermalAnnihilation &
         + rhoDotThermalAnnihilation 

#ifdef DEBUG
  if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0 &
      .and. ((debug_e == el .and. debug_i == ip)&
             .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0 )) then
    write(6,'(a,/,4(12x,12(e12.5,1x),/))')  '<< CONST >> dislocation multiplication', &
                                            rhoDotMultiplication(1:ns,1:4) * timestep
    write(6,'(a,/,8(12x,12(e12.5,1x),/))')  '<< CONST >> dislocation flux', &
                                            rhoDotFlux(1:ns,1:8) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> dipole formation by glide', &
                                            rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> athermal dipole annihilation', &
                                            rhoDotAthermalAnnihilation * timestep
    write(6,'(a,/,2(12x,12(e12.5,1x),/))')  '<< CONST >> thermally activated dipole annihilation', &
                                            rhoDotThermalAnnihilation(1:ns,9:10) * timestep
    write(6,'(a,/,10(12x,12(e12.5,1x),/))') '<< CONST >> total density change', &
                                            rhoDot * timestep
    write(6,'(a,/,10(12x,12(f12.5,1x),/))') '<< CONST >> relative density change', &
                                            rhoDot(1:ns,1:8)  * timestep / (abs(stt%rho(:,sgl))+1.0e-10), &
                                            rhoDot(1:ns,9:10) * timestep / (stt%rho(:,dip)+1.0e-10)
    write(6,*)
  endif
#endif


  if (    any(rho(:,mob) + rhoDot(1:ns,1:4)  * timestep < -prm%aTolRho) &
     .or. any(rho(:,dip) + rhoDot(1:ns,9:10) * timestep < -prm%aTolRho)) then
#ifdef DEBUG
    if (iand(debug_level(debug_constitutive),debug_levelExtensive) /= 0) then 
      write(6,'(a,i5,a,i2)') '<< CONST >> evolution rate leads to negative density at el ',el,' ip ',ip
      write(6,'(a)') '<< CONST >> enforcing cutback !!!'
    endif
#endif
    plasticState(p)%dotState = IEEE_value(1.0_pReal,IEEE_quiet_NaN)
  else
    dot%rho(:,o) = pack(rhoDot,.true.)
    forall (s = 1:ns) &
      dot%accumulatedshear(s,o) = sum(gdot(s,1:4))
  endif

  end associate

end subroutine plastic_nonlocal_dotState

!--------------------------------------------------------------------------------------------------
!> @brief Compatibility update
!> @detail Compatibility is defined as normalized product of signed cosine of the angle between the slip 
! plane normals and signed cosine of the angle between the slip directions. Only the largest values
! that sum up to a total of 1 are considered, all others are set to zero.
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_updateCompatibility(orientation,i,e) 
  
  integer, intent(in) :: &
    i, &
    e
  type(rotation), dimension(1,discretization_nIP,discretization_nElem), intent(in) :: &
    orientation                                                                                     ! crystal orientation in quaternions
                                              
  integer :: &
    Nneighbors, &                                                                                   ! number of neighbors
    n, &                                                                                            ! neighbor index 
    neighbor_e, &                                                                                   ! element index of my neighbor
    neighbor_i, &                                                                                   ! integration point index of my neighbor
    ph, &
    neighbor_phase, &
    textureID, &
    neighbor_textureID, &
    instance, &                                                                                     ! instance of plasticity
    ns, &                                                                                           ! number of active slip systems
    s1, &                                                                                           ! slip system index (me)
    s2                                                                                              ! slip system index (my neighbor)
  real(pReal), dimension(4) :: &
    absoluteMisorientation                                                                          ! absolute misorientation (without symmetry) between me and my neighbor
  real(pReal), dimension(2,totalNslip(phase_plasticityInstance(material_phaseAt(1,e))),&
                           totalNslip(phase_plasticityInstance(material_phaseAt(1,e))),&
                           nIPneighbors) :: &  
    my_compatibility                                                                                ! my_compatibility for current element and ip  
  real(pReal) :: &
    my_compatibilitySum, &
    thresholdValue, &
    nThresholdValues
  logical, dimension(totalNslip(phase_plasticityInstance(material_phaseAt(1,e)))) :: & 
    belowThreshold
  type(rotation) :: rot

  Nneighbors = nIPneighbors
  ph = material_phaseAt(1,e)
  textureID = material_texture(1,i,e)
  instance = phase_plasticityInstance(ph)
  ns = totalNslip(instance)
  associate(prm => param(instance))
  
  !*** start out fully compatible
  my_compatibility = 0.0_pReal
  
  forall(s1 = 1:ns) my_compatibility(1:2,s1,s1,1:Nneighbors) = 1.0_pReal 
  
  !*** Loop thrugh neighbors and check whether there is any compatibility.
  
  neighbors: do n = 1,Nneighbors
    neighbor_e = IPneighborhood(1,n,i,e)
    neighbor_i = IPneighborhood(2,n,i,e)
    
    
    !* FREE SURFACE
    !* Set surface transmissivity to the value specified in the material.config
    
    if (neighbor_e <= 0 .or. neighbor_i <= 0) then
      forall(s1 = 1:ns) my_compatibility(1:2,s1,s1,n) = sqrt(prm%surfaceTransmissivity)
      cycle
    endif
    
    
    !* PHASE BOUNDARY
    !* If we encounter a different nonlocal phase at the neighbor, 
    !* we consider this to be a real "physical" phase boundary, so completely incompatible.
    !* If one of the two phases has a local plasticity law, 
    !* we do not consider this to be a phase boundary, so completely compatible.
    neighbor_phase = material_phaseAt(1,neighbor_e)
    if (neighbor_phase /= ph) then
      if (.not. phase_localPlasticity(neighbor_phase) .and. .not. phase_localPlasticity(ph))&
        forall(s1 = 1:ns) my_compatibility(1:2,s1,s1,n) = 0.0_pReal
      cycle
    endif
  
      
    !* GRAIN BOUNDARY !
    !* fixed transmissivity for adjacent ips with different texture (only if explicitly given in material.config)
    if (prm%grainboundaryTransmissivity >= 0.0_pReal) then
      neighbor_textureID = material_texture(1,neighbor_i,neighbor_e)
      if (neighbor_textureID /= textureID) then
        if (.not. phase_localPlasticity(neighbor_phase)) then
          forall(s1 = 1:ns) &
            my_compatibility(1:2,s1,s1,n) = sqrt(prm%grainboundaryTransmissivity)
        endif
        cycle
      endif
    
  
    !* GRAIN BOUNDARY ?
    !* Compatibility defined by relative orientation of slip systems:
    !* The my_compatibility value is defined as the product of the slip normal projection and the slip direction projection.
    !* Its sign is always positive for screws, for edges it has the same sign as the slip normal projection. 
    !* Since the sum for each slip system can easily exceed one (which would result in a transmissivity larger than one), 
    !* only values above or equal to a certain threshold value are considered. This threshold value is chosen, such that
    !* the number of compatible slip systems is minimized with the sum of the original compatibility values exceeding one. 
    !* Finally the smallest compatibility value is decreased until the sum is exactly equal to one. 
    !* All values below the threshold are set to zero.
    else
      rot = orientation(1,i,e)%misorientation(orientation(1,neighbor_i,neighbor_e))
      absoluteMisorientation = rot%asQuaternion()
      mySlipSystems: do s1 = 1,ns
        neighborSlipSystems: do s2 = 1,ns
          my_compatibility(1,s2,s1,n) =  math_inner(prm%slip_normal(1:3,s1), &
          math_qRot(absoluteMisorientation, prm%slip_normal(1:3,s2))) &
                                  * abs(math_inner(prm%slip_direction(1:3,s1), &
                                  math_qRot(absoluteMisorientation, prm%slip_direction(1:3,s2))))
          my_compatibility(2,s2,s1,n) = abs(math_inner(prm%slip_normal(1:3,s1), &
          math_qRot(absoluteMisorientation, prm%slip_normal(1:3,s2)))) &
                                  * abs(math_inner(prm%slip_direction(1:3,s1), &
                                  math_qRot(absoluteMisorientation, prm%slip_direction(1:3,s2))))
        enddo neighborSlipSystems
        
        my_compatibilitySum = 0.0_pReal
        belowThreshold = .true.
        do while (my_compatibilitySum < 1.0_pReal .and. any(belowThreshold(1:ns)))
          thresholdValue = maxval(my_compatibility(2,1:ns,s1,n), belowThreshold(1:ns))              ! screws always positive
          nThresholdValues = real(count(my_compatibility(2,1:ns,s1,n) >= thresholdValue),pReal)
          where (my_compatibility(2,1:ns,s1,n) >= thresholdValue) &
            belowThreshold(1:ns) = .false.
          if (my_compatibilitySum + thresholdValue * nThresholdValues > 1.0_pReal) &
            where (abs(my_compatibility(1:2,1:ns,s1,n)) >= thresholdValue) &                          ! MD: rather check below threshold?
              my_compatibility(1:2,1:ns,s1,n) = sign((1.0_pReal - my_compatibilitySum) &
                                                   / nThresholdValues, my_compatibility(1:2,1:ns,s1,n))
          my_compatibilitySum = my_compatibilitySum + nThresholdValues * thresholdValue
        enddo
        where (belowThreshold(1:ns)) my_compatibility(1,1:ns,s1,n) = 0.0_pReal
        where (belowThreshold(1:ns)) my_compatibility(2,1:ns,s1,n) = 0.0_pReal
      enddo mySlipSystems
    endif
  
  enddo neighbors
  
  compatibility(1:2,1:ns,1:ns,1:Nneighbors,i,e) = my_compatibility
  
  end associate
  
end subroutine plastic_nonlocal_updateCompatibility

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_nonlocal_postResults(ph,instance,of) result(postResults)

 integer, intent(in) :: &
   ph, &
   instance, &
   of
 
 real(pReal),   dimension(sum(plastic_nonlocal_sizePostResult(:,instance))) :: &
   postResults
 
 integer :: &
   ns, &                                                                                            !< short notation for the total number of active slip systems
   c, &                                                                                             !< character of dislocation
   cs, &                                                                                            !< constitutive result index
   o, &                                                                                             !< index of current output
   t, &                                                                                             !< type of dislocation
   s                                                                                                !< index of my current slip system

 real(pReal), dimension(param(instance)%totalNslip,8) :: &
   rhoSgl, &
   rhoDotSgl                                                                                        !< evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
 real(pReal), dimension(param(instance)%totalNslip,4) :: &
   v                                                                                                !< velocities
 real(pReal), dimension(param(instance)%totalNslip,2) :: &
   rhoDotDip                                                                                        !< evolution rate of dipole dislocation densities (screw and edge dipoles)

   
ns = param(instance)%totalNslip

cs = 0

associate(prm => param(instance),dst => microstructure(instance),stt=>state(instance),dot => dotState(instance))

forall (s = 1:ns, t = 1:4)
  rhoSgl(s,t+4)    = plasticState(ph)%State(iRhoB(s,t,instance),of)
  v(s,t)           = plasticState(ph)%State(iV(s,t,instance),of)
  rhoDotSgl(s,t+4) = plasticState(ph)%dotState(iRhoB(s,t,instance),of)
endforall
forall (s = 1:ns, c = 1:2)
  rhoDotDip(s,c) = plasticState(ph)%dotState(iRhoD(s,c,instance),of)
endforall

outputsLoop: do o = 1,size(param(instance)%outputID)
  select case(param(instance)%outputID(o))
      
    case (rho_sgl_mob_edg_pos_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_mob_edg_pos(:,of)
      cs = cs + ns
      
    case (rho_sgl_imm_edg_pos_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_imm_edg_pos(:,of)
      cs = cs + ns
      
    case (rho_sgl_mob_edg_neg_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_mob_edg_neg(:,of)
      cs = cs + ns
      
    case (rho_sgl_imm_edg_neg_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_imm_edg_neg(:,of)
      cs = cs + ns
      
    case (rho_dip_edg_ID)
      postResults(cs+1:cs+ns) = stt%rho_dip_edg(:,of)
      cs = cs + ns
      
    case (rho_sgl_mob_scr_pos_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_mob_scr_pos(:,of)
      cs = cs + ns
      
    case (rho_sgl_imm_scr_pos_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_imm_scr_pos(:,of)
      cs = cs + ns

    case (rho_sgl_mob_scr_neg_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_mob_scr_neg(:,of)
      cs = cs + ns

    case (rho_sgl_imm_scr_neg_ID)
      postResults(cs+1:cs+ns) = stt%rho_sgl_imm_scr_neg(:,of)
      cs = cs + ns

    case (rho_dip_scr_ID)
      postResults(cs+1:cs+ns) = stt%rho_dip_scr(:,of)
      cs = cs + ns
      
    case (rho_forest_ID)
      postResults(cs+1:cs+ns) = stt%rho_forest(:,of)
      cs = cs + ns
      
    case (resolvedstress_back_ID)
      postResults(cs+1:cs+ns) = dst%tau_back(:,of)
      cs = cs + ns
      
    case (resistance_ID)
      postResults(cs+1:cs+ns) = dst%tau_Threshold(:,of)
      cs = cs + ns
      
    case (rho_dot_sgl_ID)
      postResults(cs+1:cs+ns) = sum(rhoDotSgl(1:ns,1:4),2) &
                              + sum(rhoDotSgl(1:ns,5:8)*sign(1.0_pReal,rhoSgl(1:ns,5:8)),2)
      cs = cs + ns
      
    case (rho_dot_sgl_mobile_ID)
      postResults(cs+1:cs+ns) = sum(rhoDotSgl(1:ns,1:4),2)
      cs = cs + ns
      
    case (rho_dot_dip_ID)
      postResults(cs+1:cs+ns) = sum(rhoDotDip,2)
      cs = cs + ns

    case (velocity_edge_pos_ID)
      postResults(cs+1:cs+ns) = v(1:ns,1)
      cs = cs + ns
    
    case (velocity_edge_neg_ID)
      postResults(cs+1:cs+ns) = v(1:ns,2)
      cs = cs + ns
    
    case (velocity_screw_pos_ID)
      postResults(cs+1:cs+ns) = v(1:ns,3)
      cs = cs + ns
    
    case (velocity_screw_neg_ID)
      postResults(cs+1:cs+ns) = v(1:ns,4)
      cs = cs + ns

    case(accumulatedshear_ID)
      postResults(cs+1:cs+ns) = stt%accumulatedshear(:,of)
      cs = cs + ns
    
  end select
enddo outputsLoop
end associate
end function plastic_nonlocal_postResults


!--------------------------------------------------------------------------------------------------
!> @brief returns copy of current dislocation densities from state
!> @details raw values is rectified
!--------------------------------------------------------------------------------------------------
function getRho(instance,of,ip,el)
  
  integer, intent(in) :: instance, of,ip,el
  real(pReal), dimension(param(instance)%totalNslip,10) :: getRho
  
  associate(prm => param(instance))

  getRho = reshape(state(instance)%rho(:,of),[prm%totalNslip,10])
  
  ! ensure positive densities (not for imm, they have a sign)
  getRho(:,mob) = max(getRho(:,mob),0.0_pReal)
  getRho(:,dip) = max(getRho(:,dip),0.0_pReal)
  
  where(abs(getRho) < max(prm%significantN/IPvolume(ip,el)**(2.0_pReal/3.0_pReal),prm%significantRho)) &
    getRho = 0.0_pReal

  end associate

end function getRho


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_nonlocal_results(instance,group)
#if defined(PETSc) || defined(DAMASK_HDF5)
  use results, only: &
    results_writeDataset

  integer, intent(in) :: instance
  character(len=*) :: group
  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
      case (rho_sgl_mob_edg_pos_ID)
        call results_writeDataset(group,stt%rho_sgl_mob_edg_pos, 'rho_sgl_mob_edg_pos', &
                                  'positive mobile edge density','1/m²')
      case (rho_sgl_imm_edg_pos_ID)
        call results_writeDataset(group,stt%rho_sgl_imm_edg_pos, 'rho_sgl_imm_edg_pos',&
                                  'positive immobile edge density','1/m²')
      case (rho_sgl_mob_edg_neg_ID)
        call results_writeDataset(group,stt%rho_sgl_mob_edg_neg, 'rho_sgl_mob_edg_neg',&
                                  'negative mobile edge density','1/m²')
      case (rho_sgl_imm_edg_neg_ID)
        call results_writeDataset(group,stt%rho_sgl_imm_edg_neg, 'rho_sgl_imm_edg_neg',&
                                  'negative immobile edge density','1/m²')
      case (rho_dip_edg_ID)
        call results_writeDataset(group,stt%rho_dip_edg, 'rho_dip_edg',&
                                  'edge dipole density','1/m²')
      case (rho_sgl_mob_scr_pos_ID)
        call results_writeDataset(group,stt%rho_sgl_mob_scr_pos, 'rho_sgl_mob_scr_pos',&
                                  'positive mobile screw density','1/m²')
      case (rho_sgl_imm_scr_pos_ID)
        call results_writeDataset(group,stt%rho_sgl_imm_scr_pos, 'rho_sgl_imm_scr_pos',&
                                  'positive immobile screw density','1/m²')
      case (rho_sgl_mob_scr_neg_ID)
        call results_writeDataset(group,stt%rho_sgl_mob_scr_neg, 'rho_sgl_mob_scr_neg',&
                                  'negative mobile screw density','1/m²')
      case (rho_sgl_imm_scr_neg_ID)
        call results_writeDataset(group,stt%rho_sgl_imm_scr_neg, 'rho_sgl_imm_scr_neg',&
                                  'negative immobile screw density','1/m²')
      case (rho_dip_scr_ID)
        call results_writeDataset(group,stt%rho_dip_scr, 'rho_dip_scr',&
                                  'screw dipole density','1/m²')
      case (rho_forest_ID)
        call results_writeDataset(group,stt%rho_forest, 'rho_forest',&
                                  'forest density','1/m²')
    end select
  enddo outputsLoop
  end associate
#else
  integer, intent(in) :: instance
  character(len=*) :: group
#endif

end subroutine plastic_nonlocal_results

end module plastic_nonlocal
