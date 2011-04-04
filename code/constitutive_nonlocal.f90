!* $Id$
!************************************ 
!*  Module: CONSTITUTIVE_NONLOCAL   *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!************************************


MODULE constitutive_nonlocal

!* Include other modules
use prec, only: pReal,pInt
implicit none


!* Definition of parameters
character (len=*), parameter :: constitutive_nonlocal_label = 'nonlocal'
character(len=22), dimension(10), parameter ::    constitutive_nonlocal_listBasicStates = (/'rhoSglEdgePosMobile   ', &
                                                                                            'rhoSglEdgeNegMobile   ', &
                                                                                            'rhoSglScrewPosMobile  ', &
                                                                                            'rhoSglScrewNegMobile  ', &
                                                                                            'rhoSglEdgePosImmobile ', &
                                                                                            'rhoSglEdgeNegImmobile ', &
                                                                                            'rhoSglScrewPosImmobile', &
                                                                                            'rhoSglScrewNegImmobile', &
                                                                                            'rhoDipEdge            ', &
                                                                                            'rhoDipScrew           ' /) ! list of "basic" microstructural state variables that are independent from other state variables
character(len=15), dimension(3), parameter :: constitutive_nonlocal_listDependentStates = (/'rhoForest      ', &
                                                                                            'tauThreshold   ', &
                                                                                            'Tdislocation_v ' /) ! list of microstructural state variables that depend on other state variables
real(pReal), parameter :: kB = 1.38e-23_pReal                                                                   ! Physical parameter, Boltzmann constant in J/Kelvin

!* Definition of global variables
integer(pInt), dimension(:), allocatable ::               constitutive_nonlocal_sizeDotState, &                 ! number of dotStates
                                                          constitutive_nonlocal_sizeState, &                    ! total number of microstructural state variables
                                                          constitutive_nonlocal_sizePostResults                 ! cumulative size of post results
integer(pInt), dimension(:,:), allocatable, target ::     constitutive_nonlocal_sizePostResult                  ! size of each post result output
character(len=64), dimension(:,:), allocatable, target :: constitutive_nonlocal_output                          ! name of each post result output

character(len=32), dimension(:), allocatable ::           constitutive_nonlocal_structureName                   ! name of the lattice structure
integer(pInt), dimension(:), allocatable ::               constitutive_nonlocal_structure, &                    ! number representing the kind of lattice structure
                                                          constitutive_nonlocal_totalNslip                      ! total number of active slip systems for each instance
integer(pInt), dimension(:,:), allocatable ::             constitutive_nonlocal_Nslip, &                        ! number of active slip systems for each family and instance
                                                          constitutive_nonlocal_slipFamily, &                   ! lookup table relating active slip system to slip family for each instance
                                                          constitutive_nonlocal_slipSystemLattice               ! lookup table relating active slip system index to lattice slip system index for each instance

real(pReal), dimension(:), allocatable ::                 constitutive_nonlocal_CoverA, &                       ! c/a ratio for hex type lattice
                                                          constitutive_nonlocal_C11, &                          ! C11 element in elasticity matrix
                                                          constitutive_nonlocal_C12, &                          ! C12 element in elasticity matrix
                                                          constitutive_nonlocal_C13, &                          ! C13 element in elasticity matrix
                                                          constitutive_nonlocal_C33, &                          ! C33 element in elasticity matrix
                                                          constitutive_nonlocal_C44, &                          ! C44 element in elasticity matrix
                                                          constitutive_nonlocal_Gmod, &                         ! shear modulus
                                                          constitutive_nonlocal_nu, &                           ! poisson's ratio
                                                          constitutive_nonlocal_atomicVolume, &                 ! atomic volume
                                                          constitutive_nonlocal_Dsd0, &                         ! prefactor for self-diffusion coefficient
                                                          constitutive_nonlocal_Qsd, &                          ! activation enthalpy for diffusion
                                                          constitutive_nonlocal_aTolRho, &                      ! absolute tolerance for dislocation density in state integration
                                                          constitutive_nonlocal_R, &                            ! cutoff radius for dislocation stress
                                                          constitutive_nonlocal_d0, &                           ! wall depth as multiple of b
                                                          constitutive_nonlocal_tauObs, &                       ! obstacle strength in Pa
                                                          constitutive_nonlocal_fattack, &                      ! attack frequency in Hz
                                                          constitutive_nonlocal_vs, &                           ! maximum dislocation velocity = velocity of sound
                                                          constitutive_nonlocal_rhoSglScatter, &                ! standard deviation of scatter in initial dislocation density
                                                          constitutive_nonlocal_surfaceTransmissivity           ! transmissivity at free surface
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_Cslip_66                        ! elasticity matrix in Mandel notation for each instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_Cslip_3333                      ! elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_nonlocal_rhoSglEdgePos0, &               ! initial edge_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglEdgeNeg0, &               ! initial edge_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewPos0, &              ! initial screw_pos dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoSglScrewNeg0, &              ! initial screw_neg dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipEdge0, &                  ! initial edge dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_rhoDipScrew0, &                 ! initial screw dipole dislocation density per slip system for each family and instance
                                                          constitutive_nonlocal_lambda0PerSlipFamily, &         ! mean free path prefactor for each family and instance
                                                          constitutive_nonlocal_lambda0PerSlipSystem, &         ! mean free path prefactor for each slip system and instance
                                                          constitutive_nonlocal_burgersPerSlipFamily, &         ! absolute length of burgers vector [m] for each family and instance
                                                          constitutive_nonlocal_burgersPerSlipSystem, &         ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_nonlocal_dLowerEdgePerSlipFamily, &      ! minimum stable edge dipole height for each family and instance
                                                          constitutive_nonlocal_dLowerEdgePerSlipSystem, &      ! minimum stable edge dipole height for each slip system and instance
                                                          constitutive_nonlocal_dLowerScrewPerSlipFamily, &     ! minimum stable screw dipole height for each family and instance
                                                          constitutive_nonlocal_dLowerScrewPerSlipSystem, &     ! minimum stable screw dipole height for each slip system and instance
                                                          constitutive_nonlocal_Qeff0, &                        ! prefactor for activation enthalpy for dislocation glide in J
                                                          constitutive_nonlocal_interactionSlipSlip             ! coefficients for slip-slip interaction for each interaction type and instance
real(pReal), dimension(:,:,:,:,:), allocatable ::         constitutive_nonlocal_v, &                            ! dislocation velocity
                                                          constitutive_nonlocal_rhoDotFlux                      ! dislocation convection term
real(pReal), dimension(:,:,:,:,:,:), allocatable ::       constitutive_nonlocal_compatibility                   ! slip system compatibility between me and my neighbors
real(pReal), dimension(:,:,:), allocatable ::             constitutive_nonlocal_forestProjectionEdge, &         ! matrix of forest projections of edge dislocations for each instance
                                                          constitutive_nonlocal_forestProjectionScrew, &        ! matrix of forest projections of screw dislocations for each instance
                                                          constitutive_nonlocal_interactionMatrixSlipSlip       ! interaction matrix of the different slip systems for each instance

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_stateInit
!* - constitutive_homogenizedC
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - constitutive_dotState
!* - constitutive_dotTemperature
!* - constitutive_postResults
!****************************************


!**************************************
!*      Module initialization         *
!**************************************
subroutine constitutive_nonlocal_init(file)

use prec,     only: pInt, pReal
use math,     only: math_Mandel3333to66, & 
                    math_Voigt66to3333, & 
                    math_mul3x3
use IO,       only: IO_lc, &
                    IO_getTag, &
                    IO_isBlank, &
                    IO_stringPos, &
                    IO_stringValue, &
                    IO_floatValue, &
                    IO_intValue, &
                    IO_error
use debug,    only: debug_verbosity
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    FE_maxNipNeighbors
use material, only: homogenization_maxNgrains, &
                    phase_constitution, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_maxNslipFamily, &
                    lattice_maxNtwinFamily, &
                    lattice_maxNslip, &
                    lattice_maxNtwin, &
                    lattice_maxNinteraction, &
                    lattice_NslipSystem, &
                    lattice_NtwinSystem, &
                    lattice_initializeStructure, &
                    lattice_Qtwin, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st, &
                    lattice_interactionSlipSlip

!*** output variables

!*** input variables
integer(pInt), intent(in) ::                file

!*** local variables
integer(pInt), parameter ::                 maxNchunks = 21
integer(pInt), dimension(1+2*maxNchunks) :: positions
integer(pInt)                               section, &
                                            maxNinstance, &
                                            maxTotalNslip, &
                                            myStructure, &
                                            f, &                ! index of my slip family
                                            i, &                ! index of my instance of this constitution
                                            j, &
                                            k, &
                                            l, &
                                            ns, &               ! short notation for total number of active slip systems for the current instance
                                            o, &                ! index of my output
                                            s, &                ! index of my slip system
                                            s1, &               ! index of my slip system
                                            s2, &               ! index of my slip system
                                            it, &               ! index of my interaction type
                                            output, &
                                            mySize
character(len=64)                           tag
character(len=1024)                         line


!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_nonlocal_label,' init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
!$OMP END CRITICAL (write2out)

maxNinstance = count(phase_constitution == constitutive_nonlocal_label)
if (maxNinstance == 0) return                                                                                                       ! we don't have to do anything if there's no instance for this constitutive law

if (debug_verbosity > 0) then
  !$OMP CRITICAL (write2out)
    write(6,'(a16,x,i5)') '# instances:',maxNinstance
  !$OMP END CRITICAL (write2out)
endif


!*** space allocation for global variables

allocate(constitutive_nonlocal_sizeDotState(maxNinstance))
allocate(constitutive_nonlocal_sizeState(maxNinstance))
allocate(constitutive_nonlocal_sizePostResults(maxNinstance))
allocate(constitutive_nonlocal_sizePostResult(maxval(phase_Noutput), maxNinstance))
allocate(constitutive_nonlocal_output(maxval(phase_Noutput), maxNinstance))
constitutive_nonlocal_sizeDotState = 0_pInt
constitutive_nonlocal_sizeState = 0_pInt
constitutive_nonlocal_sizePostResults = 0_pInt
constitutive_nonlocal_sizePostResult = 0_pInt
constitutive_nonlocal_output = ''

allocate(constitutive_nonlocal_structureName(maxNinstance))
allocate(constitutive_nonlocal_structure(maxNinstance))
allocate(constitutive_nonlocal_Nslip(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_slipFamily(lattice_maxNslip, maxNinstance))
allocate(constitutive_nonlocal_slipSystemLattice(lattice_maxNslip, maxNinstance))
allocate(constitutive_nonlocal_totalNslip(maxNinstance))
constitutive_nonlocal_structureName = ''
constitutive_nonlocal_structure = 0_pInt
constitutive_nonlocal_Nslip = 0_pInt
constitutive_nonlocal_slipFamily = 0_pInt
constitutive_nonlocal_slipSystemLattice = 0_pInt
constitutive_nonlocal_totalNslip = 0_pInt

allocate(constitutive_nonlocal_CoverA(maxNinstance))
allocate(constitutive_nonlocal_C11(maxNinstance))
allocate(constitutive_nonlocal_C12(maxNinstance))
allocate(constitutive_nonlocal_C13(maxNinstance))
allocate(constitutive_nonlocal_C33(maxNinstance))
allocate(constitutive_nonlocal_C44(maxNinstance))
allocate(constitutive_nonlocal_Gmod(maxNinstance))
allocate(constitutive_nonlocal_nu(maxNinstance))
allocate(constitutive_nonlocal_atomicVolume(maxNinstance))
allocate(constitutive_nonlocal_Dsd0(maxNinstance))
allocate(constitutive_nonlocal_Qsd(maxNinstance))
allocate(constitutive_nonlocal_aTolRho(maxNinstance))
allocate(constitutive_nonlocal_Cslip_66(6,6,maxNinstance))
allocate(constitutive_nonlocal_Cslip_3333(3,3,3,3,maxNinstance))
allocate(constitutive_nonlocal_R(maxNinstance))
allocate(constitutive_nonlocal_d0(maxNinstance))
allocate(constitutive_nonlocal_tauObs(maxNinstance))
allocate(constitutive_nonlocal_vs(maxNinstance))
allocate(constitutive_nonlocal_fattack(maxNinstance))
allocate(constitutive_nonlocal_rhoSglScatter(maxNinstance))
allocate(constitutive_nonlocal_surfaceTransmissivity(maxNinstance))
constitutive_nonlocal_CoverA = 0.0_pReal 
constitutive_nonlocal_C11 = 0.0_pReal
constitutive_nonlocal_C12 = 0.0_pReal
constitutive_nonlocal_C13 = 0.0_pReal
constitutive_nonlocal_C33 = 0.0_pReal
constitutive_nonlocal_C44 = 0.0_pReal
constitutive_nonlocal_Gmod = 0.0_pReal
constitutive_nonlocal_atomicVolume = 0.0_pReal
constitutive_nonlocal_Dsd0 = 0.0_pReal
constitutive_nonlocal_Qsd = 0.0_pReal
constitutive_nonlocal_aTolRho = 0.0_pReal
constitutive_nonlocal_nu = 0.0_pReal
constitutive_nonlocal_Cslip_66 = 0.0_pReal
constitutive_nonlocal_Cslip_3333 = 0.0_pReal
constitutive_nonlocal_R = -1.0_pReal
constitutive_nonlocal_d0 = 0.0_pReal
constitutive_nonlocal_tauObs = 0.0_pReal
constitutive_nonlocal_vs = 0.0_pReal
constitutive_nonlocal_fattack = 0.0_pReal
constitutive_nonlocal_rhoSglScatter = 0.0_pReal
constitutive_nonlocal_surfaceTransmissivity = 1.0_pReal

allocate(constitutive_nonlocal_rhoSglEdgePos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglEdgeNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewPos0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoSglScrewNeg0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoDipEdge0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_rhoDipScrew0(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_burgersPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_Lambda0PerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_interactionSlipSlip(lattice_maxNinteraction, maxNinstance))
allocate(constitutive_nonlocal_dLowerEdgePerSlipFamily(lattice_maxNslipFamily, maxNinstance))
allocate(constitutive_nonlocal_dLowerScrewPerSlipFamily(lattice_maxNslipFamily, maxNinstance))
constitutive_nonlocal_rhoSglEdgePos0 = -1.0_pReal
constitutive_nonlocal_rhoSglEdgeNeg0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewPos0 = -1.0_pReal
constitutive_nonlocal_rhoSglScrewNeg0 = -1.0_pReal
constitutive_nonlocal_rhoDipEdge0 = -1.0_pReal
constitutive_nonlocal_rhoDipScrew0 = -1.0_pReal
constitutive_nonlocal_burgersPerSlipFamily = 0.0_pReal
constitutive_nonlocal_lambda0PerSlipFamily = 0.0_pReal
constitutive_nonlocal_interactionSlipSlip = 0.0_pReal
constitutive_nonlocal_dLowerEdgePerSlipFamily = 0.0_pReal
constitutive_nonlocal_dLowerScrewPerSlipFamily = 0.0_pReal


!*** readout data from material.config file

rewind(file)
line = ''
section = 0

do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')                                                                                ! wind forward to <phase>
  read(file,'(a1024)',END=100) line
enddo

do                                                                                                                                  ! read thru sections of phase part
  read(file,'(a1024)',END=100) line
  if (IO_isBlank(line)) cycle                                                                                                       ! skip empty lines
  if (IO_getTag(line,'<','>') /= '') exit                                                                                           ! stop at next part
  if (IO_getTag(line,'[',']') /= '') then                                                                                           ! next section
    section = section + 1
    output = 0                                                                                                                      ! reset output counter
    cycle
  endif
  if (section > 0 .and. phase_constitution(section) == constitutive_nonlocal_label) then                                            ! one of my sections
    i = phase_constitutionInstance(section)                                                                                         ! which instance of my constitution is present phase
    positions = IO_stringPos(line,maxNchunks)
    tag = IO_lc(IO_stringValue(line,positions,1))                                                                                   ! extract key
    select case(tag)
      case('constitution','/nonlocal/')
        cycle
      case ('(output)')
        output = output + 1
        constitutive_nonlocal_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
      case ('lattice_structure')
        constitutive_nonlocal_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
      case ('covera_ratio')
        constitutive_nonlocal_CoverA(i) = IO_floatValue(line,positions,2)
      case ('c11')
        constitutive_nonlocal_C11(i) = IO_floatValue(line,positions,2)
      case ('c12')
        constitutive_nonlocal_C12(i) = IO_floatValue(line,positions,2)
      case ('c13')
        constitutive_nonlocal_C13(i) = IO_floatValue(line,positions,2)
      case ('c33')
        constitutive_nonlocal_C33(i) = IO_floatValue(line,positions,2)
      case ('c44')
        constitutive_nonlocal_C44(i) = IO_floatValue(line,positions,2)
      case ('nslip')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_Nslip(f,i) = IO_intValue(line,positions,1+f)
      case ('rhosgledgepos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglEdgePos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosgledgeneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglEdgeNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosglscrewpos0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglScrewPos0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhosglscrewneg0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoSglScrewNeg0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhodipedge0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoDipEdge0(f,i) = IO_floatValue(line,positions,1+f)
      case ('rhodipscrew0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_rhoDipScrew0(f,i) = IO_floatValue(line,positions,1+f)
      case ('lambda0')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_lambda0PerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case ('burgers')
        forall (f = 1:lattice_maxNslipFamily) constitutive_nonlocal_burgersPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('r')
        constitutive_nonlocal_R(i) = IO_floatValue(line,positions,2)
      case('ddipminedge')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('ddipminscrew')
        forall (f = 1:lattice_maxNslipFamily) & 
            constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i) = IO_floatValue(line,positions,1+f)
      case('atomicvolume')
        constitutive_nonlocal_atomicVolume(i) = IO_floatValue(line,positions,2)
      case('dsd0')
        constitutive_nonlocal_Dsd0(i) = IO_floatValue(line,positions,2)
      case('qsd')
        constitutive_nonlocal_Qsd(i) = IO_floatValue(line,positions,2)
      case('atol_rho')
        constitutive_nonlocal_aTolRho(i) = IO_floatValue(line,positions,2)
      case ('interaction_slipslip')
        forall (it = 1:lattice_maxNinteraction) constitutive_nonlocal_interactionSlipSlip(it,i) = IO_floatValue(line,positions,1+it)
      case('d0')
        constitutive_nonlocal_d0(i) = IO_floatValue(line,positions,2)
      case('tauobs')
        constitutive_nonlocal_tauObs(i) = IO_floatValue(line,positions,2)
      case('vs')
        constitutive_nonlocal_vs(i) = IO_floatValue(line,positions,2)
      case('fattack')
        constitutive_nonlocal_fattack(i) = IO_floatValue(line,positions,2)
      case('rhosglscatter')
        constitutive_nonlocal_rhoSglScatter(i) = IO_floatValue(line,positions,2)
      case('surfacetransmissivity')
        constitutive_nonlocal_surfaceTransmissivity(i) = IO_floatValue(line,positions,2)
      case default
        call IO_error(236,ext_msg=tag)
    end select
  endif
enddo


100 do i = 1,maxNinstance

  constitutive_nonlocal_structure(i) = &
    lattice_initializeStructure(constitutive_nonlocal_structureName(i), constitutive_nonlocal_CoverA(i))                            ! our lattice structure is defined in the material.config file by the structureName (and the c/a ratio)
  myStructure = constitutive_nonlocal_structure(i)
  
  
  !*** sanity checks
  
  if (myStructure < 1 .or. myStructure > 3)                                       call IO_error(205)
  if (sum(constitutive_nonlocal_Nslip(:,i)) <= 0_pInt)                            call IO_error(235,ext_msg='Nslip')
  do o = 1,maxval(phase_Noutput)
    if(len(constitutive_nonlocal_output(o,i)) > 64)                               call IO_error(666)
  enddo
  do f = 1,lattice_maxNslipFamily
    if (constitutive_nonlocal_Nslip(f,i) > 0_pInt) then
      if (constitutive_nonlocal_rhoSglEdgePos0(f,i) < 0.0_pReal)                  call IO_error(235,ext_msg='rhoSglEdgePos0')
      if (constitutive_nonlocal_rhoSglEdgeNeg0(f,i) < 0.0_pReal)                  call IO_error(235,ext_msg='rhoSglEdgeNeg0')
      if (constitutive_nonlocal_rhoSglScrewPos0(f,i) < 0.0_pReal)                 call IO_error(235,ext_msg='rhoSglScrewPos0')
      if (constitutive_nonlocal_rhoSglScrewNeg0(f,i) < 0.0_pReal)                 call IO_error(235,ext_msg='rhoSglScrewNeg0')
      if (constitutive_nonlocal_rhoDipEdge0(f,i) < 0.0_pReal)                     call IO_error(235,ext_msg='rhoDipEdge0')
      if (constitutive_nonlocal_rhoDipScrew0(f,i) < 0.0_pReal)                    call IO_error(235,ext_msg='rhoDipScrew0')
      if (constitutive_nonlocal_burgersPerSlipFamily(f,i) <= 0.0_pReal)           call IO_error(235,ext_msg='burgers')
      if (constitutive_nonlocal_lambda0PerSlipFamily(f,i) <= 0.0_pReal)           call IO_error(235,ext_msg='lambda0')
      if (constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i) <= 0.0_pReal)        call IO_error(235,ext_msg='dDipMinEdge')
      if (constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i) <= 0.0_pReal)       call IO_error(235,ext_msg='dDipMinScrew')
    endif
  enddo
  if (any(constitutive_nonlocal_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 0.0_pReal)) &
                                                                                  call IO_error(235,ext_msg='interaction_SlipSlip')
  if (constitutive_nonlocal_R(i) < 0.0_pReal)                                     call IO_error(235,ext_msg='r')
  if (constitutive_nonlocal_atomicVolume(i) <= 0.0_pReal)                         call IO_error(235,ext_msg='atomicVolume')
  if (constitutive_nonlocal_Dsd0(i) <= 0.0_pReal)                                 call IO_error(235,ext_msg='Dsd0')
  if (constitutive_nonlocal_Qsd(i) <= 0.0_pReal)                                  call IO_error(235,ext_msg='Qsd')
  if (constitutive_nonlocal_aTolRho(i) <= 0.0_pReal)                              call IO_error(235,ext_msg='aTol_rho')
  if (constitutive_nonlocal_d0(i) <= 0.0_pReal)                                   call IO_error(235,ext_msg='d0')
  if (constitutive_nonlocal_tauObs(i) <= 0.0_pReal)                               call IO_error(235,ext_msg='tauObs')
  if (constitutive_nonlocal_vs(i) <= 0.0_pReal)                                   call IO_error(235,ext_msg='vs')
  if (constitutive_nonlocal_fattack(i) <= 0.0_pReal)                              call IO_error(235,ext_msg='fAttack')
  if (constitutive_nonlocal_rhoSglScatter(i) < 0.0_pReal)                         call IO_error(235,ext_msg='rhoSglScatter')
  if (constitutive_nonlocal_surfaceTransmissivity(i) < 0.0_pReal &
      .or. constitutive_nonlocal_surfaceTransmissivity(i) > 1.0_pReal)            call IO_error(235,ext_msg='surfaceTransmissivity')
  
  
  !*** determine total number of active slip systems
  
  constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i) = min( lattice_NslipSystem(1:lattice_maxNslipFamily, myStructure), &
                                                                constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i) )           ! we can't use more slip systems per family than specified in lattice 
  constitutive_nonlocal_totalNslip(i) = sum(constitutive_nonlocal_Nslip(1:lattice_maxNslipFamily,i))

enddo


!*** allocation of variables whose size depends on the total number of active slip systems

maxTotalNslip = maxval(constitutive_nonlocal_totalNslip)

allocate(constitutive_nonlocal_burgersPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_burgersPerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_lambda0PerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_lambda0PerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_dLowerEdgePerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dLowerEdgePerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_dLowerScrewPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_nonlocal_dLowerScrewPerSlipSystem = 0.0_pReal

allocate(constitutive_nonlocal_Qeff0(maxTotalNslip, maxNinstance))
constitutive_nonlocal_Qeff0 = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionEdge(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionEdge = 0.0_pReal

allocate(constitutive_nonlocal_forestProjectionScrew(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_forestProjectionScrew = 0.0_pReal

allocate(constitutive_nonlocal_interactionMatrixSlipSlip(maxTotalNslip, maxTotalNslip, maxNinstance))
constitutive_nonlocal_interactionMatrixSlipSlip = 0.0_pReal

allocate(constitutive_nonlocal_v(maxTotalNslip, 4, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_v = 0.0_pReal

allocate(constitutive_nonlocal_rhoDotFlux(maxTotalNslip, 10, homogenization_maxNgrains, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_rhoDotFlux = 0.0_pReal

allocate(constitutive_nonlocal_compatibility(2,maxTotalNslip, maxTotalNslip, FE_maxNipNeighbors, mesh_maxNips, mesh_NcpElems))
constitutive_nonlocal_compatibility = 0.0_pReal

do i = 1,maxNinstance
  
  myStructure = constitutive_nonlocal_structure(i)                                                                                  ! lattice structure of this instance
    

  !*** Inverse lookup of my slip system family and the slip system in lattice
  
  l = 0_pInt
  do f = 1,lattice_maxNslipFamily
    do s = 1,constitutive_nonlocal_Nslip(f,i)
      l = l + 1
      constitutive_nonlocal_slipFamily(l,i) = f
      constitutive_nonlocal_slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1, myStructure)) + s
  enddo; enddo
  
  
  !*** determine size of state array
  
  ns = constitutive_nonlocal_totalNslip(i)
  constitutive_nonlocal_sizeState(i) = size(constitutive_nonlocal_listBasicStates) * ns & 
                                       + (size(constitutive_nonlocal_listDependentStates) - 1_pInt) * ns + 6_pInt
  constitutive_nonlocal_sizeDotState(i) = size(constitutive_nonlocal_listBasicStates) * ns
  
  
  !*** determine size of postResults array
  
  do o = 1,maxval(phase_Noutput)
    select case(constitutive_nonlocal_output(o,i))
      case( 'rho', &
            'delta', &
            'rho_edge', &
            'rho_screw', &
            'rho_sgl', &
            'delta_sgl', &
            'rho_sgl_edge', &
            'rho_sgl_edge_pos', &
            'rho_sgl_edge_neg', &
            'rho_sgl_screw', &
            'rho_sgl_screw_pos', &
            'rho_sgl_screw_neg', &
            'rho_sgl_mobile', &
            'rho_sgl_edge_mobile', &
            'rho_sgl_edge_pos_mobile', &
            'rho_sgl_edge_neg_mobile', &
            'rho_sgl_screw_mobile', &
            'rho_sgl_screw_pos_mobile', &
            'rho_sgl_screw_neg_mobile', &
            'rho_sgl_immobile', &
            'rho_sgl_edge_immobile', &
            'rho_sgl_edge_pos_immobile', &
            'rho_sgl_edge_neg_immobile', &
            'rho_sgl_screw_immobile', &
            'rho_sgl_screw_pos_immobile', &
            'rho_sgl_screw_neg_immobile', &
            'rho_dip', &
            'delta_dip', &
            'rho_dip_edge', &
            'rho_dip_screw', &
            'excess_rho', &
            'excess_rho_edge', &
            'excess_rho_screw', &
            'rho_forest', &
            'shearrate', &
            'resolvedstress', &
            'resolvedstress_internal', &
            'resolvedstress_external', &
            'resistance', &
            'rho_dot', &
            'rho_dot_sgl', &
            'rho_dot_dip', &
            'rho_dot_gen', &
            'rho_dot_gen_edge', &
            'rho_dot_gen_screw', &
            'rho_dot_sgl2dip', &
            'rho_dot_ann_ath', &
            'rho_dot_ann_the', &
            'rho_dot_flux', &
            'rho_dot_flux_edge', &
            'rho_dot_flux_screw', &
            'dislocationvelocity', &
            'fluxdensity_edge_pos_x', &
            'fluxdensity_edge_pos_y', &
            'fluxdensity_edge_pos_z', &
            'fluxdensity_edge_neg_x', &
            'fluxdensity_edge_neg_y', &
            'fluxdensity_edge_neg_z', &
            'fluxdensity_screw_pos_x', &
            'fluxdensity_screw_pos_y', &
            'fluxdensity_screw_pos_z', &
            'fluxdensity_screw_neg_x', &
            'fluxdensity_screw_neg_y', &
            'fluxdensity_screw_neg_z', &
            'd_upper_edge', &
            'd_upper_screw' )
        mySize = constitutive_nonlocal_totalNslip(i)
      case default
        mySize = 0_pInt
    end select

    if (mySize > 0_pInt) then                                                                                                       ! any meaningful output found                               
      constitutive_nonlocal_sizePostResult(o,i) = mySize
      constitutive_nonlocal_sizePostResults(i)  = constitutive_nonlocal_sizePostResults(i) + mySize
    endif
  enddo
  
  
  !*** elasticity matrix and shear modulus according to material.config
  
  select case (myStructure)
    case(1:2)                                                                                                                       ! cubic(s)
      forall(k=1:3)
        forall(j=1:3) constitutive_nonlocal_Cslip_66(k,j,i) = constitutive_nonlocal_C12(i)
        constitutive_nonlocal_Cslip_66(k,k,i) = constitutive_nonlocal_C11(i)
        constitutive_nonlocal_Cslip_66(k+3,k+3,i) = constitutive_nonlocal_C44(i)
      end forall
    case(3:)                                                                                                                        ! all hex
      constitutive_nonlocal_Cslip_66(1,1,i) = constitutive_nonlocal_C11(i)
      constitutive_nonlocal_Cslip_66(2,2,i) = constitutive_nonlocal_C11(i)
      constitutive_nonlocal_Cslip_66(3,3,i) = constitutive_nonlocal_C33(i)
      constitutive_nonlocal_Cslip_66(1,2,i) = constitutive_nonlocal_C12(i)
      constitutive_nonlocal_Cslip_66(2,1,i) = constitutive_nonlocal_C12(i)
      constitutive_nonlocal_Cslip_66(1,3,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(3,1,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(2,3,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(3,2,i) = constitutive_nonlocal_C13(i)
      constitutive_nonlocal_Cslip_66(4,4,i) = constitutive_nonlocal_C44(i)
      constitutive_nonlocal_Cslip_66(5,5,i) = constitutive_nonlocal_C44(i)
      constitutive_nonlocal_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_nonlocal_C11(i)- constitutive_nonlocal_C12(i))
  end select
  constitutive_nonlocal_Cslip_66(1:6,1:6,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_nonlocal_Cslip_66(1:6,1:6,i)))
  constitutive_nonlocal_Cslip_3333(1:3,1:3,1:3,1:3,i) = math_Voigt66to3333(constitutive_nonlocal_Cslip_66(1:6,1:6,i))

  constitutive_nonlocal_Gmod(i) = 0.2_pReal * ( constitutive_nonlocal_C11(i) - constitutive_nonlocal_C12(i) &
                                                + 3.0_pReal*constitutive_nonlocal_C44(i) )                                          ! (C11iso-C12iso)/2 with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  constitutive_nonlocal_nu(i) =   ( constitutive_nonlocal_C11(i) + 4.0_pReal*constitutive_nonlocal_C12(i) &
                                    - 2.0_pReal*constitutive_nonlocal_C44(i) ) &
                                / ( 4.0_pReal*constitutive_nonlocal_C11(i) + 6.0_pReal*constitutive_nonlocal_C12(i) &
                                    + 2.0_pReal*constitutive_nonlocal_C44(i) )                                                      ! C12iso/(C11iso+C12iso) with C11iso=(3*C11+2*C12+4*C44)/5 and C12iso=(C11+4*C12-2*C44)/5
  
  do s1 = 1,ns      
    f = constitutive_nonlocal_slipFamily(s1,i)
    
    !*** burgers vector, mean free path prefactor and minimum dipole distance for each slip system
  
    constitutive_nonlocal_burgersPerSlipSystem(s1,i) = constitutive_nonlocal_burgersPerSlipFamily(f,i)
    constitutive_nonlocal_lambda0PerSlipSystem(s1,i) = constitutive_nonlocal_lambda0PerSlipFamily(f,i)
    constitutive_nonlocal_dLowerEdgePerSlipSystem(s1,i) = constitutive_nonlocal_dLowerEdgePerSlipFamily(f,i)
    constitutive_nonlocal_dLowerScrewPerSlipSystem(s1,i) = constitutive_nonlocal_dLowerScrewPerSlipFamily(f,i)

    do s2 = 1,ns
      
      !*** calculation of forest projections for edge and screw dislocations. s2 acts as forest for s1

      constitutive_nonlocal_forestProjectionEdge(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,constitutive_nonlocal_slipSystemLattice(s1,i),myStructure), &
                            lattice_st(1:3,constitutive_nonlocal_slipSystemLattice(s2,i),myStructure)))                             ! forest projection of edge dislocations is the projection of (t = b x n) onto the slip normal of the respective slip plane
      
      constitutive_nonlocal_forestProjectionScrew(s1,s2,i) &
          = abs(math_mul3x3(lattice_sn(1:3,constitutive_nonlocal_slipSystemLattice(s1,i),myStructure), &
                            lattice_sd(1:3,constitutive_nonlocal_slipSystemLattice(s2,i),myStructure)))                             ! forest projection of screw dislocations is the projection of b onto the slip normal of the respective splip plane
  
      !*** calculation of interaction matrices

      constitutive_nonlocal_interactionMatrixSlipSlip(s1,s2,i) &
          = constitutive_nonlocal_interactionSlipSlip(lattice_interactionSlipSlip(constitutive_nonlocal_slipSystemLattice(s1,i), &
                                                                                  constitutive_nonlocal_slipSystemLattice(s2,i), &
                                                                                  myStructure), &
                                                      i)
  
    enddo
  enddo
  
  
  !*** calculation of prefactor for activation enthalpy for dislocation glide

  constitutive_nonlocal_Qeff0(1:ns,i) = constitutive_nonlocal_burgersPerSlipSystem(1:ns,i) ** 3.0_pReal &
                                      * sqrt(0.5_pReal * constitutive_nonlocal_d0(i) ** 3.0_pReal &
                                                       * constitutive_nonlocal_Gmod(i) * constitutive_nonlocal_tauObs(i))

enddo

endsubroutine



!*********************************************************************
!* initial microstructural state (just the "basic" states)           *
!*********************************************************************
function constitutive_nonlocal_stateInit(myInstance)

use prec,     only: pReal, &
                    pInt
use lattice,  only: lattice_maxNslipFamily
use math,     only: math_sampleGaussVar

implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_stateInit

!*** local variables
real(pReal), dimension(constitutive_nonlocal_totalNslip(myInstance)) :: &              
                              rhoSglEdgePos, &                ! positive edge dislocation density
                              rhoSglEdgeNeg, &                ! negative edge dislocation density
                              rhoSglScrewPos, &               ! positive screw dislocation density
                              rhoSglScrewNeg, &               ! negative screw dislocation density
                              rhoSglEdgePosUsed, &            ! used positive edge dislocation density
                              rhoSglEdgeNegUsed, &            ! used negative edge dislocation density
                              rhoSglScrewPosUsed, &           ! used positive screw dislocation density
                              rhoSglScrewNegUsed, &           ! used negative screw dislocation density
                              rhoDipEdge, &                   ! edge dipole dislocation density
                              rhoDipScrew, &                  ! screw dipole dislocation density
                              rhoForest, &                    ! forest dislocation density
                              tauSlipThreshold                ! threshold shear stress for slip
integer(pInt)                 ns, &                           ! short notation for total number of active slip systems 
                              f, &                            ! index of lattice family
                              from, &
                              upto, &
                              s, &                            ! index of slip system
                              i
real(pReal), dimension(2) ::  noise

constitutive_nonlocal_stateInit = 0.0_pReal
ns = constitutive_nonlocal_totalNslip(myInstance)

!*** set the basic state variables

do f = 1,lattice_maxNslipFamily
  from = 1 + sum(constitutive_nonlocal_Nslip(1:f-1,myInstance))
  upto = sum(constitutive_nonlocal_Nslip(1:f,myInstance))
  do s = from,upto
    do i = 1,2
      noise(i) = math_sampleGaussVar(0.0_pReal, constitutive_nonlocal_rhoSglScatter(myInstance))
    enddo
    rhoSglEdgePos(s) = constitutive_nonlocal_rhoSglEdgePos0(f, myInstance) + noise(1)
    rhoSglEdgeNeg(s) = constitutive_nonlocal_rhoSglEdgeNeg0(f, myInstance) + noise(1)
    rhoSglScrewPos(s) = constitutive_nonlocal_rhoSglScrewPos0(f, myInstance) + noise(2)
    rhoSglScrewNeg(s) = constitutive_nonlocal_rhoSglScrewNeg0(f, myInstance) + noise(2)
  enddo 
  rhoSglEdgePosUsed(from:upto)  = 0.0_pReal
  rhoSglEdgeNegUsed(from:upto)  = 0.0_pReal
  rhoSglScrewPosUsed(from:upto) = 0.0_pReal
  rhoSglScrewNegUsed(from:upto) = 0.0_pReal
  rhoDipEdge(from:upto)  = constitutive_nonlocal_rhoDipEdge0(f, myInstance)
  rhoDipScrew(from:upto) = constitutive_nonlocal_rhoDipScrew0(f, myInstance)
enddo


!*** put everything together and in right order

constitutive_nonlocal_stateInit(      1:   ns) = rhoSglEdgePos
constitutive_nonlocal_stateInit(   ns+1: 2*ns) = rhoSglEdgeNeg
constitutive_nonlocal_stateInit( 2*ns+1: 3*ns) = rhoSglScrewPos
constitutive_nonlocal_stateInit( 3*ns+1: 4*ns) = rhoSglScrewNeg
constitutive_nonlocal_stateInit( 4*ns+1: 5*ns) = rhoSglEdgePosUsed
constitutive_nonlocal_stateInit( 5*ns+1: 6*ns) = rhoSglEdgeNegUsed
constitutive_nonlocal_stateInit( 6*ns+1: 7*ns) = rhoSglScrewPosUsed
constitutive_nonlocal_stateInit( 7*ns+1: 8*ns) = rhoSglScrewNegUsed
constitutive_nonlocal_stateInit( 8*ns+1: 9*ns) = rhoDipEdge
constitutive_nonlocal_stateInit( 9*ns+1:10*ns) = rhoDipScrew

endfunction



!*********************************************************************
!* absolute state tolerance                                          *
!*********************************************************************
pure function constitutive_nonlocal_aTolState(myInstance)

use prec,     only: pReal, &
                    pInt
implicit none

!*** input variables
integer(pInt), intent(in) ::  myInstance                      ! number specifying the current instance of the constitution

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizeState(myInstance)) :: &
                              constitutive_nonlocal_aTolState ! absolute state tolerance for the current instance of this constitution

!*** local variables

constitutive_nonlocal_aTolState = constitutive_nonlocal_aTolRho(myInstance)

endfunction



!*********************************************************************
!* calculates homogenized elacticity matrix                          *
!*********************************************************************
pure function constitutive_nonlocal_homogenizedC(state,g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
implicit none

!*** input variables
integer(pInt), intent(in) ::    g, &                                ! current grain ID
                                ip, &                               ! current integration point
                                el                                  ! current element
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state ! microstructural state

!*** output variables
real(pReal), dimension(6,6) ::  constitutive_nonlocal_homogenizedC  ! homogenized elasticity matrix

!*** local variables
integer(pInt)                   myInstance                          ! current instance of this constitution

myInstance = phase_constitutionInstance(material_phase(g,ip,el))

constitutive_nonlocal_homogenizedC = constitutive_nonlocal_Cslip_66(1:6,1:6,myInstance)
 
endfunction



!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*********************************************************************
subroutine constitutive_nonlocal_microstructure(state, Temperature, Fe, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use IO,       only: IO_error
use math,     only: math_Mandel33to6, &
                    math_mul33x33, &
                    math_mul33x3, &
                    math_inv3x3, &
                    math_transpose3x3, &
                    pi
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element, &
                    mesh_node, &
                    FE_Nips, &
                    mesh_ipCenterOfGravity, &
                    mesh_ipVolume, &
                    mesh_periodicSurface
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_localConstitution, &
                    phase_constitutionInstance
use lattice,  only: lattice_sd, &
                    lattice_sn, &
                    lattice_st

implicit none

!*** input variables
integer(pInt), intent(in) ::    g, &                          ! current grain ID
                                ip, &                         ! current integration point
                                el                            ! current element
real(pReal), intent(in) ::      Temperature                   ! temperature
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                Fe                            ! elastic deformation gradient

  !*** input/output variables
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
                                state                         ! microstructural state

  !*** output variables

  !*** local variables
integer(pInt)                   neighboring_el, &             ! element number of neighboring material point
                                neighboring_ip, &             ! integration point of neighboring material point
                                instance, &                   ! my instance of this constitution
                                neighboring_instance, &       ! instance of this constitution of neighboring material point
                                latticeStruct, &              ! my lattice structure
                                neighboring_latticeStruct, &  ! lattice structure of neighboring material point
                                phase, &
                                neighboring_phase, &
                                ns, &                         ! total number of active slip systems at my material point
                                neighboring_ns, &             ! total number of active slip systems at neighboring material point
                                c, &                          ! index of dilsocation character (edge, screw)
                                s, &                          ! slip system index
                                s2, &                         ! slip system index according to ordering in "lattice.f90"
                                t, &                          ! index of dilsocation type (e+, e-, s+, s-, used e+, used e-, used s+, used s-)
                                dir, &
                                deltaX, deltaY, deltaZ, &
                                side
integer(pInt), dimension(2,3) :: periodicImages
real(pReal)                     nu, &                         ! poisson's ratio
                                x, y, z, &                    ! coordinates of connection vector in neighboring lattice frame
                                xsquare, ysquare, zsquare, &  ! squares of respective coordinates
                                distance, &                   ! length of connection vector
                                segmentLength, &              ! segment length of dislocations
                                neighboring_Nexcess, &        ! excess number of dislocation segments at neighboring material point for specific slip system and dislocation character
                                lambda, &
                                R, Rsquare, Rcube, &
                                denominator, &
                                flipSign
real(pReal), dimension(3) ::    connection, &                 ! connection vector between me and my neighbor in the deformed configuration
                                connection_neighboringLattice, & ! connection vector between me and my neighbor in the lattice configuration of my neighbor
                                connection_neighboringSlip, & ! connection vector between me and my neighbor in the slip system frame of my neighbor
                                maxCoord, minCoord, &
                                meshSize, &
                                ipCoords, &
                                neighboring_ipCoords
real(pReal), dimension(3,3) ::  sigma, &                      ! dislocation stress for one slip system in neighboring material point's slip system frame
                                neighboringLattice2neighboringSlip, & ! orthogonal transformation matrix from lattice coordinate system to slip coordinate system (passive rotation  !  !  !)
                                Tdislo_neighboringLattice, &  ! dislocation stress as 2nd Piola-Kirchhoff stress at neighboring material point
                                Tdislo, &                     ! dislocation stress as 2nd Piola-Kirchhoff stress at my material point
                                invFe, &                      ! inverse of my elastic deformation gradient
                                neighboringLattice2myLattice  ! mapping from neighboring MPs lattice configuration to my lattice configuration
real(pReal), dimension(2,maxval(constitutive_nonlocal_totalNslip)) :: &
                                neighboring_rhoExcess         ! excess density at neighboring material point
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                rhoSgl                        ! single dislocation density (edge+, edge-, screw+, screw-, used edge+, used edge-, used screw+, used screw-)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                rhoDip                        ! dipole dislocation density (edge, screw)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                rhoForest, &                  ! forest dislocation density
                                tauThreshold, &               ! threshold shear stress
                                tau                           ! resolved shear stress

phase = material_phase(g,ip,el)
instance = phase_constitutionInstance(phase)
latticeStruct = constitutive_nonlocal_structure(instance)
ns = constitutive_nonlocal_totalNslip(instance)



!*** get basic states

forall (t = 1:4) rhoSgl(1:ns,t) = max(state(g,ip,el)%p((t-1)*ns+1:t*ns), 0.0_pReal)                       ! ensure positive single mobile densities
forall (t = 5:8) rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(1:ns,c) = max(state(g,ip,el)%p((c+7)*ns+1:(c+8)*ns), 0.0_pReal)                   ! ensure positive dipole densities
where(rhoSgl(1:ns,1:4) < min(0.1, 0.01*constitutive_nonlocal_aTolRho(instance))) &
  rhoSgl(1:ns,1:4) = 0.0_pReal                                                                            ! delete non-significant single density



!*** calculate the forest dislocation density
!*** (= projection of screw and edge dislocations)

forall (s = 1:ns) &
  rhoForest(s) = dot_product((sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2) + rhoDip(1:ns,1)), &
                              constitutive_nonlocal_forestProjectionEdge(s,1:ns,instance)) & 
                + dot_product((sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2) + rhoDip(1:ns,2)), &
                              constitutive_nonlocal_forestProjectionScrew(s,1:ns,instance))



!*** calculate the threshold shear stress for dislocation slip 

forall (s = 1:ns) &
  tauThreshold(s) =   constitutive_nonlocal_Gmod(instance) & 
                    * constitutive_nonlocal_burgersPerSlipSystem(s,instance) &
                    * sqrt(dot_product((sum(abs(rhoSgl),2) + sum(abs(rhoDip),2)), &
                                        constitutive_nonlocal_interactionMatrixSlipSlip(s,1:ns,instance)))



!*** calculate the dislocation stress of the neighboring excess dislocation densities
!*** zero for material points of local constitution

Tdislo = 0.0_pReal

if (.not. phase_localConstitution(phase)) then
  invFe = math_inv3x3(Fe(1:3,1:3,1,ip,el))


  !* in case of periodic surfaces we have to find out how many periodic images in each direction we need
  
  do dir = 1,3
    maxCoord(dir) = maxval(mesh_node(dir,:))
    minCoord(dir) = minval(mesh_node(dir,:))
  enddo
  meshSize = maxCoord - minCoord
  ipCoords = mesh_ipCenterOfGravity(1:3,ip,el)
  periodicImages = 0_pInt
  do dir = 1,3
    if (mesh_periodicSurface(dir)) then
      periodicImages(1,dir) =   floor((ipCoords(dir) - constitutive_nonlocal_R(instance) - minCoord(dir)) / meshSize(dir), pInt)
      periodicImages(2,dir) = ceiling((ipCoords(dir) + constitutive_nonlocal_R(instance) - maxCoord(dir)) / meshSize(dir), pInt)
    endif
  enddo

      
  !* loop through all material points (also through their periodic images if present), 
  !* but only consider nonlocal neighbors within a certain cutoff radius R
  
  do neighboring_el = 1,mesh_NcpElems
    do neighboring_ip = 1,FE_Nips(mesh_element(2,neighboring_el))
      neighboring_phase = material_phase(g,neighboring_ip,neighboring_el)
      neighboring_instance = phase_constitutionInstance(neighboring_phase)
      neighboring_latticeStruct = constitutive_nonlocal_structure(neighboring_instance)
      neighboring_ns = constitutive_nonlocal_totalNslip(neighboring_instance)
      nu = constitutive_nonlocal_nu(neighboring_instance)
      Tdislo_neighboringLattice = 0.0_pReal
      do deltaX = periodicImages(1,1),periodicImages(2,1)
        do deltaY = periodicImages(1,2),periodicImages(2,2)
          do deltaZ = periodicImages(1,3),periodicImages(2,3)
            if (neighboring_el == el .and. neighboring_ip == ip &
                .and. deltaX == 0 .and. deltaY == 0 .and. deltaZ == 0) then
              cycle                                                                                 ! this is myself
            endif
            neighboring_ipCoords = mesh_ipCenterOfGravity(1:3,neighboring_ip,neighboring_el) &
                                 + (/real(deltaX,pReal), real(deltaY,pReal), real(deltaZ,pReal)/) * meshSize
            connection = neighboring_ipCoords - ipCoords
            distance = sqrt(sum(connection ** 2.0_pReal))
            if (.not. phase_localConstitution(neighboring_phase) &
                .and. distance <= constitutive_nonlocal_R(instance)) then                   
              

              !* determine the effective number of dislocations
              !* the segment length is the minimum of the third root of the control volume and the ip distance
              !* this ensures, that the central MP never sits on a neighboring dislocation segment
              
              connection_neighboringLattice = math_mul33x3(math_inv3x3(Fe(1:3,1:3,1,neighboring_ip,neighboring_el)), connection)
              forall (s = 1:neighboring_ns, c = 1:2) &
                neighboring_rhoExcess(c,s) =      state(g,neighboring_ip,neighboring_el)%p((2*c-2)*neighboring_ns+s) &
                                            + abs(state(g,neighboring_ip,neighboring_el)%p((2*c+2)*neighboring_ns+s)) &
                                            -     state(g,neighboring_ip,neighboring_el)%p((2*c-1)*neighboring_ns+s) &
                                            - abs(state(g,neighboring_ip,neighboring_el)%p((2*c+3)*neighboring_ns+s))
              segmentLength = min(mesh_ipVolume(neighboring_ip,neighboring_el)**(1.0_pReal/3.0_pReal), distance)
      

              !* loop through all slip systems of the neighboring material point
              !* and add up the stress contributions from egde and screw excess on these slip systems
      
              do s = 1,neighboring_ns
                if (all(abs(neighboring_rhoExcess(:,s)) < 1.0_pReal)) then
                  cycle                                                                             ! not significant
                endif
                
                
                !* map the connection vector from the lattice into the slip system frame
                
                s2 = constitutive_nonlocal_slipSystemLattice(s,neighboring_instance)
                neighboringLattice2neighboringSlip = math_transpose3x3( &
                                      reshape((/ lattice_sd(1:3, s2, neighboring_latticeStruct), &
                                                -lattice_st(1:3, s2, neighboring_latticeStruct), &
                                                 lattice_sn(1:3, s2, neighboring_latticeStruct)/), (/3,3/)))
                connection_neighboringSlip = math_mul33x3(neighboringLattice2neighboringSlip, connection_neighboringLattice)
                x = connection_neighboringSlip(1)
                y = connection_neighboringSlip(2)
                z = connection_neighboringSlip(3)
                xsquare = x ** 2.0_pReal
                ysquare = y ** 2.0_pReal
                zsquare = z ** 2.0_pReal
                
                
                !* edge contribution to stress
                
                sigma = 0.0_pReal
                neighboring_Nexcess = neighboring_rhoExcess(1,s) * mesh_ipVolume(neighboring_ip,neighboring_el) / segmentLength
                flipSign = sign(1.0_pReal, -y)
                do side = 1,-1,-2
                  lambda = real(side,pReal) * 0.5_pReal * segmentLength - y
                  R = sqrt(xsquare + zsquare + lambda**2.0_pReal)
                  Rsquare = R ** 2.0_pReal
                  Rcube = R**3.0_pReal
                  denominator = R * (R + flipSign * lambda)
                  if (denominator == 0.0_pReal) then
                    call IO_error(237,el,ip,g)
                  endif
                    
                  sigma(1,1) = sigma(1,1) - real(side,pReal) * flipSign * z / denominator &
                                                             * (1.0_pReal + xsquare / Rsquare + xsquare / denominator) &
                                                             * neighboring_Nexcess
                  sigma(2,2) = sigma(2,2) - real(side,pReal) * (flipSign * 2.0_pReal * nu * z / denominator + z * lambda / Rcube) &
                                                             * neighboring_Nexcess
                  sigma(3,3) = sigma(3,3) + real(side,pReal) * flipSign * z / denominator &
                                                             * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                                             * neighboring_Nexcess
                  sigma(1,2) = sigma(1,2) + real(side,pReal) * x * z / Rcube * neighboring_Nexcess
                  sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * x / denominator &
                                                             * (1.0_pReal - zsquare / Rsquare - zsquare / denominator) &
                                                             * neighboring_Nexcess
                  sigma(2,3) = sigma(2,3) - real(side,pReal) * (nu / R - zsquare / Rcube) * neighboring_Nexcess
                enddo
                
                
                !* screw contribution to stress
                
                neighboring_Nexcess = neighboring_rhoExcess(2,s) * mesh_ipVolume(neighboring_ip,neighboring_el) / segmentLength
                flipSign = sign(1.0_pReal, x)
                do side = 1,-1,-2
                  lambda = x + real(side,pReal) * 0.5_pReal * segmentLength
                  R = sqrt(ysquare + zsquare + lambda**2.0_pReal)
                  Rsquare = R ** 2.0_pReal
                  Rcube = R**3.0_pReal
                  denominator = R * (R + flipSign * lambda)
                  if (denominator == 0.0_pReal) then
                    call IO_error(237,el,ip,g)
                  endif
                  
                  sigma(1,2) = sigma(1,2) - real(side,pReal) * flipSign * z * (1.0_pReal - nu) / denominator * neighboring_Nexcess
                  sigma(1,3) = sigma(1,3) + real(side,pReal) * flipSign * y * (1.0_pReal - nu) / denominator * neighboring_Nexcess
                enddo
                

                !* copy symmetric parts
                
                sigma(2,1) = sigma(1,2)
                sigma(3,1) = sigma(1,3)
                sigma(3,2) = sigma(2,3)

                
                !* scale stresses and map them into the neighboring material point's lattice configuration
                
                sigma = sigma * constitutive_nonlocal_Gmod(neighboring_instance) &
                              * constitutive_nonlocal_burgersPerSlipSystem(s,neighboring_instance) &
                              / (4.0_pReal * pi * (1.0_pReal - nu))
                Tdislo_neighboringLattice = Tdislo_neighboringLattice &
                                          + math_mul33x33(math_transpose3x3(neighboringLattice2neighboringSlip), &
                                            math_mul33x33(sigma, neighboringLattice2neighboringSlip))
                                            
              enddo ! slip system loop
            endif
          enddo ! deltaZ loop
        enddo ! deltaY loop
      enddo ! deltaX loop


      !* map the stress from the neighboring MP's lattice configuration into the deformed configuration 
      !* and back into my lattice configuration

      neighboringLattice2myLattice = math_mul33x33(invFe, Fe(1:3,1:3,1,neighboring_ip,neighboring_el))
      Tdislo = Tdislo + math_mul33x33(neighboringLattice2myLattice, &
                        math_mul33x33(Tdislo_neighboringLattice, math_transpose3x3(neighboringLattice2myLattice)))
                        
    enddo ! ip loop
  enddo ! element loop
    
endif


!*** set states

state(g,ip,el)%p(1:8*ns) = reshape(rhoSgl,(/8*ns/))                                                       ! ensure positive single mobile densities
state(g,ip,el)%p(8*ns+1:10*ns) = reshape(rhoDip,(/2*ns/))                                                 ! ensure positive dipole densities
state(g,ip,el)%p(10*ns+1:11*ns) = rhoForest
state(g,ip,el)%p(11*ns+1:12*ns) = tauThreshold
state(g,ip,el)%p(12*ns+1:12*ns+6) = math_Mandel33to6(Tdislo)


#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i5,x,i2,x,i1)') '<< CONST >> nonlocal_microstructure at el ip g',el,ip,g
    write(6,*)
    write(6,'(a,/,12(x),12(e10.3,x))') '<< CONST >> rhoForest', rhoForest
    write(6,'(a,/,12(x),12(f10.5,x))') '<< CONST >> tauThreshold / MPa', tauThreshold/1e6
    write(6,'(a,/,3(12(x),3(f10.5,x),/))') '<< CONST >> Tdislocation / MPa', Tdislo/1e6
  endif
#endif

endsubroutine



!*********************************************************************
!* calculates kinetics                                               *
!*********************************************************************
subroutine constitutive_nonlocal_kinetics(Tstar_v, Temperature, state, g, ip, el, dv_dtau)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                        ! current grain number
                                            ip, &                       ! current integration point
                                            el                          ! current element number
real(pReal), intent(in) ::                  Temperature                 ! temperature
type(p_vec), intent(in) ::                  state                       ! microstructural state
real(pReal), dimension(6), intent(in) ::    Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation

!*** output variables
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))), &
                   intent(out), optional :: dv_dtau                     ! velocity derivative with respect to resolved shear stress

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            t, &                        ! dislocation type
                                            s                           ! index of my current slip system
real(pReal), dimension(6) ::                Tdislocation_v              ! dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            tau, &                      ! resolved shear stress
                                            rhoForest                   ! forest dislocation density
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            v                           ! velocity for the current element and ip
real(pReal)                                 boltzmannProbability, &
                                            tauRel, &                   ! relative thermally active resolved shear stress
                                            wallFunc, &                 ! functions reflecting the shape of the obstacle wall (see PhD thesis Mohles p.53)
                                            timeRatio                   ! ratio of travel to dwell time


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

rhoForest = state%p(10*ns+1:11*ns)
tauThreshold = state%p(11*ns+1:12*ns)
Tdislocation_v = state%p(12*ns+1:12*ns+6)

tau = 0.0_pReal
v = 0.0_pReal
if (present(dv_dtau)) dv_dtau = 0.0_pReal


if (Temperature > 0.0_pReal) then
  do s = 1,ns
  
    tau(s) = math_mul6x6(Tstar_v + Tdislocation_v, &
                         lattice_Sslip_v(:,constitutive_nonlocal_slipSystemLattice(s,myInstance),myStructure))

    !*** Only if the resolved shear stress exceeds the threshold stress, dislocations are able to cut the dislocation forest.
    !*** In contrast to small atomic obstacles the forest can't be overcome by thermal activation.
    !*** 
    !***                                                     mean travel distance
    !*** The mean dislocation velocity is calculated as:  --------------------------
    !***                                                   dwell time + travel time
    !***
    !*** with :   mean travel distance = inverse of the root of forest density
    !***          dwell time = inverse of attack frequency times probability of success
    !***          travel time = mean travel distance over velocity of sound
    
    tauRel = (abs(tau(s)) - tauThreshold(s)) / constitutive_nonlocal_tauObs(myInstance)
    if (tauRel > 0.0_pReal .and. tauRel < 1.0_pReal) then
      wallFunc = 4.0_pReal * sqrt(2.0_pReal) / 3.0_pReal * sqrt(1.0_pReal - tauRel) / tauRel
      boltzmannProbability = exp(- constitutive_nonlocal_Qeff0(s,myInstance) * wallFunc / (kB * Temperature))
      timeRatio = boltzmannProbability * constitutive_nonlocal_fattack(myInstance) &
                / (constitutive_nonlocal_vs(myInstance) * sqrt(rhoForest(s)))
        
      v(s,1:4) = sign(constitutive_nonlocal_vs(myInstance),tau(s)) * timeRatio / (1.0_pReal + timeRatio)
      
      if (present(dv_dtau)) then
        dv_dtau(s) = abs(v(s,1)) * constitutive_nonlocal_Qeff0(s,myInstance) / (kB * Temperature * (1.0_pReal + timeRatio)) &
                   * 0.5_pReal * wallFunc * (2.0_pReal - tauRel) / ((1.0_pReal - tauRel) * (abs(tau(s)) - tauThreshold(s)))
      endif
    
    !*** If resolved stress exceeds threshold plus obstacle stress, the probability for thermal activation is 1.
    !*** The tangent is zero, since no dependency of tau.
    
    elseif (tauRel >= 1.0_pReal) then
      v(s,1:4) = sign(constitutive_nonlocal_vs(myInstance), tau(s)) * constitutive_nonlocal_fattack(myInstance) &
               / (constitutive_nonlocal_vs(myInstance) * sqrt(rhoForest(s)) + constitutive_nonlocal_fattack(myInstance))
    endif
  enddo
endif

constitutive_nonlocal_v(1:ns,1:4,g,ip,el) = v
!$OMP FLUSH(constitutive_nonlocal_v)

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i5,x,i2,x,i1)') '<< CONST >> nonlocal_kinetics at el ip g',el,ip,g
    write(6,*)
    write(6,'(a,/,12(x),12(f12.5,x))') '<< CONST >> tau / MPa', tau/1e6_pReal
    write(6,'(a,/,4(12(x),12(f12.5,x),/))') '<< CONST >> v / 1e-3m/s', constitutive_nonlocal_v(:,:,g,ip,el)*1e3
  endif
#endif

endsubroutine



!*********************************************************************
!* calculates plastic velocity gradient and its tangent              *
!*********************************************************************
subroutine constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar99, Tstar_v, Temperature, state, g, ip, el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_Plain3333to99, &
                    math_mul6x6, &
                    math_Mandel6to33
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                        ! current grain number
                                            ip, &                       ! current integration point
                                            el                          ! current element number
real(pReal), intent(in) ::                  Temperature                 ! temperature
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state                       ! microstructural state
real(pReal), dimension(6), intent(in) ::    Tstar_v                     ! 2nd Piola-Kirchhoff stress in Mandel notation

!*** output variables
real(pReal), dimension(3,3), intent(out) :: Lp                          ! plastic velocity gradient
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar99                ! derivative of Lp with respect to Tstar (9x9 matrix)

!*** local variables
integer(pInt)                               myInstance, &               ! current instance of this constitution
                                            myStructure, &              ! current lattice structure
                                            ns, &                       ! short notation for the total number of active slip systems
                                            i, &
                                            j, &
                                            k, &
                                            l, &
                                            t, &                        ! dislocation type
                                            s, &                        ! index of my current slip system
                                            sLattice                    ! index of my current slip system according to lattice order
real(pReal), dimension(3,3,3,3) ::          dLp_dTstar3333              ! derivative of Lp with respect to Tstar (3x3x3x3 matrix)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                      ! single dislocation densities (including used) 
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            gdot                        ! shear rate per dislocation type 
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            tauThreshold, &             ! threshold shear stress
                                            gdotTotal, &                ! shear rate
                                            dv_dtau, &                  ! velocity derivative with respect to the shear stress
                                            dgdotTotal_dtau, &          ! derivative of the shear rate with respect to the shear stress
                                            rhoForest                   ! forest dislocation density


!*** initialize local variables

gdot = 0.0_pReal
Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)


!*** update dislocation velocity

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state(g,ip,el), g, ip, el, dv_dtau)


!*** shortcut to state variables 

forall (t = 1:8) &
  rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (s = 1:ns, t = 5:8, rhoSgl(s,t) * constitutive_nonlocal_v(s,t-4,g,ip,el) < 0.0_pReal) &                                      ! contribution of used rho for changing sign of v
  rhoSgl(s,t-4) = rhoSgl(s,t-4) + abs(rhoSgl(s,t))
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)


!*** Calculation of gdot and its tangent

forall (t = 1:4) &
  gdot(1:ns,t) =  rhoSgl(1:ns,t) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                 * constitutive_nonlocal_v(1:ns,t,g,ip,el)
gdotTotal = sum(gdot,2)
dgdotTotal_dtau = sum(rhoSgl,2) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) * dv_dtau


!*** Calculation of Lp and its tangent

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)  
  Lp = Lp + gdotTotal(s) * lattice_Sslip(1:3,1:3,sLattice,myStructure)
  forall (i=1:3,j=1:3,k=1:3,l=1:3) &
    dLp_dTstar3333(i,j,k,l) = dLp_dTstar3333(i,j,k,l) + dgdotTotal_dtau(s) * lattice_Sslip(i,j, sLattice,myStructure) &
                                                                           * lattice_Sslip(k,l, sLattice,myStructure) 
enddo
dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i5,x,i2,x,i1)') '<< CONST >> nonlocal_LpandItsTangent at el ip g ',el,ip,g
    write(6,*)
    write(6,'(a,/,4(12(x),12(f12.5,x)),/)') '<< CONST >> gdot / 1e-3',gdot*1e3_pReal
    write(6,'(a,/,12(x),12(f12.5,x))') '<< CONST >> gdot total / 1e-3',gdotTotal*1e3_pReal
    write(6,'(a,/,3(12(x),3(f12.7,x),/))') '<< CONST >> Lp',Lp
  endif
#endif

endsubroutine



!*********************************************************************
!* rate of change of microstructure                                  *
!*********************************************************************
subroutine constitutive_nonlocal_dotState(dotState, Tstar_v, Fe, Fp, Temperature, state, aTolState, timestep, orientation, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use numerics, only: numerics_integrationMode
use IO,       only: IO_error
use debug,    only: debug_verbosity, &
                    debug_selectiveDebugger, &
                    debug_g, &
                    debug_i, &
                    debug_e
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv3x3, &
                    math_det3x3, &
                    math_Mandel6to33, &
                    math_QuaternionDisorientation, &
                    math_qRot, &
                    pi, &
                    NaN
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors, &
                    mesh_element, &
                    FE_NipNeighbors, &
                    mesh_ipNeighborhood, &
                    mesh_ipVolume, &
                    mesh_ipArea, &
                    mesh_ipAreaNormal, &
                    mesh_ipCenterOfGravity
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_localConstitution
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_sn, &
                    lattice_st, &
                    lattice_maxNslipFamily, &
                    lattice_NslipSystem 
use FEsolving, only:theInc, &
                    FEsolving_execElem, & 
                    FEsolving_execIP

implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature, &            ! temperature
                                            timestep                  ! substepped crystallite time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            Fe, &                     ! elastic deformation gradient
                                            Fp                        ! plastic deformation gradient
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            orientation               ! crystal lattice orientation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            aTolState                 ! absolute state tolerance

!*** input/output variables
type(p_vec), intent(inout) ::               dotState                  ! evolution of state variables / microstructure
 
!*** output variables
 
!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            n, &                      ! index of my current neighbor
                                            neighboring_el, &         ! element number of my neighbor
                                            neighboring_ip, &         ! integration point of my neighbor
                                            neighboring_n, &          ! neighbor index pointing to me when looking from my neighbor
                                            opposite_n, &             ! index of my opposite neighbor
                                            opposite_ip, &            ! ip of my opposite neighbor
                                            opposite_el, &            ! element index of my opposite neighbor
                                            t, &                      ! type of dislocation
                                            topp, &                   ! type of dislocation with opposite sign to t
                                            s, &                      ! index of my current slip system
                                            sLattice, &               ! index of my current slip system according to lattice order
                                            i
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),10) :: &
                                            rhoDot, &                     ! density evolution
                                            rhoDotRemobilization, &       ! density evolution by remobilization
                                            rhoDotMultiplication, &       ! density evolution by multiplication
                                            rhoDotFlux, &                 ! density evolution by flux
                                            neighboring_rhoDotFlux, &     ! density evolution by flux at neighbor
                                            rhoDotSingle2DipoleGlide, &   ! density evolution by dipole formation (by glide)
                                            rhoDotAthermalAnnihilation, & ! density evolution by athermal annihilation
                                            rhoDotThermalAnnihilation     ! density evolution by thermal annihilation
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl                        ! current single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            fluxdensity, &                ! flux density at central material point
                                            neighboring_fluxdensity, &    ! flux density at neighboring material point
                                            gdot                          ! shear rates
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &                  ! forest dislocation density
                                            tauThreshold, &               ! threshold shear stress
                                            tau, &                        ! current resolved shear stress
                                            invLambda, &                  ! inverse of mean free path for dislocations
                                            vClimb                        ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                     ! current dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                     ! minimum stable dipole distance for edges and screws
                                            dUpper                        ! current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            m                             ! direction of dislocation motion
real(pReal), dimension(3,3) ::              my_F, &                       ! my total deformation gradient
                                            neighboring_F, &              ! total deformation gradient of my neighbor
                                            my_Fe, &                      ! my elastic deformation gradient
                                            neighboring_Fe, &             ! elastic deformation gradient of my neighbor
                                            Favg                          ! average total deformation gradient of me and my neighbor
real(pReal), dimension(6) ::                Tdislocation_v                ! current dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal), dimension(3) ::                normal_neighbor2me, &         ! interface normal pointing from my neighbor to me in neighbor's lattice configuration
                                            normal_neighbor2me_defConf, & ! interface normal pointing from my neighbor to me in shared deformed configuration
                                            normal_me2neighbor, &         ! interface normal pointing from me to my neighbor in my lattice configuration
                                            normal_me2neighbor_defConf    ! interface normal pointing from me to my neighbor in shared deformed configuration
real(pReal)                                 area, &                       ! area of the current interface
                                            transmissivity, &             ! overall transmissivity of dislocation flux to neighboring material point
                                            lineLength, &                 ! dislocation line length leaving the current interface
                                            D, &                          ! self diffusion
                                            correction
logical                                     considerEnteringFlux, &
                                            considerLeavingFlux

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,*)
    write(6,'(a,i5,x,i2,x,i1)') '<< CONST >> nonlocal_dotState at el ip g ',el,ip,g
    write(6,*)
  endif
#endif

select case(mesh_element(2,el))
  case (1,6,7,8,9)
    ! all fine
  case default
    call IO_error(-1,el,ip,g,'element type not supported for nonlocal constitution')
end select

myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

tau = 0.0_pReal
gdot = 0.0_pReal
dLower = 0.0_pReal
dUpper = 0.0_pReal


!*** shortcut to state variables 

forall (t = 1:8) rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(1:ns,c) = state(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
Tdislocation_v = state(g,ip,el)%p(12*ns+1:12*ns+6)


!*** sanity check for timestep

if (timestep <= 0.0_pReal) then                                                                                                     ! if illegal timestep...
  dotState%p = 0.0_pReal                                                                                                            ! ...return without doing anything (-> zero dotState)
  return
endif



!****************************************************************************
!*** Calculate shear rate

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state(g,ip,el), g, ip, el)                                                ! velocities

forall (t = 1:4) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                * constitutive_nonlocal_v(1:ns,t,g,ip,el)
forall (s = 1:ns, t = 1:4, rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) &                                      ! contribution of used rho for changing sign of v
  gdot(s,t) = gdot(s,t) + abs(rhoSgl(s,t+4)) * constitutive_nonlocal_burgersPerSlipSystem(s,myInstance) &
                                             * constitutive_nonlocal_v(s,t,g,ip,el)

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> rho / 1/m^2', rhoSgl, rhoDip
    write(6,'(a,/,4(12(x),12(e12.5,x),/))') '<< CONST >> gdot / 1/s',gdot
  endif
#endif



!****************************************************************************
!*** calculate limits for stable dipole height

do s = 1,ns   ! loop over slip systems
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)  
  tau(s) = math_mul6x6( Tstar_v + Tdislocation_v, lattice_Sslip_v(1:6,sLattice,myStructure) )
enddo

dLower(1:ns,1) = constitutive_nonlocal_dLowerEdgePerSlipSystem(1:ns,myInstance)
dLower(1:ns,2) = constitutive_nonlocal_dLowerScrewPerSlipSystem(1:ns,myInstance)
dUpper(1:ns,2) = min( 1.0_pReal / sqrt( sum(abs(rhoSgl),2)+sum(rhoDip,2) ), &
                      constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                             / ( 8.0_pReal * pi * abs(tau) ) )
dUpper(1:ns,1) = dUpper(1:ns,2) / ( 1.0_pReal - constitutive_nonlocal_nu(myInstance) )


!****************************************************************************
!*** dislocation remobilization (bauschinger effect)

rhoDotRemobilization = 0.0_pReal
if (timestep > 0.0_pReal) then
  do t = 1,4
    do s = 1,ns
      if (rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) then
        rhoDotRemobilization(s,t) = abs(rhoSgl(s,t+4)) / timestep
        rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4))
        rhoDotRemobilization(s,t+4) = - rhoSgl(s,t+4) / timestep
        rhoSgl(s,t+4) = 0.0_pReal
      endif
    enddo
  enddo
endif



!****************************************************************************
!*** calculate dislocation multiplication

rhoDotMultiplication = 0.0_pReal
where (rhoSgl(1:ns,3:4) > 0.0_pReal) &
  rhoDotMultiplication(1:ns,1:2) = spread(0.5_pReal * sum(abs(gdot(1:ns,3:4)),2) * sqrt(rhoForest)  &
                                                    / constitutive_nonlocal_lambda0PerSlipSystem(1:ns,myInstance) &
                                                    / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance), 2, 2)
where (rhoSgl(1:ns,1:2) > 0.0_pReal) &
  rhoDotMultiplication(1:ns,3:4) = spread(0.5_pReal * sum(abs(gdot(1:ns,1:2)),2) * sqrt(rhoForest)  &
                                                    / constitutive_nonlocal_lambda0PerSlipSystem(1:ns,myInstance) &
                                                    / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance), 2, 2)



!****************************************************************************
!*** calculate dislocation fluxes (only for nonlocal constitution)

rhoDotFlux = 0.0_pReal

if (.not. phase_localConstitution(material_phase(g,ip,el))) then                                                                    ! only for nonlocal constitution
  
  !*** take care of the definition of lattice_st = lattice_sd x lattice_sn !!!
  !*** opposite sign to our p vector in the (s,p,n) triplet !!!
  
  m(1:3,1:ns,1) =  lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,2) = -lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,3) = -lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  m(1:3,1:ns,4) =  lattice_st(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,myInstance), myStructure)
  
  my_Fe = Fe(1:3,1:3,g,ip,el)
  my_F = math_mul33x33(my_Fe, Fp(1:3,1:3,g,ip,el))
  
  fluxdensity = rhoSgl(1:ns,1:4) * constitutive_nonlocal_v(1:ns,1:4,g,ip,el)
    
  do n = 1,FE_NipNeighbors(mesh_element(2,el))                                                                                      ! loop through my neighbors
    neighboring_el = mesh_ipNeighborhood(1,n,ip,el)
    neighboring_ip = mesh_ipNeighborhood(2,n,ip,el)
    if (neighboring_el > 0_pInt .and. neighboring_ip > 0_pInt) then                                                                 ! if neighbor exists ...
      do neighboring_n = 1,FE_NipNeighbors(mesh_element(2,neighboring_el))                                                          ! find neighboring index that points from my neighbor to myself
        if (      el == mesh_ipNeighborhood(1,neighboring_n,neighboring_ip,neighboring_el) &
            .and. ip == mesh_ipNeighborhood(2,neighboring_n,neighboring_ip,neighboring_el)) then                                    ! possible candidate
          if (math_mul3x3(mesh_ipAreaNormal(1:3,n,ip,el),&
                          mesh_ipAreaNormal(1:3,neighboring_n,neighboring_ip,neighboring_el)) < 0.0_pReal) then                     ! area normals have opposite orientation (we have to check that because of special case for single element with two ips and periodicity. In this case the neighbor is identical in two different directions.)
            exit
          endif
        endif
      enddo
    endif
  
    opposite_n = n + mod(n,2) - mod(n+1,2)
    opposite_el = mesh_ipNeighborhood(1,opposite_n,ip,el)
    opposite_ip = mesh_ipNeighborhood(2,opposite_n,ip,el)
  
    if (neighboring_el > 0_pInt .and. neighboring_ip > 0_pInt) then                                                                 ! if neighbor exists, average deformation gradient
      neighboring_Fe = Fe(1:3,1:3,g,neighboring_ip,neighboring_el)
      neighboring_F = math_mul33x33(neighboring_Fe, Fp(1:3,1:3,g,neighboring_ip,neighboring_el))
      Favg = 0.5_pReal * (my_F + neighboring_F)
    else                                                                                                                            ! if no neighbor, take my value as average
      Favg = my_F
    endif
    

    !* FLUX FROM ME TO MY NEIGHBOR
    !* This is not considered, if my opposite neighbor has a local constitution. 
    !* Then, we assume, that the opposite(!) neighbor sends an equal amount of dislocations to me.
    !* So the net flux in the direction of my neighbor is equal to zero:
    !*    leaving flux to neighbor == entering flux from opposite neighbor
    !* In case of reduced transmissivity, part of the leaving flux is stored as dead dislocation density.
    !* That means for an interface of zero transmissivity the leaving flux is fully converted to dead dislocations.
    
    considerLeavingFlux = .true.
    if (opposite_el > 0 .and. opposite_ip > 0) then
      if (phase_localConstitution(material_phase(1,opposite_ip,opposite_el))) &
        considerLeavingFlux = .false.
    endif

    if (considerLeavingFlux) then
      normal_me2neighbor_defConf = math_det3x3(Favg) * math_mul33x3(math_inv3x3(transpose(Favg)), mesh_ipAreaNormal(1:3,n,ip,el))   ! calculate the normal of the interface in (average) deformed configuration (pointing from me to my neighbor!!!)
      normal_me2neighbor = math_mul33x3(transpose(my_Fe), normal_me2neighbor_defConf) / math_det3x3(my_Fe)                          ! interface normal in my lattice configuration
      area = mesh_ipArea(n,ip,el) * math_norm3(normal_me2neighbor)
      normal_me2neighbor = normal_me2neighbor / math_norm3(normal_me2neighbor)                                                      ! normalize the surface normal to unit length    
      do s = 1,ns
        do t = 1,4
          c = (t + 1) / 2        
          if (fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_me2neighbor) > 0.0_pReal ) then                                     ! flux from me to my neighbor == leaving flux for me (might also be a pure flux from my mobile density to dead density if interface not at all transmissive)
            lineLength = fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_me2neighbor) * area                                      ! positive line length that wants to leave through this interface
            transmissivity = sum(constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el)**2.0_pReal)                                  ! overall transmissivity from this slip system to my neighbor
            rhoDotFlux(s,t) = rhoDotFlux(s,t) - lineLength / mesh_ipVolume(ip,el)                                                   ! subtract dislocation flux from current mobile type
            rhoDotFlux(s,t+4) = rhoDotFlux(s,t+4) + lineLength / mesh_ipVolume(ip,el) * (1.0_pReal - transmissivity) &
                                                               * sign(1.0_pReal, fluxdensity(s,t))                                  ! dislocation flux that is not able to leave through interface (because of low transmissivity) will remain as immobile single density at the material point
          endif
        enddo
      enddo
    endif    
    
    
    !* FLUX FROM MY NEIGHBOR TO ME
    !* This is only considered, if I have a neighbor of nonlocal constitution that is at least a little bit compatible.
    !* If it's not at all compatible, no flux is arriving, because everything is dammed in front of my neighbor's interface.
    !* The entering flux from my neighbor will be distributed on my slip systems according to the compatibility
    
    considerEnteringFlux = .false.
    if (neighboring_el > 0_pInt .or. neighboring_ip > 0_pInt) then
      if (.not. phase_localConstitution(material_phase(1,neighboring_ip,neighboring_el)) &
          .and. any(constitutive_nonlocal_compatibility(:,:,:,n,ip,el) > 0.0_pReal)) &
      considerEnteringFlux = .true.
    endif
    
    if (considerEnteringFlux) then
      forall (t = 1:4) &
        neighboring_fluxdensity(1:ns,t) = state(g,neighboring_ip,neighboring_el)%p((t-1)*ns+1:t*ns) &
                                        * constitutive_nonlocal_v(1:ns,t,g,neighboring_ip,neighboring_el)
      normal_neighbor2me_defConf = math_det3x3(Favg) &
                  * math_mul33x3(math_inv3x3(transpose(Favg)), mesh_ipAreaNormal(1:3,neighboring_n,neighboring_ip,neighboring_el))  ! calculate the normal of the interface in (average) deformed configuration (now pointing from my neighbor to me!!!)
      normal_neighbor2me = math_mul33x3(transpose(neighboring_Fe), normal_neighbor2me_defConf) / math_det3x3(neighboring_Fe)        ! interface normal in the lattice configuration of my neighbor
      area = mesh_ipArea(neighboring_n,neighboring_ip,neighboring_el) * math_norm3(normal_neighbor2me)
      normal_neighbor2me = normal_neighbor2me / math_norm3(normal_neighbor2me)                                                      ! normalize the surface normal to unit length
      do s = 1,ns
        do t = 1,4
          c = (t + 1) / 2
          topp = t + mod(t,2) - mod(t+1,2)
          if (neighboring_fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_neighbor2me) > 0.0_pReal) then                          ! flux from my neighbor to me == entering flux for me
            lineLength = neighboring_fluxdensity(s,t) * math_mul3x3(m(1:3,s,t), normal_neighbor2me) * area                          ! positive line length that wants to enter through this interface
            where (constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) > 0.0_pReal) &                                             ! positive compatibility...
              rhoDotFlux(1:ns,t) = rhoDotFlux(1:ns,t) + lineLength / mesh_ipVolume(ip,el) &                                         ! ... transferring to equally signed dislocation type
                                                      * constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
            where (constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) < 0.0_pReal) &                                             ! ..negative compatibility...
              rhoDotFlux(1:ns,topp) = rhoDotFlux(1:ns,topp) + lineLength / mesh_ipVolume(ip,el) &                                   ! ... transferring to opposite signed dislocation type
                                                      * constitutive_nonlocal_compatibility(c,1:ns,s,n,ip,el) ** 2.0_pReal
          endif
        enddo
      enddo
    endif
    
  enddo ! neighbor loop  
endif

if (numerics_integrationMode == 1_pInt) then
  constitutive_nonlocal_rhoDotFlux(1:ns,1:10,g,ip,el) = rhoDotFlux(1:ns,1:10)                                                       ! save flux calculation for output (if in central integration mode)
endif



!****************************************************************************
!*** calculate dipole formation and annihilation

!*** formation by glide

do c = 1,2

  rhoDotSingle2DipoleGlide(1:ns,2*c-1) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                    * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                   ! negative mobile --> positive mobile
                                                       + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                 ! positive mobile --> negative mobile
                                                       + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1)))                           ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                  * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &                                     ! negative mobile --> positive mobile
                                                     + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1)) &                                   ! positive mobile --> negative mobile
                                                     + abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)))                               ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+3) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+3) * abs(gdot(1:ns,2*c))                                      ! negative mobile --> positive immobile

  rhoDotSingle2DipoleGlide(1:ns,2*c+4) = -2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                    * rhoSgl(1:ns,2*c+4) * abs(gdot(1:ns,2*c-1))                                    ! positive mobile --> negative immobile

  rhoDotSingle2DipoleGlide(1:ns,c+8) = - rhoDotSingle2DipoleGlide(1:ns,2*c-1) - rhoDotSingle2DipoleGlide(1:ns,2*c) &
                                       + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+3)) + abs(rhoDotSingle2DipoleGlide(1:ns,2*c+4))
enddo


!*** athermal annihilation

rhoDotAthermalAnnihilation = 0.0_pReal

forall (c=1:2) &  
  rhoDotAthermalAnnihilation(1:ns,c+8) = -2.0_pReal * dLower(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
               * (  2.0_pReal * (rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &             ! was single hitting single
                  + 2.0_pReal * (abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) & ! was single hitting immobile single or was immobile single hit by single
                  + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                                 ! single knocks dipole constituent
  
  
!*** thermally activated annihilation of dipoles

rhoDotThermalAnnihilation = 0.0_pReal

D = constitutive_nonlocal_Dsd0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / ( kB * Temperature ) &
          * constitutive_nonlocal_Gmod(myInstance) / ( 2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance)) ) &
          * 2.0_pReal / ( dUpper(1:ns,1) + dLower(1:ns,1) )
          
rhoDotThermalAnnihilation(1:ns,9) = - 4.0_pReal * rhoDip(1:ns,1) * vClimb / (dUpper(1:ns,1) - dLower(1:ns,1))                       ! edge climb
rhoDotThermalAnnihilation(1:ns,10) = 0.0_pReal                                                                                      !!! cross slipping still has to be implemented !!!


!****************************************************************************
!*** assign the rates of dislocation densities to my dotState

rhoDot = 0.0_pReal
forall (t = 1:10) &
  rhoDot(1:ns,t) = rhoDotFlux(1:ns,t) &
              + rhoDotMultiplication(1:ns,t) &
              + rhoDotRemobilization(1:ns,t) &
              + rhoDotSingle2DipoleGlide(1:ns,t) &
              + rhoDotAthermalAnnihilation(1:ns,t) &
              + rhoDotThermalAnnihilation(1:ns,t) 

dotState%p(1:10*ns) = dotState%p(1:10*ns) + reshape(rhoDot,(/10*ns/))

#ifndef _OPENMP
  if (debug_verbosity > 6 .and. ((debug_e == el .and. debug_i == ip .and. debug_g == g) .or. .not. debug_selectiveDebugger)) then
    write(6,'(a,/,8(12(x),12(e12.5,x),/))') '<< CONST >> dislocation remobilization', rhoDotRemobilization(1:ns,1:8) * timestep
    write(6,'(a,/,4(12(x),12(e12.5,x),/))') '<< CONST >> dislocation multiplication', rhoDotMultiplication(1:ns,1:4) * timestep
    write(6,'(a,/,8(12(x),12(e12.5,x),/))') '<< CONST >> dislocation flux', rhoDotFlux(1:ns,1:8) * timestep
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> dipole formation by glide', rhoDotSingle2DipoleGlide * timestep
    write(6,'(a,/,2(12(x),12(e12.5,x),/))') '<< CONST >> athermal dipole annihilation', &
                                            rhoDotAthermalAnnihilation(1:ns,1:2) * timestep
    write(6,'(a,/,2(12(x),12(e12.5,x),/))') '<< CONST >> thermally activated dipole annihilation', &
                                            rhoDotThermalAnnihilation(1:ns,9:10) * timestep
    write(6,'(a,/,10(12(x),12(e12.5,x),/))') '<< CONST >> total density change', rhoDot * timestep
    write(6,'(a,/,10(12(x),12(f12.7,x),/))') '<< CONST >> relative density change', &
                                            rhoDot(1:ns,1:8) * timestep / (abs(rhoSgl)+1.0e-10), &
                                            rhoDot(1:ns,9:10) * timestep / (rhoDip+1.0e-10)
  endif
#endif

endsubroutine



!*********************************************************************
!* COMPATIBILITY UPDATE                                              *
!* Compatibility is defined as normalized product of signed cosine   *
!* of the angle between the slip plane normals and signed cosine of  *
!* the angle between the slip directions. Only the largest values    *
!* that sum up to a total of 1 are considered, all others are set to *
!* zero.                                                             *
!*********************************************************************
subroutine constitutive_nonlocal_updateCompatibility(orientation,i,e)

use prec,     only:   pReal, &
                      pInt
use math, only:       math_QuaternionDisorientation, &
                      math_mul3x3, &
                      math_qRot
use material, only:   material_phase, &
                      phase_constitution, &
                      phase_localConstitution, &
                      phase_constitutionInstance, &
                      homogenization_maxNgrains
use mesh, only:       mesh_element, &
                      mesh_ipNeighborhood, &
                      FE_NipNeighbors, &
                      FE_maxNipNeighbors, &
                      mesh_maxNips, &
                      mesh_NcpElems
use lattice, only:    lattice_sn, &
                      lattice_sd, &
                      lattice_st

implicit none

!* input variables
integer(pInt), intent(in) ::                    i, &                          ! ip index
                                                e                             ! element index
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                                orientation                   ! crystal orientation in quaternions
                                            
!* output variables

!* local variables
integer(pInt)                                   Nneighbors, &                 ! number of neighbors
                                                n, &                          ! neighbor index 
                                                neighboring_e, &              ! element index of my neighbor
                                                neighboring_i, &              ! integration point index of my neighbor
                                                my_phase, &
                                                neighboring_phase, &
                                                my_structure, &               ! lattice structure
                                                my_instance, &                ! instance of constitution
                                                ns, &                         ! number of active slip systems
                                                s1, &                         ! slip system index (me)
                                                s2                            ! slip system index (my neighbor)
real(pReal), dimension(4) ::                    absoluteMisorientation        ! absolute misorientation (without symmetry) between me and my neighbor
real(pReal), dimension(2,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e))),&
                         constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e))),&
                         FE_NipNeighbors(mesh_element(2,e))) :: &  
                                                compatibility                 ! compatibility for current element and ip
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e)))) :: &  
                                                slipNormal, &
                                                slipDirection
real(pReal)                                     compatibilitySum, &
                                                thresholdValue, &
                                                nThresholdValues
logical, dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(1,i,e)))) :: & 
                                                belowThreshold


Nneighbors = FE_NipNeighbors(mesh_element(2,e))
my_phase = material_phase(1,i,e)
my_instance = phase_constitutionInstance(my_phase)
my_structure = constitutive_nonlocal_structure(my_instance)
ns = constitutive_nonlocal_totalNslip(my_instance)
slipNormal(1:3,1:ns) =    lattice_sn(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,my_instance), my_structure)
slipDirection(1:3,1:ns) = lattice_sd(1:3, constitutive_nonlocal_slipSystemLattice(1:ns,my_instance), my_structure)


!*** start out fully compatible

compatibility = 0.0_pReal
forall(s1 = 1:ns) &
  compatibility(1:2,s1,s1,1:Nneighbors) = 1.0_pReal


!*** Loop thrugh neighbors and check whether there is any compatibility.

do n = 1,Nneighbors
  neighboring_e = mesh_ipNeighborhood(1,n,i,e)
  neighboring_i = mesh_ipNeighborhood(2,n,i,e)
  
  
  !* FREE SURFACE
  !* Set surface transmissivity to the value specified in the material.config
  
  if (neighboring_e <= 0 .or. neighboring_i <= 0) then
    forall(s1 = 1:ns) &
      compatibility(1:2,s1,s1,n) = sqrt(constitutive_nonlocal_surfaceTransmissivity(my_instance))
    cycle
  endif
  
  
  !* PHASE BOUNDARY
  !* If we encounter a different nonlocal "cpfem" phase at the neighbor, 
  !* we consider this to be a real "physical" phase boundary, so completely incompatible.
  !* If the neighboring "cpfem" phase has a local constitution, 
  !* we do not consider this to be a phase boundary, so completely compatible.
  
  neighboring_phase = material_phase(1,neighboring_i,neighboring_e)
  if (neighboring_phase /= my_phase) then
    if (.not. phase_localConstitution(neighboring_phase)) then
      compatibility(1:2,1:ns,1:ns,n) = 0.0_pReal
    endif
    cycle
  endif

    
  !* GRAIN BOUNDARY ?
  !* The compatibility value is defined as the product of the slip normal projection and the slip direction projection.
  !* Its sign is always positive for screws, for edges it has the same sign as the slip normal projection. 
  !* Since the sum for each slip system can easily exceed one (which would result in a transmissivity larger than one), 
  !* only values above or equal to a certain threshold value are considered. This threshold value is chosen, such that
  !* the number of compatible slip systems is minimized with the sum of the original compatibility values exceeding one. 
  !* Finally the smallest compatibility value is decreased until the sum is exactly equal to one. 
  !* All values below the threshold are set to zero.
  
  absoluteMisorientation = math_QuaternionDisorientation(orientation(1:4,1,i,e), &
                                                         orientation(1:4,1,neighboring_i,neighboring_e), &
                                                         0_pInt)      ! no symmetry
                                                         
  do s1 = 1,ns    ! my slip systems
    do s2 = 1,ns  ! my neighbor's slip systems
      compatibility(1,s2,s1,n) =     math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2))) &
                               * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
      compatibility(2,s2,s1,n) = abs(math_mul3x3(slipNormal(1:3,s1), math_qRot(absoluteMisorientation, slipNormal(1:3,s2)))) &
                               * abs(math_mul3x3(slipDirection(1:3,s1), math_qRot(absoluteMisorientation, slipDirection(1:3,s2))))
    enddo
    
    compatibilitySum = 0.0_pReal
    belowThreshold = .true.
    do while (compatibilitySum < 1.0_pReal .and. any(belowThreshold(1:ns)))
      thresholdValue = maxval(compatibility(2,1:ns,s1,n), belowThreshold(1:ns))              ! screws always positive
      nThresholdValues = real(count(compatibility(2,1:ns,s1,n) == thresholdValue),pReal)
      where (compatibility(2,1:ns,s1,n) >= thresholdValue) &
        belowThreshold(1:ns) = .false.
      if (compatibilitySum + thresholdValue * nThresholdValues > 1.0_pReal) &
        where (abs(compatibility(1:2,1:ns,s1,n)) == thresholdValue) &
          compatibility(1:2,1:ns,s1,n) = sign((1.0_pReal - compatibilitySum) / nThresholdValues, compatibility(1:2,1:ns,s1,n))
      compatibilitySum = compatibilitySum + nThresholdValues * thresholdValue
    enddo
    where (belowThreshold(1:ns)) compatibility(1,1:ns,s1,n) = 0.0_pReal
    where (belowThreshold(1:ns)) compatibility(2,1:ns,s1,n) = 0.0_pReal
  enddo ! my slip systems cycle
enddo   ! neighbor cycle

constitutive_nonlocal_compatibility(1:2,1:ns,1:ns,1:Nneighbors,i,e) = compatibility

endsubroutine 



!*********************************************************************
!* rate of change of temperature                                     *
!*********************************************************************
pure function constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,state,g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips
use material, only: homogenization_maxNgrains
implicit none

!* input variables
integer(pInt), intent(in) ::              g, &              ! current grain ID
                                          ip, &             ! current integration point
                                          el                ! current element
real(pReal), intent(in) ::                Temperature       ! temperature
real(pReal), dimension(6), intent(in) ::  Tstar_v           ! 2nd Piola-Kirchhoff stress in Mandel notation
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: & 
                                          state             ! microstructural state

!* output variables
real(pReal) constitutive_nonlocal_dotTemperature            ! evolution of Temperature

!* local variables
   
constitutive_nonlocal_dotTemperature = 0.0_pReal

endfunction



!*********************************************************************
!* return array of constitutive results                              *
!*********************************************************************
function constitutive_nonlocal_postResults(Tstar_v, Fe, Temperature, dt, state, dotState, g,ip,el)

use prec,     only: pReal, &
                    pInt, &
                    p_vec
use math,     only: math_norm3, &
                    math_mul6x6, &
                    math_mul3x3, &
                    math_mul33x3, &
                    math_mul33x33, &
                    math_inv3x3, &
                    math_det3x3, &
                    math_Mandel6to33, &
                    pi
use mesh,     only: mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_element
use material, only: homogenization_maxNgrains, &
                    material_phase, &
                    phase_constitutionInstance, &
                    phase_Noutput
use lattice,  only: lattice_Sslip, &
                    lattice_Sslip_v, &
                    lattice_sd, &
                    lattice_st, &
                    lattice_maxNslipFamily, &
                    lattice_NslipSystem
implicit none

!*** input variables
integer(pInt), intent(in) ::                g, &                      ! current grain number
                                            ip, &                     ! current integration point
                                            el                        ! current element number
real(pReal), intent(in) ::                  Temperature, &            ! temperature
                                            dt                        ! time increment
real(pReal), dimension(6), intent(in) ::    Tstar_v                   ! current 2nd Piola-Kirchhoff stress in Mandel notation
real(pReal), dimension(3,3), intent(in) ::  Fe                        ! elastic deformation gradient
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                            state, &                  ! current microstructural state
                                            dotState                  ! evolution rate of microstructural state

!*** output variables
real(pReal), dimension(constitutive_nonlocal_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            constitutive_nonlocal_postResults

!*** local variables
integer(pInt)                               myInstance, &             ! current instance of this constitution
                                            myStructure, &            ! current lattice structure
                                            ns, &                     ! short notation for the total number of active slip systems
                                            c, &                      ! character of dislocation
                                            cs, &                     ! constitutive result index
                                            o, &                      ! index of current output
                                            t, &                      ! type of dislocation
                                            s, &                      ! index of my current slip system
                                            sLattice                  ! index of my current slip system according to lattice order
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),8) :: &
                                            rhoSgl, &                 ! current single dislocation densities (positive/negative screw and edge without dipoles)
                                            rhoDotSgl                 ! evolution rate of single dislocation densities (positive/negative screw and edge without dipoles)
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),4) :: &
                                            gdot                      ! shear rates
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
                                            rhoForest, &              ! forest dislocation density
                                            tauThreshold, &           ! threshold shear stress
                                            tau, &                    ! current resolved shear stress
                                            vClimb                    ! climb velocity of edge dipoles
real(pReal), dimension(constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            rhoDip, &                 ! current dipole dislocation densities (screw and edge dipoles)
                                            rhoDotDip, &              ! evolution rate of dipole dislocation densities (screw and edge dipoles)
                                            dLower, &                 ! minimum stable dipole distance for edges and screws
                                            dUpper                    ! current maximum stable dipole distance for edges and screws
real(pReal), dimension(3,constitutive_nonlocal_totalNslip(phase_constitutionInstance(material_phase(g,ip,el))),2) :: &
                                            m, &                      ! direction of dislocation motion for edge and screw (unit vector)
                                            m_currentconf             ! direction of dislocation motion for edge and screw (unit vector) in current configuration
real(pReal), dimension(6) ::                Tdislocation_v            ! current dislocation stress (resulting from the neighboring excess dislocation densities) as 2nd Piola-Kirchhoff stress
real(pReal)                                 D                         ! self diffusion


myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_nonlocal_structure(myInstance) 
ns = constitutive_nonlocal_totalNslip(myInstance)

cs = 0_pInt
constitutive_nonlocal_postResults = 0.0_pReal


!* short hand notations for state variables

forall (t = 1:8) rhoSgl(1:ns,t) = state(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDip(1:ns,c) = state(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)
rhoForest = state(g,ip,el)%p(10*ns+1:11*ns)
tauThreshold = state(g,ip,el)%p(11*ns+1:12*ns)
Tdislocation_v = state(g,ip,el)%p(12*ns+1:12*ns+6)
forall (t = 1:8) rhoDotSgl(1:ns,t) = dotState(g,ip,el)%p((t-1)*ns+1:t*ns)
forall (c = 1:2) rhoDotDip(1:ns,c) = dotState(g,ip,el)%p((7+c)*ns+1:(8+c)*ns)


!* Calculate shear rate

call constitutive_nonlocal_kinetics(Tstar_v, Temperature, state(g,ip,el), g, ip, el)                                  ! need to calculate dislocation velocity again, because it was overwritten during stiffness calculation

do t = 1,4
  do s = 1,ns
    if (rhoSgl(s,t+4) * constitutive_nonlocal_v(s,t,g,ip,el) < 0.0_pReal) then
      rhoSgl(s,t) = rhoSgl(s,t) + abs(rhoSgl(s,t+4))                                                                  ! remobilization of immobile singles for changing sign of v (bauschinger effect)
      rhoSgl(s,t+4) = 0.0_pReal                                                                                       ! remobilization of immobile singles for changing sign of v (bauschinger effect)
    endif
  enddo
enddo

forall (t = 1:4) &
  gdot(1:ns,t) = rhoSgl(1:ns,t) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                * constitutive_nonlocal_v(1:ns,t,g,ip,el)
  
!* calculate limits for stable dipole height

do s = 1,ns
  sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
  tau(s) = math_mul6x6(Tstar_v + Tdislocation_v, lattice_Sslip_v(1:6,sLattice,myStructure))
enddo

dLower(1:ns,1) = constitutive_nonlocal_dLowerEdgePerSlipSystem(1:ns,myInstance)
dLower(1:ns,2) = constitutive_nonlocal_dLowerScrewPerSlipSystem(1:ns,myInstance)
dUpper(1:ns,2) = min( constitutive_nonlocal_Gmod(myInstance) * constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                                                             / (8.0_pReal * pi * abs(tau)), &
                      1.0_pReal / sqrt(sum(abs(rhoSgl),2)+sum(rhoDip,2)) )
dUpper(1:ns,1) = dUpper(1:ns,2) / (1.0_pReal - constitutive_nonlocal_nu(myInstance))


!*** dislocation motion

m(1:3,1:ns,1) = lattice_sd(1:3,constitutive_nonlocal_slipSystemLattice(1:ns,myInstance),myStructure)
m(1:3,1:ns,2) = -lattice_st(1:3,constitutive_nonlocal_slipSystemLattice(1:ns,myInstance),myStructure)
forall (c = 1:2, s = 1:ns) &
  m_currentconf(1:3,s,c) = math_mul33x3(Fe, m(1:3,s,c))


do o = 1,phase_Noutput(material_phase(g,ip,el))
  select case(constitutive_nonlocal_output(o,myInstance))
    
    case ('rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2) + sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl),2)
      cs = cs + ns
      
    case ('rho_sgl_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,1:4)),2)
      cs = cs + ns
      
    case ('rho_sgl_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,5:8)),2)
      cs = cs + ns
      
    case ('rho_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDip,2)
      cs = cs + ns
      
    case ('rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2) + rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/1,2,5,6/))),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,1:2),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,5:6)),2)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1)
      cs = cs + ns
      
    case ('rho_sgl_edge_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(1:ns,5))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_edge_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(1:ns,6))
      cs = cs + ns
      
    case ('rho_dip_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(1:ns,1)
      cs = cs + ns
      
    case ('rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2) + rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('rho_sgl_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,(/3,4,7,8/))),2)
      cs = cs + ns
            
    case ('rho_sgl_screw_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoSgl(1:ns,3:4),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(abs(rhoSgl(1:ns,7:8)),2)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3)
      cs = cs + ns
      
    case ('rho_sgl_screw_pos_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(1:ns,7))
      cs = cs + ns
      
    case ('rho_sgl_screw_neg')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8))
      cs = cs + ns

    case ('rho_sgl_screw_neg_mobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,4)
      cs = cs + ns

    case ('rho_sgl_screw_neg_immobile')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = abs(rhoSgl(1:ns,8))
      cs = cs + ns

    case ('rho_dip_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoDip(1:ns,2)
      cs = cs + ns
      
    case ('excess_rho')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                    - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6))) &
                                                    + (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                    - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('excess_rho_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,1) + abs(rhoSgl(1:ns,5))) &
                                                    - (rhoSgl(1:ns,2) + abs(rhoSgl(1:ns,6)))
      cs = cs + ns
      
    case ('excess_rho_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = (rhoSgl(1:ns,3) + abs(rhoSgl(1:ns,7))) &
                                                    - (rhoSgl(1:ns,4) + abs(rhoSgl(1:ns,8)))
      cs = cs + ns
      
    case ('rho_forest')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoForest
      cs = cs + ns
    
    case ('delta')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2) + sum(rhoDip,2))
      cs = cs + ns
      
    case ('delta_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(abs(rhoSgl),2))
      cs = cs + ns
      
    case ('delta_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 1.0_pReal / sqrt(sum(rhoDip,2))
      cs = cs + ns
      
    case ('shearrate')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(gdot,2)
      cs = cs + ns
      
    case ('resolvedstress')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v + Tdislocation_v, lattice_Sslip_v(1:6,sLattice,myStructure))
      enddo
      cs = cs + ns
      
    case ('resolvedstress_internal')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tdislocation_v, lattice_Sslip_v(1:6,sLattice,myStructure))
      enddo
      cs = cs + ns
      
    case ('resolvedstress_external')
      do s = 1,ns  
        sLattice = constitutive_nonlocal_slipSystemLattice(s,myInstance)
        constitutive_nonlocal_postResults(cs+s) = math_mul6x6(Tstar_v, lattice_Sslip_v(1:6,sLattice,myStructure))
      enddo
      cs = cs + ns
      
    case ('resistance')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = tauThreshold
      cs = cs + ns
    
    case ('rho_dot')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotSgl,2) + sum(rhoDotDip,2)
      cs = cs + ns
      
    case ('rho_dot_sgl')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotSgl,2)
      cs = cs + ns
      
    case ('rho_dot_dip')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(rhoDotDip,2)
      cs = cs + ns
    
    case ('rho_dot_gen')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(1:ns,3:4)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance)
      cs = cs + ns

    case ('rho_dot_gen_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) =   sum(abs(gdot(1:ns,1:2)),2) * sqrt(rhoForest)  &
                                                      / constitutive_nonlocal_lambda0PerSlipSystem(1:ns,myInstance) &
                                                      / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance)
      cs = cs + ns
      
    case ('rho_dot_sgl2dip')
      do c=1,2                                                                                                                      ! dipole formation by glide
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dUpper(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                      * (  2.0_pReal * (  rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &
                                        + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &                                               ! was single hitting single
                         + 2.0_pReal * (  abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) &
                                        + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))))                                         ! was single hitting immobile/used single
      enddo
      cs = cs + ns
    
    case ('rho_dot_ann_ath')
      do c=1,2
        constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_postResults(cs+1:cs+ns) + &
            2.0_pReal * dLower(1:ns,c) / constitutive_nonlocal_burgersPerSlipSystem(1:ns,myInstance) &
                      * (  2.0_pReal * (  rhoSgl(1:ns,2*c-1) * abs(gdot(1:ns,2*c)) &
                                        + rhoSgl(1:ns,2*c) * abs(gdot(1:ns,2*c-1))) &                                               ! was single hitting single
                         + 2.0_pReal * (  abs(rhoSgl(1:ns,2*c+3)) * abs(gdot(1:ns,2*c)) &
                                        + abs(rhoSgl(1:ns,2*c+4)) * abs(gdot(1:ns,2*c-1))) &                                        ! was single hitting immobile/used single
                         + rhoDip(1:ns,c) * (abs(gdot(1:ns,2*c-1)) + abs(gdot(1:ns,2*c))))                                          ! single knocks dipole constituent
      enddo
      cs = cs + ns
      
    case ('rho_dot_ann_the') 
      D = constitutive_nonlocal_Dsd0(myInstance) * exp(-constitutive_nonlocal_Qsd(myInstance) / (kB * Temperature))

      vClimb =  constitutive_nonlocal_atomicVolume(myInstance) * D / (kB * Temperature) &
          * constitutive_nonlocal_Gmod(myInstance) / (2.0_pReal * pi * (1.0_pReal-constitutive_nonlocal_nu(myInstance))) &
          * 2.0_pReal / (dUpper(1:ns,1) + dLower(1:ns,1))
          
      constitutive_nonlocal_postResults(cs+1:cs+ns) = 4.0_pReal * rhoDip(1:ns,1) * vClimb / (dUpper(1:ns,1) - dLower(1:ns,1))
      ! !!! cross-slip of screws missing !!!
      cs = cs + ns

    case ('rho_dot_flux')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,1:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,5:8,g,ip,el)),2)
      cs = cs + ns
    
    case ('rho_dot_flux_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,1:2,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,5:6,g,ip,el)),2)
      cs = cs + ns
      
    case ('rho_dot_flux_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = sum(constitutive_nonlocal_rhoDotFlux(1:ns,3:4,g,ip,el),2) &
                                                      + sum(abs(constitutive_nonlocal_rhoDotFlux(1:ns,7:8,g,ip,el)),2)
      cs = cs + ns
            
    case ('dislocationvelocity')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = constitutive_nonlocal_v(1:ns,1,g,ip,el)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * constitutive_nonlocal_v(1:ns,1,g,ip,el) &
                                                                     * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * constitutive_nonlocal_v(1:ns,1,g,ip,el) &
                                                                     * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,1) * constitutive_nonlocal_v(1:ns,1,g,ip,el) &
                                                                     * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * constitutive_nonlocal_v(1:ns,2,g,ip,el) &
                                                                       * m_currentconf(1,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * constitutive_nonlocal_v(1:ns,2,g,ip,el) &
                                                                       * m_currentconf(2,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_edge_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,2) * constitutive_nonlocal_v(1:ns,2,g,ip,el) &
                                                                       * m_currentconf(3,1:ns,1)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * constitutive_nonlocal_v(1:ns,3,g,ip,el) &
                                                                     * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * constitutive_nonlocal_v(1:ns,3,g,ip,el) &
                                                                     * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_pos_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = rhoSgl(1:ns,3) * constitutive_nonlocal_v(1:ns,3,g,ip,el) &
                                                                     * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_x')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * constitutive_nonlocal_v(1:ns,4,g,ip,el) &
                                                                       * m_currentconf(1,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_y')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * constitutive_nonlocal_v(1:ns,4,g,ip,el) &
                                                                       * m_currentconf(2,1:ns,2)
      cs = cs + ns
    
    case ('fluxdensity_screw_neg_z')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = - rhoSgl(1:ns,4) * constitutive_nonlocal_v(1:ns,4,g,ip,el) &
                                                                     * m_currentconf(3,1:ns,2)
      cs = cs + ns
    
    case ('d_upper_edge')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(1:ns,1)
      cs = cs + ns
      
    case ('d_upper_screw')
      constitutive_nonlocal_postResults(cs+1:cs+ns) = dUpper(1:ns,2)
      cs = cs + ns
            
 end select
enddo

endfunction

END MODULE
