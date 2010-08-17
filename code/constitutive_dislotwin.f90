!* $Id$
!************************************
!*      Module: CONSTITUTIVE        *
!************************************

MODULE constitutive_dislotwin

!* Include other modules
use prec, only: pReal,pInt
implicit none

!* Lists of states and physical parameters
character(len=*), parameter :: constitutive_dislotwin_label = 'dislotwin'
character(len=18), dimension(2), parameter:: constitutive_dislotwin_listBasicSlipStates = (/'rhoEdge   ', &
                                                                                            'rhoEdgeDip'/)
character(len=18), dimension(1), parameter:: constitutive_dislotwin_listBasicTwinStates = (/'twinFraction'/)
character(len=18), dimension(4), parameter:: constitutive_dislotwin_listDependentSlipStates =(/'invLambdaSlip    ', &
                                                                                               'invLambdaSlipTwin', &
                                                                                               'meanFreePathSlip ', &
                                                                                               'tauSlipThreshold '/)
character(len=18), dimension(4), parameter:: constitutive_dislotwin_listDependentTwinStates =(/'invLambdaTwin   ', &
                                                                                               'meanFreePathTwin', &
                                                                                               'tauTwinThreshold', &
                                                                                               'twinVolume      '/)
real(pReal), parameter :: kB = 1.38e-23_pReal ! Boltzmann constant in J/Kelvin

!* Definition of global variables
integer(pInt), dimension(:), allocatable ::               constitutive_dislotwin_sizeDotState, &                ! number of dotStates
                                                          constitutive_dislotwin_sizeState, &                   ! total number of microstructural state variables
                                                          constitutive_dislotwin_sizePostResults                ! cumulative size of post results
integer(pInt), dimension(:,:), allocatable, target ::     constitutive_dislotwin_sizePostResult                 ! size of each post result output
character(len=64), dimension(:,:), allocatable, target :: constitutive_dislotwin_output                         ! name of each post result output
character(len=32), dimension(:), allocatable ::           constitutive_dislotwin_structureName                  ! name of the lattice structure
integer(pInt), dimension(:), allocatable ::               constitutive_dislotwin_structure, &                   ! number representing the kind of lattice structure
                                                          constitutive_dislotwin_totalNslip, &                  ! total number of active slip systems for each instance
                                                          constitutive_dislotwin_totalNtwin                     ! total number of active twin systems for each instance
integer(pInt), dimension(:,:), allocatable ::             constitutive_dislotwin_Nslip, &                       ! number of active slip systems for each family and instance
                                                          constitutive_dislotwin_Ntwin, &                       ! number of active twin systems for each family and instance
                                                          constitutive_dislotwin_slipFamily, &                  ! lookup table relating active slip system to slip family for each instance
                                                          constitutive_dislotwin_twinFamily, &                  ! lookup table relating active twin system to twin family for each instance
                                                          constitutive_dislotwin_slipSystemLattice, &           ! lookup table relating active slip system index to lattice slip system index for each instance
                                                          constitutive_dislotwin_twinSystemLattice              ! lookup table relating active twin system index to lattice twin system index for each instance
real(pReal), dimension(:), allocatable ::                 constitutive_dislotwin_CoverA, &                      ! c/a ratio for hex type lattice
                                                          constitutive_dislotwin_C11, &                         ! C11 element in elasticity matrix
                                                          constitutive_dislotwin_C12, &                         ! C12 element in elasticity matrix
                                                          constitutive_dislotwin_C13, &                         ! C13 element in elasticity matrix
                                                          constitutive_dislotwin_C33, &                         ! C33 element in elasticity matrix
                                                          constitutive_dislotwin_C44, &                         ! C44 element in elasticity matrix
                                                          constitutive_dislotwin_Gmod, &                        ! shear modulus
                                                          constitutive_dislotwin_CAtomicVolume, &               ! atomic volume in Bugers vector unit
                                                          constitutive_dislotwin_D0, &                          ! prefactor for self-diffusion coefficient
                                                          constitutive_dislotwin_Qsd, &                         ! activation energy for dislocation climb
                                                          constitutive_dislotwin_GrainSize, &                   ! grain size
                                                          constitutive_dislotwin_p, &                           ! p-exponent in glide velocity
                                                          constitutive_dislotwin_q, &                           ! q-exponent in glide velocity
                                                          constitutive_dislotwin_MaxTwinFraction, &             ! maximum allowed total twin volume fraction
                                                          constitutive_dislotwin_r, &                           ! r-exponent in twin nucleation rate
                                                          constitutive_dislotwin_CEdgeDipMinDistance, &         !
                                                          constitutive_dislotwin_Cmfptwin, &                    !
                                                          constitutive_dislotwin_Cthresholdtwin, &              !
                                                          constitutive_dislotwin_SolidSolutionStrength, &       ! Strength due to elements in solid solution
                                                          constitutive_dislotwin_L0, &                          ! Length of twin nuclei in Burgers vectors
                                                          constitutive_dislotwin_relevantRho                    ! dislocation density considered relevant
real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislotwin_Cslip_66                       ! elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:),     allocatable :: constitutive_dislotwin_Ctwin_66                       ! twin elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:,:),   allocatable :: constitutive_dislotwin_Cslip_3333                     ! elasticity matrix for each instance
real(pReal),       dimension(:,:,:,:,:,:), allocatable :: constitutive_dislotwin_Ctwin_3333                     ! twin elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_dislotwin_rhoEdge0, &                    ! initial edge dislocation density per slip system for each family and instance
                                                          constitutive_dislotwin_rhoEdgeDip0, &                 ! initial edge dipole density per slip system for each family and instance
                                                          constitutive_dislotwin_burgersPerSlipFamily, &        ! absolute length of burgers vector [m] for each slip family and instance
                                                          constitutive_dislotwin_burgersPerSlipSystem, &        ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_dislotwin_burgersPerTwinFamily, &        ! absolute length of burgers vector [m] for each twin family and instance
                                                          constitutive_dislotwin_burgersPerTwinSystem, &        ! absolute length of burgers vector [m] for each twin system and instance
                                                          constitutive_dislotwin_QedgePerSlipFamily, &          ! activation energy for glide [J] for each slip family and instance
                                                          constitutive_dislotwin_QedgePerSlipSystem, &          ! activation energy for glide [J] for each slip system and instance
                                                          constitutive_dislotwin_v0PerSlipFamily, &             ! dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_dislotwin_v0PerSlipSystem, &             ! dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_dislotwin_Ndot0PerTwinFamily, &          ! twin nucleation rate [1/m³s] for each twin family and instance
                                                          constitutive_dislotwin_Ndot0PerTwinSystem, &          ! twin nucleation rate [1/m³s] for each twin system and instance
                                                          constitutive_dislotwin_twinsizePerTwinFamily, &       ! twin thickness [m] for each twin family and instance
                                                          constitutive_dislotwin_twinsizePerTwinSystem, &       ! twin thickness [m] for each twin system and instance
                                                          constitutive_dislotwin_CLambdaSlipPerSlipFamily, &    ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_dislotwin_CLambdaSlipPerSlipSystem, &    ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_dislotwin_interactionSlipSlip, &         ! coefficients for slip-slip interaction for each interaction type and instance
                                                          constitutive_dislotwin_interactionSlipTwin, &         ! coefficients for slip-twin interaction for each interaction type and instance
                                                          constitutive_dislotwin_interactionTwinSlip, &         ! coefficients for twin-slip interaction for each interaction type and instance
                                                          constitutive_dislotwin_interactionTwinTwin            ! coefficients for twin-twin interaction for each interaction type and instance
real(pReal),       dimension(:,:,:),       allocatable :: constitutive_dislotwin_interactionMatrixSlipSlip, &   ! interaction matrix of the different slip systems for each instance
                                                          constitutive_dislotwin_interactionMatrixSlipTwin, &   ! interaction matrix of slip systems with twin systems for each instance
                                                          constitutive_dislotwin_interactionMatrixTwinSlip, &   ! interaction matrix of twin systems with slip systems for each instance
                                                          constitutive_dislotwin_interactionMatrixTwinTwin, &   ! interaction matrix of the different twin systems for each instance
                                                          constitutive_dislotwin_forestProjectionEdge           ! matrix of forest projections of edge dislocations for each instance
CONTAINS
!****************************************
!* - constitutive_dislotwin_init
!* - constitutive_dislotwin_stateInit
!* - constitutive_dislotwin_relevantState
!* - constitutive_dislotwin_homogenizedC
!* - constitutive_dislotwin_microstructure
!* - constitutive_dislotwin_LpAndItsTangent
!* - consistutive_dislotwin_dotState
!* - constitutive_dislotwin_dotTemperature
!* - consistutive_dislotwin_postResults
!****************************************

subroutine constitutive_dislotwin_init(file)
!**************************************
!*      Module initialization         *
!**************************************
use prec,    only: pInt,pReal
use math,    only: math_Mandel3333to66,math_Voigt66to3333,math_mul3x3
use IO
use material
use lattice

!* Input variables
integer(pInt), intent(in) :: file
!* Local variables
integer(pInt), parameter :: maxNchunks = 21
integer(pInt), dimension(1+2*maxNchunks) :: positions
integer(pInt) section,maxNinstance,f,i,j,k,l,m,n,o,p,q,r,s,s1,s2,t1,t2,ns,nt,output,mySize,myStructure,maxTotalNslip,maxTotalNtwin
character(len=64) tag
character(len=1024) line

!write(6,*)
!write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_dislotwin_label,' init  -+>>>'
!write(6,*) '$Id$'
!write(6,*)

maxNinstance = count(phase_constitution == constitutive_dislotwin_label)
if (maxNinstance == 0) return

!* Space allocation for global variables
allocate(constitutive_dislotwin_sizeDotState(maxNinstance))
allocate(constitutive_dislotwin_sizeState(maxNinstance))
allocate(constitutive_dislotwin_sizePostResults(maxNinstance))
allocate(constitutive_dislotwin_sizePostResult(maxval(phase_Noutput),maxNinstance))
allocate(constitutive_dislotwin_output(maxval(phase_Noutput),maxNinstance))
constitutive_dislotwin_sizeDotState    = 0_pInt
constitutive_dislotwin_sizeState       = 0_pInt
constitutive_dislotwin_sizePostResults = 0_pInt
constitutive_dislotwin_sizePostResult  = 0_pInt
constitutive_dislotwin_output          = ''

allocate(constitutive_dislotwin_structureName(maxNinstance))
allocate(constitutive_dislotwin_structure(maxNinstance))
allocate(constitutive_dislotwin_Nslip(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_Ntwin(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_dislotwin_slipFamily(lattice_maxNslip,maxNinstance))
allocate(constitutive_dislotwin_twinFamily(lattice_maxNtwin,maxNinstance))
allocate(constitutive_dislotwin_slipSystemLattice(lattice_maxNslip,maxNinstance))
allocate(constitutive_dislotwin_twinSystemLattice(lattice_maxNtwin,maxNinstance))
allocate(constitutive_dislotwin_totalNslip(maxNinstance))
allocate(constitutive_dislotwin_totalNtwin(maxNinstance))
constitutive_dislotwin_structureName     = ''
constitutive_dislotwin_structure         = 0_pInt
constitutive_dislotwin_Nslip             = 0_pInt
constitutive_dislotwin_Ntwin             = 0_pInt
constitutive_dislotwin_slipFamily        = 0_pInt
constitutive_dislotwin_twinFamily        = 0_pInt
constitutive_dislotwin_slipSystemLattice = 0.0_pReal
constitutive_dislotwin_twinSystemLattice = 0.0_pReal
constitutive_dislotwin_totalNslip        = 0_pInt
constitutive_dislotwin_totalNtwin        = 0_pInt
allocate(constitutive_dislotwin_CoverA(maxNinstance))
allocate(constitutive_dislotwin_C11(maxNinstance))
allocate(constitutive_dislotwin_C12(maxNinstance))
allocate(constitutive_dislotwin_C13(maxNinstance))
allocate(constitutive_dislotwin_C33(maxNinstance))
allocate(constitutive_dislotwin_C44(maxNinstance))
allocate(constitutive_dislotwin_Gmod(maxNinstance))
allocate(constitutive_dislotwin_CAtomicVolume(maxNinstance))
allocate(constitutive_dislotwin_D0(maxNinstance))
allocate(constitutive_dislotwin_Qsd(maxNinstance))
allocate(constitutive_dislotwin_GrainSize(maxNinstance))
allocate(constitutive_dislotwin_p(maxNinstance))
allocate(constitutive_dislotwin_q(maxNinstance))
allocate(constitutive_dislotwin_MaxTwinFraction(maxNinstance))
allocate(constitutive_dislotwin_r(maxNinstance))
allocate(constitutive_dislotwin_CEdgeDipMinDistance(maxNinstance))
allocate(constitutive_dislotwin_Cmfptwin(maxNinstance))
allocate(constitutive_dislotwin_Cthresholdtwin(maxNinstance))
allocate(constitutive_dislotwin_SolidSolutionStrength(maxNinstance))
allocate(constitutive_dislotwin_L0(maxNinstance))
allocate(constitutive_dislotwin_relevantRho(maxNinstance))
allocate(constitutive_dislotwin_Cslip_66(6,6,maxNinstance))
allocate(constitutive_dislotwin_Cslip_3333(3,3,3,3,maxNinstance))
constitutive_dislotwin_CoverA               = 0.0_pReal
constitutive_dislotwin_C11                  = 0.0_pReal
constitutive_dislotwin_C12                  = 0.0_pReal
constitutive_dislotwin_C13                  = 0.0_pReal
constitutive_dislotwin_C33                  = 0.0_pReal
constitutive_dislotwin_C44                  = 0.0_pReal
constitutive_dislotwin_Gmod                 = 0.0_pReal
constitutive_dislotwin_CAtomicVolume        = 0.0_pReal
constitutive_dislotwin_D0                   = 0.0_pReal
constitutive_dislotwin_Qsd                  = 0.0_pReal
constitutive_dislotwin_GrainSize            = 0.0_pReal
constitutive_dislotwin_p                    = 0.0_pReal
constitutive_dislotwin_q                    = 0.0_pReal
constitutive_dislotwin_MaxTwinFraction      = 0.0_pReal
constitutive_dislotwin_r                    = 0.0_pReal
constitutive_dislotwin_CEdgeDipMinDistance  = 0.0_pReal
constitutive_dislotwin_Cmfptwin             = 0.0_pReal
constitutive_dislotwin_Cthresholdtwin       = 0.0_pReal
constitutive_dislotwin_SolidSolutionStrength= 0.0_pReal
constitutive_dislotwin_L0                   = 0.0_pReal
constitutive_dislotwin_relevantRho          = 0.0_pReal
constitutive_dislotwin_Cslip_66             = 0.0_pReal
constitutive_dislotwin_Cslip_3333           = 0.0_pReal
allocate(constitutive_dislotwin_rhoEdge0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_rhoEdgeDip0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_burgersPerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_dislotwin_QedgePerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_v0PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_dislotwin_Ndot0PerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_dislotwin_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_dislotwin_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
constitutive_dislotwin_rhoEdge0                 = 0.0_pReal
constitutive_dislotwin_rhoEdgeDip0              = 0.0_pReal
constitutive_dislotwin_burgersPerSlipFamily     = 0.0_pReal
constitutive_dislotwin_burgersPerTwinFamily     = 0.0_pReal
constitutive_dislotwin_QedgePerSlipFamily       = 0.0_pReal
constitutive_dislotwin_v0PerSlipFamily          = 0.0_pReal
constitutive_dislotwin_Ndot0PerTwinFamily       = 0.0_pReal
constitutive_dislotwin_twinsizePerTwinFamily    = 0.0_pReal
constitutive_dislotwin_CLambdaSlipPerSlipFamily = 0.0_pReal
allocate(constitutive_dislotwin_interactionSlipSlip(lattice_maxNinteraction,maxNinstance))
allocate(constitutive_dislotwin_interactionSlipTwin(lattice_maxNinteraction,maxNinstance))
allocate(constitutive_dislotwin_interactionTwinSlip(lattice_maxNinteraction,maxNinstance))
allocate(constitutive_dislotwin_interactionTwinTwin(lattice_maxNinteraction,maxNinstance))
constitutive_dislotwin_interactionSlipSlip = 0.0_pReal
constitutive_dislotwin_interactionSlipTwin = 0.0_pReal
constitutive_dislotwin_interactionTwinSlip = 0.0_pReal
constitutive_dislotwin_interactionTwinTwin = 0.0_pReal

!* Readout data from material.config file
rewind(file)
line = ''
section = 0

do while (IO_lc(IO_getTag(line,'<','>')) /= 'phase')     ! wind forward to <phase>
   read(file,'(a1024)',END=100) line
enddo

do                                                       ! read thru sections of phase part
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1
     output = 0                                           ! reset output counter
   endif
   if (section > 0 .and. phase_constitution(section) == constitutive_dislotwin_label) then  ! one of my sections
     i = phase_constitutionInstance(section)     ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_dislotwin_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
       case ('lattice_structure')
              constitutive_dislotwin_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
       case ('covera_ratio')
              constitutive_dislotwin_CoverA(i) = IO_floatValue(line,positions,2)
       case ('c11')
              constitutive_dislotwin_C11(i) = IO_floatValue(line,positions,2)
       case ('c12')
              constitutive_dislotwin_C12(i) = IO_floatValue(line,positions,2)
       case ('c13')
              constitutive_dislotwin_C13(i) = IO_floatValue(line,positions,2)
       case ('c33')
              constitutive_dislotwin_C33(i) = IO_floatValue(line,positions,2)
       case ('c44')
              constitutive_dislotwin_C44(i) = IO_floatValue(line,positions,2)
       case ('nslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_Nslip(j,i) = IO_intValue(line,positions,1+j)
       case ('ntwin')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_dislotwin_Ntwin(j,i) = IO_intValue(line,positions,1+j)
       case ('rhoedge0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_rhoEdge0(j,i) = IO_floatValue(line,positions,1+j)
       case ('rhoedgedip0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_rhoEdgeDip0(j,i) = IO_floatValue(line,positions,1+j)
       case ('slipburgers')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_burgersPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('twinburgers')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_dislotwin_burgersPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('qedge')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_QedgePerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('v0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_v0PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('ndot0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_dislotwin_Ndot0PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('twinsize')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_dislotwin_twinsizePerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('clambdaslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_dislotwin_CLambdaSlipPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
       case ('grainsize')
              constitutive_dislotwin_GrainSize(i) = IO_floatValue(line,positions,2)
       case ('maxtwinfraction')
              constitutive_dislotwin_MaxTwinFraction(i) = IO_floatValue(line,positions,2)
       case ('pexponent')
              constitutive_dislotwin_p(i) = IO_floatValue(line,positions,2)
       case ('qexponent')
              constitutive_dislotwin_q(i) = IO_floatValue(line,positions,2)
       case ('rexponent')
              constitutive_dislotwin_r(i) = IO_floatValue(line,positions,2)
       case ('d0')
              constitutive_dislotwin_D0(i) = IO_floatValue(line,positions,2)
       case ('qsd')
              constitutive_dislotwin_Qsd(i) = IO_floatValue(line,positions,2)
       case ('relevantrho')
              constitutive_dislotwin_relevantRho(i) = IO_floatValue(line,positions,2)
       case ('cmfptwin')
              constitutive_dislotwin_Cmfptwin(i) = IO_floatValue(line,positions,2)
       case ('cthresholdtwin')
              constitutive_dislotwin_Cthresholdtwin(i) = IO_floatValue(line,positions,2)
       case ('SolidSolutionStrength')
              constitutive_dislotwin_SolidSolutionStrength(i) = IO_floatValue(line,positions,2)
       case ('L0')
              constitutive_dislotwin_L0(i) = IO_floatValue(line,positions,2)
       case ('cedgedipmindistance')
              constitutive_dislotwin_CEdgeDipMinDistance(i) = IO_floatValue(line,positions,2)
       case ('catomicvolume')
              constitutive_dislotwin_CAtomicVolume(i) = IO_floatValue(line,positions,2)
       case ('interactionslipslip')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislotwin_interactionSlipSlip(j,i) = IO_floatValue(line,positions,1+j)
       case ('interactionsliptwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislotwin_interactionSlipTwin(j,i) = IO_floatValue(line,positions,1+j)
       case ('interactiontwinslip')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislotwin_interactionTwinSlip(j,i) = IO_floatValue(line,positions,1+j)
       case ('interactiontwintwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_dislotwin_interactionTwinTwin(j,i) = IO_floatValue(line,positions,1+j)
     end select
   endif
enddo

100 do i = 1,maxNinstance
   constitutive_dislotwin_structure(i) = &
   lattice_initializeStructure(constitutive_dislotwin_structureName(i),constitutive_dislotwin_CoverA(i))
   myStructure = constitutive_dislotwin_structure(i)

   !* Sanity checks
   if (myStructure < 1 .or. myStructure > 3)                                call IO_error(205)
   if (sum(constitutive_dislotwin_Nslip(:,i)) <= 0_pInt)                    call IO_error(225)
   if (sum(constitutive_dislotwin_Ntwin(:,i)) < 0_pInt)                     call IO_error(225) !***
   do f = 1,lattice_maxNslipFamily
     if (constitutive_dislotwin_Nslip(f,i) > 0_pInt) then
       if (constitutive_dislotwin_rhoEdge0(f,i) < 0.0_pReal)                 call IO_error(220)
       if (constitutive_dislotwin_rhoEdgeDip0(f,i) < 0.0_pReal)              call IO_error(220)
       if (constitutive_dislotwin_burgersPerSlipFamily(f,i) <= 0.0_pReal)    call IO_error(221)
       if (constitutive_dislotwin_v0PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(226)
     endif
   enddo
   do f = 1,lattice_maxNtwinFamily
     if (constitutive_dislotwin_Nslip(f,i) > 0_pInt) then
       if (constitutive_dislotwin_burgersPerTwinFamily(f,i) <= 0.0_pReal)    call IO_error(221) !***
       if (constitutive_dislotwin_Ndot0PerTwinFamily(f,i) < 0.0_pReal)       call IO_error(226) !***
     endif
   enddo
!   if (any(constitutive_dislotwin_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 1.0_pReal)) call IO_error(229)
   if (constitutive_dislotwin_CAtomicVolume(i) <= 0.0_pReal)                 call IO_error(230)
   if (constitutive_dislotwin_D0(i) <= 0.0_pReal)                            call IO_error(231)
   if (constitutive_dislotwin_Qsd(i) <= 0.0_pReal)                           call IO_error(232)
   if (constitutive_dislotwin_relevantRho(i) <= 0.0_pReal)                   call IO_error(233)

   !* Determine total number of active slip or twin systems
   constitutive_dislotwin_Nslip(:,i) = min(lattice_NslipSystem(:,myStructure),constitutive_dislotwin_Nslip(:,i))
   constitutive_dislotwin_Ntwin(:,i) = min(lattice_NtwinSystem(:,myStructure),constitutive_dislotwin_Ntwin(:,i))
   constitutive_dislotwin_totalNslip(i) = sum(constitutive_dislotwin_Nslip(:,i))
   constitutive_dislotwin_totalNtwin(i) = sum(constitutive_dislotwin_Ntwin(:,i))

enddo

!* Allocation of variables whose size depends on the total number of active slip systems
maxTotalNslip = maxval(constitutive_dislotwin_totalNslip)
maxTotalNtwin = maxval(constitutive_dislotwin_totalNtwin)

allocate(constitutive_dislotwin_burgersPerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_dislotwin_burgersPerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_dislotwin_QedgePerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_dislotwin_v0PerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_dislotwin_Ndot0PerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_dislotwin_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_dislotwin_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance))
constitutive_dislotwin_burgersPerSlipSystem     = 0.0_pReal
constitutive_dislotwin_burgersPerTwinSystem     = 0.0_pReal
constitutive_dislotwin_QedgePerSlipSystem       = 0.0_pReal
constitutive_dislotwin_v0PerSlipSystem          = 0.0_pReal
constitutive_dislotwin_Ndot0PerTwinSystem       = 0.0_pReal
constitutive_dislotwin_twinsizePerTwinSystem    = 0.0_pReal
constitutive_dislotwin_CLambdaSlipPerSlipSystem = 0.0_pReal

allocate(constitutive_dislotwin_interactionMatrixSlipSlip(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_dislotwin_interactionMatrixSlipTwin(maxTotalNslip,maxTotalNtwin,maxNinstance))
allocate(constitutive_dislotwin_interactionMatrixTwinSlip(maxTotalNtwin,maxTotalNslip,maxNinstance))
allocate(constitutive_dislotwin_interactionMatrixTwinTwin(maxTotalNtwin,maxTotalNtwin,maxNinstance))
allocate(constitutive_dislotwin_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance))
constitutive_dislotwin_interactionMatrixSlipSlip = 0.0_pReal
constitutive_dislotwin_interactionMatrixSlipTwin = 0.0_pReal
constitutive_dislotwin_interactionMatrixTwinSlip = 0.0_pReal
constitutive_dislotwin_interactionMatrixTwinTwin = 0.0_pReal
constitutive_dislotwin_forestProjectionEdge      = 0.0_pReal

allocate(constitutive_dislotwin_Ctwin_66(6,6,maxTotalNtwin,maxNinstance))
allocate(constitutive_dislotwin_Ctwin_3333(3,3,3,3,maxTotalNtwin,maxNinstance))
constitutive_dislotwin_Ctwin_66 = 0.0_pReal
constitutive_dislotwin_Ctwin_3333 = 0.0_pReal

do i = 1,maxNinstance
   myStructure = constitutive_dislotwin_structure(i)

   !* Inverse lookup of my slip system family
   l = 0_pInt
   do f = 1,lattice_maxNslipFamily
      do k = 1,constitutive_dislotwin_Nslip(f,i)
         l = l + 1
         constitutive_dislotwin_slipFamily(l,i) = f
         constitutive_dislotwin_slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1,myStructure)) + k
   enddo; enddo

   !* Inverse lookup of my twin system family
   l = 0_pInt
   do f = 1,lattice_maxNtwinFamily
      do k = 1,constitutive_dislotwin_Ntwin(f,i)
         l = l + 1
         constitutive_dislotwin_twinFamily(l,i) = f
         constitutive_dislotwin_twinSystemLattice(l,i) = sum(lattice_NtwinSystem(1:f-1,myStructure)) + k
   enddo; enddo

   !* Determine size of state array
   ns = constitutive_dislotwin_totalNslip(i)
   nt = constitutive_dislotwin_totalNtwin(i)
   constitutive_dislotwin_sizeDotState(i) = &
   size(constitutive_dislotwin_listBasicSlipStates)*ns+size(constitutive_dislotwin_listBasicTwinStates)*nt
   constitutive_dislotwin_sizeState(i) = &
   constitutive_dislotwin_sizeDotState(i)+ &
   size(constitutive_dislotwin_listDependentSlipStates)*ns+size(constitutive_dislotwin_listDependentTwinStates)*nt

   !* Determine size of postResults array
   do o = 1,maxval(phase_Noutput)
      select case(constitutive_dislotwin_output(o,i))
        case('edge_density', &
             'dipole_density', &
             'shear_rate_slip', &
             'mfp_slip', &
             'resolved_stress_slip', &
             'threshold_stress_slip', &
             'edge_dipole_distance', &
             'stress_exponent' &
             )
           mySize = constitutive_dislotwin_totalNslip(i)
        case('twin_fraction', &
             'shear_rate_twin', &
             'mfp_twin', &
             'resolved_stress_twin', &
             'threshold_stress_twin' &
             )
           mySize = constitutive_dislotwin_totalNtwin(i)
        case default
           mySize = 0_pInt
      end select

       if (mySize > 0_pInt) then  ! any meaningful output found
          constitutive_dislotwin_sizePostResult(o,i) = mySize
          constitutive_dislotwin_sizePostResults(i)  = constitutive_dislotwin_sizePostResults(i) + mySize
       endif
   enddo

   !* Elasticity matrix and shear modulus according to material.config
   select case (myStructure)
   case(1:2) ! cubic(s)
     forall(k=1:3)
       forall(j=1:3) &
         constitutive_dislotwin_Cslip_66(k,j,i)     = constitutive_dislotwin_C12(i)
         constitutive_dislotwin_Cslip_66(k,k,i)     = constitutive_dislotwin_C11(i)
         constitutive_dislotwin_Cslip_66(k+3,k+3,i) = constitutive_dislotwin_C44(i)
     end forall
   case(3:)   ! all hex
     constitutive_dislotwin_Cslip_66(1,1,i) = constitutive_dislotwin_C11(i)
     constitutive_dislotwin_Cslip_66(2,2,i) = constitutive_dislotwin_C11(i)
     constitutive_dislotwin_Cslip_66(3,3,i) = constitutive_dislotwin_C33(i)
     constitutive_dislotwin_Cslip_66(1,2,i) = constitutive_dislotwin_C12(i)
     constitutive_dislotwin_Cslip_66(2,1,i) = constitutive_dislotwin_C12(i)
     constitutive_dislotwin_Cslip_66(1,3,i) = constitutive_dislotwin_C13(i)
     constitutive_dislotwin_Cslip_66(3,1,i) = constitutive_dislotwin_C13(i)
     constitutive_dislotwin_Cslip_66(2,3,i) = constitutive_dislotwin_C13(i)
     constitutive_dislotwin_Cslip_66(3,2,i) = constitutive_dislotwin_C13(i)
     constitutive_dislotwin_Cslip_66(4,4,i) = constitutive_dislotwin_C44(i)
     constitutive_dislotwin_Cslip_66(5,5,i) = constitutive_dislotwin_C44(i)
     constitutive_dislotwin_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_dislotwin_C11(i)-constitutive_dislotwin_C12(i))
   end select
   constitutive_dislotwin_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_dislotwin_Cslip_66(:,:,i)))
   constitutive_dislotwin_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(constitutive_dislotwin_Cslip_66(:,:,i))
   constitutive_dislotwin_Gmod(i) = &
   0.2_pReal*(constitutive_dislotwin_C11(i)-constitutive_dislotwin_C12(i))+0.3_pReal*constitutive_dislotwin_C44(i)

   !* Construction of the twin elasticity matrices
   do j=1,lattice_maxNtwinFamily
      do k=1,constitutive_dislotwin_Ntwin(j,i)
         do l=1,3 ; do m=1,3 ; do n=1,3 ; do o=1,3 ; do p=1,3 ; do q=1,3 ; do r=1,3 ; do s=1,3
           constitutive_dislotwin_Ctwin_3333(l,m,n,o,sum(constitutive_dislotwin_Nslip(1:j-1,i))+k,i) = &
             constitutive_dislotwin_Ctwin_3333(l,m,n,o,sum(constitutive_dislotwin_Nslip(1:j-1,i))+k,i) + &
             constitutive_dislotwin_Cslip_3333(p,q,r,s,i)*&
             lattice_Qtwin(l,p,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(m,q,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(n,r,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(o,s,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)
           enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
         constitutive_dislotwin_Ctwin_66(:,:,k,i) = math_Mandel3333to66(constitutive_dislotwin_Ctwin_3333(:,:,:,:,k,i))
        enddo
   enddo

   !* Burgers vector, dislocation velocity prefactor, mean free path prefactor and minimum dipole distance for each slip system
   do s = 1,constitutive_dislotwin_totalNslip(i)
      f = constitutive_dislotwin_slipFamily(s,i)
      constitutive_dislotwin_burgersPerSlipSystem(s,i)     = constitutive_dislotwin_burgersPerSlipFamily(f,i)
      constitutive_dislotwin_QedgePerSlipSystem(s,i)       = constitutive_dislotwin_QedgePerSlipFamily(f,i)
      constitutive_dislotwin_v0PerSlipSystem(s,i)          = constitutive_dislotwin_v0PerSlipFamily(f,i)
      constitutive_dislotwin_CLambdaSlipPerSlipSystem(s,i) = constitutive_dislotwin_CLambdaSlipPerSlipFamily(f,i)
   enddo

   !* Burgers vector, nucleation rate prefactor and twin size for each twin system
   do s = 1,constitutive_dislotwin_totalNtwin(i)
      f = constitutive_dislotwin_twinFamily(s,i)
      constitutive_dislotwin_burgersPerTwinSystem(s,i)  = constitutive_dislotwin_burgersPerTwinFamily(f,i)
      constitutive_dislotwin_Ndot0PerTwinSystem(s,i)    = constitutive_dislotwin_Ndot0PerTwinFamily(f,i)
      constitutive_dislotwin_twinsizePerTwinSystem(s,i) = constitutive_dislotwin_twinsizePerTwinFamily(f,i)
   enddo

   !* Construction of interaction matrices
   do s1 = 1,constitutive_dislotwin_totalNslip(i)
      do s2 = 1,constitutive_dislotwin_totalNslip(i)
         constitutive_dislotwin_interactionMatrixSlipSlip(s1,s2,i) = &
         constitutive_dislotwin_interactionSlipSlip(lattice_interactionSlipSlip(constitutive_dislotwin_slipSystemLattice(s1,i), &
                                                                                constitutive_dislotwin_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do s1 = 1,constitutive_dislotwin_totalNslip(i)
      do t2 = 1,constitutive_dislotwin_totalNtwin(i)
         constitutive_dislotwin_interactionMatrixSlipTwin(s1,t2,i) = &
         constitutive_dislotwin_interactionSlipTwin(lattice_interactionSlipTwin(constitutive_dislotwin_slipSystemLattice(s1,i), &
                                                                                constitutive_dislotwin_twinSystemLattice(t2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do t1 = 1,constitutive_dislotwin_totalNtwin(i)
      do s2 = 1,constitutive_dislotwin_totalNslip(i)
         constitutive_dislotwin_interactionMatrixTwinSlip(t1,s2,i) = &
         constitutive_dislotwin_interactionTwinSlip(lattice_interactionTwinSlip(constitutive_dislotwin_twinSystemLattice(t1,i), &
                                                                                constitutive_dislotwin_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do t1 = 1,constitutive_dislotwin_totalNtwin(i)
      do t2 = 1,constitutive_dislotwin_totalNtwin(i)
         constitutive_dislotwin_interactionMatrixTwinTwin(t1,t2,i) = &
         constitutive_dislotwin_interactionTwinTwin(lattice_interactionTwinTwin(constitutive_dislotwin_twinSystemLattice(t1,i), &
                                                                                constitutive_dislotwin_twinSystemLattice(t2,i), &
                                                                                myStructure),i)
   enddo; enddo

   !* Calculation of forest projections for edge dislocations
   do s1 = 1,constitutive_dislotwin_totalNslip(i)
      do s2 = 1,constitutive_dislotwin_totalNslip(i)
         constitutive_dislotwin_forestProjectionEdge(s1,s2,i) = &
         abs(math_mul3x3(lattice_sn(:,constitutive_dislotwin_slipSystemLattice(s1,i),myStructure), &
                         lattice_st(:,constitutive_dislotwin_slipSystemLattice(s2,i),myStructure)))
   enddo; enddo

enddo

return
end subroutine


function constitutive_dislotwin_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
use prec,    only: pReal,pInt
use math,    only: pi
use lattice, only: lattice_maxNslipFamily,lattice_maxNtwinFamily
implicit none

!* Input-Output variables
integer(pInt) :: myInstance
real(pReal), dimension(constitutive_dislotwin_sizeState(myInstance))  :: constitutive_dislotwin_stateInit
!* Local variables
integer(pInt) s0,s1,s,t,f,ns,nt
real(pReal), dimension(constitutive_dislotwin_totalNslip(myInstance)) :: rhoEdge0, &
                                                                         rhoEdgeDip0, &
                                                                         invLambdaSlip0, &
                                                                         MeanFreePathSlip0, &
                                                                         tauSlipThreshold0
real(pReal), dimension(constitutive_dislotwin_totalNtwin(myInstance)) :: MeanFreePathTwin0,TwinVolume0

ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)
constitutive_dislotwin_stateInit = 0.0_pReal

!* Initialize basic slip state variables
s1 = 0_pInt
do f = 1,lattice_maxNslipFamily
   s0 = s1 + 1_pInt
   s1 = s0 + constitutive_dislotwin_Nslip(f,myInstance) - 1_pInt
   do s = s0,s1
      rhoEdge0(s)    = constitutive_dislotwin_rhoEdge0(f,myInstance)
      rhoEdgeDip0(s) = constitutive_dislotwin_rhoEdgeDip0(f,myInstance)
   enddo
enddo
constitutive_dislotwin_stateInit(1:ns)      = rhoEdge0
constitutive_dislotwin_stateInit(ns+1:2*ns) = rhoEdgeDip0

!* Initialize dependent slip microstructural variables
forall (s = 1:ns) &
invLambdaSlip0(s) = sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislotwin_forestProjectionEdge(1:ns,s,myInstance)))/ &
constitutive_dislotwin_CLambdaSlipPerSlipSystem(s,myInstance)
constitutive_dislotwin_stateInit(2*ns+nt+1:3*ns+nt) = invLambdaSlip0

forall (s = 1:ns) &
MeanFreePathSlip0(s) = &
constitutive_dislotwin_GrainSize(myInstance)/(1.0_pReal+invLambdaSlip0(s)*constitutive_dislotwin_GrainSize(myInstance))
constitutive_dislotwin_stateInit(4*ns+2*nt+1:5*ns+2*nt) = MeanFreePathSlip0

forall (s = 1:ns) &
tauSlipThreshold0(s) = constitutive_dislotwin_SolidSolutionStrength(myInstance)+ &
constitutive_dislotwin_Gmod(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(s,myInstance)* &
sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislotwin_interactionMatrixSlipSlip(1:ns,s,myInstance)))
constitutive_dislotwin_stateInit(5*ns+3*nt+1:6*ns+3*nt) = tauSlipThreshold0

!* Initialize dependent twin microstructural variables
forall (t = 1:nt) &
MeanFreePathTwin0(t) = constitutive_dislotwin_GrainSize(myInstance)
constitutive_dislotwin_stateInit(5*ns+2*nt+1:5*ns+3*nt) = MeanFreePathTwin0

forall (t = 1:nt) &
TwinVolume0(t) = &
(pi/6.0_pReal)*constitutive_dislotwin_twinsizePerTwinSystem(t,myInstance)*MeanFreePathTwin0(t)**(2.0_pReal)
constitutive_dislotwin_stateInit(6*ns+4*nt+1:6*ns+5*nt) = TwinVolume0

!write(6,*) '#STATEINIT#'
!write(6,*)
!write(6,'(a,/,4(3(f30.20,x)/))') 'RhoEdge',rhoEdge0
!write(6,'(a,/,4(3(f30.20,x)/))') 'RhoEdgedip',rhoEdgeDip0
!write(6,'(a,/,4(3(f30.20,x)/))') 'invLambdaSlip',invLambdaSlip0
!write(6,'(a,/,4(3(f30.20,x)/))') 'MeanFreePathSlip',MeanFreePathSlip0
!write(6,'(a,/,4(3(f30.20,x)/))') 'tauSlipThreshold', tauSlipThreshold0
!write(6,'(a,/,4(3(f30.20,x)/))') 'MeanFreePathTwin', MeanFreePathTwin0
!write(6,'(a,/,4(3(f30.20,x)/))') 'TwinVolume', TwinVolume0

return
end function


pure function constitutive_dislotwin_relevantState(myInstance)
!*********************************************************************
!* relevant microstructural state                                    *
!*********************************************************************
use prec,     only: pReal, pInt
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: myInstance
real(pReal), dimension(constitutive_dislotwin_sizeState(myInstance)) :: constitutive_dislotwin_relevantState

constitutive_dislotwin_relevantState = constitutive_dislotwin_relevantRho(myInstance)

return
endfunction


pure function constitutive_dislotwin_homogenizedC(state,g,ip,el)
!*********************************************************************
!* calculates homogenized elacticity matrix                          *
!*  - state           : microstructure quantities                    *
!*  - g               : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
real(pReal), dimension(6,6) :: constitutive_dislotwin_homogenizedC
!* Local variables
integer(pInt) myInstance,ns,nt,i
real(pReal) sumf

!* Shortened notation
myInstance = phase_constitutionInstance(material_phase(g,ip,el))
ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)

!* Total twin volume fraction
sumf = sum(state(g,ip,el)%p((2*ns+1):(2*ns+nt))) ! safe for nt == 0

!* Homogenized elasticity matrix
constitutive_dislotwin_homogenizedC = (1.0_pReal-sumf)*constitutive_dislotwin_Cslip_66(:,:,myInstance)
do i=1,nt
   constitutive_dislotwin_homogenizedC = &
   constitutive_dislotwin_homogenizedC + state(g,ip,el)%p(2*ns+i)*constitutive_dislotwin_Ctwin_66(:,:,i,myInstance)
enddo

return
end function


subroutine constitutive_dislotwin_microstructure(Temperature,state,g,ip,el)
!*********************************************************************
!* calculates quantities characterizing the microstructure           *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec
use math,     only: pi
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
use lattice,  only: lattice_interactionSlipTwin,lattice_interactionTwinTwin
!use debug,    only: debugger
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: Temperature
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: state
!* Local variables
integer(pInt) myInstance,myStructure,ns,nt,s,t
real(pReal) sumf,sfe
real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: fOverStacksize

!* Shortened notation
myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_dislotwin_structure(myInstance)
ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)
!* State: 1           :  ns         rho_edge
!* State: ns+1        :  2*ns       rho_dipole
!* State: 2*ns+1      :  2*ns+nt    f
!* State: 2*ns+nt+1   :  3*ns+nt    1/lambda_slip
!* State: 3*ns+nt+1   :  4*ns+nt    1/lambda_sliptwin
!* State: 4*ns+nt+1   :  4*ns+2*nt  1/lambda_twin
!* State: 4*ns+2*nt+1 :  5*ns+2*nt  mfp_slip
!* State: 5*ns+2*nt+1 :  5*ns+3*nt  mfp_twin
!* State: 5*ns+3*nt+1 :  6*ns+3*nt  threshold_stress_slip
!* State: 6*ns+3*nt+1 :  6*ns+4*nt  threshold_stress_twin
!* State: 6*ns+4*nt+1 :  6*ns+5*nt  twin volume

!* Total twin volume fraction
sumf = sum(state(g,ip,el)%p((2*ns+1):(2*ns+nt))) ! safe for nt == 0

!* Stacking fault energy
sfe = 0.0002_pReal*Temperature-0.0396_pReal

!* rescaled twin volume fraction for topology
forall (t = 1:nt) &
  fOverStacksize(t) = &
    state(g,ip,el)%p(2*ns+t)/constitutive_dislotwin_twinsizePerTwinSystem(t,myInstance)

!* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
forall (s = 1:ns) &
  state(g,ip,el)%p(2*ns+nt+s) = &
    sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)),&
                     constitutive_dislotwin_forestProjectionEdge(1:ns,s,myInstance)))/ &
    constitutive_dislotwin_CLambdaSlipPerSlipSystem(s,myInstance)

!* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
!$OMP CRITICAL (evilmatmul)
state(g,ip,el)%p((3*ns+nt+1):(4*ns+nt)) = 0.0_pReal
if (nt > 0_pInt) &
  state(g,ip,el)%p((3*ns+nt+1):(4*ns+nt)) = &
    matmul(constitutive_dislotwin_interactionMatrixSlipTwin(1:ns,1:nt,myInstance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
!$OMP END CRITICAL (evilmatmul)

!* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
!$OMP CRITICAL (evilmatmul)
if (nt > 0_pInt) &
  state(g,ip,el)%p((4*ns+nt+1):(4*ns+2*nt)) = &
    matmul(constitutive_dislotwin_interactionMatrixTwinTwin(1:nt,1:nt,myInstance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
!$OMP END CRITICAL (evilmatmul)

!* mean free path between 2 obstacles seen by a moving dislocation
do s = 1,ns
   if (nt > 0_pInt) then
      state(g,ip,el)%p(4*ns+2*nt+s) = &
        constitutive_dislotwin_GrainSize(myInstance)/(1.0_pReal+constitutive_dislotwin_GrainSize(myInstance)*&
        (state(g,ip,el)%p(2*ns+nt+s)+state(g,ip,el)%p(3*ns+nt+s)))
   else
      state(g,ip,el)%p(4*ns+s) = &
        constitutive_dislotwin_GrainSize(myInstance)/&
        (1.0_pReal+constitutive_dislotwin_GrainSize(myInstance)*(state(g,ip,el)%p(2*ns+s)))
   endif
enddo

!* mean free path between 2 obstacles seen by a growing twin
forall (t = 1:nt) &
  state(g,ip,el)%p(5*ns+2*nt+t) = &
    (constitutive_dislotwin_Cmfptwin(myInstance)*constitutive_dislotwin_GrainSize(myInstance))/&
    (1.0_pReal+constitutive_dislotwin_GrainSize(myInstance)*state(g,ip,el)%p(4*ns+nt+t))

!* threshold stress for dislocation motion
forall (s = 1:ns) &
  state(g,ip,el)%p(5*ns+3*nt+s) = constitutive_dislotwin_SolidSolutionStrength(myInstance)+ &
    constitutive_dislotwin_Gmod(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(s,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)),&
                     constitutive_dislotwin_interactionMatrixSlipSlip(1:ns,s,myInstance)))

!* threshold stress for growing twin
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+3*nt+t) = &
    constitutive_dislotwin_Cthresholdtwin(myInstance)*&
    (sfe/(3.0_pReal*constitutive_dislotwin_burgersPerTwinSystem(t,myInstance))+&
    3.0_pReal*constitutive_dislotwin_burgersPerTwinSystem(t,myInstance)*constitutive_dislotwin_Gmod(myInstance)/&
    (constitutive_dislotwin_L0(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(t,myInstance)))

!* final twin volume after growth
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+4*nt+t) = &
    (pi/6.0_pReal)*constitutive_dislotwin_twinsizePerTwinSystem(t,myInstance)*state(g,ip,el)%p(5*ns+2*nt+t)**(2.0_pReal)

!if ((ip==1).and.(el==1)) then
!   write(6,*) '#MICROSTRUCTURE#'
! write(6,*)
! write(6,'(a,/,4(3(f10.4,x)/))') 'rhoEdge',state(g,ip,el)%p(1:ns)/1e9
! write(6,'(a,/,4(3(f10.4,x)/))') 'rhoEdgeDip',state(g,ip,el)%p(ns+1:2*ns)/1e9
! write(6,'(a,/,4(3(f10.4,x)/))') 'Fraction',state(g,ip,el)%p(2*ns+1:2*ns+nt)
!endif


return
end subroutine


subroutine constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,g,ip,el)
!*********************************************************************
!* calculates plastic velocity gradient and its tangent              *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-rank tensor)           *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec
use math,     only: math_Plain3333to99
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance
use lattice,  only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: state
real(pReal), dimension(3,3), intent(out) :: Lp
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar
!* Local variables
integer(pInt) myInstance,myStructure,ns,nt,f,i,j,k,l,m,n,index_myFamily
real(pReal) sumf,StressRatio_p,StressRatio_pminus1,StressRatio_r,BoltzmannRatio,DotGamma0
real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
real(pReal), dimension(constitutive_dislotwin_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip
real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_dislotwin_structure(myInstance)
ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)

!* Total twin volume fraction
sumf = sum(state(g,ip,el)%p((2*ns+1):(2*ns+nt))) ! safe for nt == 0

Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal
dLp_dTstar = 0.0_pReal

!* Dislocation glide part
gdot_slip = 0.0_pReal
dgdot_dtauslip = 0.0_pReal
j = 0_pInt
do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
   do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))

      !* Stress ratios
      StressRatio_p = (abs(tau_slip(j))/state(g,ip,el)%p(5*ns+3*nt+j))**constitutive_dislotwin_p(myInstance)
      StressRatio_pminus1 = (abs(tau_slip(j))/state(g,ip,el)%p(5*ns+3*nt+j))**(constitutive_dislotwin_p(myInstance)-1.0_pReal)
      !* Boltzmann ratio
      BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(f,myInstance)/(kB*Temperature)
      !* Initial shear rates
      DotGamma0 = &
        state(g,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)*&
        constitutive_dislotwin_v0PerSlipSystem(f,myInstance)

      !* Shear rates due to slip
      gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1-StressRatio_p)**constitutive_dislotwin_q(myInstance))*&
                     sign(1.0_pReal,tau_slip(j))

      !* Derivatives of shear rates
      dgdot_dtauslip(j) = &
        ((abs(gdot_slip(j))*BoltzmannRatio*&
        constitutive_dislotwin_p(myInstance)*constitutive_dislotwin_q(myInstance))/state(g,ip,el)%p(5*ns+3*nt+j))*&
        StressRatio_pminus1*(1-StressRatio_p)**(constitutive_dislotwin_q(myInstance)-1.0_pReal)

      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + (1.0_pReal - sumf)*gdot_slip(j)*lattice_Sslip(:,:,index_myFamily+i,myStructure)

      !* Calculation of the tangent of Lp
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*&
                                  lattice_Sslip(k,l,index_myFamily+i,myStructure)*&
                                  lattice_Sslip(m,n,index_myFamily+i,myStructure)
   enddo
enddo

!* Mechanical twinning part
gdot_twin = 0.0_pReal
dgdot_dtautwin = 0.0_pReal
j = 0_pInt
do f = 1,lattice_maxNtwinFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1,myStructure)) ! at which index starts my family
   do i = 1,constitutive_dislotwin_Ntwin(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))

      !* Stress ratios
      StressRatio_r = (state(g,ip,el)%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_dislotwin_r(myInstance)

      !* Shear rates and their derivatives due to twin
      if ( tau_twin(j) > 0.0_pReal ) then
        gdot_twin(j) = &
          (constitutive_dislotwin_MaxTwinFraction(myInstance)-sumf)*lattice_shearTwin(index_myFamily+i,myStructure)*&
          state(g,ip,el)%p(6*ns+4*nt+j)*constitutive_dislotwin_Ndot0PerTwinSystem(f,myInstance)*exp(-StressRatio_r)
        dgdot_dtautwin(j) = ((gdot_twin(j)*constitutive_dislotwin_r(myInstance))/tau_twin(j))*StressRatio_r
      endif

      !* Plastic velocity gradient for mechanical twinning
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,myStructure)

      !* Calculation of the tangent of Lp
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,myStructure)*&
                                  lattice_Stwin(m,n,index_myFamily+i,myStructure)
   enddo
enddo

dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)

!if ((ip==1).and.(el==1)) then
!   write(6,*) '#LP/TANGENT#'
!   write(6,*)
!   write(6,*) 'Tstar_v', Tstar_v
!   write(6,*) 'tau_slip', tau_slip
!   write(6,'(a10,/,4(3(e20.8,x),/))') 'state',state(1,1,1)%p
!   write(6,'(a,/,3(3(f10.4,x)/))') 'Lp',Lp
!   write(6,'(a,/,9(9(f10.4,x)/))') 'dLp_dTstar',dLp_dTstar
!endif

return
end subroutine


function constitutive_dislotwin_dotState(Tstar_v,Temperature,state,g,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - state           : microstructure quantities                    *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec

use math,     only: pi
use mesh,     only: mesh_NcpElems, mesh_maxNips
use material, only: homogenization_maxNgrains, material_phase, phase_constitutionInstance
use lattice,  only: lattice_Sslip,lattice_Sslip_v, &
                    lattice_Stwin,lattice_Stwin_v, &
                    lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem, &
                    lattice_shearTwin
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
real(pReal), dimension(constitutive_dislotwin_sizeDotState(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
constitutive_dislotwin_dotState
!* Local variables
integer(pInt) MyInstance,MyStructure,ns,nt,f,i,j,index_myFamily
real(pReal) sumf,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,&
            EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r
real(pReal), dimension(constitutive_dislotwin_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
gdot_slip,tau_slip,DotRhoMultiplication,EdgeDipDistance,DotRhoEdgeEdgeAnnihilation,DotRhoEdgeDipAnnihilation,&

ClimbVelocity,DotRhoEdgeDipClimb,DotRhoDipFormation
real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
             tau_twin

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
MyStructure = constitutive_dislotwin_structure(myInstance)
ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)

!* Total twin volume fraction
sumf = sum(state(g,ip,el)%p((2*ns+1):(2*ns+nt))) ! safe for nt == 0

constitutive_dislotwin_dotState = 0.0_pReal

!* Dislocation density evolution
gdot_slip = 0.0_pReal
j = 0_pInt
do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,MyStructure)) ! at which index starts my family
   do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt


      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))
      !* Stress ratios
      StressRatio_p = (abs(tau_slip(j))/state(g,ip,el)%p(5*ns+3*nt+j))**constitutive_dislotwin_p(myInstance)
      StressRatio_pminus1 = (abs(tau_slip(j))/state(g,ip,el)%p(5*ns+3*nt+j))**(constitutive_dislotwin_p(myInstance)-1.0_pReal)
      !* Boltzmann ratio
      BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(f,myInstance)/(kB*Temperature)
      !* Initial shear rates
      DotGamma0 = &
        state(g,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)*&
        constitutive_dislotwin_v0PerSlipSystem(f,myInstance)

      !* Shear rates due to slip
      gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1-StressRatio_p)**constitutive_dislotwin_q(myInstance))*&
                     sign(1.0_pReal,tau_slip(j))

      !* Multiplication
      DotRhoMultiplication(j) = abs(gdot_slip(j))/&
                                (constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)*state(g,ip,el)%p(4*ns+2*nt+j))

      !* Dipole formation
      EdgeDipMinDistance = &
        constitutive_dislotwin_CEdgeDipMinDistance(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)
      if (tau_slip(j) == 0.0_pReal) then
        DotRhoDipFormation(j) = 0.0_pReal
      else
        EdgeDipDistance(j) = &
          (3.0_pReal*constitutive_dislotwin_Gmod(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance))/&
          (16.0_pReal*pi*abs(tau_slip(j)))
      if (EdgeDipDistance(j)>state(g,ip,el)%p(4*ns+2*nt+j)) EdgeDipDistance(j)=state(g,ip,el)%p(4*ns+2*nt+j)
      if (EdgeDipDistance(j)<EdgeDipMinDistance) EdgeDipDistance(j)=EdgeDipMinDistance
        DotRhoDipFormation(j) = &
          ((2.0_pReal*EdgeDipDistance(j))/constitutive_dislotwin_burgersPerSlipSystem(f,myInstance))*&
          state(g,ip,el)%p(j)*abs(gdot_slip(j))
      endif

      !* Spontaneous annihilation of 2 single edge dislocations
      DotRhoEdgeEdgeAnnihilation(j) = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislotwin_burgersPerSlipSystem(f,myInstance))*&
        state(g,ip,el)%p(j)*abs(gdot_slip(j))

      !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
      DotRhoEdgeDipAnnihilation(j) = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislotwin_burgersPerSlipSystem(f,myInstance))*&
        state(g,ip,el)%p(ns+j)*abs(gdot_slip(j))

      !* Dislocation dipole climb
      AtomicVolume = &
        constitutive_dislotwin_CAtomicVolume(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)**(3.0_pReal)
      VacancyDiffusion = &
        constitutive_dislotwin_D0(myInstance)*exp(-constitutive_dislotwin_Qsd(myInstance)/(kB*Temperature))
      if (tau_slip(j) == 0.0_pReal) then
        DotRhoEdgeDipClimb(j) = 0.0_pReal
      else
        ClimbVelocity(j) = &
          ((3.0_pReal*constitutive_dislotwin_Gmod(myInstance)*VacancyDiffusion*AtomicVolume)/(2.0_pReal*pi*kB*Temperature))*&
          (1/(EdgeDipDistance(j)+EdgeDipMinDistance))
        DotRhoEdgeDipClimb(j) = &
          (4.0_pReal*ClimbVelocity(j)*state(g,ip,el)%p(ns+j))/(EdgeDipDistance(j)-EdgeDipMinDistance)
      endif

      !* Edge dislocation density rate of change
      constitutive_dislotwin_dotState(j) = &
        DotRhoMultiplication(j)-DotRhoDipFormation(j)-DotRhoEdgeEdgeAnnihilation(j)

      !* Edge dislocation dipole density rate of change
      constitutive_dislotwin_dotState(ns+j) = &
        DotRhoDipFormation(j)-DotRhoEdgeDipAnnihilation(j)-DotRhoEdgeDipClimb(j)

   enddo
enddo

!* Twin volume fraction evolution
j = 0_pInt
do f = 1,lattice_maxNtwinFamily                                 ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1,MyStructure)) ! at which index starts my family
   do i = 1,constitutive_dislotwin_Ntwin(f,myInstance)          ! process each (active) twin system in family
      j = j+1_pInt

      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))
      !* Stress ratios
      StressRatio_r = (state(g,ip,el)%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_dislotwin_r(myInstance)

      !* Shear rates and their derivatives due to twin
      if ( tau_twin(j) > 0.0_pReal ) then
        constitutive_dislotwin_dotState(2*ns+j) = &
          (constitutive_dislotwin_MaxTwinFraction(myInstance)-sumf)*&
          state(g,ip,el)%p(6*ns+4*nt+j)*constitutive_dislotwin_Ndot0PerTwinSystem(f,myInstance)*exp(-StressRatio_r)
      endif

   enddo
enddo

!write(6,*) '#DOTSTATE#'
!write(6,*)
!write(6,'(a,/,4(3(f30.20,x)/))') 'tau slip',tau_slip
!write(6,'(a,/,4(3(f30.20,x)/))') 'gamma slip',gdot_slip
!write(6,'(a,/,4(3(f30.20,x)/))') 'RhoEdge',state(g,ip,el)%p(1:ns)
!write(6,'(a,/,4(3(f30.20,x)/))') 'Threshold Slip', state(g,ip,el)%p(5*ns+3*nt+1:6*ns+3*nt)
!write(6,'(a,/,4(3(f30.20,x)/))') 'Multiplication',DotRhoMultiplication
!write(6,'(a,/,4(3(f30.20,x)/))') 'DipFormation',DotRhoDipFormation
!write(6,'(a,/,4(3(f30.20,x)/))') 'SingleSingle',DotRhoEdgeEdgeAnnihilation
!write(6,'(a,/,4(3(f30.20,x)/))') 'SingleDipole',DotRhoEdgeDipAnnihilation
!write(6,'(a,/,4(3(f30.20,x)/))') 'DipClimb',DotRhoEdgeDipClimb

return
end function


pure function constitutive_dislotwin_dotTemperature(Tstar_v,Temperature,state,g,ip,el)
!*********************************************************************
!* rate of change of microstructure                                  *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotTemperature : evolution of Temperature         *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
real(pReal) constitutive_dislotwin_dotTemperature

constitutive_dislotwin_dotTemperature = 0.0_pReal

return
end function


pure function constitutive_dislotwin_postResults(Tstar_v,Temperature,dt,state,g,ip,el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Temperature     : temperature                                  *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID at current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
use prec,     only: pReal,pInt,p_vec
use math,     only: pi
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_constitutionInstance,phase_Noutput
use lattice,  only: lattice_Sslip_v,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin
implicit none

!* Definition of variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: dt,Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
integer(pInt) myInstance,myStructure,ns,nt,f,o,i,c,j,index_myFamily
real(pReal) sumf,tau,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r,gdot_slip,dgdot_dtauslip
real(pReal), dimension(constitutive_dislotwin_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
constitutive_dislotwin_postResults

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_dislotwin_structure(myInstance)
ns = constitutive_dislotwin_totalNslip(myInstance)
nt = constitutive_dislotwin_totalNtwin(myInstance)

!* Total twin volume fraction
sumf = sum(state(g,ip,el)%p((2*ns+1):(2*ns+nt))) ! safe for nt == 0

!* Required output
c = 0_pInt
constitutive_dislotwin_postResults = 0.0_pReal

do o = 1,phase_Noutput(material_phase(g,ip,el))
   select case(constitutive_dislotwin_output(o,myInstance))

     case ('edge_density')
       constitutive_dislotwin_postResults(c+1:c+ns) = state(g,ip,el)%p(1:ns)
       c = c + ns
     case ('dipole_density')
       constitutive_dislotwin_postResults(c+1:c+ns) = state(g,ip,el)%p(ns+1:2*ns)
       c = c + ns
     case ('shear_rate_slip')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
          do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
             j = j + 1_pInt

             !* Resolved shear stress on slip system
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))
             !* Stress ratios
             StressRatio_p = (abs(tau)/state(g,ip,el)%p(5*ns+3*nt+j))**constitutive_dislotwin_p(myInstance)
             StressRatio_pminus1 = (abs(tau)/state(g,ip,el)%p(5*ns+3*nt+j))**(constitutive_dislotwin_p(myInstance)-1.0_pReal)
             !* Boltzmann ratio
             BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(f,myInstance)/(kB*Temperature)
             !* Initial shear rates
             DotGamma0 = &
               state(g,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)* &
               constitutive_dislotwin_v0PerSlipSystem(f,myInstance)

             !* Shear rates due to slip
             constitutive_dislotwin_postResults(c+j) = &
               DotGamma0*exp(-BoltzmannRatio*(1-StressRatio_p)**constitutive_dislotwin_q(myInstance))*sign(1.0_pReal,tau)
       enddo ; enddo
       c = c + ns
     case ('mfp_slip')
       constitutive_dislotwin_postResults(c+1:c+ns) = state(g,ip,el)%p((4*ns+2*nt+1):(5*ns+2*nt))
       c = c + ns
     case ('resolved_stress_slip')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
          do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
             j = j + 1_pInt
             constitutive_dislotwin_postResults(c+j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))
       enddo; enddo
       c = c + ns
     case ('threshold_stress_slip')
       constitutive_dislotwin_postResults(c+1:c+ns) = state(g,ip,el)%p((5*ns+3*nt+1):(6*ns+3*nt))
       c = c + ns
     case ('edge_dipole_distance')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
          do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
             j = j + 1_pInt
             constitutive_dislotwin_postResults(c+j) = &
               (3.0_pReal*constitutive_dislotwin_Gmod(myInstance)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance))/&
               (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))))
             constitutive_dislotwin_postResults(c+j) = min(constitutive_dislotwin_postResults(c+j),state(g,ip,el)%p(4*ns+2*nt+j))
!            constitutive_dislotwin_postResults(c+j) = max(constitutive_dislotwin_postResults(c+j),state(g,ip,el)%p(4*ns+2*nt+j))
       enddo; enddo
       c = c + ns
     case ('twin_fraction')
       constitutive_dislotwin_postResults(c+1:c+nt) = state(g,ip,el)%p((2*ns+1):(2*ns+nt))
       c = c + nt
     case ('shear_rate_twin')
       if (nt > 0_pInt) then
         j = 0_pInt
         do f = 1,lattice_maxNtwinFamily                                 ! loop over all twin families
           index_myFamily = sum(lattice_NtwinSystem(1:f-1,myStructure)) ! at which index starts my family
           do i = 1,constitutive_dislotwin_Ntwin(f,myInstance)          ! process each (active) twin system in family
             j = j + 1_pInt

             !* Resolved shear stress on twin system
             tau = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))
             !* Stress ratios
             StressRatio_r = (state(g,ip,el)%p(6*ns+3*nt+j)/tau)**constitutive_dislotwin_r(myInstance)

             !* Shear rates and their derivatives due to twin
             if ( tau > 0.0_pReal ) then
               constitutive_dislotwin_postResults(c+j) = &
                 (constitutive_dislotwin_MaxTwinFraction(myInstance)-sumf)*&
                 state(g,ip,el)%p(6*ns+4*nt+j)*constitutive_dislotwin_Ndot0PerTwinSystem(f,myInstance)*exp(-StressRatio_r)
             endif

         enddo ; enddo
       endif
       c = c + nt
     case ('mfp_twin')
       constitutive_dislotwin_postResults(c+1:c+nt) = state(g,ip,el)%p((5*ns+2*nt+1):(5*ns+3*nt))
       c = c + nt
     case ('resolved_stress_twin')
       if (nt > 0_pInt) then
         j = 0_pInt
         do f = 1,lattice_maxNtwinFamily                                 ! loop over all slip families
           index_myFamily = sum(lattice_NtwinSystem(1:f-1,myStructure))  ! at which index starts my family
           do i = 1,constitutive_dislotwin_Ntwin(f,myInstance)           ! process each (active) slip system in family
             j = j + 1_pInt
             constitutive_dislotwin_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))
         enddo; enddo
       endif
       c = c + nt
     case ('threshold_stress_twin')
       constitutive_dislotwin_postResults(c+1:c+nt) = state(g,ip,el)%p((6*ns+3*nt+1):(6*ns+4*nt))
       c = c + nt
     case ('stress_exponent')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
          do i = 1,constitutive_dislotwin_Nslip(f,myInstance)          ! process each (active) slip system in family
             j = j + 1_pInt

             !* Resolved shear stress on slip system
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure))
             !* Stress ratios
             StressRatio_p = (abs(tau)/state(g,ip,el)%p(5*ns+3*nt+j))**constitutive_dislotwin_p(myInstance)
             StressRatio_pminus1 = (abs(tau)/state(g,ip,el)%p(5*ns+3*nt+j))**(constitutive_dislotwin_p(myInstance)-1.0_pReal)
             !* Boltzmann ratio
             BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(f,myInstance)/(kB*Temperature)
             !* Initial shear rates
             DotGamma0 = &
               state(g,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(f,myInstance)* &
               constitutive_dislotwin_v0PerSlipSystem(f,myInstance)

             !* Shear rates due to slip
             gdot_slip = &
               DotGamma0*exp(-BoltzmannRatio*(1-StressRatio_p)**constitutive_dislotwin_q(myInstance))*sign(1.0_pReal,tau)

             !* Derivatives of shear rates
             dgdot_dtauslip = &
               ((abs(gdot_slip)*BoltzmannRatio*&
               constitutive_dislotwin_p(myInstance)*constitutive_dislotwin_q(myInstance))/state(g,ip,el)%p(5*ns+3*nt+j))*&
               StressRatio_pminus1*(1-StressRatio_p)**(constitutive_dislotwin_q(myInstance)-1.0_pReal)

             !* Stress exponent
             if (gdot_slip==0.0_pReal) then
               constitutive_dislotwin_postResults(c+j) = 0.0_pReal
             else
               constitutive_dislotwin_postResults(c+j) = (tau/gdot_slip)*dgdot_dtauslip
             endif
       enddo ; enddo
       c = c + ns

   end select
enddo

return
end function

END MODULE
