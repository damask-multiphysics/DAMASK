!* $Id: constitutive_titanmod.f90 519 2010-03-24 08:17:27Z MPIE\f.roters $
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
! Parameters for titanium
! [titanmod]
! constitution            titanmod
! (output)                rhoedge
! (output)                rhoscrew
! (output)                velocity_edge
! (output)                velocity_screw
! (output)                rss_slip
! (output)                gamma_dot
! (output)                                resistance_edge
! (output)                                resistance_screw
! (output)                                segment_edge
! (output)                                segment_screw
! (output)                                total_generation
! (output)                                total_annihilation
! (output)                                total_density
! lattice_structure       hex
! covera_ratio            1.587
! c11                     162.2e9
! c12                     91.8e9
! c13                     68.8e9
! c33                     180.5e9
! c44                     46.7e9
! nslip                   3 3 6 0               # per family
! ntwin                   0 0 0 0                # per family
! rho_edge0                                1.66e9 1.66e9 1.66e9 1.66e9
! rho_screw0                                1.66e9 1.66e9 1.66e9 1.66e9
! slipburgers             2.86e-10 2.86e-10 2.86e-10 2.86e-10 # per family
! twinburgers                                2.86e-10 2.86e-10 2.86e-10 2.86e-10
! f0                                                3.0e-19 3.0e-19 3.0e-19 3.0e-19
! v0e                                                1e-1 10e0 1e0 1e0
! v0s                                                1e-1 1e0 1e0 1e0
! ndot0                                        1 1 1 1
! twinsize                                1e-05 5e-05 1e-05 1e-05
! celambdaslip                        1.0 1.0 1.0 1.0
! cslambdaslip                        1.0 1.0 1.0 1.0
! rlengthscrew                        1e-6 1e-6 1e-6 1e-6
! grainsize                                10e-6
! maxtwinfraction                        0.7
! pe                                                0.61 0.45 0.61 0.61
! ps                                                0.61 0.45 0.61 0.61
! qe                                                1.40 1.60 1.40 1.40
! qs                                                1.40 1.60 1.40 1.40
! rexponent                                0.2 0.2 0.2 0.2
! d0                                                0.2e-5 0.2e-5 0.2e-5 0.2e-5
! qsd                                                3.0e-19 3.0e-19 3.0e-19 3.0e-19
! cmfptwin                                0.5
! cthresholdtwin                        0.1
! cedgedipmindistance                50
! catomicvolume                        10e-15
! tau0e                   20e6 24e6 47e6 47e6               # per family
! tau0s                   20e6 24e6 50e6 50e6               # per family
! capre                   1.6e-9 5.6e-9 10.6e-9 10.6e-9
! caprs                   10.6e-9 110.6e-9 10.6e-9 10.6e-9
! #                                                1     2    3          4         5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20                
! interactionslipslip     0.15 0.07 0.10 0.10 0.25 0.15 0.15 0.15 0.65 0.15 0.65 0.65 0.65 0.65 0.65 0.15 0.65 0.65 0.65 0.65
! interactionsliptwin     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! interactiontwinslip     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! interactiontwintwin     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! relevantRho             1.e1

! Parameters for aluminum
! [alumod]
! constitution            titanmod
! (output)                rhoedge
! (output)                rhoscrew
! (output)                                velocity_edge
! (output)                                velocity_screw
! (output)                                rss_slip
! (output)                                gamma_dot
! (output)                                dgdotdtau
! (output)                                resistance_edge
! (output)                                resistance_screw
! (output)                                segment_edge
! (output)                                segment_screw
! (output)                                total_generation
! (output)                                total_annihilation
! (output)                                total_density
! (output)                                stressratio_edgep
! (output)                                stressratio_screwp
! lattice_structure       fcc
! covera_ratio            1.587
! c11                     106.75e9
! c12                     60.41e9
! c44                     28.34e9
! c13                     68.8e9
! c33                     180.5e9
! nslip                   12 0 0 0               # per family
! ntwin                   0 0 0 0                # per family
! rho_edge0                                5.56e8 0 0 0
! rho_screw0                                5.56e8 0 0 0
! slipburgers             2.86e-10 0 0 0 # per family
! twinburgers                                2.86e-10 0 0 0
! f0                                                3.2e-19 0 0 0
! v0e                                                1e-1 0 0 0
! v0s                                                1e-1 0 0 0
! ndot0                                        1 0 0 0
! twinsize                                1e-05 0 0 0
! celambdaslip                        1.0 0 0 0
! cslambdaslip                        1.0 0 0 0
! rlengthscrew                        1e-6 1e-6 1e-6 1e-6
! grainsize                                10e-6
! maxtwinfraction                        0.7
! pe                                                0.11 0 0 0
! ps                                                0.11 0 0 0
! qe                                                1.41 0 0 0
! qs                                                1.41 0 0 0
! rexponent                                0.2 0 0 0
! d0                                                0.2e-5 0 0 0
! qsd                                                3.0e-19 0 0 0
! cmfptwin                                0.5
! cthresholdtwin                        0.1
! cedgedipmindistance                50
! catomicvolume                        10e-15
! tau0e                   5e6 0 0 0               # per family
! tau0s                   5e6 0 0 0               # per family
! capre                   10.6e-9 0 0 0
! caprs                   100.0e-9 0 0 0
! #interactionslipslip        0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
! interactionslipslip     0.10 0.22 0.30 0.16 0.38 0.45
! # Devincre : 0.122 0.122 0.625 0.07 0.137 0.122
! # Arsenlis : G0=0.10, G1=0.22, G2=0.30, G3=0.38, G4=0.16, G5=0.40
! #G3=glissile, G4=Hirth lock, G5=Lomer-Cottrell lock, G2=cross-slip interaction, G0=same Burgers vector and parallel plane (self interaction)
! #G1=different Burgers vectors but parallel slip planes
! #! Interaction types
! #! 1 --- self interaction, G0
! #! 2 --- coplanar interaction, G1
! #! 3 --- collinear interaction, G2
! #! 4 --- Hirth locks, G3
! #! 5 --- glissile junctions, G4
! #! 6 --- Lomer locks, G5
! interactionsliptwin     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! interactiontwinslip     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! interactiontwintwin     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
! relevantRho             1.e1

MODULE constitutive_titanmod

!* Include other modules
use prec, only: pReal,pInt
implicit none

!* Lists of states and physical parameters
character(len=*), parameter :: constitutive_titanmod_label = 'titanmod'
character(len=18), dimension(2), parameter:: constitutive_titanmod_listBasicSlipStates = (/'rho_edge', &
                                                                                           'rho_screw'/)
character(len=18), dimension(3), parameter:: constitutive_titanmod_listBasicTwinStates = (/'twinrho_edge', & 
                                                                                           'twinrho_screw', &
                                                                                           'gdot_twin'/)
                                                                                            
character(len=18), dimension(8), parameter:: constitutive_titanmod_listDependentSlipStates =(/'segment_edge', &
                                                                                              'segment_screw', &
                                                                                              'resistance_edge', &
                                                                                              'resistance_screw', &
                                                                                              'tau_slip', &
                                                                                              'gdot_slip', &
                                                                                              'velocity_edge', &
                                                                                              'velocity_screw' &
                                                                                              /)

character(len=18), dimension(8), parameter:: constitutive_titanmod_listDependentTwinStates =(/'twin_fraction', &
                                                                                              'tau_twin', &
                                                                                              'twinsegment_edge', &
                                                                                              'twinsegment_screw', &
                                                                                              'twinresistance_edge', &
                                                                                              'twinresistance_screw', &
                                                                                              'twinvelocity_edge', &
                                                                                              'twinvelocity_screw' &
                                                                                             /)
real(pReal), parameter :: kB = 1.38e-23_pReal ! Boltzmann constant in J/Kelvin

!* Definition of global variables
integer(pInt), dimension(:), allocatable ::               constitutive_titanmod_sizeDotState, &                ! number of dotStates
                                                          constitutive_titanmod_sizeState, &                   ! total number of microstructural state variables
                                                          constitutive_titanmod_sizePostResults                ! cumulative size of post results
integer(pInt), dimension(:,:), allocatable, target ::     constitutive_titanmod_sizePostResult                 ! size of each post result output
character(len=64), dimension(:,:), allocatable, target :: constitutive_titanmod_output                         ! name of each post result output 
character(len=32), dimension(:), allocatable ::           constitutive_titanmod_structureName                  ! name of the lattice structure
integer(pInt), dimension(:), allocatable ::               constitutive_titanmod_structure, &                   ! number representing the kind of lattice structure
                                                          constitutive_titanmod_totalNslip, &                  ! total number of active slip systems for each instance
                                                          constitutive_titanmod_totalNtwin                     ! total number of active twin systems for each instance
integer(pInt), dimension(:,:), allocatable ::             constitutive_titanmod_Nslip, &                       ! number of active slip systems for each family and instance
                                                          constitutive_titanmod_Ntwin, &                       ! number of active twin systems for each family and instance
                                                          constitutive_titanmod_slipFamily, &                  ! lookup table relating active slip system to slip family for each instance
                                                          constitutive_titanmod_twinFamily, &                  ! lookup table relating active twin system to twin family for each instance
                                                          constitutive_titanmod_slipSystemLattice, &           ! lookup table relating active slip system index to lattice slip system index for each instance
                                                          constitutive_titanmod_twinSystemLattice              ! lookup table relating active twin system index to lattice twin system index for each instance
real(pReal), dimension(:), allocatable ::                 constitutive_titanmod_CoverA, &                      ! c/a ratio for hex type lattice
                                                          constitutive_titanmod_C11, &                         ! C11 element in elasticity matrix
                                                          constitutive_titanmod_C12, &                         ! C12 element in elasticity matrix
                                                          constitutive_titanmod_C13, &                         ! C13 element in elasticity matrix
                                                          constitutive_titanmod_C33, &                         ! C33 element in elasticity matrix
                                                          constitutive_titanmod_C44, &                         ! C44 element in elasticity matrix
                                                          constitutive_titanmod_Gmod, &                        ! shear modulus
                                                          constitutive_titanmod_CAtomicVolume, &               ! atomic volume in Bugers vector unit
                                                          constitutive_titanmod_dc, &                          ! prefactor for self-diffusion coefficient
                                                          constitutive_titanmod_twinhpconstant, &                         ! activation energy for dislocation climb
                                                          constitutive_titanmod_GrainSize, &                   ! grain size - Not being used
                                                          constitutive_titanmod_MaxTwinFraction, &             ! maximum allowed total twin volume fraction
                                                          constitutive_titanmod_r, &                           ! r-exponent in twin nucleation rate
                                                          constitutive_titanmod_CEdgeDipMinDistance, &         ! Not being used
                                                          constitutive_titanmod_Cmfptwin, &                    ! Not being used
                                                          constitutive_titanmod_Cthresholdtwin, &              ! Not being used
                                                          constitutive_titanmod_relevantRho                    ! dislocation density considered relevant                                                                                                  
real(pReal),       dimension(:,:,:),       allocatable :: constitutive_titanmod_Cslip_66                       ! elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:),     allocatable :: constitutive_titanmod_Ctwin_66                       ! twin elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:,:),   allocatable :: constitutive_titanmod_Cslip_3333                     ! elasticity matrix for each instance
real(pReal),       dimension(:,:,:,:,:,:), allocatable :: constitutive_titanmod_Ctwin_3333                     ! twin elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_titanmod_rho_edge0, &                   ! initial edge dislocation density per slip system for each family and instance
                                                          constitutive_titanmod_rho_screw0, &                  ! initial screw dislocation density per slip system for each family and instance
                                                          constitutive_titanmod_twinrho_edge0, &                   ! initial edge dislocation density per twin system for each family and instance
                                                          constitutive_titanmod_twinrho_screw0, &                  ! initial screw dislocation density per twin system for each family and instance
                                                          constitutive_titanmod_burgersPerSlipFamily, &        ! absolute length of burgers vector [m] for each slip family and instance
                                                          constitutive_titanmod_burgersPerSlipSystem, &        ! absolute length of burgers vector [m] for each slip system and instance
                                                          constitutive_titanmod_burgersPerTwinFamily, &        ! absolute length of burgers vector [m] for each twin family and instance
                                                          constitutive_titanmod_burgersPerTwinSystem, &        ! absolute length of burgers vector [m] for each twin system and instance
                                                          constitutive_titanmod_f0_PerSlipFamily, &            ! activation energy for glide [J] for each slip family and instance
                                                          constitutive_titanmod_f0_PerSlipSystem, &            ! activation energy for glide [J] for each slip system and instance
                                                          constitutive_titanmod_twinf0_PerTwinFamily, &        ! activation energy for glide [J] for each twin family and instance
                                                          constitutive_titanmod_twinf0_PerTwinSystem, &        ! activation energy for glide [J] for each twin system and instance
                                                          constitutive_titanmod_twinshearconstant_PerTwinFamily, &        ! activation energy for glide [J] for each twin family and instance
                                                          constitutive_titanmod_twinshearconstant_PerTwinSystem, &        ! activation energy for glide [J] for each twin system and instance
                                                          constitutive_titanmod_tau0e_PerSlipFamily, &         ! Initial yield stress for edge dislocations per slip family
                                                          constitutive_titanmod_tau0e_PerSlipSystem, &         ! Initial yield stress for edge dislocations per slip system
                                                          constitutive_titanmod_tau0s_PerSlipFamily, &         ! Initial yield stress for screw dislocations per slip family
                                                          constitutive_titanmod_tau0s_PerSlipSystem, &         ! Initial yield stress for screw dislocations per slip system
                                                          constitutive_titanmod_twintau0e_PerTwinFamily, &         ! Initial yield stress for edge dislocations per twin family
                                                          constitutive_titanmod_twintau0e_PerTwinSystem, &         ! Initial yield stress for edge dislocations per twin system
                                                          constitutive_titanmod_twintau0s_PerTwinFamily, &         ! Initial yield stress for screw dislocations per twin family
                                                          constitutive_titanmod_twintau0s_PerTwinSystem, &         ! Initial yield stress for screw dislocations per twin system
                                                          constitutive_titanmod_capre_PerSlipFamily, &         ! Capture radii for edge dislocations per slip family
                                                          constitutive_titanmod_capre_PerSlipSystem, &         ! Capture radii for edge dislocations per slip system
                                                          constitutive_titanmod_caprs_PerSlipFamily, &         ! Capture radii for screw dislocations per slip family
                                                          constitutive_titanmod_caprs_PerSlipSystem, &         ! Capture radii for screw dislocations per slip system
                                                          constitutive_titanmod_twincapre_PerTwinFamily, &         ! Capture radii for edge dislocations per twin family
                                                          constitutive_titanmod_twincapre_PerTwinSystem, &         ! Capture radii for edge dislocations per twin system
                                                          constitutive_titanmod_twincaprs_PerTwinFamily, &         ! Capture radii for screw dislocations per twin family
                                                          constitutive_titanmod_twincaprs_PerTwinSystem, &         ! Capture radii for screw dislocations per twin system
                                                          constitutive_titanmod_pe_PerSlipFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_ps_PerSlipFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_qe_PerSlipFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_qs_PerSlipFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_pe_PerSlipSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_ps_PerSlipSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_qe_PerSlipSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_qs_PerSlipSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinpe_PerTwinFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinps_PerTwinFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinqe_PerTwinFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinqs_PerTwinFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinpe_PerTwinSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinps_PerTwinSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinqe_PerTwinSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinqs_PerTwinSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_v0e_PerSlipFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_v0e_PerSlipSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_v0s_PerSlipFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_v0s_PerSlipSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_twinv0e_PerTwinFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_twinv0e_PerTwinSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_twinv0s_PerTwinFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_twinv0s_PerTwinSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_rlengthscrew_PerSlipFamily, &  ! screw dislocation mobility prefactor for kink-pairs per slip family
                                                          constitutive_titanmod_rlengthscrew_PerSlipSystem, &  ! screw dislocation mobility prefactor for kink-pairs per slip system
                                                          constitutive_titanmod_Ndot0PerTwinFamily, &          ! twin nucleation rate [1/m³s] for each twin family and instance
                                                          constitutive_titanmod_Ndot0PerTwinSystem, &          ! twin nucleation rate [1/m³s] for each twin system and instance
                                                          constitutive_titanmod_twinsizePerTwinFamily, &       ! twin thickness [m] for each twin family and instance
                                                          constitutive_titanmod_twinsizePerTwinSystem, &       ! twin thickness [m] for each twin system and instance
                                                          constitutive_titanmod_CeLambdaSlipPerSlipFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_CeLambdaSlipPerSlipSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_CsLambdaSlipPerSlipFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_CsLambdaSlipPerSlipSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_twinCeLambdaSlipPerTwinFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_twinCeLambdaSlipPerTwinSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_twinCsLambdaSlipPerTwinFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_twinCsLambdaSlipPerTwinSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_interactionSlipSlip, &         ! coefficients for slip-slip interaction for each interaction type and instance
                                                          constitutive_titanmod_interaction_ee, &                    ! coefficients for e-e interaction for each interaction type and instance
                                                          constitutive_titanmod_interaction_ss, &               ! coefficients for s-s interaction for each interaction type and instance
                                                          constitutive_titanmod_interaction_es, &               ! coefficients for e-s-twin interaction for each interaction type and instance
                                                          constitutive_titanmod_interactionSlipTwin, &         ! coefficients for twin-slip interaction for each interaction type and instance
                                                          constitutive_titanmod_interactionTwinSlip, &         ! coefficients for twin-slip interaction for each interaction type and instance
                                                          constitutive_titanmod_interactionTwinTwin            ! coefficients for twin-twin interaction for each interaction type and instance
real(pReal),       dimension(:,:,:),       allocatable :: constitutive_titanmod_interactionMatrixSlipSlip, &   ! interaction matrix of the different slip systems for each instance
                                                          constitutive_titanmod_interactionMatrix_ee, &         ! interaction matrix of e-e for each instance
                                                          constitutive_titanmod_interactionMatrix_ss, &         ! interaction matrix of s-s for each instance
                                                          constitutive_titanmod_interactionMatrix_es, &         ! interaction matrix of e-s for each instance
                                                          constitutive_titanmod_interactionMatrixSlipTwin, &   ! interaction matrix of slip systems with twin systems for each instance
                                                          constitutive_titanmod_interactionMatrixTwinSlip, &   ! interaction matrix of twin systems with slip systems for each instance
                                                          constitutive_titanmod_interactionMatrixTwinTwin, &   ! interaction matrix of the different twin systems for each instance                                                          
                                                          constitutive_titanmod_forestProjectionEdge, &           ! matrix of forest projections of edge dislocations for each instance  
                                                          constitutive_titanmod_forestProjectionScrew, &           ! matrix of forest projections of screw dislocations for each instance  
                                                          constitutive_titanmod_TwinforestProjectionEdge, &     ! matrix of forest projections of edge dislocations in twin system for each instance  
                                                          constitutive_titanmod_TwinforestProjectionScrew       ! matrix of forest projections of screw dislocations in twin system for each instance  
CONTAINS
!****************************************
!* - constitutive_titanmod_init
!* - constitutive_titanmod_stateInit
!* - constitutive_titanmod_relevantState
!* - constitutive_titanmod_homogenizedC
!* - constitutive_titanmod_microstructure
!* - constitutive_titanmod_LpAndItsTangent
!* - consistutive_titanmod_dotState
!* - constitutive_titanmod_dotTemperature
!* - consistutive_titanmod_postResults
!****************************************


subroutine constitutive_titanmod_init(file)
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
integer(pInt) section,maxNinstance,f,i,j,k,l,m,n,o,p,q,r,s,s1,s2,t1,t2,ns,nt,output,mySize,myStructure,maxTotalNslip, &
maxTotalNtwin
character(len=64) tag
character(len=1024) line

!write(6,*)
!write(6,'(a20,a20,a12)') '<<<+-  constitutive_',constitutive_titanmod_label,' init  -+>>>'
!write(6,*) '$Id: constitutive_titanmod.f90 519 2010-03-24 08:17:27Z MPIE\f.roters $'
!write(6,*)

maxNinstance = count(phase_constitution == constitutive_titanmod_label)
if (maxNinstance == 0) return

!* Space allocation for global variables
allocate(constitutive_titanmod_sizeDotState(maxNinstance))       
allocate(constitutive_titanmod_sizeState(maxNinstance)) 
allocate(constitutive_titanmod_sizePostResults(maxNinstance)) 
allocate(constitutive_titanmod_sizePostResult(maxval(phase_Noutput),maxNinstance))
allocate(constitutive_titanmod_output(maxval(phase_Noutput),maxNinstance))
constitutive_titanmod_sizeDotState    = 0_pInt
constitutive_titanmod_sizeState       = 0_pInt
constitutive_titanmod_sizePostResults = 0_pInt
constitutive_titanmod_sizePostResult  = 0_pInt
constitutive_titanmod_output          = ''

allocate(constitutive_titanmod_structureName(maxNinstance)) 
allocate(constitutive_titanmod_structure(maxNinstance))
allocate(constitutive_titanmod_Nslip(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_Ntwin(lattice_maxNtwinFamily,maxNinstance))  
allocate(constitutive_titanmod_slipFamily(lattice_maxNslip,maxNinstance))
allocate(constitutive_titanmod_twinFamily(lattice_maxNtwin,maxNinstance))
allocate(constitutive_titanmod_slipSystemLattice(lattice_maxNslip,maxNinstance))
allocate(constitutive_titanmod_twinSystemLattice(lattice_maxNtwin,maxNinstance)) 
allocate(constitutive_titanmod_totalNslip(maxNinstance))   
allocate(constitutive_titanmod_totalNtwin(maxNinstance))   
constitutive_titanmod_structureName     = ''
constitutive_titanmod_structure         = 0_pInt
constitutive_titanmod_Nslip             = 0_pInt
constitutive_titanmod_Ntwin             = 0_pInt
constitutive_titanmod_slipFamily        = 0_pInt
constitutive_titanmod_twinFamily        = 0_pInt
constitutive_titanmod_slipSystemLattice = 0.0_pReal
constitutive_titanmod_twinSystemLattice = 0.0_pReal
constitutive_titanmod_totalNslip        = 0_pInt
constitutive_titanmod_totalNtwin        = 0_pInt
allocate(constitutive_titanmod_CoverA(maxNinstance))
allocate(constitutive_titanmod_C11(maxNinstance))
allocate(constitutive_titanmod_C12(maxNinstance))
allocate(constitutive_titanmod_C13(maxNinstance))
allocate(constitutive_titanmod_C33(maxNinstance))
allocate(constitutive_titanmod_C44(maxNinstance))
allocate(constitutive_titanmod_Gmod(maxNinstance))
allocate(constitutive_titanmod_CAtomicVolume(maxNinstance))
allocate(constitutive_titanmod_dc(maxNinstance))
allocate(constitutive_titanmod_twinhpconstant(maxNinstance))
allocate(constitutive_titanmod_GrainSize(maxNinstance))
allocate(constitutive_titanmod_MaxTwinFraction(maxNinstance))
allocate(constitutive_titanmod_r(maxNinstance))
allocate(constitutive_titanmod_CEdgeDipMinDistance(maxNinstance))
allocate(constitutive_titanmod_Cmfptwin(maxNinstance))
allocate(constitutive_titanmod_Cthresholdtwin(maxNinstance))
allocate(constitutive_titanmod_relevantRho(maxNinstance))
allocate(constitutive_titanmod_Cslip_66(6,6,maxNinstance))
allocate(constitutive_titanmod_Cslip_3333(3,3,3,3,maxNinstance))
constitutive_titanmod_CoverA              = 0.0_pReal 
constitutive_titanmod_C11                 = 0.0_pReal
constitutive_titanmod_C12                 = 0.0_pReal
constitutive_titanmod_C13                 = 0.0_pReal
constitutive_titanmod_C33                 = 0.0_pReal
constitutive_titanmod_C44                 = 0.0_pReal
constitutive_titanmod_Gmod                = 0.0_pReal
constitutive_titanmod_CAtomicVolume       = 0.0_pReal
constitutive_titanmod_dc                  = 0.0_pReal
constitutive_titanmod_twinhpconstant      = 0.0_pReal
constitutive_titanmod_GrainSize           = 0.0_pReal
constitutive_titanmod_MaxTwinFraction     = 0.0_pReal
constitutive_titanmod_r                   = 0.0_pReal
constitutive_titanmod_CEdgeDipMinDistance = 0.0_pReal
constitutive_titanmod_Cmfptwin            = 0.0_pReal
constitutive_titanmod_Cthresholdtwin      = 0.0_pReal
constitutive_titanmod_relevantRho         = 0.0_pReal
constitutive_titanmod_Cslip_66            = 0.0_pReal
constitutive_titanmod_Cslip_3333          = 0.0_pReal
allocate(constitutive_titanmod_rho_edge0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_rho_screw0(lattice_maxNslipFamily,maxNinstance)) 
allocate(constitutive_titanmod_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_burgersPerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_titanmod_f0_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_tau0e_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_tau0s_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_capre_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_caprs_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_pe_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_ps_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_qe_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_qs_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_v0e_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_v0s_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_rlengthscrew_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_Ndot0PerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_titanmod_CeLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_CsLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance))

allocate(constitutive_titanmod_twinrho_edge0(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinrho_screw0(lattice_maxNtwinFamily,maxNinstance)) 
allocate(constitutive_titanmod_twinf0_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinshearconstant_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twintau0e_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twintau0s_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twincapre_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twincaprs_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinpe_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinps_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinqe_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinqs_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinv0e_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinv0s_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinCeLambdaSlipPerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinCsLambdaSlipPerTwinFamily(lattice_maxNTwinFamily,maxNinstance))

constitutive_titanmod_rho_edge0                 = 0.0_pReal
constitutive_titanmod_rho_screw0              = 0.0_pReal
constitutive_titanmod_burgersPerSlipFamily     = 0.0_pReal
constitutive_titanmod_burgersPerTwinFamily     = 0.0_pReal
constitutive_titanmod_f0_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_tau0e_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_tau0s_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_capre_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_caprs_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_v0e_PerSlipFamily          = 0.0_pReal
constitutive_titanmod_v0s_PerSlipFamily          = 0.0_pReal
constitutive_titanmod_rlengthscrew_PerSlipFamily = 0.0_pReal
constitutive_titanmod_Ndot0PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twinsizePerTwinFamily    = 0.0_pReal
constitutive_titanmod_CeLambdaSlipPerSlipFamily = 0.0_pReal
constitutive_titanmod_CsLambdaSlipPerSlipFamily = 0.0_pReal
constitutive_titanmod_pe_PerSlipFamily = 0.0_pReal
constitutive_titanmod_ps_PerSlipFamily = 0.0_pReal
constitutive_titanmod_qe_PerSlipFamily = 0.0_pReal
constitutive_titanmod_qs_PerSlipFamily = 0.0_pReal

constitutive_titanmod_twinrho_edge0                 = 0.0_pReal
constitutive_titanmod_twinrho_screw0              = 0.0_pReal
constitutive_titanmod_twinf0_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twinshearconstant_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twintau0e_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twintau0s_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twincapre_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twincaprs_PerTwinFamily       = 0.0_pReal
constitutive_titanmod_twinv0e_PerTwinFamily          = 0.0_pReal
constitutive_titanmod_twinv0s_PerTwinFamily          = 0.0_pReal
constitutive_titanmod_twinCeLambdaSlipPerTwinFamily = 0.0_pReal
constitutive_titanmod_twinCsLambdaSlipPerTwinFamily = 0.0_pReal
constitutive_titanmod_twinpe_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinps_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinqe_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinqs_PerTwinFamily = 0.0_pReal

allocate(constitutive_titanmod_interactionSlipSlip(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interaction_ee(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interaction_ss(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interaction_es(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interactionSlipTwin(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interactionTwinSlip(lattice_maxNinteraction,maxNinstance)) 
allocate(constitutive_titanmod_interactionTwinTwin(lattice_maxNinteraction,maxNinstance)) 
constitutive_titanmod_interactionSlipSlip = 0.0_pReal
constitutive_titanmod_interaction_ee = 0.0_pReal
constitutive_titanmod_interaction_ss = 0.0_pReal
constitutive_titanmod_interaction_ss = 0.0_pReal
constitutive_titanmod_interactionSlipTwin = 0.0_pReal
constitutive_titanmod_interactionTwinSlip = 0.0_pReal
constitutive_titanmod_interactionTwinTwin = 0.0_pReal

!* Readout data from material.config file
rewind(file)
line = ''
section = 0

write(6,*) 'Reading material parameters from material config file'

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
   if (section > 0 .and. phase_constitution(section) == constitutive_titanmod_label) then  ! one of my sections
     i = phase_constitutionInstance(section)     ! which instance of my constitution is present phase
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1
         constitutive_titanmod_output(output,i) = IO_lc(IO_stringValue(line,positions,2))
                write(6,*) tag
       case ('lattice_structure')
              constitutive_titanmod_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
                write(6,*) tag
       case ('covera_ratio')
              constitutive_titanmod_CoverA(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('c11')
              constitutive_titanmod_C11(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_C11(i)
       case ('c12')
              constitutive_titanmod_C12(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_C12(i)
       case ('c13')
              constitutive_titanmod_C13(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_C13(i)
       case ('c33')
              constitutive_titanmod_C33(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_C33(i)
       case ('c44')
              constitutive_titanmod_C44(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_C44(i)
       case ('nslip')
            forall (j = 1:lattice_maxNslipFamily) &
            constitutive_titanmod_Nslip(j,i) = IO_intValue(line,positions,1+j)
            write(6,*) tag,constitutive_titanmod_Nslip(1,i),constitutive_titanmod_Nslip(2,i),constitutive_titanmod_Nslip(3,i), &
                    constitutive_titanmod_Nslip(4,i)
       case ('ntwin')
            forall (j = 1:lattice_maxNtwinFamily) &
             constitutive_titanmod_Ntwin(j,i) = IO_intValue(line,positions,1+j)
            write(6,*) tag,constitutive_titanmod_Ntwin(1,i),constitutive_titanmod_Ntwin(2,i),constitutive_titanmod_Ntwin(3,i), &
                        constitutive_titanmod_Ntwin(4,i)
       case ('rho_edge0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_rho_edge0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_rho_edge0(1,i),constitutive_titanmod_rho_edge0(2,i), &
                        constitutive_titanmod_rho_edge0(3,i), constitutive_titanmod_rho_edge0(4,i)
       case ('rho_screw0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_rho_screw0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_rho_screw0(1,i),constitutive_titanmod_rho_screw0(2,i), &
                        constitutive_titanmod_rho_screw0(3,i), constitutive_titanmod_rho_screw0(4,i)
       case ('twinrho_edge0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinrho_edge0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinrho_edge0(1,i),constitutive_titanmod_twinrho_edge0(2,i), &
                        constitutive_titanmod_twinrho_edge0(3,i), constitutive_titanmod_twinrho_edge0(4,i)
       case ('twinrho_screw0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinrho_screw0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinrho_screw0(1,i),constitutive_titanmod_twinrho_screw0(2,i), &
                        constitutive_titanmod_twinrho_screw0(3,i), constitutive_titanmod_twinrho_screw0(4,i)
       case ('slipburgers')
             forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_burgersPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
              write(6,*) tag,constitutive_titanmod_burgersPerSlipFamily(1,i),constitutive_titanmod_burgersPerSlipFamily(2,i), &
                        constitutive_titanmod_burgersPerSlipFamily(3,i), constitutive_titanmod_burgersPerSlipFamily(4,i)
       case ('twinburgers')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_burgersPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('f0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_f0_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_f0_PerSlipFamily(1,i),constitutive_titanmod_f0_PerSlipFamily(2,i), &
                        constitutive_titanmod_f0_PerSlipFamily(3,i), constitutive_titanmod_f0_PerSlipFamily(4,i)
       case ('twinf0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinf0_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinf0_PerTwinFamily(1,i),constitutive_titanmod_twinf0_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinf0_PerTwinFamily(3,i), constitutive_titanmod_twinf0_PerTwinFamily(4,i)
       case ('tau0e')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_tau0e_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_tau0e_PerSlipFamily(1,i),constitutive_titanmod_tau0e_PerSlipFamily(2,i), &
                        constitutive_titanmod_tau0e_PerSlipFamily(3,i), constitutive_titanmod_tau0e_PerSlipFamily(4,i)
       case ('twintau0e')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twintau0e_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twintau0e_PerTwinFamily(1,i),constitutive_titanmod_twintau0e_PerTwinFamily(2,i), &
                        constitutive_titanmod_twintau0e_PerTwinFamily(3,i), constitutive_titanmod_twintau0e_PerTwinFamily(4,i)
       case ('tau0s')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_tau0s_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_tau0s_PerSlipFamily(1,i),constitutive_titanmod_tau0s_PerSlipFamily(2,i), &
                        constitutive_titanmod_tau0s_PerSlipFamily(3,i), constitutive_titanmod_tau0s_PerSlipFamily(4,i)
       case ('twintau0s')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twintau0s_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twintau0s_PerTwinFamily(1,i),constitutive_titanmod_twintau0s_PerTwinFamily(2,i), &
                        constitutive_titanmod_twintau0s_PerTwinFamily(3,i), constitutive_titanmod_twintau0s_PerTwinFamily(4,i)
       case ('capre')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_capre_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_capre_PerSlipFamily(1,i),constitutive_titanmod_capre_PerSlipFamily(2,i), &
                        constitutive_titanmod_capre_PerSlipFamily(3,i), constitutive_titanmod_capre_PerSlipFamily(4,i)
       case ('twincapre')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twincapre_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twincapre_PerTwinFamily(1,i),constitutive_titanmod_twincapre_PerTwinFamily(2,i), &
                        constitutive_titanmod_twincapre_PerTwinFamily(3,i), constitutive_titanmod_twincapre_PerTwinFamily(4,i)
       case ('caprs')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_caprs_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_caprs_PerSlipFamily(1,i),constitutive_titanmod_caprs_PerSlipFamily(2,i), &
                        constitutive_titanmod_caprs_PerSlipFamily(3,i), constitutive_titanmod_caprs_PerSlipFamily(4,i)
       case ('twincaprs')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twincaprs_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twincaprs_PerTwinFamily(1,i),constitutive_titanmod_twincaprs_PerTwinFamily(2,i), &
                        constitutive_titanmod_twincaprs_PerTwinFamily(3,i), constitutive_titanmod_twincaprs_PerTwinFamily(4,i)
       case ('v0e')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_v0e_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_v0e_PerSlipFamily(1,i),constitutive_titanmod_v0e_PerSlipFamily(2,i), &
                        constitutive_titanmod_v0e_PerSlipFamily(3,i), constitutive_titanmod_v0e_PerSlipFamily(4,i)
       case ('twinv0e')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinv0e_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinv0e_PerTwinFamily(1,i),constitutive_titanmod_twinv0e_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinv0e_PerTwinFamily(3,i), constitutive_titanmod_twinv0e_PerTwinFamily(4,i)
       case ('v0s')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_v0s_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_v0s_PerSlipFamily(1,i),constitutive_titanmod_v0s_PerSlipFamily(2,i), &
                        constitutive_titanmod_v0s_PerSlipFamily(3,i), constitutive_titanmod_v0s_PerSlipFamily(4,i)
       case ('twinv0s')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinv0s_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinv0s_PerTwinFamily(1,i),constitutive_titanmod_twinv0s_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinv0s_PerTwinFamily(3,i), constitutive_titanmod_twinv0s_PerTwinFamily(4,i)
       case ('rlengthscrew')
       forall (j = 1:lattice_maxNslipFamily) &
       constitutive_titanmod_rlengthscrew_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
       write(6,*) tag,constitutive_titanmod_rlengthscrew_PerSlipFamily(1,i), &
         constitutive_titanmod_rlengthscrew_PerSlipFamily(2,i), &
            constitutive_titanmod_rlengthscrew_PerSlipFamily(3,i), constitutive_titanmod_rlengthscrew_PerSlipFamily(4,i)
       case ('ndot0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_Ndot0PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('twinsize')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinsizePerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('celambdaslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_CeLambdaSlipPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('twincelambdaslip')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twincelambdaslipPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twincelambdaslipPerTwinFamily(1,i),constitutive_titanmod_twincelambdaslipPerTwinFamily(2,i), &
                        constitutive_titanmod_twincelambdaslipPerTwinFamily(3,i), constitutive_titanmod_twincelambdaslipPerTwinFamily(4,i)
       case ('cslambdaslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_CsLambdaSlipPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('twincslambdaslip')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twincslambdaslipPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twincslambdaslipPerTwinFamily(1,i),constitutive_titanmod_twincslambdaslipPerTwinFamily(2,i), &
                        constitutive_titanmod_twincslambdaslipPerTwinFamily(3,i), constitutive_titanmod_twincslambdaslipPerTwinFamily(4,i)
       case ('grainsize')
              constitutive_titanmod_GrainSize(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('maxtwinfraction')
              constitutive_titanmod_MaxTwinFraction(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('pe')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_pe_PerSlipFamily(j,i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_pe_PerSlipFamily(1,i),constitutive_titanmod_pe_PerSlipFamily(2,i), &
                        constitutive_titanmod_pe_PerSlipFamily(3,i), constitutive_titanmod_pe_PerSlipFamily(4,i),i
       case ('twinpe')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinpe_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinpe_PerTwinFamily(1,i),constitutive_titanmod_twinpe_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinpe_PerTwinFamily(3,i), constitutive_titanmod_twinpe_PerTwinFamily(4,i)
       case ('ps')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_ps_PerSlipFamily(j,i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_ps_PerSlipFamily(1,i),constitutive_titanmod_ps_PerSlipFamily(2,i), &
                        constitutive_titanmod_ps_PerSlipFamily(3,i), constitutive_titanmod_ps_PerSlipFamily(4,i),i
       case ('twinps')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinps_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinps_PerTwinFamily(1,i),constitutive_titanmod_twinps_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinps_PerTwinFamily(3,i), constitutive_titanmod_twinps_PerTwinFamily(4,i)
       case ('qe')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_qe_PerSlipFamily(j,i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_qe_PerSlipFamily(1,i),constitutive_titanmod_qe_PerSlipFamily(2,i), &
                        constitutive_titanmod_qe_PerSlipFamily(3,i), constitutive_titanmod_qe_PerSlipFamily(4,i),i
       case ('twinqe')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinqe_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinqe_PerTwinFamily(1,i),constitutive_titanmod_twinqe_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinqe_PerTwinFamily(3,i), constitutive_titanmod_twinqe_PerTwinFamily(4,i)
       case ('qs')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_qs_PerSlipFamily(j,i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_qs_PerSlipFamily(1,i),constitutive_titanmod_qs_PerSlipFamily(2,i), &
                        constitutive_titanmod_qs_PerSlipFamily(3,i), constitutive_titanmod_qs_PerSlipFamily(4,i),i
       case ('twinqs')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinqs_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinqs_PerTwinFamily(1,i),constitutive_titanmod_twinqs_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinqs_PerTwinFamily(3,i), constitutive_titanmod_twinqs_PerTwinFamily(4,i)
       case ('twinshearconstant')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinshearconstant_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinshearconstant_PerTwinFamily(1,i), &
                 constitutive_titanmod_twinshearconstant_PerTwinFamily(2,i), &
                        constitutive_titanmod_twinshearconstant_PerTwinFamily(3,i), &
                           constitutive_titanmod_twinshearconstant_PerTwinFamily(4,i)
       case ('dc')
              constitutive_titanmod_dc(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('twinhpconstant')
              constitutive_titanmod_twinhpconstant(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('relevantrho')
              constitutive_titanmod_relevantRho(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('interactionslipslip')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interactionSlipSlip(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactionee')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interaction_ee(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactionss')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interaction_ss(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactiones')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interaction_es(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactionsliptwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interactionSlipTwin(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactiontwinslip')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interactionTwinSlip(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('interactiontwintwin')
              forall (j = 1:lattice_maxNinteraction) &
                constitutive_titanmod_interactionTwinTwin(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
     end select
   endif
enddo

write(6,*) 'Material Property reading done'
 
100 do i = 1,maxNinstance
   constitutive_titanmod_structure(i) = &
   lattice_initializeStructure(constitutive_titanmod_structureName(i),constitutive_titanmod_CoverA(i))
   myStructure = constitutive_titanmod_structure(i)

   !* Sanity checks
   if (myStructure < 1 .or. myStructure > 3)                                call IO_error(205)
   if (sum(constitutive_titanmod_Nslip(:,i)) <= 0_pInt)                    call IO_error(207)
   if (sum(constitutive_titanmod_Ntwin(:,i)) < 0_pInt)                     call IO_error(208) !***
   do f = 1,lattice_maxNslipFamily
     if (constitutive_titanmod_Nslip(f,i) > 0_pInt) then   
       if (constitutive_titanmod_rho_edge0(f,i) < 0.0_pReal)                 call IO_error(209)
       if (constitutive_titanmod_rho_screw0(f,i) < 0.0_pReal)              call IO_error(210)
       if (constitutive_titanmod_burgersPerSlipFamily(f,i) <= 0.0_pReal)    call IO_error(211)
       if (constitutive_titanmod_f0_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(212)
       if (constitutive_titanmod_tau0e_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(229)
       if (constitutive_titanmod_tau0s_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(233)
       if (constitutive_titanmod_capre_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(234)
       if (constitutive_titanmod_caprs_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(235)
       if (constitutive_titanmod_v0e_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(226)
       if (constitutive_titanmod_v0s_PerSlipFamily(f,i) <= 0.0_pReal)         call IO_error(226)
       if (constitutive_titanmod_rlengthscrew_PerSlipFamily(f,i) <= 0.0_pReal) call IO_error(238)
     endif
   enddo
   do f = 1,lattice_maxNtwinFamily
     if (constitutive_titanmod_Ntwin(f,i) > 0_pInt) then   
       if (constitutive_titanmod_twinrho_edge0(f,i) < 0.0_pReal)                 call IO_error(209)
       if (constitutive_titanmod_twinrho_screw0(f,i) < 0.0_pReal)              call IO_error(210)
       if (constitutive_titanmod_burgersPerTwinFamily(f,i) <= 0.0_pReal)    call IO_error(221) !***
       if (constitutive_titanmod_Ndot0PerTwinFamily(f,i) < 0.0_pReal)       call IO_error(226) !***
       if (constitutive_titanmod_twinf0_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(228)
       if (constitutive_titanmod_twinshearconstant_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(228)
       if (constitutive_titanmod_twintau0e_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(229)
       if (constitutive_titanmod_twintau0s_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(233)
       if (constitutive_titanmod_twincapre_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(234)
       if (constitutive_titanmod_twincaprs_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(235)
       if (constitutive_titanmod_twinv0e_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(226)
       if (constitutive_titanmod_twinv0s_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(226)
     endif
   enddo
!   if (any(constitutive_titanmod_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 1.0_pReal)) call IO_error(229)
!   if (constitutive_titanmod_CAtomicVolume(i) <= 0.0_pReal)                 call IO_error(230)
   if (constitutive_titanmod_dc(i) <= 0.0_pReal)                            call IO_error(231)
   if (constitutive_titanmod_twinhpconstant(i) <= 0.0_pReal)                call IO_error(232)
   if (constitutive_titanmod_relevantRho(i) <= 0.0_pReal)                   call IO_error(233)
   
   !* Determine total number of active slip or twin systems
   constitutive_titanmod_Nslip(:,i) = min(lattice_NslipSystem(:,myStructure),constitutive_titanmod_Nslip(:,i))
   constitutive_titanmod_Ntwin(:,i) = min(lattice_NtwinSystem(:,myStructure),constitutive_titanmod_Ntwin(:,i))
   constitutive_titanmod_totalNslip(i) = sum(constitutive_titanmod_Nslip(:,i))
   constitutive_titanmod_totalNtwin(i) = sum(constitutive_titanmod_Ntwin(:,i))
  write(6,*) 'Sanity Checks done !'
enddo   
   
!* Allocation of variables whose size depends on the total number of active slip systems
maxTotalNslip = maxval(constitutive_titanmod_totalNslip)
maxTotalNtwin = maxval(constitutive_titanmod_totalNtwin)      
write(6,*) 'maxTotalNslip',maxTotalNslip
write(6,*) 'maxTotalNtwin',maxTotalNtwin
allocate(constitutive_titanmod_burgersPerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_titanmod_burgersPerTwinSystem(maxTotalNtwin, maxNinstance))

allocate(constitutive_titanmod_f0_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_tau0e_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_tau0s_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_capre_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_caprs_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_pe_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_ps_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_qe_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_qs_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_v0e_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_v0s_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_rlengthscrew_PerSlipSystem(maxTotalNslip,maxNinstance))

allocate(constitutive_titanmod_twinf0_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinshearconstant_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twintau0e_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twintau0s_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twincapre_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twincaprs_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinpe_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinps_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinqe_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinqs_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinv0e_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinv0s_PerTwinSystem(maxTotalNTwin,maxNinstance))

allocate(constitutive_titanmod_Ndot0PerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_titanmod_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_titanmod_CeLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_titanmod_CsLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance))

allocate(constitutive_titanmod_twinCeLambdaSlipPerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_titanmod_twinCsLambdaSlipPerTwinSystem(maxTotalNtwin, maxNinstance))

constitutive_titanmod_burgersPerSlipSystem     = 0.0_pReal
constitutive_titanmod_burgersPerTwinSystem     = 0.0_pReal
constitutive_titanmod_f0_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_tau0e_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_tau0s_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_capre_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_caprs_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_v0e_PerSlipSystem          = 0.0_pReal
constitutive_titanmod_v0s_PerSlipSystem          = 0.0_pReal
constitutive_titanmod_rlengthscrew_PerSlipSystem = 0.0_pReal
constitutive_titanmod_pe_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_ps_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_qe_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_qs_PerSlipSystem        = 0.0_pReal

constitutive_titanmod_twinf0_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twinshearconstant_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twintau0e_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twintau0s_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twincapre_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twincaprs_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twinv0e_PerTwinSystem          = 0.0_pReal
constitutive_titanmod_twinv0s_PerTwinSystem          = 0.0_pReal
constitutive_titanmod_twinpe_PerTwinSystem        = 0.0_pReal
constitutive_titanmod_twinps_PerTwinSystem        = 0.0_pReal
constitutive_titanmod_twinqe_PerTwinSystem        = 0.0_pReal
constitutive_titanmod_twinqs_PerTwinSystem        = 0.0_pReal

constitutive_titanmod_Ndot0PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twinsizePerTwinSystem    = 0.0_pReal
constitutive_titanmod_CeLambdaSlipPerSlipSystem = 0.0_pReal
constitutive_titanmod_CsLambdaSlipPerSlipSystem = 0.0_pReal
constitutive_titanmod_twinCeLambdaSlipPerTwinSystem = 0.0_pReal
constitutive_titanmod_twinCsLambdaSlipPerTwinSystem = 0.0_pReal

allocate(constitutive_titanmod_interactionMatrixSlipSlip(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_interactionMatrix_ee(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_interactionMatrix_ss(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_interactionMatrix_es(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_interactionMatrixSlipTwin(maxTotalNslip,maxTotalNtwin,maxNinstance))
allocate(constitutive_titanmod_interactionMatrixTwinSlip(maxTotalNtwin,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_interactionMatrixTwinTwin(maxTotalNtwin,maxTotalNtwin,maxNinstance))
allocate(constitutive_titanmod_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_forestProjectionScrew(maxTotalNslip,maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_TwinforestProjectionEdge(maxTotalNtwin,maxTotalNtwin,maxNinstance))
allocate(constitutive_titanmod_TwinforestProjectionScrew(maxTotalNtwin,maxTotalNtwin,maxNinstance))
constitutive_titanmod_interactionMatrixSlipSlip = 0.0_pReal
constitutive_titanmod_interactionMatrix_ee = 0.0_pReal
constitutive_titanmod_interactionMatrix_ss = 0.0_pReal
constitutive_titanmod_interactionMatrix_es = 0.0_pReal
constitutive_titanmod_interactionMatrixSlipTwin = 0.0_pReal
constitutive_titanmod_interactionMatrixTwinSlip = 0.0_pReal
constitutive_titanmod_interactionMatrixTwinTwin = 0.0_pReal
constitutive_titanmod_forestProjectionEdge      = 0.0_pReal
constitutive_titanmod_forestProjectionScrew      = 0.0_pReal
constitutive_titanmod_TwinforestProjectionEdge      = 0.0_pReal
constitutive_titanmod_TwinforestProjectionScrew      = 0.0_pReal

allocate(constitutive_titanmod_Ctwin_66(6,6,maxTotalNtwin,maxNinstance))
allocate(constitutive_titanmod_Ctwin_3333(3,3,3,3,maxTotalNtwin,maxNinstance))
constitutive_titanmod_Ctwin_66 = 0.0_pReal
constitutive_titanmod_Ctwin_3333 = 0.0_pReal

write(6,*) 'Allocated slip system variables'
do i = 1,maxNinstance 
   myStructure = constitutive_titanmod_structure(i)

   !* Inverse lookup of slip system family
   l = 0_pInt
   do f = 1,lattice_maxNslipFamily
      do k = 1,constitutive_titanmod_Nslip(f,i)
         l = l + 1
         constitutive_titanmod_slipFamily(l,i) = f
         constitutive_titanmod_slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1,myStructure)) + k
   enddo; enddo
   
   !* Inverse lookup of twin system family
   l = 0_pInt
   do f = 1,lattice_maxNtwinFamily
      do k = 1,constitutive_titanmod_Ntwin(f,i)
         l = l + 1
         constitutive_titanmod_twinFamily(l,i) = f
         constitutive_titanmod_twinSystemLattice(l,i) = sum(lattice_NtwinSystem(1:f-1,myStructure)) + k
   enddo; enddo
   
   !* Determine size of state array  
   ns = constitutive_titanmod_totalNslip(i)
   nt = constitutive_titanmod_totalNtwin(i)
   constitutive_titanmod_sizeDotState(i) = &
   size(constitutive_titanmod_listBasicSlipStates)*ns+size(constitutive_titanmod_listBasicTwinStates)*nt
   constitutive_titanmod_sizeState(i) = &
   constitutive_titanmod_sizeDotState(i)+ &
   size(constitutive_titanmod_listDependentSlipStates)*ns+size(constitutive_titanmod_listDependentTwinStates)*nt
  write(6,*) 'Determined size of state and dot state' 
   !* Determine size of postResults array   
   do o = 1,maxval(phase_Noutput)
      select case(constitutive_titanmod_output(o,i))
        case('rhoedge', &
             'rhoscrew', &
             'segment_edge', &
             'segment_screw', &
             'resistance_edge', &
             'resistance_screw', &
             'tau_slip', &
             'gdot_slip', &
             'velocity_edge', &
             'velocity_screw', &
             'total_density' &
             )
           mySize = constitutive_titanmod_totalNslip(i)
        case('twinrhoedge', &
             'twinrhoscrew', &
             'twin_fraction', &
             'gdot_twin', &
             'tau_twin', &
             'twinsegment_edge', &
             'twinsegment_screw', &
             'twinresistance_edge', &
             'twinresistance_screw', &
             'twinvelocity_edge', &
             'twinvelocity_screw', &
             'twintotal_density' &             
             )
           mySize = constitutive_titanmod_totalNtwin(i)
        case default
           mySize = 0_pInt
      end select

       if (mySize > 0_pInt) then  ! any meaningful output found                               
          constitutive_titanmod_sizePostResult(o,i) = mySize
          constitutive_titanmod_sizePostResults(i)  = constitutive_titanmod_sizePostResults(i) + mySize
       endif
   enddo
   
write(6,*) 'Determining elasticity matrix'

   !* Elasticity matrix and shear modulus according to material.config
   select case (myStructure)
   case(1:2) ! cubic(s)
     forall(k=1:3)
       forall(j=1:3) &
         constitutive_titanmod_Cslip_66(k,j,i)     = constitutive_titanmod_C12(i)
         constitutive_titanmod_Cslip_66(k,k,i)     = constitutive_titanmod_C11(i)
         constitutive_titanmod_Cslip_66(k+3,k+3,i) = constitutive_titanmod_C44(i)
     end forall
   case(3:)   ! all hex
     constitutive_titanmod_Cslip_66(1,1,i) = constitutive_titanmod_C11(i)
     constitutive_titanmod_Cslip_66(2,2,i) = constitutive_titanmod_C11(i)
     constitutive_titanmod_Cslip_66(3,3,i) = constitutive_titanmod_C33(i)
     constitutive_titanmod_Cslip_66(1,2,i) = constitutive_titanmod_C12(i)
     constitutive_titanmod_Cslip_66(2,1,i) = constitutive_titanmod_C12(i)
     constitutive_titanmod_Cslip_66(1,3,i) = constitutive_titanmod_C13(i)
     constitutive_titanmod_Cslip_66(3,1,i) = constitutive_titanmod_C13(i)
     constitutive_titanmod_Cslip_66(2,3,i) = constitutive_titanmod_C13(i)
     constitutive_titanmod_Cslip_66(3,2,i) = constitutive_titanmod_C13(i)
     constitutive_titanmod_Cslip_66(4,4,i) = constitutive_titanmod_C44(i)
     constitutive_titanmod_Cslip_66(5,5,i) = constitutive_titanmod_C44(i)
     constitutive_titanmod_Cslip_66(6,6,i) = 0.5_pReal*(constitutive_titanmod_C11(i)-constitutive_titanmod_C12(i))
   end select
   constitutive_titanmod_Cslip_66(:,:,i) = math_Mandel3333to66(math_Voigt66to3333(constitutive_titanmod_Cslip_66(:,:,i)))
   constitutive_titanmod_Cslip_3333(:,:,:,:,i) = math_Voigt66to3333(constitutive_titanmod_Cslip_66(:,:,i))
   constitutive_titanmod_Gmod(i) = &
   0.2_pReal*(constitutive_titanmod_C11(i)-constitutive_titanmod_C12(i))+0.3_pReal*constitutive_titanmod_C44(i)
   
   !* Construction of the twin elasticity matrices
   do j=1,lattice_maxNtwinFamily
      do k=1,constitutive_titanmod_Ntwin(j,i)           
         do l=1,3 ; do m=1,3 ; do n=1,3 ; do o=1,3 ; do p=1,3 ; do q=1,3 ; do r=1,3 ; do s=1,3
           constitutive_titanmod_Ctwin_3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1,i))+k,i) = &
             constitutive_titanmod_Ctwin_3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1,i))+k,i) + &
             constitutive_titanmod_Cslip_3333(p,q,r,s,i)*&
             lattice_Qtwin(l,p,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(m,q,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(n,r,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)* &
             lattice_Qtwin(o,s,sum(lattice_NslipSystem(1:j-1,myStructure))+k,myStructure)
           enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo ; enddo
         constitutive_titanmod_Ctwin_66(:,:,k,i) = math_Mandel3333to66(constitutive_titanmod_Ctwin_3333(:,:,:,:,k,i))
        enddo
   enddo

   !* Burgers vector, dislocation velocity prefactor for each slip system 
   do s = 1,constitutive_titanmod_totalNslip(i)   
      f = constitutive_titanmod_slipFamily(s,i)    
      constitutive_titanmod_burgersPerSlipSystem(s,i)     = constitutive_titanmod_burgersPerSlipFamily(f,i)
      constitutive_titanmod_f0_PerSlipSystem(s,i)       = constitutive_titanmod_f0_PerSlipFamily(f,i)
      constitutive_titanmod_tau0e_PerSlipSystem(s,i)       = constitutive_titanmod_tau0e_PerSlipFamily(f,i)
      constitutive_titanmod_tau0s_PerSlipSystem(s,i)       = constitutive_titanmod_tau0s_PerSlipFamily(f,i)
      constitutive_titanmod_capre_PerSlipSystem(s,i)       = constitutive_titanmod_capre_PerSlipFamily(f,i)
      constitutive_titanmod_caprs_PerSlipSystem(s,i)       = constitutive_titanmod_caprs_PerSlipFamily(f,i)
      constitutive_titanmod_v0e_PerSlipSystem(s,i)          = constitutive_titanmod_v0e_PerSlipFamily(f,i)
      constitutive_titanmod_v0s_PerSlipSystem(s,i)          = constitutive_titanmod_v0s_PerSlipFamily(f,i)
      constitutive_titanmod_rlengthscrew_PerSlipSystem(s,i) = constitutive_titanmod_rlengthscrew_PerSlipFamily(f,i)
      constitutive_titanmod_pe_PerSlipSystem(s,i)          = constitutive_titanmod_pe_PerSlipFamily(f,i)
      constitutive_titanmod_ps_PerSlipSystem(s,i)          = constitutive_titanmod_ps_PerSlipFamily(f,i)
      constitutive_titanmod_qe_PerSlipSystem(s,i)          = constitutive_titanmod_qe_PerSlipFamily(f,i)
      constitutive_titanmod_qs_PerSlipSystem(s,i)          = constitutive_titanmod_qs_PerSlipFamily(f,i)
      constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,i) = constitutive_titanmod_CeLambdaSlipPerSlipFamily(f,i)
      constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,i) = constitutive_titanmod_CsLambdaSlipPerSlipFamily(f,i)
   enddo   
   
   !* Burgers vector, nucleation rate prefactor and twin size for each twin system 
   do s = 1,constitutive_titanmod_totalNtwin(i)   
      f = constitutive_titanmod_twinFamily(s,i)    
      constitutive_titanmod_burgersPerTwinSystem(s,i)  = constitutive_titanmod_burgersPerTwinFamily(f,i)
      constitutive_titanmod_Ndot0PerTwinSystem(s,i)    = constitutive_titanmod_Ndot0PerTwinFamily(f,i)
      constitutive_titanmod_twinsizePerTwinSystem(s,i) = constitutive_titanmod_twinsizePerTwinFamily(f,i)

      constitutive_titanmod_twinf0_PerTwinSystem(s,i)       = constitutive_titanmod_twinf0_PerTwinFamily(f,i)
      constitutive_titanmod_twinshearconstant_PerTwinSystem(s,i) = constitutive_titanmod_twinshearconstant_PerTwinFamily(f,i)
      constitutive_titanmod_twintau0e_PerTwinSystem(s,i)       = constitutive_titanmod_twintau0e_PerTwinFamily(f,i)
      constitutive_titanmod_twintau0s_PerTwinSystem(s,i)       = constitutive_titanmod_twintau0s_PerTwinFamily(f,i)
      constitutive_titanmod_twincapre_PerTwinSystem(s,i)       = constitutive_titanmod_twincapre_PerTwinFamily(f,i)
      constitutive_titanmod_twincaprs_PerTwinSystem(s,i)       = constitutive_titanmod_twincaprs_PerTwinFamily(f,i)
      constitutive_titanmod_twinv0e_PerTwinSystem(s,i)          = constitutive_titanmod_twinv0e_PerTwinFamily(f,i)
      constitutive_titanmod_twinv0s_PerTwinSystem(s,i)          = constitutive_titanmod_twinv0s_PerTwinFamily(f,i)
      constitutive_titanmod_twinpe_PerTwinSystem(s,i)          = constitutive_titanmod_twinpe_PerTwinFamily(f,i)
      constitutive_titanmod_twinps_PerTwinSystem(s,i)          = constitutive_titanmod_twinps_PerTwinFamily(f,i)
      constitutive_titanmod_twinqe_PerTwinSystem(s,i)          = constitutive_titanmod_twinqe_PerTwinFamily(f,i)
      constitutive_titanmod_twinqs_PerTwinSystem(s,i)          = constitutive_titanmod_twinqs_PerTwinFamily(f,i)
      constitutive_titanmod_twinCeLambdaSlipPerTwinSystem(s,i) = constitutive_titanmod_twinCeLambdaSlipPerTwinFamily(f,i)
      constitutive_titanmod_twinCsLambdaSlipPerTwinSystem(s,i) = constitutive_titanmod_twinCsLambdaSlipPerTwinFamily(f,i)

    enddo   
     
   !* Construction of interaction matrices
   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrixSlipSlip(s1,s2,i) = &
         constitutive_titanmod_interactionSlipSlip(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                                constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrix_ee(s1,s2,i) = &
         constitutive_titanmod_interaction_ee(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                                constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrix_ss(s1,s2,i) = &
         constitutive_titanmod_interaction_ss(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                                constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo

   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrix_es(s1,s2,i) = &
         constitutive_titanmod_interaction_es(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                                constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                                myStructure),i)
   enddo; enddo
   
   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do t2 = 1,constitutive_titanmod_totalNtwin(i)     
         constitutive_titanmod_interactionMatrixSlipTwin(s1,t2,i) = &
         constitutive_titanmod_interactionSlipTwin(lattice_interactionSlipTwin(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                                constitutive_titanmod_twinSystemLattice(t2,i), &
                                                                                myStructure),i)         
   enddo; enddo
   
   do t1 = 1,constitutive_titanmod_totalNtwin(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrixTwinSlip(t1,s2,i) = &
         constitutive_titanmod_interactionTwinSlip(lattice_interactionTwinSlip(constitutive_titanmod_twinSystemLattice(t1,i), &
                                                                                constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                                myStructure),i)         
   enddo; enddo

   do t1 = 1,constitutive_titanmod_totalNtwin(i)
      do t2 = 1,constitutive_titanmod_totalNtwin(i)     
         constitutive_titanmod_interactionMatrixTwinTwin(t1,t2,i) = &
         constitutive_titanmod_interactionTwinTwin(lattice_interactionTwinTwin(constitutive_titanmod_twinSystemLattice(t1,i), &
                                                                                constitutive_titanmod_twinSystemLattice(t2,i), &
                                                                                myStructure),i)         
   enddo; enddo
   
   !* Calculation of forest projections for edge dislocations 
   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)      
         constitutive_titanmod_forestProjectionEdge(s1,s2,i) = &
         abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,i),myStructure), &
                         lattice_st(:,constitutive_titanmod_slipSystemLattice(s2,i),myStructure))) 
   enddo; enddo

   !* Calculation of forest projections for screw dislocations 
   do s1 = 1,constitutive_titanmod_totalNslip(i)
      do s2 = 1,constitutive_titanmod_totalNslip(i)      
         constitutive_titanmod_forestProjectionScrew(s1,s2,i) = &
         abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,i),myStructure), &
                         lattice_sd(:,constitutive_titanmod_slipSystemLattice(s2,i),myStructure))) 
   enddo; enddo
  
! Same framework for twins as for slip.

   !* Calculation of forest projections for edge dislocations in twin system
   do t1 = 1,constitutive_titanmod_totalNtwin(i)
      do t2 = 1,constitutive_titanmod_totalNtwin(i)      
         constitutive_titanmod_TwinforestProjectionEdge(t1,t2,i) = &
         abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,i),myStructure), &
                         lattice_tt(:,constitutive_titanmod_twinSystemLattice(t2,i),myStructure))) 
   enddo; enddo

   !* Calculation of forest projections for screw dislocations in twin system
   do t1 = 1,constitutive_titanmod_totalNtwin(i)
      do t2 = 1,constitutive_titanmod_totalNtwin(i)      
         constitutive_titanmod_TwinforestProjectionScrew(t1,t2,i) = &
         abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,i),myStructure), &
                         lattice_td(:,constitutive_titanmod_twinSystemLattice(t2,i),myStructure))) 
   enddo; enddo

   
  enddo
write(6,*) 'Init All done'
return
end subroutine


function constitutive_titanmod_stateInit(myInstance)
!*********************************************************************
!* initial microstructural state                                     *
!*********************************************************************
use prec,    only: pReal,pInt
use math,    only: pi
use lattice, only: lattice_maxNslipFamily,lattice_maxNtwinFamily
implicit none

!* Input-Output variables
integer(pInt) :: myInstance
real(pReal), dimension(constitutive_titanmod_sizeState(myInstance))  :: constitutive_titanmod_stateInit
!* Local variables
integer(pInt) s0,s1,s,t,f,ns,nt,ts0,ts1,tf,ts
real(pReal), dimension(constitutive_titanmod_totalNslip(myInstance)) ::  rho_edge0, &
                                                                         rho_screw0, &
                                                                         segment_edge0, &
                                                                         segment_screw0, &
                                                                         resistance_edge0, &
                                                                         resistance_screw0
real(pReal), dimension(constitutive_titanmod_totalNtwin(myInstance)) ::  twinrho_edge0, &
                                                                         twinrho_screw0, &
                                                                         twinsegment_edge0, &
                                                                         twinsegment_screw0, &
                                                                         twinresistance_edge0, &
                                                                         twinresistance_screw0, &
                                                                         twingamma_dot0

ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)
constitutive_titanmod_stateInit = 0.0_pReal

!* Initialize basic slip state variables
! For slip
s1 = 0_pInt
do f = 1,lattice_maxNslipFamily
   s0 = s1 + 1_pInt
   s1 = s0 + constitutive_titanmod_Nslip(f,myInstance) - 1_pInt 
   do s = s0,s1
      rho_edge0(s)    = constitutive_titanmod_rho_edge0(f,myInstance)
      rho_screw0(s) = constitutive_titanmod_rho_screw0(f,myInstance)
   enddo 
enddo

!* Initialize basic slip state variables
! For twin
ts1 = 0_pInt
do tf = 1,lattice_maxNtwinFamily
   ts0 = ts1 + 1_pInt
   ts1 = ts0 + constitutive_titanmod_Ntwin(tf,myInstance) - 1_pInt 
   do ts = ts0,ts1
      twinrho_edge0(ts)    = constitutive_titanmod_twinrho_edge0(tf,myInstance)
      twinrho_screw0(ts) = constitutive_titanmod_twinrho_screw0(tf,myInstance)
      twingamma_dot0(ts)=0.0_pReal
   enddo 
enddo

constitutive_titanmod_stateInit(1:ns)      = rho_edge0
constitutive_titanmod_stateInit(ns+1:2*ns) = rho_screw0
constitutive_titanmod_stateInit(2*ns+1:2*ns+nt)      = twinrho_edge0
constitutive_titanmod_stateInit(2*ns+nt+1:2*ns+2*nt) = twinrho_screw0
constitutive_titanmod_stateInit(2*ns+2*nt+1:2*ns+3*nt)=twingamma_dot0

!* Initialize dependent slip microstructural variables
forall (s = 1:ns) &
segment_edge0(s) = constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((rho_edge0+rho_screw0),constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance)))
 
constitutive_titanmod_stateInit(2*ns+2*nt+1:3*ns+3*nt) = segment_edge0

forall (s = 1:ns) &
segment_screw0(s) = constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((rho_edge0+rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
  
constitutive_titanmod_stateInit(3*ns+2*nt+1:4*ns+3*nt) = segment_screw0

forall (s = 1:ns) &
resistance_edge0(s) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)* &
sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+dot_product((rho_screw0), &
        constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))
constitutive_titanmod_stateInit(4*ns+2*nt+1:5*ns+3*nt) = resistance_edge0

forall (s = 1:ns) &
resistance_screw0(s) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)* &
sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+dot_product((rho_screw0), &
constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))
constitutive_titanmod_stateInit(5*ns+2*nt+1:6*ns+3*nt) = resistance_screw0

!* Initialize dependent twin microstructural variables
forall (t = 1:nt) &
twinsegment_edge0(t) = constitutive_titanmod_twinCeLambdaSlipPertwinSystem(t,myInstance)/ &
        sqrt(dot_product((twinrho_edge0+twinrho_screw0),constitutive_titanmod_twinforestProjectionEdge(1:nt,t,myInstance)))
 
constitutive_titanmod_stateInit(6*ns+2*nt+1:6*ns+4*nt) = twinsegment_edge0

forall (t = 1:nt) &
twinsegment_screw0(t) = constitutive_titanmod_twinCsLambdaSlipPertwinSystem(t,myInstance)/ &
        sqrt(dot_product((twinrho_edge0+twinrho_screw0),constitutive_titanmod_twinforestProjectionScrew(1:nt,t,myInstance)))
 
constitutive_titanmod_stateInit(6*ns+3*nt+1:6*ns+5*nt) = twinsegment_edge0

forall (t = 1:nt) &
twinresistance_edge0(t) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSystem(t,myInstance)* &
sqrt(dot_product((twinrho_edge0),constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance))+ &
        dot_product((twinrho_screw0),constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance)))

constitutive_titanmod_stateInit(6*ns+4*nt+1:6*ns+6*nt) = twinresistance_edge0

forall (t = 1:nt) &
twinresistance_screw0(t) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSystem(t,myInstance)* &
sqrt(dot_product((twinrho_edge0),constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance))+ &
        dot_product((twinrho_screw0),constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance)))

constitutive_titanmod_stateInit(6*ns+5*nt+1:6*ns+7*nt) = twinresistance_screw0

!forall (t = 1:nt) &
!MeanFreePathTwin0(t) = constitutive_titanmod_GrainSize(myInstance)
!constitutive_titanmod_stateInit(5*ns+2*nt+1:5*ns+3*nt) = MeanFreePathTwin0

!forall (t = 1:nt) &
!TwinVolume0(t) = & 
!(pi/6.0_pReal)*constitutive_titanmod_twinsizePerTwinSystem(t,myInstance)*MeanFreePathTwin0(t)**(2.0_pReal)
!constitutive_titanmod_stateInit(6*ns+4*nt+1:6*ns+5*nt) = TwinVolume0

!write(6,*) '#STATEINIT#'
!write(6,*)
!write(6,'(a,/,4(3(f30.20,x)/))') 'rho_edge',rho_edge0
!write(6,'(a,/,4(3(f30.20,x)/))') 'rho_screw',rho_screw0
!write(6,'(a,/,4(3(f30.20,x)/))') 'segment_edge',segment_edge0
!write(6,'(a,/,4(3(f30.20,x)/))') 'segment_screw',segment_screw0
!write(6,'(a,/,4(3(f30.20,x)/))') 'tauSlipThreshold', tauSlipThreshold0
!write(6,'(a,/,4(3(f30.20,x)/))') 'MeanFreePathTwin', MeanFreePathTwin0
!write(6,'(a,/,4(3(f30.20,x)/))') 'TwinVolume', TwinVolume0

return
end function

pure function constitutive_titanmod_relevantState(myInstance)
!*********************************************************************
!* relevant microstructural state                                    *
!*********************************************************************
use prec,     only: pReal, pInt
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: myInstance
real(pReal), dimension(constitutive_titanmod_sizeState(myInstance)) :: constitutive_titanmod_relevantState

constitutive_titanmod_relevantState = constitutive_titanmod_relevantRho(myInstance)

return
endfunction

pure function constitutive_titanmod_homogenizedC(state,g,ip,el)
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
real(pReal), dimension(6,6) :: constitutive_titanmod_homogenizedC
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   volumefraction_pertwinsystem
!* Local variables 
integer(pInt) myInstance,ns,nt,i
real(pReal) sumf
 
!* Shortened notation
myInstance = phase_constitutionInstance(material_phase(g,ip,el))
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

!* Total twin volume fraction
do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(2*ns+2*nt+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)
enddo
!sumf = sum(state(g,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0
sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0

!* Homogenized elasticity matrix
constitutive_titanmod_homogenizedC = (1.0_pReal-sumf)*constitutive_titanmod_Cslip_66(:,:,myInstance)
do i=1,nt
   constitutive_titanmod_homogenizedC = &
!   constitutive_titanmod_homogenizedC + state(g,ip,el)%p(6*ns+7*nt+i)*constitutive_titanmod_Ctwin_66(:,:,i,myInstance)
   constitutive_titanmod_homogenizedC + volumefraction_pertwinsystem(i)*constitutive_titanmod_Ctwin_66(:,:,i,myInstance)

enddo 

return
end function


subroutine constitutive_titanmod_microstructure(Temperature,state,g,ip,el)
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
integer(pInt) myInstance,myStructure,ns,nt,s,t,i
real(pReal) sumf,sfe
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
            fOverStacksize, volumefraction_pertwinsystem
 
!* Shortened notation
myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_titanmod_structure(myInstance)
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

! Need to update this list
!* State: 1           :  ns         rho_edge
!* State: ns+1        :  2*ns       rho_screw
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
do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(2*ns+2*nt+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

!sumf = sum(state(g,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0

!* Stacking fault energy
sfe = 0.0002_pReal*Temperature-0.0396_pReal

!* rescaled twin volume fraction for topology
!forall (t = 1:nt) &
!  fOverStacksize(t) = &
!    state(g,ip,el)%p(2*ns+t)/constitutive_titanmod_twinsizePerTwinSystem(t,myInstance)
        
! average segment length for edge dislocations in matrix
forall (s = 1:ns) &
  state(g,ip,el)%p(2*ns+3*nt+s) = constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)), &
        constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance)))
     
! average segment length for edge dislocations in matrix
forall (s = 1:ns) &
  state(g,ip,el)%p(3*ns+3*nt+s) = constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)), &
        constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance)))

!* Average segment length for screw dislocations in matrix
!do s = 1,ns
!   if (nt > 0_pInt) then
!      state(g,ip,el)%p(4*ns+2*nt+s) = &
!                constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance) / &
!        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns:2*ns)), &
!        constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
!   else
!      state(g,ip,el)%p(4*ns+s) = &
!                constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance) / &
!        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns:2*ns)), &
!        constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
!   endif
!enddo

!* threshold stress or slip resistance for edge dislocation motion
forall (s = 1:ns) &
  state(g,ip,el)%p(4*ns+3*nt+s) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(1:ns)),&
                         constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+ &
                      dot_product((state(g,ip,el)%p(ns+1:2*ns)),&
                         constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))

!* threshold stress or slip resistance for screw dislocation motion
forall (s = 1:ns) &
  state(g,ip,el)%p(5*ns+3*nt+s) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(1:ns)),&
                         constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+ &
                        dot_product((state(g,ip,el)%p(ns+1:2*ns)),&
                         constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))

! average segment length for edge dislocations in twin
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+3*nt+t) = constitutive_titanmod_twinCeLambdaSlipPerTwinSystem(t,myInstance)/ &
        sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)+state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)), &
        constitutive_titanmod_twinforestProjectionEdge(1:nt,t,myInstance)))
     
! average segment length for screw dislocations in twin
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+4*nt+t) = constitutive_titanmod_twinCeLambdaSlipPerTwinSystem(t,myInstance)/ &
        sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)+state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)), &
        constitutive_titanmod_twinforestProjectionScrew(1:nt,t,myInstance)))

!* threshold stress or slip resistance for edge dislocation motion in twin
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+5*nt+t) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSystem(t,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)),&
                         constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance))+ &
                      dot_product((state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)),&
                         constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance)))

!* threshold stress or slip resistance for screw dislocation motion in twin
forall (t = 1:nt) &
  state(g,ip,el)%p(6*ns+6*nt+t) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSystem(t,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)),&
                         constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance))+ &
                      dot_product((state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)),&
                         constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance)))


return
end subroutine


subroutine constitutive_titanmod_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,g,ip,el)
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
real(pReal) sumf,StressRatio_edge_p,minusStressRatio_edge_p,StressRatio_edge_pminus1,StressRatio_screw_p, &
        StressRatio_screw_pminus1, StressRatio_r,BoltzmannRatio,DotGamma0, minusStressRatio_screw_p,gdotTotal, &
        screwvelocity_kink_prefactor,twinStressRatio_edge_p,twinminusStressRatio_edge_p,twinStressRatio_edge_pminus1, &
   twinStressRatio_screw_p, twinStressRatio_screw_pminus1, twinStressRatio_r, twinDotGamma0, &
   twinminusStressRatio_screw_p
real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
real(pReal), dimension(constitutive_titanmod_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip, edge_velocity, screw_velocity
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin, twinedge_velocity, twinscrew_velocity,volumefraction_pertwinsystem

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(2*ns+2*nt+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

!sumf = sum(state(g,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0


Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal
dLp_dTstar = 0.0_pReal

!* Dislocation glide part
gdot_slip = 0.0_pReal
dgdot_dtauslip = 0.0_pReal
j = 0_pInt
do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
   do i = 1,constitutive_titanmod_Nslip(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure)) 
!                state(g,ip,el)%p(9*ns+3*nt+j)=tau_slip(j)
!*************************************************
 
      if(myStructure==3.and.j>3) then ! only for hex and for all the non-basal slip systems
      screwvelocity_kink_prefactor=state(g,ip,el)%p(3*ns+3*nt+j)/constitutive_titanmod_rlengthscrew_PerSlipSystem(j,myInstance)
        else
        screwvelocity_kink_prefactor=1.0_pReal
        endif
  
!       state(g,ip,el)%p(14*ns+3*nt+j)=screwvelocity_kink_prefactor
   
     !* Stress ratio for edge
         StressRatio_edge_p = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0e_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(4*ns+3*nt+j)) &
        )**constitutive_titanmod_pe_PerSlipSystem(j,myInstance)
        
     !* Stress ratio for screw
         StressRatio_screw_p = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(5*ns+3*nt+j)) &
        )**constitutive_titanmod_ps_PerSlipSystem(j,myInstance)
                
!        state(g,ip,el)%p(10*ns+3*nt+j)=StressRatio_edge_p
!        state(g,ip,el)%p(11*ns+3*nt+j)=StressRatio_screw_p
        
        if((1.0_pReal-StressRatio_edge_p)>0.001_pReal) then
        minusStressRatio_edge_p=1.0_pReal-StressRatio_edge_p
        else
        minusStressRatio_edge_p=0.001_pReal
        endif
        
        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif
        
      StressRatio_edge_pminus1 = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0e_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(4*ns+3*nt+j)) &
        )**(constitutive_titanmod_pe_PerSlipSystem(j,myInstance)-1.0_pReal)

      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(5*ns+3*nt+j)) &
        )**(constitutive_titanmod_ps_PerSlipSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_f0_PerSlipSystem(j,myInstance)/(kB*Temperature)

      !* Initial shear rates
      DotGamma0 = &
        constitutive_titanmod_burgersPerSlipSystem(j,myInstance)*(state(g,ip,el)%p(j)*&
        + constitutive_titanmod_v0e_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(ns+j)* &
                constitutive_titanmod_v0s_PerSlipSystem(j,myInstance))

         edge_velocity(j) =constitutive_titanmod_v0e_PerSlipSystem(j,myInstance)*exp(-BoltzmannRatio* &
        (minusStressRatio_edge_p)** &
        constitutive_titanmod_qe_PerSlipSystem(j,myInstance))

        screw_velocity(j) =screwvelocity_kink_prefactor * constitutive_titanmod_v0s_PerSlipSystem(j,myInstance)* &
                exp(-BoltzmannRatio*(minusStressRatio_screw_p)** &
        constitutive_titanmod_qs_PerSlipSystem(j,myInstance))

                !* Shear rates due to slip
       gdot_slip(j) = constitutive_titanmod_burgersPerSlipSystem(j,myInstance)*(state(g,ip,el)%p(j)* &
                edge_velocity(j)+state(g,ip,el)%p(ns+j) * screw_velocity(j))* sign(1.0_pReal,tau_slip(j))
                
      !* Derivatives of shear rates
      dgdot_dtauslip(j) = ( &
                ( &
                ( &
                ( &
                (abs(gdot_slip(j))) * &
                BoltzmannRatio*&
        constitutive_titanmod_pe_PerSlipSystem(j,myInstance)* &
                constitutive_titanmod_qe_PerSlipSystem(j,myInstance) &
                )/ &
                constitutive_titanmod_tau0e_PerSlipSystem(j,myInstance) &
                )*&
        StressRatio_edge_pminus1*(minusStressRatio_edge_p)** &
                (constitutive_titanmod_qe_PerSlipSystem(j,myInstance)-1.0_pReal) &
                ) + &
                ( &
                ( &
                ( &
                (abs(gdot_slip(j))) * &
                BoltzmannRatio* screwvelocity_kink_prefactor *&
        constitutive_titanmod_ps_PerSlipSystem(j,myInstance)* &
                constitutive_titanmod_qs_PerSlipSystem(j,myInstance) &
                )/ &
                constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance) &
                )*&
        StressRatio_screw_pminus1*(minusStressRatio_screw_p)**(constitutive_titanmod_qs_PerSlipSystem(j,myInstance)-1.0_pReal) &
                ) &
                ) !* sign(1.0_pReal,tau_slip(j))
                
!                state(g,ip,el)%p(13*ns+3*nt+j)=dgdot_dtauslip(j)
                
!*************************************************                
!sumf=0.0_pReal
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
   do i = 1,constitutive_titanmod_Ntwin(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))        
     
!**************************************************************************************
      !* Stress ratios
!      StressRatio_r = (state(g,ip,el)%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_titanmod_r(myInstance)      
      
          !* Shear rates and their derivatives due to twin
!      if ( tau_twin(j) > 0.0_pReal ) !then          
!        gdot_twin(j) =  0.0_pReal!&
!          (constitutive_titanmod_MaxTwinFraction(myInstance)-sumf)*lattice_shearTwin(index_myFamily+i,myStructure)*&
!          state(g,ip,el)%p(6*ns+4*nt+j)*constitutive_titanmod_Ndot0PerTwinSystem(f,myInstance)*exp(-StressRatio_r) 
!        dgdot_dtautwin(j) = ((gdot_twin(j)*constitutive_titanmod_r(myInstance))/tau_twin(j))*StressRatio_r
!      endif
!**************************************************************************************
   
     !* Stress ratio for edge
         twinStressRatio_edge_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0e_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+5*nt+j)) &
        )**constitutive_titanmod_twinpe_PerTwinSystem(j,myInstance)
        
     !* Stress ratio for screw
         twinStressRatio_screw_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0s_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+6*nt+j)) &
        )**constitutive_titanmod_twinps_PerTwinSystem(j,myInstance)
                
!        state(g,ip,el)%p(10*ns+3*nt+j)=twinStressRatio_edge_p
!        state(g,ip,el)%p(11*ns+3*nt+j)=twinStressRatio_screw_p
        
        if((1.0_pReal-twinStressRatio_edge_p)>0.001_pReal) then
        twinminusStressRatio_edge_p=1.0_pReal-twinStressRatio_edge_p
        else
        twinminusStressRatio_edge_p=0.001_pReal
        endif
        
        if((1.0_pReal-twinStressRatio_screw_p)>0.001_pReal) then
        twinminusStressRatio_screw_p=1.0_pReal-twinStressRatio_screw_p
        else
        twinminusStressRatio_screw_p=0.001_pReal
        endif
        
      twinStressRatio_edge_pminus1 = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0e_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+5*nt+j)) &
        )**(constitutive_titanmod_twinpe_PerTwinSystem(j,myInstance)-1.0_pReal)

      twinStressRatio_screw_pminus1 = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0s_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+6*nt+j)) &
        )**(constitutive_titanmod_twinps_PerTwinSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_twinf0_PerTwinSystem(j,myInstance)/(kB*Temperature)

      !* Initial shear rates
      TwinDotGamma0 = &
        constitutive_titanmod_burgersPerTwinSystem(j,myInstance)*(state(g,ip,el)%p(2*ns+j)*&
         constitutive_titanmod_twinv0e_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(2*ns+nt+j)* &
                constitutive_titanmod_twinv0s_PerTwinSystem(j,myInstance))

         twinedge_velocity(j) =constitutive_titanmod_twinv0e_PerTwinSystem(j,myInstance)*exp(-BoltzmannRatio* &
        (twinminusStressRatio_edge_p)** &
        constitutive_titanmod_twinqe_PerTwinSystem(j,myInstance))

        twinscrew_velocity(j) =constitutive_titanmod_twinv0s_PerTwinSystem(j,myInstance)* &
                exp(-BoltzmannRatio*(twinminusStressRatio_screw_p)** &
        constitutive_titanmod_twinqs_PerTwinSystem(j,myInstance))
        
        state(g,ip,el)%p(6*ns+8*nt+j)=twinscrew_velocity(j)
        
        
                !* Shear rates due to twin
       gdot_twin(j) = constitutive_titanmod_burgersPerTwinSystem(j,myInstance)*(state(g,ip,el)%p(2*ns+j)* &
                twinedge_velocity(j)+state(g,ip,el)%p(2*ns+nt+j) * twinscrew_velocity(j))* sign(1.0_pReal,tau_twin(j))

!                forall (s = 1:ns) &
!                  state(g,ip,el)%p(7*ns+3*nt+j)= twinedge_velocity(j)
!        forall (s = 1:ns) &
!                state(g,ip,el)%p(8*ns+3*nt+j)= twinscrew_velocity(j)
!                state(g,ip,el)%p(12*ns+3*nt+j)=gdot_twin(j)
                
      !* Derivatives of shear rates in twin
      dgdot_dtautwin(j) = ( &
                ( &
                ( &
                ( &
                (abs(gdot_twin(j))) * &
                BoltzmannRatio*&
        constitutive_titanmod_twinpe_PerTwinSystem(j,myInstance)* &
                constitutive_titanmod_twinqe_PerTwinSystem(j,myInstance) &
                )/ &
                constitutive_titanmod_twintau0e_PerTwinSystem(j,myInstance) &
                )*&
        twinStressRatio_edge_pminus1*(twinminusStressRatio_edge_p)** &
                (constitutive_titanmod_twinqe_PerTwinSystem(j,myInstance)-1.0_pReal) &
                ) + &
                ( &
                ( &
                ( &
                (abs(gdot_twin(j))) * &
                BoltzmannRatio* &
        constitutive_titanmod_twinps_PerTwinSystem(j,myInstance)* &
                constitutive_titanmod_twinqs_PerTwinSystem(j,myInstance) &
                )/ &
                constitutive_titanmod_twintau0s_PerTwinSystem(j,myInstance) &
                )*&
        twinStressRatio_screw_pminus1*(twinminusStressRatio_screw_p)** &
               (constitutive_titanmod_twinqs_PerTwinSystem(j,myInstance)-1.0_pReal) &
                ) &
                ) !* sign(1.0_pReal,tau_slip(j))
                
!                state(g,ip,el)%p(13*ns+3*nt+j)=dgdot_dtautwin(j)

      !* Plastic velocity gradient for mechanical twinning                                                      
!      Lp = Lp + sumf*gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,myStructure)
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


function constitutive_titanmod_dotState(Tstar_v,Temperature,state,g,ip,el)
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
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase, phase_constitutionInstance
use lattice,  only: lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v,lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                     lattice_NslipSystem,lattice_NtwinSystem,lattice_shearTwin   
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: g,ip,el
real(pReal), intent(in) :: Temperature
real(pReal), dimension(6), intent(in) :: Tstar_v
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
real(pReal), dimension(constitutive_titanmod_sizeDotState(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
constitutive_titanmod_dotState
!* Local variables
integer(pInt) MyInstance,MyStructure,ns,nt,f,i,j,k,index_myFamily,s,t
real(pReal) sumf,StressRatio_edge_p,minusStressRatio_edge_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,&
            EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r,StressRatio_screw_p,minusStressRatio_screw_p, &
            screwvelocity_kink_prefactor,twinStressRatio_edge_p,twinminusStressRatio_edge_p,twinStressRatio_pminus1, &
            twinDotGamma0,twinStressRatio_screw_p, &
            twinminusStressRatio_screw_p
real(pReal), dimension(constitutive_titanmod_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
gdot_slip,tau_slip,DotRhoEdgeGeneration,EdgeDipDistance,DotRhoEdgeAnnihilation,DotRhoScrewAnnihilation,&
ClimbVelocity,DotRhoScrewGeneration, edge_segment, screw_segment,edge_velocity,screw_velocity
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: gdot_twin, &
tau_twin,twinedge_segment,twinscrew_segment,twinedge_velocity,twinscrew_velocity,TwinDotRhoEdgeGeneration, &
TwinDotRhoEdgeAnnihilation,TwinDotRhoScrewGeneration,TwinDotRhoScrewAnnihilation,volumefraction_pertwinsystem
   
!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
MyStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(2*ns+2*nt+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

!sumf = sum(state(g,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0

constitutive_titanmod_dotState = 0.0_pReal

!* average segment length for edge dislocations in matrix
forall (s = 1:ns) &
  edge_segment(s) = &
  (constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance))/ &
        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)), &
                constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance)))

!* average segment length for screw dislocations in matrix
forall (s = 1:ns) &
  screw_segment(s) = &
  (constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance))/ &
        sqrt(dot_product((state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)), &
                constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
    
!* average segment length for edge dislocations in twin
forall (t = 1:nt) &
  twinedge_segment(t) = &
  (constitutive_titanmod_twinCeLambdaSlipPerTwinSystem(t,myInstance))/ &
        sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)+state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)), &
                constitutive_titanmod_TwinforestProjectionEdge(1:nt,t,myInstance)))

!* average segment length for screw dislocations in twin
forall (t = 1:nt) &
  twinscrew_segment(t) = &
  (constitutive_titanmod_twinCsLambdaSlipPerTwinSystem(t,myInstance))/ &
        sqrt(dot_product((state(g,ip,el)%p(2*ns+1:2*ns+nt)+state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)), &
                constitutive_titanmod_TwinforestProjectionScrew(1:nt,t,myInstance)))


    j = 0_pInt
 do f = 1,lattice_maxNslipFamily                                             ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure))                 ! at which index starts my family
   do i = 1,constitutive_titanmod_Nslip(f,myInstance)                        ! process each (active) slip system in family
     j = j+1_pInt

! Resolved shear stress
     tau_slip(j)  = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure)) 


        if(myStructure==3.and.j>3) then ! only for hex and for all the non-basal slip systems
        screwvelocity_kink_prefactor=state(g,ip,el)%p(3*ns+2*nt+j)/constitutive_titanmod_rlengthscrew_PerSlipSystem(j,myInstance)
        else
        screwvelocity_kink_prefactor=1.0_pReal
        endif

     !* Stress ratio for edge
         StressRatio_edge_p = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0e_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(4*ns+3*nt+j)) &
        )**(constitutive_titanmod_pe_PerSlipSystem(j,myInstance))
        
     !* Stress ratio for screw
         StressRatio_screw_p = ((abs(tau_slip(j)))/ &
         ( constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(5*ns+3*nt+j)) &
        )**(constitutive_titanmod_ps_PerSlipSystem(j,myInstance))

        if((1.0_pReal-StressRatio_edge_p)>0.001_pReal) then
        minusStressRatio_edge_p=1.0_pReal-StressRatio_edge_p
        else
        minusStressRatio_edge_p=0.001_pReal
        endif
        
        if((1-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif

        !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_f0_PerSlipSystem(j,myInstance)/(kB*Temperature)

!         if (tau_slip(j) == 0.0_pReal) then
!             edge_velocity(j) = 0.0_pReal
!             screw_velocity(j) = 0.0_pReal
!         else          
            edge_velocity(j) =constitutive_titanmod_v0e_PerSlipSystem(j,myInstance)*exp(-BoltzmannRatio* &
              (minusStressRatio_edge_p)** &
            constitutive_titanmod_qe_PerSlipSystem(j,myInstance))
            screw_velocity(j) =screwvelocity_kink_prefactor* constitutive_titanmod_v0s_PerSlipSystem(j,myInstance)* &
                exp(-BoltzmannRatio*(minusStressRatio_screw_p)** &
            constitutive_titanmod_qs_PerSlipSystem(j,myInstance))
!         endif

      !* Multiplication of edge dislocations
      DotRhoEdgeGeneration(j) = 2.0_pReal*(state(g,ip,el)%p(ns+j)*screw_velocity(j)/screw_segment(j))
      !* Multiplication of screw dislocations
      DotRhoScrewGeneration(j) = 2.0_pReal*(state(g,ip,el)%p(j)*edge_velocity(j)/edge_segment(j))

      !* Annihilation of edge dislocations
      DotRhoEdgeAnnihilation(j) = -((state(g,ip,el)%p(j))**2)* &
                constitutive_titanmod_capre_PerSlipSystem(j,myInstance)*edge_velocity(j)

      !* Annihilation of screw dislocations
      DotRhoScrewAnnihilation(j) = -((state(g,ip,el)%p(ns+j))**2)* &
                constitutive_titanmod_caprs_PerSlipSystem(j,myInstance)*screw_velocity(j)
       
      !* Edge dislocation density rate of change
      constitutive_titanmod_dotState(j) = &
        DotRhoEdgeGeneration(j)+DotRhoEdgeAnnihilation(j)

      !* Screw dislocation density rate of change
      constitutive_titanmod_dotState(ns+j) = &
        DotRhoScrewGeneration(j)+DotRhoScrewAnnihilation(j)

                                  
    enddo
  enddo
  
!* Twin fraction evolution
j = 0_pInt
do f = 1,lattice_maxNtwinFamily                                 ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1,MyStructure)) ! at which index starts my family
   do i = 1,constitutive_titanmod_Ntwin(f,myInstance)          ! process each (active) twin system in family
      j = j+1_pInt

!*************************************************************************
!This was in dislotwin - keeping it for safety
!*************************************************************************
!      !* Resolved shear stress on twin system
!      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))
!      !* Stress ratios
!      StressRatio_r = (state(g,ip,el)%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_titanmod_r(myInstance)
!      
!      !* Shear rates and their derivatives due to twin
!      if ( tau_twin(j) > 0.0_pReal ) then
!        constitutive_titanmod_dotState(2*ns+j) = &
!          (constitutive_titanmod_MaxTwinFraction(myInstance)-sumf)*&
!          state(g,ip,el)%p(6*ns+4*nt+j)*constitutive_titanmod_Ndot0PerTwinSystem(f,myInstance)*exp(-StressRatio_r) 
!      endif
!*************************************************************************

      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))

     !* Stress ratio for edge
         twinStressRatio_edge_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0e_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+5*nt+j)) &
        )**(constitutive_titanmod_twinpe_PerTwinSystem(j,myInstance))
        
     !* Stress ratio for screw
         twinStressRatio_screw_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0s_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(6*ns+6*nt+j)) &
        )**(constitutive_titanmod_twinps_PerTwinSystem(j,myInstance))

        if((1.0_pReal-twinStressRatio_edge_p)>0.001_pReal) then
        twinminusStressRatio_edge_p=1.0_pReal-twinStressRatio_edge_p
        else
        twinminusStressRatio_edge_p=0.001_pReal
        endif
        
        if((1-twinStressRatio_screw_p)>0.001_pReal) then
        twinminusStressRatio_screw_p=1.0_pReal-twinStressRatio_screw_p
        else
        twinminusStressRatio_screw_p=0.001_pReal
        endif

        !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_twinf0_PerTwinSystem(j,myInstance)/(kB*Temperature)

!         if (tau_slip(j) == 0.0_pReal) then
!             edge_velocity(j) = 0.0_pReal
!             screw_velocity(j) = 0.0_pReal
!         else          
            twinedge_velocity(j) =constitutive_titanmod_twinv0e_PerTwinSystem(j,myInstance)*exp(-BoltzmannRatio* &
              (twinminusStressRatio_edge_p)** &
            constitutive_titanmod_twinqe_PerTwinSystem(j,myInstance))
            twinscrew_velocity(j) =constitutive_titanmod_twinv0s_PerTwinSystem(j,myInstance)* &
                exp(-BoltzmannRatio*(twinminusStressRatio_screw_p)** &
            constitutive_titanmod_twinqs_PerTwinSystem(j,myInstance))
!         endif

      !* Multiplication of edge dislocations
      TwinDotRhoEdgeGeneration(j) = 2.0_pReal*(state(g,ip,el)%p(2*ns+nt+j)*twinscrew_velocity(j)/twinscrew_segment(j))
      !* Multiplication of screw dislocations
      TwinDotRhoScrewGeneration(j) = 2.0_pReal*(state(g,ip,el)%p(2*ns+j)*twinedge_velocity(j)/twinedge_segment(j))

      !* Annihilation of edge dislocations
      TwinDotRhoEdgeAnnihilation(j) = -((state(g,ip,el)%p(2*ns+j))**2)* &
                constitutive_titanmod_twincapre_PerTwinSystem(j,myInstance)*twinedge_velocity(j)

      !* Annihilation of screw dislocations
      TwinDotRhoScrewAnnihilation(j) = -((state(g,ip,el)%p(2*ns+nt+j))**2)* &
                constitutive_titanmod_twincaprs_PerTwinSystem(j,myInstance)*twinscrew_velocity(j)
       
      !* Edge dislocation density rate of change
      constitutive_titanmod_dotState(2*ns+j) = &
        TwinDotRhoEdgeGeneration(j)+TwinDotRhoEdgeAnnihilation(j)

      !* Screw dislocation density rate of change
      constitutive_titanmod_dotState(2*ns+nt+j) = &
        TwinDotRhoScrewGeneration(j)+TwinDotRhoScrewAnnihilation(j)

      gdot_twin(j) = constitutive_titanmod_burgersPerTwinSystem(j,myInstance)*(state(g,ip,el)%p(2*ns+j)* &
       twinedge_velocity(j)+state(g,ip,el)%p(2*ns+nt+j) * twinscrew_velocity(j))* sign(1.0_pReal,tau_twin(j))
      
      constitutive_titanmod_dotState(2*ns+2*nt+j)=gdot_twin(j)

    enddo
enddo

!write(6,*) '#DOTSTATE#'
!write(6,*)
!write(6,'(a,/,4(3(f30.20,x)/))') 'tau slip',tau_slip
!write(6,'(a,/,4(3(f30.20,x)/))') 'gamma slip',gdot_slip
!write(6,'(a,/,4(3(f30.20,x)/))') 'rho_edge',state(g,ip,el)%p(1:ns)
!write(6,'(a,/,4(3(f30.20,x)/))') 'Threshold Slip Edge', state(g,ip,el)%p(5*ns+3*nt+1:6*ns+3*nt)
!write(6,'(a,/,4(3(f30.20,x)/))') 'Threshold Slip Screw', state(g,ip,el)%p(6*ns+3*nt+1:7*ns+3*nt)
!write(6,'(a,/,4(3(f30.20,x)/))') 'EdgeGeneration',DotRhoEdgeGeneration
!write(6,'(a,/,4(3(f30.20,x)/))') 'ScrewGeneration',DotRhoScrewGeneration
!write(6,'(a,/,4(3(f30.20,x)/))') 'EdgeAnnihilation',DotRhoEdgeAnnihilation
!write(6,'(a,/,4(3(f30.20,x)/))') 'ScrewAnnihilation',DotRhoScrewAnnihilation
!write(6,'(a,/,4(3(f30.20,x)/))') 'DipClimb',DotRhoEdgeDipClimb 

return
end function


pure function constitutive_titanmod_dotTemperature(Tstar_v,Temperature,state,g,ip,el)
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
real(pReal) constitutive_titanmod_dotTemperature

constitutive_titanmod_dotTemperature = 0.0_pReal
    
return
end function


pure function constitutive_titanmod_postResults(Tstar_v,Temperature,dt,state,g,ip,el)
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
real(pReal) sumf,tau,StressRatio_edge_p,StressRatio_screw_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r, &
                gdot_slip,dgdot_dtauslip
real(pReal), dimension(constitutive_titanmod_sizePostResults(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
constitutive_titanmod_postResults
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
        volumefraction_pertwinsystem

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(2*ns+2*nt+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

!sumf = sum(state(g,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0


!* Required output 
c = 0_pInt
constitutive_titanmod_postResults = 0.0_pReal

do o = 1,phase_Noutput(material_phase(g,ip,el))
   select case(constitutive_titanmod_output(o,myInstance))

     case ('rhoedge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p(1:ns)
       c = c + ns
     case ('rhoscrew')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p(ns+1:2*ns)
       c = c + ns
     case ('twinrhoedge')
       constitutive_titanmod_postResults(c+1:c+nt) = state(g,ip,el)%p(2*ns+1:2*ns+nt)
       c = c + nt
     case ('twinrhoscrew')
       constitutive_titanmod_postResults(c+1:c+nt) = state(g,ip,el)%p(2*ns+nt+1:2*ns+2*nt)
       c = c + nt
     case ('gdot_slip')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((12*ns+3*nt+1):(13*ns+3*nt))
       c = c + ns
     case ('gdot_twin')
       constitutive_titanmod_postResults(c+1:c+nt) = state(g,ip,el)%p((2*ns+2*nt+1):(2*ns+3*nt))
       c = c + nt
     case ('dgdotdtau')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((13*ns+3*nt+1):(14*ns+3*nt))
       c = c + ns
     case ('velocity_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((7*ns+3*nt+1):(8*ns+3*nt))
       c = c + ns
     case ('twinvelocity_edge')
       constitutive_titanmod_postResults(c+1:c+nt) = state(g,ip,el)%p((6*ns+8*nt+1):(6*ns+9*nt))
       c = c + nt
     case ('velocity_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((8*ns+3*nt+1):(9*ns+3*nt))
       c = c + ns
     case ('segment_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((2*ns+nt+1):(3*ns+nt))
       c = c + ns
     case ('segment_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((4*ns+2*nt+1):(5*ns+2*nt))
       c = c + ns
     case ('resistance_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((5*ns+3*nt+1):(6*ns+3*nt))
       c = c + ns
         case ('resistance_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((6*ns+3*nt+1):(7*ns+3*nt))
       c = c + ns
         case ('tau_slip')
             constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p(9*ns+3*nt+1:10*ns+3*nt)
           c=c + ns
     case('edge_generation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
       constitutive_titanmod_postResults(c+j) = 2.0_pReal*state(g,ip,el)%p(ns+j)* &
           state(g,ip,el)%p(8*ns+3*nt+j)/state(g,ip,el)%p(4*ns+2*nt+j)
       enddo; enddo
       c = c + ns
     case('screw_generation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
       constitutive_titanmod_postResults(c+j) = 2.0_pReal*state(g,ip,el)%p(j)* &
           state(g,ip,el)%p(7*ns+3*nt+j)/state(g,ip,el)%p(2*ns+nt+j)
       enddo; enddo
       
       c = c + ns
     case('edge_annihilation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
       constitutive_titanmod_postResults(c+j) = -((state(g,ip,el)%p(j))**2)* &
                constitutive_titanmod_capre_PerSlipSystem(j,myInstance)*state(g,ip,el)%p(7*ns+3*nt+j)
       enddo; enddo
           
       c = c + ns
     case('screw_annihilation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
       constitutive_titanmod_postResults(c+j) = -((state(g,ip,el)%p(ns+j))**2)* &
                constitutive_titanmod_caprs_PerSlipSystem(j,myInstance)*state(g,ip,el)%p(8*ns+3*nt+j)
       enddo; enddo
           
       c = c + ns
     case('total_generation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
       constitutive_titanmod_postResults(c+j) = 2.0_pReal*state(g,ip,el)%p(ns+j)* &
           state(g,ip,el)%p(8*ns+3*nt+j)/state(g,ip,el)%p(4*ns+2*nt+j) + 2.0_pReal*state(g,ip,el)%p(j)* &
           state(g,ip,el)%p(7*ns+3*nt+j)/state(g,ip,el)%p(2*ns+nt+j)
       enddo; enddo
           
       c = c + ns
     case('total_annihilation')
       j = 0_pInt
       do f = 1,lattice_maxNslipFamily                                 
          index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) 
          do i = 1,constitutive_titanmod_Nslip(f,myInstance)          
             j = j + 1_pInt
        constitutive_titanmod_postResults(c+j) = -((state(g,ip,el)%p(j))**2)* &
                constitutive_titanmod_capre_PerSlipSystem(j,myInstance)*state(g,ip,el)%p(7*ns+3*nt+j)- &
                ((state(g,ip,el)%p(ns+j))**2)*constitutive_titanmod_caprs_PerSlipSystem(j,myInstance)* &
                state(g,ip,el)%p(8*ns+3*nt+j)
       enddo; enddo
           
       c = c + ns
     case('total_density')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p(1:ns)+state(g,ip,el)%p(ns+1:2*ns)
       c = c + ns
     case ('twin_fraction')
       constitutive_titanmod_postResults(c+1:c+nt) = volumefraction_pertwinsystem(1:nt)
       c = c + nt

   end select
enddo

return
end function

END MODULE
