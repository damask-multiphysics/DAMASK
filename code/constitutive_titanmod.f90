! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$

! states for titanmod
! Basic states

! rhoedge 
! rhoscrew 
! shear_system

! Dependent states
! segment_edge 
! segment_screw 
! resistance_edge 
! resistance_screw 
! tau_slip 
! gdot_slip 
! velocity_edge 
! velocity_screw 
! gdot_slip_edge 
! gdot_slip_screw 
! stressratio_edge_p 
! stressratio_screw_p
! shear_basal
! shear_prism
! shear_pyra
! shear_pyrca 

MODULE constitutive_titanmod

!* Include other modules
use prec, only: pReal,pInt
implicit none

!* Lists of states and physical parameters
character(len=*), parameter :: constitutive_titanmod_label = 'titanmod'
character(len=18), dimension(3), parameter:: constitutive_titanmod_listBasicSlipStates = (/'rho_edge    ', &
                                                                                           'rho_screw   ', &
                                                                                           'shear_system'/)

character(len=18), dimension(1), parameter:: constitutive_titanmod_listBasicTwinStates = (/'gdot_twin'/)
                                                                                            
character(len=19), dimension(11), parameter:: constitutive_titanmod_listDependentSlipStates =(/'segment_edge       ', &
                                                                                               'segment_screw      ', &
                                                                                               'resistance_edge    ', &
                                                                                               'resistance_screw   ', &
                                                                                               'tau_slip           ', &
                                                                                               'velocity_edge      ', &
                                                                                               'velocity_screw     ', &
                                                                                               'gdot_slip_edge     ', &
                                                                                               'gdot_slip_screw    ', &
                                                                                               'stressratio_edge_p ', &
                                                                                               'stressratio_screw_p' &
                                                                                              /)

character(len=18), dimension(2), parameter:: constitutive_titanmod_listDependentTwinStates =(/'twin_fraction', &
                                                                                              'tau_twin     ' &
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
                                                          constitutive_titanmod_debyefrequency, &              !Debye frequency
                                                          constitutive_titanmod_kinkf0, &                      !Debye frequency
                                                          constitutive_titanmod_Gmod, &                        ! shear modulus
                                                          constitutive_titanmod_CAtomicVolume, &               ! atomic volume in Bugers vector unit
                                                          constitutive_titanmod_dc, &                          ! prefactor for self-diffusion coefficient
                                                          constitutive_titanmod_twinhpconstant, &              ! activation energy for dislocation climb
                                                          constitutive_titanmod_GrainSize, &                   ! grain size - Not being used
                                                          constitutive_titanmod_MaxTwinFraction, &             ! maximum allowed total twin volume fraction
                                                          constitutive_titanmod_r, &                           ! r-exponent in twin nucleation rate
                                                          constitutive_titanmod_CEdgeDipMinDistance, &         ! Not being used
                                                          constitutive_titanmod_Cmfptwin, &                    ! Not being used
                                                          constitutive_titanmod_Cthresholdtwin, &              ! Not being used
                                                          constitutive_titanmod_aTolRho                        ! absolute tolerance for integration of dislocation density
real(pReal),       dimension(:,:,:),       allocatable :: constitutive_titanmod_Cslip_66                       ! elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:),     allocatable :: constitutive_titanmod_Ctwin_66                       ! twin elasticity matrix in Mandel notation for each instance
real(pReal),       dimension(:,:,:,:,:),   allocatable :: constitutive_titanmod_Cslip_3333                     ! elasticity matrix for each instance
real(pReal),       dimension(:,:,:,:,:,:), allocatable :: constitutive_titanmod_Ctwin_3333                     ! twin elasticity matrix for each instance
real(pReal), dimension(:,:), allocatable ::               constitutive_titanmod_rho_edge0, &                   ! initial edge dislocation density per slip system for each family and instance
                                                          constitutive_titanmod_rho_screw0, &                  ! initial screw dislocation density per slip system for each family and instance
                                                          constitutive_titanmod_shear_system0, &               ! accumulated shear on each system
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
                                                          constitutive_titanmod_twintau0_PerTwinFamily, &         ! Initial yield stress for edge dislocations per twin family
                                                          constitutive_titanmod_twintau0_PerTwinSystem, &         ! Initial yield stress for edge dislocations per twin system
                                                          constitutive_titanmod_capre_PerSlipFamily, &         ! Capture radii for edge dislocations per slip family
                                                          constitutive_titanmod_capre_PerSlipSystem, &         ! Capture radii for edge dislocations per slip system
                                                          constitutive_titanmod_caprs_PerSlipFamily, &         ! Capture radii for screw dislocations per slip family
                                                          constitutive_titanmod_caprs_PerSlipSystem, &         ! Capture radii for screw dislocations per slip system
                                                          constitutive_titanmod_pe_PerSlipFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_ps_PerSlipFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_qe_PerSlipFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_qs_PerSlipFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_pe_PerSlipSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_ps_PerSlipSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_qe_PerSlipSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_qs_PerSlipSystem, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinp_PerTwinFamily, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinq_PerTwinFamily, &            ! q-exponent in glide velocity
                                                          constitutive_titanmod_twinp_PerTwinSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_twinq_PerTwinSystem, &            ! p-exponent in glide velocity
                                                          constitutive_titanmod_v0e_PerSlipFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_v0e_PerSlipSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_v0s_PerSlipFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_v0s_PerSlipSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_twingamma0_PerTwinFamily, &           ! edge dislocation velocity prefactor [m/s] for each family and instance
                                                          constitutive_titanmod_twingamma0_PerTwinSystem, &           ! screw dislocation velocity prefactor [m/s] for each slip system and instance
                                                          constitutive_titanmod_kinkcriticallength_PerSlipFamily, &  ! screw dislocation mobility prefactor for kink-pairs per slip family
                                                          constitutive_titanmod_kinkcriticallength_PerSlipSystem, &  ! screw dislocation mobility prefactor for kink-pairs per slip system
                                                          constitutive_titanmod_twinsizePerTwinFamily, &       ! twin thickness [m] for each twin family and instance
                                                          constitutive_titanmod_twinsizePerTwinSystem, &       ! twin thickness [m] for each twin system and instance
                                                          constitutive_titanmod_CeLambdaSlipPerSlipFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_CeLambdaSlipPerSlipSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_CsLambdaSlipPerSlipFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_CsLambdaSlipPerSlipSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
                                                          constitutive_titanmod_twinLambdaSlipPerTwinFamily, &   ! Adj. parameter for distance between 2 forest dislocations for each slip family and instance
                                                          constitutive_titanmod_twinLambdaSlipPerTwinSystem, &   ! Adj. parameter for distance between 2 forest dislocations for each slip system and instance
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
integer(pInt) section,maxNinstance,f,i,j,k,l,m,n,o,p,q,r,s,s1,s2,t,t1,t2,ns,nt,output,mySize,myStructure,maxTotalNslip, &
maxTotalNtwin
character(len=64) tag
character(len=1024) line

write(6,*)
write(6,*) '<<<+-  constitutive_',trim(constitutive_titanmod_label),' init  -+>>>'
write(6,*) '$Id$'
#include "compilation_info.f90"

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
constitutive_titanmod_slipSystemLattice = 0_pInt
constitutive_titanmod_twinSystemLattice = 0_pInt
constitutive_titanmod_totalNslip        = 0_pInt
constitutive_titanmod_totalNtwin        = 0_pInt
allocate(constitutive_titanmod_CoverA(maxNinstance))
allocate(constitutive_titanmod_C11(maxNinstance))
allocate(constitutive_titanmod_C12(maxNinstance))
allocate(constitutive_titanmod_C13(maxNinstance))
allocate(constitutive_titanmod_C33(maxNinstance))
allocate(constitutive_titanmod_C44(maxNinstance))
allocate(constitutive_titanmod_debyefrequency(maxNinstance))
allocate(constitutive_titanmod_kinkf0(maxNinstance))
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
allocate(constitutive_titanmod_aTolRho(maxNinstance))
allocate(constitutive_titanmod_Cslip_66(6,6,maxNinstance))
allocate(constitutive_titanmod_Cslip_3333(3,3,3,3,maxNinstance))
constitutive_titanmod_CoverA              = 0.0_pReal 
constitutive_titanmod_C11                 = 0.0_pReal
constitutive_titanmod_C12                 = 0.0_pReal
constitutive_titanmod_C13                 = 0.0_pReal
constitutive_titanmod_C33                 = 0.0_pReal
constitutive_titanmod_C44                 = 0.0_pReal
constitutive_titanmod_debyefrequency      = 0.0_pReal
constitutive_titanmod_kinkf0              = 0.0_pReal
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
constitutive_titanmod_aTolRho             = 0.0_pReal
constitutive_titanmod_Cslip_66            = 0.0_pReal
constitutive_titanmod_Cslip_3333          = 0.0_pReal
allocate(constitutive_titanmod_rho_edge0(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_rho_screw0(lattice_maxNslipFamily,maxNinstance)) 
allocate(constitutive_titanmod_shear_system0(lattice_maxNslipFamily,maxNinstance)) 
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
allocate(constitutive_titanmod_kinkcriticallength_PerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance))
allocate(constitutive_titanmod_CeLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance))
allocate(constitutive_titanmod_CsLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance))

allocate(constitutive_titanmod_twinf0_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinshearconstant_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twintau0_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinp_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinq_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twingamma0_PerTwinFamily(lattice_maxNTwinFamily,maxNinstance))
allocate(constitutive_titanmod_twinLambdaSlipPerTwinFamily(lattice_maxNTwinFamily,maxNinstance))

constitutive_titanmod_rho_edge0                 = 0.0_pReal
constitutive_titanmod_rho_screw0              = 0.0_pReal
constitutive_titanmod_shear_system0              = 0.0_pReal
constitutive_titanmod_burgersPerSlipFamily     = 0.0_pReal
constitutive_titanmod_burgersPerTwinFamily     = 0.0_pReal
constitutive_titanmod_f0_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_tau0e_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_tau0s_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_capre_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_caprs_PerSlipFamily       = 0.0_pReal
constitutive_titanmod_v0e_PerSlipFamily          = 0.0_pReal
constitutive_titanmod_v0s_PerSlipFamily          = 0.0_pReal
constitutive_titanmod_kinkcriticallength_PerSlipFamily = 0.0_pReal
constitutive_titanmod_twinsizePerTwinFamily    = 0.0_pReal
constitutive_titanmod_CeLambdaSlipPerSlipFamily = 0.0_pReal
constitutive_titanmod_CsLambdaSlipPerSlipFamily = 0.0_pReal
constitutive_titanmod_pe_PerSlipFamily = 0.0_pReal
constitutive_titanmod_ps_PerSlipFamily = 0.0_pReal
constitutive_titanmod_qe_PerSlipFamily = 0.0_pReal
constitutive_titanmod_qs_PerSlipFamily = 0.0_pReal

constitutive_titanmod_twinf0_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinshearconstant_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twintau0_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twingamma0_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinLambdaSlipPerTwinFamily = 0.0_pReal
constitutive_titanmod_twinp_PerTwinFamily = 0.0_pReal
constitutive_titanmod_twinq_PerTwinFamily = 0.0_pReal

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

write(6,*) 'titanmod: Reading material parameters from material config file'

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
                write(6,*) tag,constitutive_titanmod_output(output,i)
       case ('lattice_structure')
              constitutive_titanmod_structureName(i) = IO_lc(IO_stringValue(line,positions,2))
                write(6,*) tag,constitutive_titanmod_structureName(i)
       case ('covera_ratio')
              constitutive_titanmod_CoverA(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_CoverA(i)
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
       case ('debyefrequency')
              constitutive_titanmod_debyefrequency(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_debyefrequency(i)
       case ('kinkf0')
              constitutive_titanmod_kinkf0(i) = IO_floatValue(line,positions,2)
                write(6,*) tag,constitutive_titanmod_kinkf0(i)
       case ('nslip')
            forall (j = 1:lattice_maxNslipFamily) &
            constitutive_titanmod_Nslip(j,i) = IO_intValue(line,positions,1+j)
            write(6,*) tag,constitutive_titanmod_Nslip(1:4,i)
       case ('ntwin')
            forall (j = 1:lattice_maxNtwinFamily) &
             constitutive_titanmod_Ntwin(j,i) = IO_intValue(line,positions,1+j)
            write(6,*) tag,constitutive_titanmod_Ntwin(1:4,i)
       case ('rho_edge0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_rho_edge0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_rho_edge0(1:4,i)
       case ('rho_screw0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_rho_screw0(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_rho_screw0(1:4,i)
       case ('slipburgers')
             forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_burgersPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
              write(6,*) tag,constitutive_titanmod_burgersPerSlipFamily(1:4,i)
       case ('twinburgers')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_burgersPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
              write(6,*) tag,constitutive_titanmod_burgersPerTwinFamily(1:4,i)
       case ('f0')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_f0_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_f0_PerSlipFamily(1:4,i)
       case ('twinf0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinf0_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinf0_PerTwinFamily(1:4,i)
       case ('tau0e')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_tau0e_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_tau0e_PerSlipFamily(1:4,i)
       case ('twintau0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twintau0_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twintau0_PerTwinFamily(1:4,i)
       case ('tau0s')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_tau0s_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_tau0s_PerSlipFamily(1:4,i)
       case ('capre')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_capre_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_capre_PerSlipFamily(1:4,i)
       case ('caprs')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_caprs_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_caprs_PerSlipFamily(1:4,i)
       case ('v0e')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_v0e_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_v0e_PerSlipFamily(1:4,i)
       case ('twingamma0')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twingamma0_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twingamma0_PerTwinFamily(1:4,i)
       case ('v0s')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_v0s_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_v0s_PerSlipFamily(1:4,i)
       case ('kinkcriticallength')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_kinkcriticallength_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
              write(6,*) tag,constitutive_titanmod_kinkcriticallength_PerSlipFamily(1:4,i)
       case ('twinsize')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinsizePerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('celambdaslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_CeLambdaSlipPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('twinlambdaslip')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinlambdaslipPerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinlambdaslipPerTwinFamily(1:4,i)
       case ('cslambdaslip')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_CsLambdaSlipPerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag
       case ('grainsize')
              constitutive_titanmod_GrainSize(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('maxtwinfraction')
              constitutive_titanmod_MaxTwinFraction(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('pe')
              forall (j = 1:lattice_maxNslipFamily) &
                constitutive_titanmod_pe_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_pe_PerSlipFamily(1:4,i)
       case ('twinp')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinp_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinp_PerTwinFamily(1:4,i)
       case ('ps')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_ps_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_ps_PerSlipFamily(1:4,i)
       case ('qe')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_qe_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_qe_PerSlipFamily(1:4,i)
       case ('twinq')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinq_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinq_PerTwinFamily(1:4,i)
       case ('qs')
                          forall (j = 1:lattice_maxNslipFamily) &
                                constitutive_titanmod_qs_PerSlipFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_qs_PerSlipFamily(1:4,i)
       case ('twinshearconstant')
              forall (j = 1:lattice_maxNtwinFamily) &
                constitutive_titanmod_twinshearconstant_PerTwinFamily(j,i) = IO_floatValue(line,positions,1+j)
                write(6,*) tag,constitutive_titanmod_twinshearconstant_PerTwinFamily(1:4,i)
       case ('dc')
              constitutive_titanmod_dc(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('twinhpconstant')
              constitutive_titanmod_twinhpconstant(i) = IO_floatValue(line,positions,2)
                write(6,*) tag
       case ('atol_rho')
              constitutive_titanmod_aTolRho(i) = IO_floatValue(line,positions,2)
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
       if (constitutive_titanmod_kinkcriticallength_PerSlipFamily(f,i) <= 0.0_pReal) call IO_error(238)
     endif
   enddo
   do f = 1,lattice_maxNtwinFamily
     if (constitutive_titanmod_Ntwin(f,i) > 0_pInt) then   
       if (constitutive_titanmod_burgersPerTwinFamily(f,i) <= 0.0_pReal)    call IO_error(221) !***
       if (constitutive_titanmod_twinf0_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(228)
       if (constitutive_titanmod_twinshearconstant_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(228)
       if (constitutive_titanmod_twintau0_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(229)
       if (constitutive_titanmod_twingamma0_PerTwinFamily(f,i) <= 0.0_pReal)         call IO_error(226)
     endif
   enddo
!   if (any(constitutive_titanmod_interactionSlipSlip(1:maxval(lattice_interactionSlipSlip(:,:,myStructure)),i) < 1.0_pReal)) call IO_error(229)
   if (constitutive_titanmod_dc(i) <= 0.0_pReal)                            call IO_error(231)
   if (constitutive_titanmod_twinhpconstant(i) <= 0.0_pReal)                call IO_error(232)
   if (constitutive_titanmod_aTolRho(i) <= 0.0_pReal)                       call IO_error(233)
   
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
allocate(constitutive_titanmod_kinkcriticallength_PerSlipSystem(maxTotalNslip,maxNinstance))
allocate(constitutive_titanmod_CeLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance))
allocate(constitutive_titanmod_CsLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance))

allocate(constitutive_titanmod_twinf0_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinshearconstant_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twintau0_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinp_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twinq_PerTwinSystem(maxTotalNTwin,maxNinstance))
allocate(constitutive_titanmod_twingamma0_PerTwinSystem(maxTotalNTwin,maxNinstance))

allocate(constitutive_titanmod_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance))
allocate(constitutive_titanmod_twinLambdaSlipPerTwinSystem(maxTotalNtwin, maxNinstance))

constitutive_titanmod_burgersPerSlipSystem     = 0.0_pReal
constitutive_titanmod_burgersPerTwinSystem     = 0.0_pReal
constitutive_titanmod_f0_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_tau0e_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_tau0s_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_capre_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_caprs_PerSlipSystem       = 0.0_pReal
constitutive_titanmod_v0e_PerSlipSystem          = 0.0_pReal
constitutive_titanmod_v0s_PerSlipSystem          = 0.0_pReal
constitutive_titanmod_kinkcriticallength_PerSlipSystem = 0.0_pReal
constitutive_titanmod_pe_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_ps_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_qe_PerSlipSystem        = 0.0_pReal
constitutive_titanmod_qs_PerSlipSystem        = 0.0_pReal

constitutive_titanmod_twinf0_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twinshearconstant_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twintau0_PerTwinSystem       = 0.0_pReal
constitutive_titanmod_twingamma0_PerTwinSystem          = 0.0_pReal
constitutive_titanmod_twinp_PerTwinSystem        = 0.0_pReal
constitutive_titanmod_twinq_PerTwinSystem        = 0.0_pReal

constitutive_titanmod_twinsizePerTwinSystem    = 0.0_pReal
constitutive_titanmod_CeLambdaSlipPerSlipSystem = 0.0_pReal
constitutive_titanmod_CsLambdaSlipPerSlipSystem = 0.0_pReal
constitutive_titanmod_twinLambdaSlipPerTwinSystem = 0.0_pReal

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
             'velocity_edge', &
             'velocity_screw', &
             'tau_slip', &
             'gdot_slip_edge', &
             'gdot_slip_screw', &
             'gdot_slip', &
             'stressratio_edge_p', &
             'stressratio_screw_p', &
             'shear_system' &
             )
           mySize = constitutive_titanmod_totalNslip(i)
        case('twin_fraction', &
             'gdot_twin', &
             'tau_twin' &
             )
           mySize = constitutive_titanmod_totalNtwin(i)
        case('shear_basal', & ! use only if all 4 slip familiies in hex are considered
             'shear_prism', & ! use only if all 4 slip familiies in hex are considered
             'shear_pyra', & ! use only if all 4 slip familiies in hex are considered
             'shear_pyrca', & ! use only if all 4 slip familiies in hex are considered
             'rhoedge_basal', &
             'rhoedge_prism', &
             'rhoedge_pyra', &
             'rhoedge_pyrca', &
             'rhoscrew_basal', &
             'rhoscrew_prism', &
             'rhoscrew_pyra', &
             'rhoscrew_pyrca', &
             'shear_total' &
             )
           mySize = 1_pInt
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
   constitutive_titanmod_Cslip_66(:,:,i) =     math_Mandel3333to66(math_Voigt66to3333(constitutive_titanmod_Cslip_66(:,:,i)))
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
      constitutive_titanmod_kinkcriticallength_PerSlipSystem(s,i) = constitutive_titanmod_kinkcriticallength_PerSlipFamily(f,i)
      constitutive_titanmod_pe_PerSlipSystem(s,i)          = constitutive_titanmod_pe_PerSlipFamily(f,i)
      constitutive_titanmod_ps_PerSlipSystem(s,i)          = constitutive_titanmod_ps_PerSlipFamily(f,i)
      constitutive_titanmod_qe_PerSlipSystem(s,i)          = constitutive_titanmod_qe_PerSlipFamily(f,i)
      constitutive_titanmod_qs_PerSlipSystem(s,i)          = constitutive_titanmod_qs_PerSlipFamily(f,i)
      constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,i) = constitutive_titanmod_CeLambdaSlipPerSlipFamily(f,i)
      constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,i) = constitutive_titanmod_CsLambdaSlipPerSlipFamily(f,i)
   enddo   
   
   !* Burgers vector, nucleation rate prefactor and twin size for each twin system 
   do t = 1,constitutive_titanmod_totalNtwin(i)   
      f = constitutive_titanmod_twinFamily(t,i)    
      constitutive_titanmod_burgersPerTwinSystem(t,i)  = constitutive_titanmod_burgersPerTwinFamily(f,i)
      constitutive_titanmod_twinsizePerTwinSystem(t,i) = constitutive_titanmod_twinsizePerTwinFamily(f,i)
      constitutive_titanmod_twinf0_PerTwinSystem(t,i)       = constitutive_titanmod_twinf0_PerTwinFamily(f,i)
      constitutive_titanmod_twinshearconstant_PerTwinSystem(t,i) = constitutive_titanmod_twinshearconstant_PerTwinFamily(f,i)
      constitutive_titanmod_twintau0_PerTwinSystem(t,i)       = constitutive_titanmod_twintau0_PerTwinFamily(f,i)
      constitutive_titanmod_twingamma0_PerTwinSystem(t,i)          = constitutive_titanmod_twingamma0_PerTwinFamily(f,i)
      constitutive_titanmod_twinp_PerTwinSystem(t,i)          = constitutive_titanmod_twinp_PerTwinFamily(f,i)
      constitutive_titanmod_twinq_PerTwinSystem(t,i)          = constitutive_titanmod_twinq_PerTwinFamily(f,i)
      constitutive_titanmod_twinLambdaSlipPerTwinSystem(t,i) = constitutive_titanmod_twinLambdaSlipPerTwinFamily(f,i)

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
                                                                         shear_system0, &
                                                                         segment_edge0, &
                                                                         segment_screw0, &
                                                                         resistance_edge0, &
                                                                         resistance_screw0
real(pReal), dimension(constitutive_titanmod_totalNtwin(myInstance)) ::  twingamma_dot0, &
                                                                         resistance_twin0

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
      shear_system0(s) = 0.0_pReal
   enddo 
enddo

!* Initialize basic slip state variables
! For twin
ts1 = 0_pInt
do tf = 1,lattice_maxNtwinFamily
   ts0 = ts1 + 1_pInt
   ts1 = ts0 + constitutive_titanmod_Ntwin(tf,myInstance) - 1_pInt 
   do ts = ts0,ts1
      twingamma_dot0(ts)=0.0_pReal
   enddo 
enddo

constitutive_titanmod_stateInit(1:ns)      = rho_edge0
constitutive_titanmod_stateInit(ns+1:2*ns) = rho_screw0
constitutive_titanmod_stateInit(2*ns+1:3*ns) = shear_system0
constitutive_titanmod_stateInit(3*ns+1:3*ns+nt)  = twingamma_dot0

!* Initialize dependent slip microstructural variables
forall (s = 1:ns) &
segment_edge0(s) = constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
        dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
 
constitutive_titanmod_stateInit(3*ns+nt+1:4*ns+nt) = segment_edge0

forall (s = 1:ns) &
segment_screw0(s) = constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
        dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
  
constitutive_titanmod_stateInit(4*ns+nt+1:5*ns+nt) = segment_screw0

forall (s = 1:ns) &
resistance_edge0(s) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)* &
sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+dot_product((rho_screw0), &
        constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))
        
constitutive_titanmod_stateInit(5*ns+nt+1:6*ns+nt) = resistance_edge0

forall (s = 1:ns) &
resistance_screw0(s) = &
constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)* &
sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+dot_product((rho_screw0), &
constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))

constitutive_titanmod_stateInit(6*ns+nt+1:7*ns+nt) = resistance_screw0

forall (t = 1:nt) &
resistance_twin0(t) = 0.0_pReal
constitutive_titanmod_stateInit(7*ns+nt+1:7*ns+2*nt)=resistance_twin0

return
end function

pure function constitutive_titanmod_aTolState(myInstance)
!*********************************************************************
!* absolute state tolerance                                          *
!*********************************************************************
use prec,     only: pReal, pInt
implicit none

!* Input-Output variables
integer(pInt), intent(in) :: myInstance
real(pReal), dimension(constitutive_titanmod_sizeState(myInstance)) :: constitutive_titanmod_aTolState

constitutive_titanmod_aTolState = constitutive_titanmod_aTolRho(myInstance)

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
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(3*ns+i)/ &
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
             volumefraction_pertwinsystem
 
!* Shortened notation
myInstance = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_titanmod_structure(myInstance)
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

! Need to update this list
!* State: 1           :  ns         rho_edge
!* State: ns+1        :  2*ns       rho_screw
!* State: 2*ns+1      :  3*ns       shear_system
!* State: 3*ns+1      :  3*ns+nt    gamma_twin
!* State: 3*ns+nt+1   :  4*ns+nt    segment_edge
!* State: 4*ns+nt+1   :  5*ns+nt    segment_screw
!* State: 5*ns+nt+1   :  6*ns+nt    resistance_edge
!* State: 6*ns+nt+1   :  7*ns+nt    resistance_screw
!* State: 7*ns+nt+1   :  7*ns+2*nt  resistance_twin
!* State: 7*ns+2*nt+1 :  8*ns+2*nt  velocity_edge
!* State: 8*ns+2*nt+1 :  9*ns+2*nt  velocity_screw
!* State: 9*ns+2*nt+1  :  10*ns+2*nt  tau_slip
!* State: 10*ns+2*nt+1 :  11*ns+2*nt  gdot_slip_edge
!* State: 11*ns+2*nt+1 :  12*ns+2*nt  gdot_slip_screw
!* State: 12*ns+2*nt+1 :  13*ns+2*nt  StressRatio_edge_p
!* State: 13*ns+2*nt+1 :  14*ns+2*nt  StressRatio_screw_p

!* Total twin volume fraction
do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(3*ns+i)/ &
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
  state(g,ip,el)%p(3*ns+nt+s) = constitutive_titanmod_CeLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product(state(g,ip,el)%p(1:ns), &
        constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
        dot_product(state(g,ip,el)%p(ns+1:2*ns), &
        constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
        
! average segment length for screw dislocations in matrix
forall (s = 1:ns) &
  state(g,ip,el)%p(4*ns+nt+s) = constitutive_titanmod_CsLambdaSlipPerSlipSystem(s,myInstance)/ &
        sqrt(dot_product(state(g,ip,el)%p(1:ns), &
        constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
        dot_product(state(g,ip,el)%p(ns+1:2*ns), &
        constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))


!* threshold stress or slip resistance for edge dislocation motion
forall (s = 1:ns) &
  state(g,ip,el)%p(5*ns+nt+s) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(1:ns)),&
                         constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+ &
                      dot_product((state(g,ip,el)%p(ns+1:2*ns)),&
                         constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))

!* threshold stress or slip resistance for screw dislocation motion
forall (s = 1:ns) &
  state(g,ip,el)%p(6*ns+nt+s) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSystem(s,myInstance)*&
    sqrt(dot_product((state(g,ip,el)%p(1:ns)),&
                         constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+ &
                        dot_product((state(g,ip,el)%p(ns+1:2*ns)),&
                         constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))

!* threshold stress or slip resistance for dislocation motion in twin
forall (t = 1:nt) &
  state(g,ip,el)%p(7*ns+nt+t) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSystem(t,myInstance)*&
    (dot_product((abs(state(g,ip,el)%p(2*ns+1:2*ns+nt))),&
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
        StressRatio_screw_pminus1, BoltzmannRatioedge, minusStressRatio_screw_p, &
        screwvelocity_prefactor,twinStressRatio_p,twinminusStressRatio_p,twinStressRatio_pminus1, &
        twinDotGamma0,BoltzmannRatioscrew,BoltzmannRatiotwin,bottomstress_edge,bottomstress_screw
real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
real(pReal), dimension(constitutive_titanmod_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_slip,dgdot_dtauslip,tau_slip, edge_velocity, screw_velocity,gdot_slip_edge,gdot_slip_screw
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
   gdot_twin,dgdot_dtautwin,tau_twin, volumefraction_pertwinsystem

!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
myStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(3*ns+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0


Lp = 0.0_pReal
dLp_dTstar3333 = 0.0_pReal
dLp_dTstar = 0.0_pReal

!* Dislocation glide part
gdot_slip = 0.0_pReal
gdot_slip_edge = 0.0_pReal
gdot_slip_screw = 0.0_pReal
dgdot_dtauslip = 0.0_pReal
j = 0_pInt
do f = 1,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure)) ! at which index starts my family
   do i = 1,constitutive_titanmod_Nslip(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,index_myFamily+i,myStructure)) 
      !*************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!      if(myStructure>=3.and.j>3) then ! for all non-basal slip systems
      if(myStructure==3) then ! only for prismatic and pyr <a> systems in hex
      screwvelocity_prefactor=constitutive_titanmod_debyefrequency(myInstance)* &
        state(g,ip,el)%p(4*ns+nt+j)*(constitutive_titanmod_burgersPerSlipSystem(j,myInstance)/ &
        constitutive_titanmod_kinkcriticallength_PerSlipSystem(j,myInstance))**2

     !* Stress ratio for screw ! No slip resistance for screw dislocations, only Peierls stress
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)
         StressRatio_screw_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**constitutive_titanmod_ps_PerSlipSystem(j,myInstance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif
        
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)
      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**(constitutive_titanmod_ps_PerSlipSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_kinkf0(myInstance)/(kB*Temperature)

        else  ! if the structure is not hex or the slip family is basal
        screwvelocity_prefactor=constitutive_titanmod_v0s_PerSlipSystem(j,myInstance)
        bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(6*ns+nt+j)
        StressRatio_screw_p = ((abs(tau_slip(j)))/( bottomstress_screw ))**constitutive_titanmod_ps_PerSlipSystem(j,myInstance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif

      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/( bottomstress_screw))** &
      (constitutive_titanmod_ps_PerSlipSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_f0_PerSlipSystem(j,myInstance)/(kB*Temperature)

        endif

        !* Stress ratio for edge
         bottomstress_edge=constitutive_titanmod_tau0e_PerSlipSystem(j,myInstance)+state(g,ip,el)%p(5*ns+nt+j)
         StressRatio_edge_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_edge) &
        )**constitutive_titanmod_pe_PerSlipSystem(j,myInstance)
                                
        if((1.0_pReal-StressRatio_edge_p)>0.001_pReal) then
        minusStressRatio_edge_p=1.0_pReal-StressRatio_edge_p
        else
        minusStressRatio_edge_p=0.001_pReal
        endif
        
      StressRatio_edge_pminus1 = ((abs(tau_slip(j)))/( bottomstress_edge))** &
      (constitutive_titanmod_pe_PerSlipSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for edge. For screws it is defined above
      BoltzmannRatioedge = constitutive_titanmod_f0_PerSlipSystem(j,myInstance)/(kB*Temperature)
            
        screw_velocity(j) =screwvelocity_prefactor * & ! there is no v0 for screw now because it is included in the prefactor
                exp(-BoltzmannRatioscrew*(minusStressRatio_screw_p)** &
        constitutive_titanmod_qs_PerSlipSystem(j,myInstance))

        edge_velocity(j) =constitutive_titanmod_v0e_PerSlipSystem(j,myInstance)*exp(-BoltzmannRatioedge* &
        (minusStressRatio_edge_p)** &
        constitutive_titanmod_qe_PerSlipSystem(j,myInstance))

                !* Shear rates due to edge slip
       gdot_slip_edge(j) = constitutive_titanmod_burgersPerSlipSystem(j,myInstance)*(state(g,ip,el)%p(j)* &
                edge_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !* Shear rates due to screw slip
       gdot_slip_screw(j) = constitutive_titanmod_burgersPerSlipSystem(j,myInstance)*(state(g,ip,el)%p(ns+j) * &
                screw_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !Total shear rate
            
       gdot_slip(j) = gdot_slip_edge(j) + gdot_slip_screw(j)
                
       state(g,ip,el)%p(7*ns+2*nt+j)=edge_velocity(j)
       state(g,ip,el)%p(8*ns+2*nt+j)=screw_velocity(j)
       state(g,ip,el)%p(9*ns+2*nt+j)=tau_slip(j)
       state(g,ip,el)%p(10*ns+2*nt+j)=gdot_slip_edge(j)
       state(g,ip,el)%p(11*ns+2*nt+j)=gdot_slip_screw(j)
       state(g,ip,el)%p(12*ns+2*nt+j)=StressRatio_edge_p
       state(g,ip,el)%p(13*ns+2*nt+j)=StressRatio_screw_p
                
      !* Derivatives of shear rates
      dgdot_dtauslip(j) = constitutive_titanmod_burgersPerSlipSystem(j,myInstance)*(( &
                ( &
                ( &
                ( &
                (edge_velocity(j)*state(g,ip,el)%p(j))) * &
                BoltzmannRatioedge*&
        constitutive_titanmod_pe_PerSlipSystem(j,myInstance)* &
                constitutive_titanmod_qe_PerSlipSystem(j,myInstance) &
                )/ &
                bottomstress_edge &
                )*&
        StressRatio_edge_pminus1*(minusStressRatio_edge_p)** &
                (constitutive_titanmod_qe_PerSlipSystem(j,myInstance)-1.0_pReal) &
                ) + &
                ( &
                ( &
                ( &
                (state(g,ip,el)%p(ns+j) * screw_velocity(j)) * &
                BoltzmannRatioscrew* &
        constitutive_titanmod_ps_PerSlipSystem(j,myInstance)* &
                constitutive_titanmod_qs_PerSlipSystem(j,myInstance) &
                )/ &
                bottomstress_screw &
                )*&
        StressRatio_screw_pminus1*(minusStressRatio_screw_p)**(constitutive_titanmod_qs_PerSlipSystem(j,myInstance)-1.0_pReal) &
                ) &
                ) !* sign(1.0_pReal,tau_slip(j))
                

                
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
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(7*ns+nt+j)) &
        )**constitutive_titanmod_twinp_PerTwinSystem(j,myInstance)
               
        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
              
      twinStressRatio_pminus1 = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(7*ns+nt+j)) &
        )**(constitutive_titanmod_twinp_PerTwinSystem(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio
      BoltzmannRatiotwin = constitutive_titanmod_twinf0_PerTwinSystem(j,myInstance)/(kB*Temperature)

      !* Initial twin shear rates
      TwinDotGamma0 = &
        constitutive_titanmod_twingamma0_PerTwinSystem(j,myInstance)

      !* Shear rates due to twin
         gdot_twin(j) =sign(1.0_pReal,tau_twin(j))*constitutive_titanmod_twingamma0_PerTwinSystem(j,myInstance)* &
         exp(-BoltzmannRatiotwin*(twinminusStressRatio_p)**constitutive_titanmod_twinq_PerTwinSystem(j,myInstance))
         
                                      
      !* Derivatives of shear rates in twin
      dgdot_dtautwin(j) = ( &
                ( &
                ( &
                (abs(gdot_twin(j))) * &
                BoltzmannRatiotwin*&
        constitutive_titanmod_twinp_PerTwinSystem(j,myInstance)* &
                constitutive_titanmod_twinq_PerTwinSystem(j,myInstance) &
                )/ &
                constitutive_titanmod_twintau0_PerTwinSystem(j,myInstance) &
                )*&
        twinStressRatio_pminus1*(twinminusStressRatio_p)** &
                (constitutive_titanmod_twinq_PerTwinSystem(j,myInstance)-1.0_pReal) &
                ) !* sign(1.0_pReal,tau_slip(j))
                
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
!   write(6,'(a10,/,4(3(e20.8,1x),/))') 'state',state(1,1,1)%p
!   write(6,'(a,/,3(3(f10.4,1x)/))') 'Lp',Lp
!   write(6,'(a,/,9(9(f10.4,1x)/))') 'dLp_dTstar',dLp_dTstar
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
integer(pInt) MyInstance,MyStructure,ns,nt,f,i,j,index_myFamily
real(pReal) sumf,BoltzmannRatio,&
            twinStressRatio_p,twinminusStressRatio_p
real(pReal), dimension(constitutive_titanmod_totalNslip(phase_constitutionInstance(material_phase(g,ip,el)))) :: &
DotRhoEdgeGeneration,DotRhoEdgeAnnihilation,DotRhoScrewAnnihilation,&
DotRhoScrewGeneration
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_constitutionInstance(material_phase(g,ip,el)))) :: gdot_twin, &
tau_twin, &
volumefraction_pertwinsystem
   
!* Shortened notation
myInstance  = phase_constitutionInstance(material_phase(g,ip,el))
MyStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1,nt
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(3*ns+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSystem(i,myInstance)

enddo

sumf = sum(abs(volumefraction_pertwinsystem(1:nt))) ! safe for nt == 0

constitutive_titanmod_dotState = 0.0_pReal

    j = 0_pInt
 do f = 1,lattice_maxNslipFamily                                             ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1,myStructure))                 ! at which index starts my family
   do i = 1,constitutive_titanmod_Nslip(f,myInstance)                        ! process each (active) slip system in family
     j = j+1_pInt

      !* Multiplication of edge dislocations
      DotRhoEdgeGeneration(j) = (state(g,ip,el)%p(ns+j)*state(g,ip,el)%p(8*ns+2*nt+j)/state(g,ip,el)%p(4*ns+nt+j))
      !* Multiplication of screw dislocations
      DotRhoScrewGeneration(j) = (state(g,ip,el)%p(j)*state(g,ip,el)%p(7*ns+2*nt+j)/state(g,ip,el)%p(3*ns+nt+j))

      !* Annihilation of edge dislocations
      DotRhoEdgeAnnihilation(j) = -((state(g,ip,el)%p(j))**2)* &
                constitutive_titanmod_capre_PerSlipSystem(j,myInstance)*state(g,ip,el)%p(7*ns+2*nt+j)/2.0_pReal

      !* Annihilation of screw dislocations
      DotRhoScrewAnnihilation(j) = -((state(g,ip,el)%p(ns+j))**2)* &
                constitutive_titanmod_caprs_PerSlipSystem(j,myInstance)*state(g,ip,el)%p(8*ns+2*nt+j)/2.0_pReal
       
      !* Edge dislocation density rate of change
      constitutive_titanmod_dotState(j) = &
        DotRhoEdgeGeneration(j)+DotRhoEdgeAnnihilation(j)

      !* Screw dislocation density rate of change
      constitutive_titanmod_dotState(ns+j) = &
        DotRhoScrewGeneration(j)+DotRhoScrewAnnihilation(j)

      constitutive_titanmod_dotState(2*ns+j) = &
                                  state(g,ip,el)%p(10*ns+2*nt+j)+state(g,ip,el)%p(11*ns+2*nt+j) ! sum of shear due to edge and screw 
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
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSystem(j,myInstance)+state(g,ip,el)%p(7*ns+nt+j)) &
        )**(constitutive_titanmod_twinp_PerTwinSystem(j,myInstance))
        

        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
        
        !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_twinf0_PerTwinSystem(j,myInstance)/(kB*Temperature)
       
            gdot_twin(j) =constitutive_titanmod_twingamma0_PerTwinSystem(j,myInstance)*exp(-BoltzmannRatio* &
              (twinminusStressRatio_p)** &
            constitutive_titanmod_twinq_PerTwinSystem(j,myInstance))*sign(1.0_pReal,tau_twin(j))
             
      constitutive_titanmod_dotState(3*ns+j)=gdot_twin(j)

    enddo
enddo

!write(6,*) '#DOTSTATE#'
!write(6,*)
!write(6,'(a,/,4(3(f30.20,1x)/))') 'EdgeGeneration',DotRhoEdgeGeneration
!write(6,'(a,/,4(3(f30.20,1x)/))') 'ScrewGeneration',DotRhoScrewGeneration
!write(6,'(a,/,4(3(f30.20,1x)/))') 'EdgeAnnihilation',DotRhoEdgeAnnihilation
!write(6,'(a,/,4(3(f30.20,1x)/))') 'ScrewAnnihilation',DotRhoScrewAnnihilation


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
integer(pInt) myInstance,myStructure,ns,nt,o,i,c
real(pReal) sumf
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
volumefraction_pertwinsystem(i)=state(g,ip,el)%p(3*ns+i)/ &
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
     case ('segment_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((3*ns+nt+1):(4*ns+nt))
       c = c + ns
     case ('segment_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((4*ns+nt+1):(5*ns+nt))
       c = c + ns
     case ('resistance_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((5*ns+nt+1):(6*ns+nt))
       c = c + ns
     case ('resistance_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((6*ns+nt+1):(7*ns+nt))
       c = c + ns
     case ('velocity_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((7*ns+2*nt+1):(8*ns+2*nt))
       c = c + ns
     case ('velocity_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = state(g,ip,el)%p((8*ns+2*nt+1):(9*ns+2*nt))
       c = c + ns
     case ('tau_slip')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((9*ns+2*nt+1):(10*ns+2*nt)))
       c = c + ns
     case ('gdot_slip_edge')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((10*ns+2*nt+1):(11*ns+2*nt)))
       c = c + ns
     case ('gdot_slip_screw')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((11*ns+2*nt+1):(12*ns+2*nt)))
       c = c + ns
     case ('gdot_slip')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((10*ns+2*nt+1):(11*ns+2*nt))) + &
                                                  abs(state(g,ip,el)%p((11*ns+2*nt+1):(12*ns+2*nt)))
       c = c + ns
     case ('stressratio_edge_p')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((12*ns+2*nt+1):(13*ns+2*nt)))
       c = c + ns
     case ('stressratio_screw_p')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((13*ns+2*nt+1):(14*ns+2*nt)))
       c = c + ns
     case ('shear_system')
       constitutive_titanmod_postResults(c+1:c+ns) = abs(state(g,ip,el)%p((2*ns+1):(3*ns)))
       c = c + ns
     case ('shear_basal')
       constitutive_titanmod_postResults(c+1:c+1) = sum(abs(state(g,ip,el)%p((2*ns+1):(2*ns+3))))
       c = c + 1
     case ('shear_prism')
       constitutive_titanmod_postResults(c+1:c+1) = sum(abs(state(g,ip,el)%p((2*ns+4):(2*ns+6))))
       c = c + 1
     case ('shear_pyra')
       constitutive_titanmod_postResults(c+1:c+1) = sum(abs(state(g,ip,el)%p((2*ns+7):(2*ns+12))))
       c = c + 1
     case ('shear_pyrca')
       constitutive_titanmod_postResults(c+1:c+1) = sum(abs(state(g,ip,el)%p((2*ns+13):(2*ns+24))))
       c = c + 1

     case ('rhoedge_basal')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((1):(3)))
       c = c + 1
     case ('rhoedge_prism')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((4):(6)))
       c = c + 1
     case ('rhoedge_pyra')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((7):(12)))
       c = c + 1
     case ('rhoedge_pyrca')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((13):(24)))
       c = c + 1

     case ('rhoscrew_basal')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((ns+1):(ns+3)))
       c = c + 1
     case ('rhoscrew_prism')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((ns+4):(ns+6)))
       c = c + 1
     case ('rhoscrew_pyra')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((ns+7):(ns+12)))
       c = c + 1
     case ('rhoscrew_pyrca')
       constitutive_titanmod_postResults(c+1:c+1) = sum(state(g,ip,el)%p((ns+13):(ns+24)))
       c = c + 1

       case ('shear_total')
       constitutive_titanmod_postResults(c+1:c+1) = sum(abs(state(g,ip,el)%p((2*ns+1):(3*ns))))
       c = c + 1
     case ('twin_fraction')
       constitutive_titanmod_postResults(c+1:c+nt) = abs(volumefraction_pertwinsystem(1:nt))
       c = c + nt
             ! 'rhoedge', &
             ! 'rhoscrew', &
             ! 'segment_edge', &
             ! 'segment_screw', &
             ! 'resistance_edge', &
             ! 'resistance_screw', &
             ! 'velocity_edge', &
             ! 'velocity_screw', &
             ! 'tau_slip', &
             ! 'gdot_slip_edge' &
             ! 'gdot_slip_screw' &
             ! 'gdot_slip', &
             ! 'StressRatio_edge_p' &
             ! 'StressRatio_screw_p'
             ! 'shear_total', &
             ! 'shear_basal' &
             ! 'shear_prism', &
             ! 'shear_pyra', &
             ! 'shear_pyrca', &
             ! 'twin_fraction', &

   end select
enddo

return
end function

END MODULE
