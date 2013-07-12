! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Alankar Alankar, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for titanium
!--------------------------------------------------------------------------------------------------
module constitutive_titanmod
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 character(len=*), parameter, public :: &
   CONSTITUTIVE_TITANMOD_label = 'titanmod'
 character(len=18), dimension(3), parameter :: &
   CONSTITUTIVE_TITANMOD_listBasicSlipStates = ['rho_edge    ', &
                                                'rho_screw   ', &
                                                'shear_system']

 character(len=18), dimension(1), parameter :: &
   CONSTITUTIVE_TITANMOD_listBasicTwinStates = ['gdot_twin']
                                                                                            
 character(len=19), dimension(11), parameter :: &
   CONSTITUTIVE_TITANMOD_listDependentSlipStates =['segment_edge       ', &
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
                                                   ]

 character(len=18), dimension(2), parameter :: &
   constitutive_titanmod_listDependentTwinStates =['twin_fraction', &
                                                   'tau_twin     ' &
                                                  ]
 real(pReal), parameter :: kB = 1.38e-23_pReal                                                       !< Boltzmann constant in J/Kelvin


 integer(pInt), dimension(:), allocatable, public, protected :: &
   constitutive_titanmod_sizeState, &                                                                !< total number of microstructural state variables
   constitutive_titanmod_sizeDotState, &                                                             !< number of dotStates
   constitutive_titanmod_sizePostResults                                                             !<  cumulative size of post results

 integer(pInt), dimension(:,:), allocatable, target, public :: &
   constitutive_titanmod_sizePostResult                                                              !<  size of each post result output

 character(len=64), dimension(:,:), allocatable, target, public :: &
   constitutive_titanmod_output                                                                      !<  name of each post result output

 integer(pInt), dimension(:), allocatable :: & 
   constitutive_titanmod_Noutput                                                                     !<  number of outputs per instance of this plasticity 

 character(len=32), dimension(:), allocatable, public, protected :: &
   constitutive_titanmod_structureName                                                               !<  name of the lattice structure

 integer(pInt), dimension(:), allocatable :: & 
   constitutive_titanmod_structure, &                                                                !<  number representing the kind of lattice structure
   constitutive_titanmod_totalNslip, &                                                               !<  total number of active slip systems for each instance
   constitutive_titanmod_totalNtwin                                                                  !<  total number of active twin systems for each instance

 integer(pInt), dimension(:,:), allocatable :: &
   constitutive_titanmod_Nslip, &                                                                    !<  number of active slip systems for each family and instance
   constitutive_titanmod_Ntwin, &                                                                    !< number of active twin systems for each family and instance
   constitutive_titanmod_slipFamily, &                                                               !< lookup table relating active slip system to slip family for each instance
   constitutive_titanmod_twinFamily, &                                                               !< lookup table relating active twin system to twin family for each instance
   constitutive_titanmod_slipSystemLattice, &                                                        !< lookup table relating active slip system index to lattice slip system index for each instance
   constitutive_titanmod_twinSystemLattice                                                           !< lookup table relating active twin system index to lattice twin system index for each instance

real(pReal), dimension(:), allocatable :: &
   constitutive_titanmod_CoverA, &                                                                   !< c/a ratio for hex type lattice
   constitutive_titanmod_debyefrequency, &                                                           !< Debye frequency
   constitutive_titanmod_kinkf0, &                                                                   !<
   constitutive_titanmod_Gmod, &                                                                     !< shear modulus
   constitutive_titanmod_CAtomicVolume, &                                                            !< atomic volume in Bugers vector unit
   constitutive_titanmod_dc, &                                                                       !< prefactor for self-diffusion coefficient
   constitutive_titanmod_twinhpconstant, &                                                           !< activation energy for dislocation climb
   constitutive_titanmod_GrainSize, &                                                                !< grain size - Not being used
   constitutive_titanmod_MaxTwinFraction, &                                                          !< maximum allowed total twin volume fraction
   constitutive_titanmod_r, &                                                                        !< r-exponent in twin nucleation rate
   constitutive_titanmod_CEdgeDipMinDistance, &                                                      !< Not being used
   constitutive_titanmod_Cmfptwin, &                                                                 !< Not being used
   constitutive_titanmod_Cthresholdtwin, &                                                           !< Not being used
   constitutive_titanmod_aTolRho                                                                     !< absolute tolerance for integration of dislocation density

real(pReal), dimension(:,:), allocatable :: &
   constitutive_titanmod_rho_edge0, &                                                                !< initial edge dislocation density per slip system for each family and instance
   constitutive_titanmod_rho_screw0, &                                                               !< initial screw dislocation density per slip system for each family and instance
   constitutive_titanmod_shear_system0, &                                                            !< accumulated shear on each system
   constitutive_titanmod_burgersPerSlipFam, &                                                     !< absolute length of burgers vector [m] for each slip family and instance
   constitutive_titanmod_burgersPerSlipSys, &                                                     !< absolute length of burgers vector [m] for each slip system and instance
   constitutive_titanmod_burgersPerTwinFam, &                                                     !< absolute length of burgers vector [m] for each twin family and instance
   constitutive_titanmod_burgersPerTwinSys, &                                                     !< absolute length of burgers vector [m] for each twin system and instance
   constitutive_titanmod_f0_PerSlipFam, &                                                         !< activation energy for glide [J] for each slip family and instance
   constitutive_titanmod_f0_PerSlipSys, &                                                         !< activation energy for glide [J] for each slip system and instance
   constitutive_titanmod_twinf0_PerTwinFam, &                                                     !< activation energy for glide [J] for each twin family and instance
   constitutive_titanmod_twinf0_PerTwinSys, &                                                     !< activation energy for glide [J] for each twin system and instance
   constitutive_titanmod_twinshearconstant_PerTwinFam, &                                          !< activation energy for glide [J] for each twin family and instance
   constitutive_titanmod_twinshearconstant_PerTwinSys, &                                          !< activation energy for glide [J] for each twin system and instance
   constitutive_titanmod_tau0e_PerSlipFam, &                                                      !< Initial yield stress for edge dislocations per slip family
   constitutive_titanmod_tau0e_PerSlipSys, &                                                      !< Initial yield stress for edge dislocations per slip system
   constitutive_titanmod_tau0s_PerSlipFam, &                                                      !< Initial yield stress for screw dislocations per slip family
   constitutive_titanmod_tau0s_PerSlipSys, &                                                      !< Initial yield stress for screw dislocations per slip system
   constitutive_titanmod_twintau0_PerTwinFam, &                                                   !< Initial yield stress for edge dislocations per twin family
   constitutive_titanmod_twintau0_PerTwinSys, &                                                   !< Initial yield stress for edge dislocations per twin system
   constitutive_titanmod_capre_PerSlipFam, &                                                      !< Capture radii for edge dislocations per slip family
   constitutive_titanmod_capre_PerSlipSys, &                                                      !< Capture radii for edge dislocations per slip system
   constitutive_titanmod_caprs_PerSlipFam, &                                                      !< Capture radii for screw dislocations per slip family
   constitutive_titanmod_caprs_PerSlipSys, &                                                      !< Capture radii for screw dislocations per slip system
   constitutive_titanmod_pe_PerSlipFam, &                                                         !< p-exponent in glide velocity
   constitutive_titanmod_ps_PerSlipFam, &                                                         !< p-exponent in glide velocity
   constitutive_titanmod_qe_PerSlipFam, &                                                         !< q-exponent in glide velocity
   constitutive_titanmod_qs_PerSlipFam, &                                                         !< q-exponent in glide velocity
   constitutive_titanmod_pe_PerSlipSys, &                                                         !< p-exponent in glide velocity
   constitutive_titanmod_ps_PerSlipSys, &                                                         !< p-exponent in glide velocity
   constitutive_titanmod_qe_PerSlipSys, &                                                         !< q-exponent in glide velocity
   constitutive_titanmod_qs_PerSlipSys, &                                                         !< q-exponent in glide velocity
   constitutive_titanmod_twinp_PerTwinFam, &                                                      !< p-exponent in glide velocity
   constitutive_titanmod_twinq_PerTwinFam, &                                                      !< q-exponent in glide velocity
   constitutive_titanmod_twinp_PerTwinSys, &                                                      !< p-exponent in glide velocity
   constitutive_titanmod_twinq_PerTwinSys, &                                                      !< p-exponent in glide velocity
   constitutive_titanmod_v0e_PerSlipFam, &                                                        !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_v0e_PerSlipSys, &                                                        !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_v0s_PerSlipFam, &                                                        !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_v0s_PerSlipSys, &                                                        !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_twingamma0_PerTwinFam, &                                                 !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_twingamma0_PerTwinSys, &                                                 !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_kinkcriticallength_PerSlipFam, &                                         !< screw dislocation mobility prefactor for kink-pairs per slip family
   constitutive_titanmod_kinkcriticallength_PerSlipSys, &                                         !< screw dislocation mobility prefactor for kink-pairs per slip system
   constitutive_titanmod_twinsizePerTwinFam, &                                                    !< twin thickness [m] for each twin family and instance
   constitutive_titanmod_twinsizePerTwinSys, &                                                    !< twin thickness [m] for each twin system and instance
   constitutive_titanmod_CeLambdaSlipPerSlipFam, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_CeLambdaSlipPerSlipSys, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_CsLambdaSlipPerSlipFam, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_CsLambdaSlipPerSlipSys, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_twinLambdaSlipPerTwinFam, &                                              !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_twinLambdaSlipPerTwinSys, &                                              !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_interactionSlipSlip, &                                                      !< coefficients for slip-slip interaction for each interaction type and instance
   constitutive_titanmod_interaction_ee, &                                                           !< coefficients for e-e interaction for each interaction type and instance
   constitutive_titanmod_interaction_ss, &                                                           !< coefficients for s-s interaction for each interaction type and instance
   constitutive_titanmod_interaction_es, &                                                           !< coefficients for e-s-twin interaction for each interaction type and instance
   constitutive_titanmod_interactionSlipTwin, &                                                      !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_titanmod_interactionTwinSlip, &                                                      !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_titanmod_interactionTwinTwin                                                         !< coefficients for twin-twin interaction for each interaction type and instance

  real(pReal), dimension(:,:,:),allocatable :: &
   constitutive_titanmod_Cslip_66, &                                                                 !< elasticity matrix in Mandel notation for each instance
   constitutive_titanmod_interactionMatrixSlipSlip, &                                                !< interaction matrix of the different slip systems for each instance
   constitutive_titanmod_interactionMatrix_ee, &                                                     !< interaction matrix of e-e for each instance
   constitutive_titanmod_interactionMatrix_ss, &                                                     !< interaction matrix of s-s for each instance
   constitutive_titanmod_interactionMatrix_es, &                                                     !< interaction matrix of e-s for each instance
   constitutive_titanmod_interactionMatrixSlipTwin, &                                                !< interaction matrix of slip systems with twin systems for each instance
   constitutive_titanmod_interactionMatrixTwinSlip, &                                                !< interaction matrix of twin systems with slip systems for each instance
   constitutive_titanmod_interactionMatrixTwinTwin, &                                                !< interaction matrix of the different twin systems for each instance                                                          
   constitutive_titanmod_forestProjectionEdge, &                                                     !< matrix of forest projections of edge dislocations for each instance  
   constitutive_titanmod_forestProjectionScrew, &                                                    !< matrix of forest projections of screw dislocations for each instance  
   constitutive_titanmod_TwinforestProjectionEdge, &                                                 !< matrix of forest projections of edge dislocations in twin system for each instance  
   constitutive_titanmod_TwinforestProjectionScrew                                                   !< matrix of forest projections of screw dislocations in twin system for each instance  

 real(pReal),       dimension(:,:,:,:),     allocatable :: &
  constitutive_titanmod_Ctwin_66                                                                    !< twin elasticity matrix in Mandel notation for each instance

 real(pReal),       dimension(:,:,:,:,:),   allocatable :: &
   constitutive_titanmod_Cslip_3333                                                                  !< elasticity matrix for each instance

 real(pReal),       dimension(:,:,:,:,:,:), allocatable :: &
   constitutive_titanmod_Ctwin_3333                                                                  !< twin elasticity matrix for each instance


 public :: &
   constitutive_titanmod_microstructure, &
   constitutive_titanmod_stateInit, &
   constitutive_titanmod_init, &
   constitutive_titanmod_LpAndItsTangent, &
   constitutive_titanmod_dotState, &
   constitutive_titanmod_deltaState, &
   constitutive_titanmod_dotTemperature, &
   constitutive_titanmod_postResults, &
   constitutive_titanmod_homogenizedC, &
   constitutive_titanmod_aTolState


 contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_titanmod_init(myFile)
 use, intrinsic :: iso_fortran_env                                        ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_Mandel3333to66,&
   math_Voigt66to3333,&
   math_mul3x3
 use IO
 use material 
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic

 use lattice
 
 implicit none
 integer(pInt), intent(in) :: myFile

 integer(pInt), parameter :: MAXNCHUNKS = 21_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt), dimension(6) :: configNchunks
 integer(pInt) :: section = 0_pInt,f,i,j,k,l,m,n,o,p,q,r,s,s1,s2,t,t1,t2,ns,nt,&
                  Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
                  Nchunks_SlipFamilies, Nchunks_TwinFamilies, &
                  mySize,myStructure,maxTotalNslip,maxTotalNtwin, &
                  maxNinstance
 character(len=65536) :: &
   tag  = '', &
   line = ''                                                                                        ! to start initialized
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//trim(CONSTITUTIVE_TITANMOD_label)//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == CONSTITUTIVE_TITANMOD_label),pInt)
 if (maxNinstance == 0) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 Nchunks_SlipFamilies = lattice_maxNslipFamily
 Nchunks_TwinFamilies = lattice_maxNtwinFamily
 Nchunks_SlipSlip =     lattice_maxNinteraction
 Nchunks_SlipTwin =     lattice_maxNinteraction
 Nchunks_TwinSlip =     lattice_maxNinteraction
 Nchunks_TwinTwin =     lattice_maxNinteraction


 allocate(constitutive_titanmod_sizeDotState(maxNinstance)) 
          constitutive_titanmod_sizeDotState = 0_pInt
 allocate(constitutive_titanmod_sizeState(maxNinstance)) 
          constitutive_titanmod_sizeState = 0_pInt
 allocate(constitutive_titanmod_sizePostResults(maxNinstance)) 
          constitutive_titanmod_sizePostResults = 0_pInt
 allocate(constitutive_titanmod_sizePostResult(maxval(phase_Noutput),maxNinstance))
          constitutive_titanmod_sizePostResult = 0_pInt
 allocate(constitutive_titanmod_output(maxval(phase_Noutput),maxNinstance))
          constitutive_titanmod_output = ''
 allocate(constitutive_titanmod_Noutput(maxNinstance))
          constitutive_titanmod_Noutput = 0_pInt
 
 allocate(constitutive_titanmod_structureName(maxNinstance)) 
          constitutive_titanmod_structureName = ''
 allocate(constitutive_titanmod_structure(maxNinstance))
          constitutive_titanmod_structure = 0_pInt
 allocate(constitutive_titanmod_Nslip(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_Nslip = 0_pInt
 allocate(constitutive_titanmod_Ntwin(lattice_maxNtwinFamily,maxNinstance))  
          constitutive_titanmod_Ntwin = 0_pInt
 allocate(constitutive_titanmod_slipFamily(lattice_maxNslip,maxNinstance))
          constitutive_titanmod_slipFamily = 0_pInt
 allocate(constitutive_titanmod_twinFamily(lattice_maxNtwin,maxNinstance))
          constitutive_titanmod_twinFamily = 0_pInt
 allocate(constitutive_titanmod_slipSystemLattice(lattice_maxNslip,maxNinstance))
          constitutive_titanmod_slipSystemLattice = 0_pInt
 allocate(constitutive_titanmod_twinSystemLattice(lattice_maxNtwin,maxNinstance)) 
          constitutive_titanmod_twinSystemLattice = 0_pInt
 allocate(constitutive_titanmod_totalNslip(maxNinstance))   
          constitutive_titanmod_totalNslip = 0_pInt
 allocate(constitutive_titanmod_totalNtwin(maxNinstance))   
          constitutive_titanmod_totalNtwin = 0_pInt
 allocate(constitutive_titanmod_CoverA(maxNinstance))
          constitutive_titanmod_CoverA = 0.0_pReal 
 allocate(constitutive_titanmod_debyefrequency(maxNinstance))
          constitutive_titanmod_debyefrequency = 0.0_pReal
 allocate(constitutive_titanmod_kinkf0(maxNinstance))
          constitutive_titanmod_kinkf0 = 0.0_pReal
 allocate(constitutive_titanmod_Gmod(maxNinstance))
          constitutive_titanmod_Gmod = 0.0_pReal
 allocate(constitutive_titanmod_CAtomicVolume(maxNinstance))
          constitutive_titanmod_CAtomicVolume = 0.0_pReal
 allocate(constitutive_titanmod_dc(maxNinstance))
          constitutive_titanmod_dc = 0.0_pReal
 allocate(constitutive_titanmod_twinhpconstant(maxNinstance))
          constitutive_titanmod_twinhpconstant = 0.0_pReal
 allocate(constitutive_titanmod_GrainSize(maxNinstance))
          constitutive_titanmod_GrainSize = 0.0_pReal
 allocate(constitutive_titanmod_MaxTwinFraction(maxNinstance))
          constitutive_titanmod_MaxTwinFraction = 0.0_pReal
 allocate(constitutive_titanmod_r(maxNinstance))
          constitutive_titanmod_r = 0.0_pReal
 allocate(constitutive_titanmod_CEdgeDipMinDistance(maxNinstance))
          constitutive_titanmod_CEdgeDipMinDistance = 0.0_pReal
 allocate(constitutive_titanmod_Cmfptwin(maxNinstance))
          constitutive_titanmod_Cmfptwin = 0.0_pReal
 allocate(constitutive_titanmod_Cthresholdtwin(maxNinstance))
          constitutive_titanmod_Cthresholdtwin = 0.0_pReal
 allocate(constitutive_titanmod_aTolRho(maxNinstance))
          constitutive_titanmod_aTolRho = 0.0_pReal
 allocate(constitutive_titanmod_Cslip_66(6,6,maxNinstance))
          constitutive_titanmod_Cslip_66 = 0.0_pReal
 allocate(constitutive_titanmod_Cslip_3333(3,3,3,3,maxNinstance))
          constitutive_titanmod_Cslip_3333 = 0.0_pReal
 allocate(constitutive_titanmod_rho_edge0(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_rho_edge0 = 0.0_pReal
 allocate(constitutive_titanmod_rho_screw0(lattice_maxNslipFamily,maxNinstance)) 
          constitutive_titanmod_rho_screw0 = 0.0_pReal
 allocate(constitutive_titanmod_shear_system0(lattice_maxNslipFamily,maxNinstance)) 
          constitutive_titanmod_shear_system0 = 0.0_pReal
 allocate(constitutive_titanmod_burgersPerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_burgersPerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_burgersPerTwinFam(lattice_maxNtwinFamily,maxNinstance))
          constitutive_titanmod_burgersPerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_f0_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_f0_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_tau0e_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_tau0e_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_tau0s_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_tau0s_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_capre_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_capre_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_caprs_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_caprs_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_pe_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_pe_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_ps_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_ps_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_qe_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_qe_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_qs_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_qs_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_v0e_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_v0e_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_v0s_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_v0s_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_kinkcriticallength_PerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_kinkcriticallength_PerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_twinsizePerTwinFam(lattice_maxNtwinFamily,maxNinstance))
          constitutive_titanmod_twinsizePerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_CeLambdaSlipPerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_CeLambdaSlipPerSlipFam = 0.0_pReal
 allocate(constitutive_titanmod_CsLambdaSlipPerSlipFam(lattice_maxNslipFamily,maxNinstance))
          constitutive_titanmod_CsLambdaSlipPerSlipFam = 0.0_pReal
 
 allocate(constitutive_titanmod_twinf0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twinf0_PerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twinshearconstant_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twinshearconstant_PerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twintau0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twintau0_PerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twinp_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twingamma0_PerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twinq_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twinLambdaSlipPerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twingamma0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twinp_PerTwinFam = 0.0_pReal
 allocate(constitutive_titanmod_twinLambdaSlipPerTwinFam(lattice_maxNTwinFamily,maxNinstance))
          constitutive_titanmod_twinq_PerTwinFam = 0.0_pReal
 
 allocate(constitutive_titanmod_interactionSlipSlip(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interactionSlipSlip = 0.0_pReal
 allocate(constitutive_titanmod_interaction_ee(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interaction_ee = 0.0_pReal
 allocate(constitutive_titanmod_interaction_ss(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interaction_ss = 0.0_pReal
 allocate(constitutive_titanmod_interaction_es(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interaction_ss = 0.0_pReal
 allocate(constitutive_titanmod_interactionSlipTwin(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interactionSlipTwin = 0.0_pReal
 allocate(constitutive_titanmod_interactionTwinSlip(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interactionTwinSlip = 0.0_pReal
 allocate(constitutive_titanmod_interactionTwinTwin(lattice_maxNinteraction,maxNinstance)) 
          constitutive_titanmod_interactionTwinTwin = 0.0_pReal
 
 rewind(myFile)
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= 'phase')                   ! wind forward to <phase>
   line = IO_read(myFile)
 enddo
 
 do while (trim(line) /= '#EOF#')                                                                   ! read through sections of phase part
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt                                                                     ! advance section counter
     cycle                                                                                          ! skip to next line
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (phase_plasticity(section) == CONSTITUTIVE_TITANMOD_label) then                             ! one of my sections
       i = phase_plasticityInstance(section)                                                        ! which instance of my plasticity is present phase
       positions = IO_stringPos(line,MAXNCHUNKS)
       tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                           ! extract key
       select case(tag)
         case ('plasticity','elasticity')
           cycle
         case ('(output)')
           constitutive_titanmod_Noutput(i) = constitutive_titanmod_Noutput(i) + 1_pInt
           constitutive_titanmod_output(constitutive_titanmod_Noutput(i),i) = IO_lc(IO_stringValue(line,positions,2_pInt))
         case ('lattice_structure')
           constitutive_titanmod_structureName(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
           configNchunks = lattice_configNchunks(constitutive_titanmod_structureName(i))
           Nchunks_SlipFamilies = configNchunks(1)
           Nchunks_TwinFamilies = configNchunks(2)
           Nchunks_SlipSlip =     configNchunks(3)
           Nchunks_SlipTwin =     configNchunks(4)
           Nchunks_TwinSlip =     configNchunks(5)
           Nchunks_TwinTwin =     configNchunks(6)
         case ('covera_ratio')
           constitutive_titanmod_CoverA(i)       = IO_floatValue(line,positions,2_pInt)
         case ('c11')
           constitutive_titanmod_Cslip_66(1,1,i) = IO_floatValue(line,positions,2_pInt)
         case ('c12')
           constitutive_titanmod_Cslip_66(1,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c13')
           constitutive_titanmod_Cslip_66(1,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c22')
           constitutive_titanmod_Cslip_66(2,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c23')
           constitutive_titanmod_Cslip_66(2,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c33')
           constitutive_titanmod_Cslip_66(3,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c44')
           constitutive_titanmod_Cslip_66(4,4,i) = IO_floatValue(line,positions,2_pInt)
         case ('c55')
           constitutive_titanmod_Cslip_66(5,5,i) = IO_floatValue(line,positions,2_pInt)
         case ('c66')
           constitutive_titanmod_Cslip_66(1,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('debyefrequency')
           constitutive_titanmod_debyefrequency(i) = IO_floatValue(line,positions,2_pInt)
         case ('kinkf0')
           constitutive_titanmod_kinkf0(i) = IO_floatValue(line,positions,2_pInt)
         case ('nslip')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_Nslip(j,i) = IO_intValue(line,positions,1_pInt+j)
           enddo
         case ('ntwin')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_Ntwin(j,i) = IO_intValue(line,positions,1_pInt+j)
           enddo
         case ('rho_edge0')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_rho_edge0(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('rho_screw0')
           do j = 1_pInt, Nchunks_SlipFamilies 
             constitutive_titanmod_rho_screw0(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('slipburgers')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_burgersPerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinburgers')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_burgersPerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('f0')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_f0_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinf0')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinf0_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('tau0e')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_tau0e_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twintau0')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twintau0_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('tau0s')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_tau0s_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('capre')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_capre_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('caprs')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_caprs_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('v0e')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_v0e_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twingamma0')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twingamma0_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('v0s')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_v0s_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('kinkcriticallength')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_kinkcriticallength_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinsize')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinsizePerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('celambdaslip')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_CeLambdaSlipPerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinlambdaslip')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinlambdaslipPerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('cslambdaslip')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_CsLambdaSlipPerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('grainsize')
           constitutive_titanmod_GrainSize(i) = IO_floatValue(line,positions,2_pInt)
         case ('maxtwinfraction')
           constitutive_titanmod_MaxTwinFraction(i) = IO_floatValue(line,positions,2_pInt)
         case ('pe')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_pe_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinp')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinp_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('ps')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_ps_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('qe')
           do j = 1_pInt, Nchunks_SlipFamilies 
             constitutive_titanmod_qe_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinq')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinq_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('qs')
           do j = 1_pInt, Nchunks_SlipFamilies
             constitutive_titanmod_qs_PerSlipFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('twinshearconstant')
           do j = 1_pInt, Nchunks_TwinFamilies
             constitutive_titanmod_twinshearconstant_PerTwinFam(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('dc')
           constitutive_titanmod_dc(i) = IO_floatValue(line,positions,2_pInt)
         case ('twinhpconstant')
           constitutive_titanmod_twinhpconstant(i) = IO_floatValue(line,positions,2_pInt)
         case ('atol_rho')
           constitutive_titanmod_aTolRho(i) = IO_floatValue(line,positions,2_pInt)
         case ('interactionslipslip')
           do j = 1_pInt, Nchunks_SlipSlip
             constitutive_titanmod_interactionSlipSlip(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactionee')
           do j = 1_pInt, lattice_maxNinteraction
             constitutive_titanmod_interaction_ee(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactionss')
           do j = 1_pInt, lattice_maxNinteraction
             constitutive_titanmod_interaction_ss(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactiones')
           do j = 1_pInt, lattice_maxNinteraction
             constitutive_titanmod_interaction_es(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactionsliptwin')
           do j = 1_pInt, Nchunks_SlipTwin
             constitutive_titanmod_interactionSlipTwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactiontwinslip')
           do j = 1_pInt, Nchunks_TwinSlip
             constitutive_titanmod_interactionTwinSlip(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case ('interactiontwintwin')
           do j = 1_pInt, Nchunks_TwinTwin
             constitutive_titanmod_interactionTwinTwin(j,i) = IO_floatValue(line,positions,1_pInt+j)
           enddo
         case default
           call IO_error(210_pInt,ext_msg=trim(tag)//' ('//CONSTITUTIVE_TITANMOD_label//')')
      end select
     endif
   endif
 enddo

 sanityChecks: do i = 1_pInt,maxNinstance
   constitutive_titanmod_structure(i) = &
   lattice_initializeStructure(constitutive_titanmod_structureName(i),constitutive_titanmod_CoverA(i))
   myStructure = constitutive_titanmod_structure(i)

   if (myStructure < 1_pInt)                                                 call IO_error(205_pInt,e=i)
   if (sum(constitutive_titanmod_Nslip(:,i)) <= 0_pInt)                      call IO_error(211_pInt,e=i,ext_msg='nslip (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
   if (sum(constitutive_titanmod_Ntwin(:,i)) < 0_pInt)                       call IO_error(211_pInt,e=i,ext_msg='ntwin (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
   do f = 1_pInt,lattice_maxNslipFamily
     if (constitutive_titanmod_Nslip(f,i) > 0_pInt) then   
       if (constitutive_titanmod_rho_edge0(f,i) < 0.0_pReal)                 call IO_error(211_pInt,e=i,ext_msg='rho_edge0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_rho_screw0(f,i) < 0.0_pReal)                call IO_error(211_pInt,e=i,ext_msg='rho_screw0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_burgersPerSlipFam(f,i) <= 0.0_pReal)     call IO_error(211_pInt,e=i,ext_msg='slipburgers (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_f0_PerSlipFam(f,i) <= 0.0_pReal)         call IO_error(211_pInt,e=i,ext_msg='f0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_tau0e_PerSlipFam(f,i) <= 0.0_pReal)      call IO_error(211_pInt,e=i,ext_msg='tau0e (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_tau0s_PerSlipFam(f,i) <= 0.0_pReal)      call IO_error(211_pInt,e=i,ext_msg='tau0s (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_capre_PerSlipFam(f,i) <= 0.0_pReal)      call IO_error(211_pInt,e=i,ext_msg='capre (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_caprs_PerSlipFam(f,i) <= 0.0_pReal)      call IO_error(211_pInt,e=i,ext_msg='caprs (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_v0e_PerSlipFam(f,i) <= 0.0_pReal)        call IO_error(211_pInt,e=i,ext_msg='v0e (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_v0s_PerSlipFam(f,i) <= 0.0_pReal)        call IO_error(211_pInt,e=i,ext_msg='v0s (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_kinkcriticallength_PerSlipFam(f,i) <= 0.0_pReal) &
                                                                      call IO_error(211_pInt,e=i,ext_msg='kinkCriticalLength (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
     endif
   enddo
   do f = 1_pInt,lattice_maxNtwinFamily
     if (constitutive_titanmod_Ntwin(f,i) > 0_pInt) then   
       if (constitutive_titanmod_burgersPerTwinFam(f,i) <= 0.0_pReal)     call IO_error(211_pInt,e=i,ext_msg='twinburgers (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_twinf0_PerTwinFam(f,i) <= 0.0_pReal)     call IO_error(211_pInt,e=i,ext_msg='twinf0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_twinshearconstant_PerTwinFam(f,i) <= 0.0_pReal) &
                                                                        call IO_error(211_pInt,e=i,ext_msg='twinshearconstant (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_twintau0_PerTwinFam(f,i) <= 0.0_pReal)   call IO_error(211_pInt,e=i,ext_msg='twintau0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
       if (constitutive_titanmod_twingamma0_PerTwinFam(f,i) <= 0.0_pReal) call IO_error(211_pInt,e=i,ext_msg='twingamma0 (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
     endif
   enddo
   if (constitutive_titanmod_dc(i) <= 0.0_pReal)                             call IO_error(211_pInt,e=i,ext_msg='dc (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
   if (constitutive_titanmod_twinhpconstant(i) <= 0.0_pReal)                 call IO_error(211_pInt,e=i,ext_msg='twinhpconstant (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
   if (constitutive_titanmod_aTolRho(i) <= 0.0_pReal)                        call IO_error(211_pInt,e=i,ext_msg='aTolRho (' &
                                                                                  //CONSTITUTIVE_TITANMOD_label//')')
   
   !* Determine total number of active slip or twin systems
   constitutive_titanmod_Nslip(:,i) = min(lattice_NslipSystem(:,myStructure),constitutive_titanmod_Nslip(:,i))
   constitutive_titanmod_Ntwin(:,i) = min(lattice_NtwinSystem(:,myStructure),constitutive_titanmod_Ntwin(:,i))
   constitutive_titanmod_totalNslip(i) = sum(constitutive_titanmod_Nslip(:,i))
   constitutive_titanmod_totalNtwin(i) = sum(constitutive_titanmod_Ntwin(:,i))
 enddo sanityChecks

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip = maxval(constitutive_titanmod_totalNslip)
 maxTotalNtwin = maxval(constitutive_titanmod_totalNtwin)
 
 allocate(constitutive_titanmod_burgersPerSlipSys(maxTotalNslip, maxNinstance))
 allocate(constitutive_titanmod_burgersPerTwinSys(maxTotalNtwin, maxNinstance))
          constitutive_titanmod_burgersPerTwinSys        = 0.0_pReal
 
 allocate(constitutive_titanmod_f0_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_f0_PerSlipSys            = 0.0_pReal
 allocate(constitutive_titanmod_tau0e_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_tau0e_PerSlipSys         = 0.0_pReal
 allocate(constitutive_titanmod_tau0s_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_tau0s_PerSlipSys         = 0.0_pReal
 allocate(constitutive_titanmod_capre_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_capre_PerSlipSys         = 0.0_pReal
 allocate(constitutive_titanmod_caprs_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_caprs_PerSlipSys         = 0.0_pReal
 allocate(constitutive_titanmod_pe_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_pe_PerSlipSys            = 0.0_pReal
 allocate(constitutive_titanmod_ps_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_ps_PerSlipSys            = 0.0_pReal
 allocate(constitutive_titanmod_qe_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_qe_PerSlipSys            = 0.0_pReal
 allocate(constitutive_titanmod_qs_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_qs_PerSlipSys            = 0.0_pReal
 allocate(constitutive_titanmod_v0e_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_v0e_PerSlipSys           = 0.0_pReal
 allocate(constitutive_titanmod_v0s_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_v0s_PerSlipSys           = 0.0_pReal
 allocate(constitutive_titanmod_kinkcriticallength_PerSlipSys(maxTotalNslip,maxNinstance))
          constitutive_titanmod_kinkcriticallength_PerSlipSys = 0.0_pReal
 allocate(constitutive_titanmod_CeLambdaSlipPerSlipSys(maxTotalNslip, maxNinstance))
          constitutive_titanmod_CeLambdaSlipPerSlipSys   = 0.0_pReal
 allocate(constitutive_titanmod_CsLambdaSlipPerSlipSys(maxTotalNslip, maxNinstance))
          constitutive_titanmod_CsLambdaSlipPerSlipSys   = 0.0_pReal
 
 allocate(constitutive_titanmod_twinf0_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twinf0_PerTwinSys        = 0.0_pReal
 allocate(constitutive_titanmod_twinshearconstant_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twinshearconstant_PerTwinSys = 0.0_pReal
 allocate(constitutive_titanmod_twintau0_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twintau0_PerTwinSys      = 0.0_pReal
 allocate(constitutive_titanmod_twinp_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twinp_PerTwinSys         = 0.0_pReal
 allocate(constitutive_titanmod_twinq_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twinq_PerTwinSys         = 0.0_pReal
 allocate(constitutive_titanmod_twingamma0_PerTwinSys(maxTotalNTwin,maxNinstance))
          constitutive_titanmod_twingamma0_PerTwinSys    = 0.0_pReal
 allocate(constitutive_titanmod_twinsizePerTwinSys(maxTotalNtwin, maxNinstance))
          constitutive_titanmod_twinsizePerTwinSys       = 0.0_pReal
 allocate(constitutive_titanmod_twinLambdaSlipPerTwinSys(maxTotalNtwin, maxNinstance))
          constitutive_titanmod_twinLambdaSlipPerTwinSys = 0.0_pReal

 allocate(constitutive_titanmod_interactionMatrixSlipSlip(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_interactionMatrixSlipSlip = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrix_ee(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_interactionMatrix_ee = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrix_ss(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_interactionMatrix_ss = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrix_es(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_interactionMatrix_es = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrixSlipTwin(maxTotalNslip,maxTotalNtwin,maxNinstance))
          constitutive_titanmod_interactionMatrixSlipTwin = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrixTwinSlip(maxTotalNtwin,maxTotalNslip,maxNinstance))
          constitutive_titanmod_interactionMatrixTwinSlip = 0.0_pReal
 allocate(constitutive_titanmod_interactionMatrixTwinTwin(maxTotalNtwin,maxTotalNtwin,maxNinstance))
          constitutive_titanmod_interactionMatrixTwinTwin = 0.0_pReal
 allocate(constitutive_titanmod_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_forestProjectionEdge      = 0.0_pReal
 allocate(constitutive_titanmod_forestProjectionScrew(maxTotalNslip,maxTotalNslip,maxNinstance))
          constitutive_titanmod_forestProjectionScrew     = 0.0_pReal
 allocate(constitutive_titanmod_TwinforestProjectionEdge(maxTotalNtwin,maxTotalNtwin,maxNinstance))
          constitutive_titanmod_TwinforestProjectionEdge  = 0.0_pReal
 allocate(constitutive_titanmod_TwinforestProjectionScrew(maxTotalNtwin,maxTotalNtwin,maxNinstance))
          constitutive_titanmod_TwinforestProjectionScrew = 0.0_pReal
 
 allocate(constitutive_titanmod_Ctwin_66(6,6,maxTotalNtwin,maxNinstance))
 allocate(constitutive_titanmod_Ctwin_3333(3,3,3,3,maxTotalNtwin,maxNinstance))
 constitutive_titanmod_Ctwin_66 = 0.0_pReal
 constitutive_titanmod_Ctwin_3333 = 0.0_pReal

 instancesLoop: do i = 1_pInt,maxNinstance
   myStructure = constitutive_titanmod_structure(i)

   !* Inverse lookup of slip system family
   l = 0_pInt
   do f = 1_pInt,lattice_maxNslipFamily
     do k = 1_pInt,constitutive_titanmod_Nslip(f,i)
       l = l + 1_pInt
       constitutive_titanmod_slipFamily(l,i) = f
       constitutive_titanmod_slipSystemLattice(l,i) = sum(lattice_NslipSystem(1:f-1_pInt,myStructure)) + k
   enddo; enddo

   !* Inverse lookup of twin system family
   l = 0_pInt
   do f = 1_pInt,lattice_maxNtwinFamily
     do k = 1_pInt,constitutive_titanmod_Ntwin(f,i)
       l = l + 1_pInt
       constitutive_titanmod_twinFamily(l,i) = f
       constitutive_titanmod_twinSystemLattice(l,i) = sum(lattice_NtwinSystem(1:f-1_pInt,myStructure)) + k
   enddo; enddo
   
   !* Determine size of state array  
   ns = constitutive_titanmod_totalNslip(i)
   nt = constitutive_titanmod_totalNtwin(i)
   constitutive_titanmod_sizeDotState(i) = &
   size(constitutive_titanmod_listBasicSlipStates)*ns+size(constitutive_titanmod_listBasicTwinStates)*nt
   constitutive_titanmod_sizeState(i) = &
   constitutive_titanmod_sizeDotState(i)+ &
   size(constitutive_titanmod_listDependentSlipStates)*ns+size(constitutive_titanmod_listDependentTwinStates)*nt

   !* Determine size of postResults array  
   
   do o = 1_pInt,constitutive_titanmod_Noutput(i)
     mySize = 0_pInt
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
            'shear_system')
          mySize = constitutive_titanmod_totalNslip(i)
        case('twin_fraction', &
             'gdot_twin', &
             'tau_twin' )
          mySize = constitutive_titanmod_totalNtwin(i)
        case('shear_basal', &                                                                       ! use only if all 4 slip families in hex are considered
             'shear_prism', &                                                                       ! use only if all 4 slip families in hex are considered
             'shear_pyra', &                                                                        ! use only if all 4 slip families in hex are considered
             'shear_pyrca', &                                                                       ! use only if all 4 slip families in hex are considered
             'rhoedge_basal', &
             'rhoedge_prism', &
             'rhoedge_pyra', &
             'rhoedge_pyrca', &
             'rhoscrew_basal', &
             'rhoscrew_prism', &
             'rhoscrew_pyra', &
             'rhoscrew_pyrca', &
             'shear_total')
          mySize = 1_pInt
        case default
          call IO_error(212_pInt,ext_msg=constitutive_titanmod_output(o,i)//' ('//CONSTITUTIVE_TITANMOD_label//')')
      end select

      if (mySize > 0_pInt) then  ! any meaningful output found
        constitutive_titanmod_sizePostResult(o,i) = mySize
        constitutive_titanmod_sizePostResults(i)  = constitutive_titanmod_sizePostResults(i) + mySize
      endif
   enddo
   
   !* Elasticity matrix and shear modulus according to material.config
   constitutive_titanmod_Cslip_66(:,:,i) = lattice_symmetrizeC66(constitutive_titanmod_structureName(i),&
                                                                 constitutive_titanmod_Cslip_66(:,:,i))  
   constitutive_titanmod_Gmod(i) = &
       0.2_pReal*(constitutive_titanmod_Cslip_66(1,1,i)-constitutive_titanmod_Cslip_66(1,2,i))&
     + 0.3_pReal*constitutive_titanmod_Cslip_66(4,4,i)
   constitutive_titanmod_Cslip_66(1:6,1:6,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_titanmod_Cslip_66(1:6,1:6,i)))
   constitutive_titanmod_Cslip_3333(1:3,1:3,1:3,1:3,i) = &
     math_Voigt66to3333(constitutive_titanmod_Cslip_66(1:6,1:6,i))
   
   !* Construction of the twin elasticity matrices
   do j=1_pInt,lattice_maxNtwinFamily
     do k=1_pInt,constitutive_titanmod_Ntwin(j,i)           
       do l=1_pInt,3_pInt ; do m=1_pInt,3_pInt ; do n=1_pInt,3_pInt ; do o=1_pInt,3_pInt
         do p=1_pInt,3_pInt ; do q=1_pInt,3_pInt ; do r=1_pInt,3_pInt ; do s=1_pInt,3_pInt
           constitutive_titanmod_Ctwin_3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1_pInt,i))+k,i) = &
             constitutive_titanmod_Ctwin_3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1_pInt,i))+k,i) + &
             constitutive_titanmod_Cslip_3333(p,q,r,s,i)*&
             lattice_Qtwin(l,p,sum(lattice_NslipSystem(1:j-1_pInt,myStructure))+k,myStructure)* &
             lattice_Qtwin(m,q,sum(lattice_NslipSystem(1:j-1_pInt,myStructure))+k,myStructure)* &
             lattice_Qtwin(n,r,sum(lattice_NslipSystem(1:j-1_pInt,myStructure))+k,myStructure)* &
             lattice_Qtwin(o,s,sum(lattice_NslipSystem(1:j-1_pInt,myStructure))+k,myStructure)
         enddo; enddo; enddo; enddo
       enddo; enddo; enddo ; enddo
       constitutive_titanmod_Ctwin_66(1:6,1:6,k,i) =  &
         math_Mandel3333to66(constitutive_titanmod_Ctwin_3333(1:3,1:3,1:3,1:3,k,i))
   enddo; enddo

   !* Burgers vector, dislocation velocity prefactor for each slip system 
   do s = 1_pInt,constitutive_titanmod_totalNslip(i)   
     f = constitutive_titanmod_slipFamily(s,i)    
     constitutive_titanmod_burgersPerSlipSys(s,i)      = constitutive_titanmod_burgersPerSlipFam(f,i)
     constitutive_titanmod_f0_PerSlipSys(s,i)          = constitutive_titanmod_f0_PerSlipFam(f,i)
     constitutive_titanmod_tau0e_PerSlipSys(s,i)       = constitutive_titanmod_tau0e_PerSlipFam(f,i)
     constitutive_titanmod_tau0s_PerSlipSys(s,i)       = constitutive_titanmod_tau0s_PerSlipFam(f,i)
     constitutive_titanmod_capre_PerSlipSys(s,i)       = constitutive_titanmod_capre_PerSlipFam(f,i)
     constitutive_titanmod_caprs_PerSlipSys(s,i)       = constitutive_titanmod_caprs_PerSlipFam(f,i)
     constitutive_titanmod_v0e_PerSlipSys(s,i)         = constitutive_titanmod_v0e_PerSlipFam(f,i)
     constitutive_titanmod_v0s_PerSlipSys(s,i)         = constitutive_titanmod_v0s_PerSlipFam(f,i)
     constitutive_titanmod_kinkcriticallength_PerSlipSys(s,i) = constitutive_titanmod_kinkcriticallength_PerSlipFam(f,i)
     constitutive_titanmod_pe_PerSlipSys(s,i)          = constitutive_titanmod_pe_PerSlipFam(f,i)
     constitutive_titanmod_ps_PerSlipSys(s,i)          = constitutive_titanmod_ps_PerSlipFam(f,i)
     constitutive_titanmod_qe_PerSlipSys(s,i)          = constitutive_titanmod_qe_PerSlipFam(f,i)
     constitutive_titanmod_qs_PerSlipSys(s,i)          = constitutive_titanmod_qs_PerSlipFam(f,i)
     constitutive_titanmod_CeLambdaSlipPerSlipSys(s,i) = constitutive_titanmod_CeLambdaSlipPerSlipFam(f,i)
     constitutive_titanmod_CsLambdaSlipPerSlipSys(s,i) = constitutive_titanmod_CsLambdaSlipPerSlipFam(f,i)
   enddo   
   
   !* Burgers vector, nucleation rate prefactor and twin size for each twin system 
   do t = 1_pInt,constitutive_titanmod_totalNtwin(i)   
     f = constitutive_titanmod_twinFamily(t,i)    
     constitutive_titanmod_burgersPerTwinSys(t,i)        = constitutive_titanmod_burgersPerTwinFam(f,i)
     constitutive_titanmod_twinsizePerTwinSys(t,i)       = constitutive_titanmod_twinsizePerTwinFam(f,i)
     constitutive_titanmod_twinf0_PerTwinSys(t,i)        = constitutive_titanmod_twinf0_PerTwinFam(f,i)
     constitutive_titanmod_twinshearconstant_PerTwinSys(t,i) = constitutive_titanmod_twinshearconstant_PerTwinFam(f,i)
     constitutive_titanmod_twintau0_PerTwinSys(t,i)      = constitutive_titanmod_twintau0_PerTwinFam(f,i)
     constitutive_titanmod_twingamma0_PerTwinSys(t,i)    = constitutive_titanmod_twingamma0_PerTwinFam(f,i)
     constitutive_titanmod_twinp_PerTwinSys(t,i)         = constitutive_titanmod_twinp_PerTwinFam(f,i)
     constitutive_titanmod_twinq_PerTwinSys(t,i)         = constitutive_titanmod_twinq_PerTwinFam(f,i)
     constitutive_titanmod_twinLambdaSlipPerTwinSys(t,i) = constitutive_titanmod_twinLambdaSlipPerTwinFam(f,i)
   enddo   
     
   !* Construction of interaction matrices
   do s1 = 1_pInt,constitutive_titanmod_totalNslip(i)
      do s2 = 1_pInt,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrixSlipSlip(s1,s2,i) = &
         constitutive_titanmod_interactionSlipSlip(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                               constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                               myStructure),i)
         constitutive_titanmod_interactionMatrix_ee(s1,s2,i) = &
         constitutive_titanmod_interaction_ee(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                          constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                          myStructure),i)
         constitutive_titanmod_interactionMatrix_ss(s1,s2,i) = &
         constitutive_titanmod_interaction_ss(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                          constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                          myStructure),i)
         constitutive_titanmod_interactionMatrix_es(s1,s2,i) = &
         constitutive_titanmod_interaction_es(lattice_interactionSlipSlip(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                          constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                          myStructure),i)
   enddo; enddo
   
   do s1 = 1_pInt,constitutive_titanmod_totalNslip(i)
      do t2 = 1_pInt,constitutive_titanmod_totalNtwin(i)     
         constitutive_titanmod_interactionMatrixSlipTwin(s1,t2,i) = &
         constitutive_titanmod_interactionSlipTwin(lattice_interactionSlipTwin(constitutive_titanmod_slipSystemLattice(s1,i), &
                                                                               constitutive_titanmod_twinSystemLattice(t2,i), &
                                                                               myStructure),i)         
   enddo; enddo
   
   do t1 = 1_pInt,constitutive_titanmod_totalNtwin(i)
      do s2 = 1_pInt,constitutive_titanmod_totalNslip(i)     
         constitutive_titanmod_interactionMatrixTwinSlip(t1,s2,i) = &
         constitutive_titanmod_interactionTwinSlip(lattice_interactionTwinSlip(constitutive_titanmod_twinSystemLattice(t1,i), &
                                                                               constitutive_titanmod_slipSystemLattice(s2,i), &
                                                                               myStructure),i)         
   enddo; enddo

   do t1 = 1_pInt,constitutive_titanmod_totalNtwin(i)
      do t2 = 1_pInt,constitutive_titanmod_totalNtwin(i)     
         constitutive_titanmod_interactionMatrixTwinTwin(t1,t2,i) = &
         constitutive_titanmod_interactionTwinTwin(lattice_interactionTwinTwin(constitutive_titanmod_twinSystemLattice(t1,i), &
                                                                               constitutive_titanmod_twinSystemLattice(t2,i), &
                                                                               myStructure),i)         
   enddo; enddo
   
   !* Calculation of forest projections for edge dislocations 
   do s1 = 1_pInt,constitutive_titanmod_totalNslip(i)
      do s2 = 1_pInt,constitutive_titanmod_totalNslip(i)      
         constitutive_titanmod_forestProjectionEdge(s1,s2,i) = &
         abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,i),myStructure), &
                         lattice_st(:,constitutive_titanmod_slipSystemLattice(s2,i),myStructure))) 
   !* Calculation of forest projections for screw dislocations 
         constitutive_titanmod_forestProjectionScrew(s1,s2,i) = &
         abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,i),myStructure), &
                         lattice_sd(:,constitutive_titanmod_slipSystemLattice(s2,i),myStructure))) 
   enddo; enddo
  

   !* Calculation of forest projections for edge dislocations in twin system
   do t1 = 1_pInt,constitutive_titanmod_totalNtwin(i)
      do t2 = 1_pInt,constitutive_titanmod_totalNtwin(i)      
         constitutive_titanmod_TwinforestProjectionEdge(t1,t2,i) = &
         abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,i),myStructure), &
                         lattice_tt(:,constitutive_titanmod_twinSystemLattice(t2,i),myStructure))) 
   !* Calculation of forest projections for screw dislocations in twin system
         constitutive_titanmod_TwinforestProjectionScrew(t1,t2,i) = &
         abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,i),myStructure), &
                         lattice_td(:,constitutive_titanmod_twinSystemLattice(t2,i),myStructure))) 
   enddo; enddo

 enddo instancesLoop

end subroutine constitutive_titanmod_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_stateInit(myInstance)
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily

 implicit none
 integer(pInt), intent(in) :: myInstance                                                            !< number specifying the instance of the plasticity

 real(pReal), dimension(constitutive_titanmod_sizeState(myInstance)) :: &
                              constitutive_titanmod_stateInit
 integer(pInt) :: s,s0,s1, &
                  t,t0,t1, & 
                  ns,nt,f
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

 
 !* Initialize basic slip state variables
 ! For slip
 s1 = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily
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
 t1 = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily
    t0 = t1 + 1_pInt
    t1 = t0 + constitutive_titanmod_Ntwin(f,myInstance) - 1_pInt 
    do t = t0,t1
       twingamma_dot0(t)=0.0_pReal
    enddo 
 enddo
 
 !* Initialize dependent slip microstructural variables
 forall (s = 1_pInt:ns)
   segment_edge0(s) = constitutive_titanmod_CeLambdaSlipPerSlipSys(s,myInstance)/ &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
     dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
   segment_screw0(s) = constitutive_titanmod_CsLambdaSlipPerSlipSys(s,myInstance)/ &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
     dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
   resistance_edge0(s) = &
     constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSys(s,myInstance)* &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+ &
     dot_product((rho_screw0),constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))
   resistance_screw0(s) = &
    constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSys(s,myInstance)* &
    sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+ &
    dot_product((rho_screw0), constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))
 end forall
 
 forall (t = 1_pInt:nt) &
   resistance_twin0(t) = 0.0_pReal

 constitutive_titanmod_stateInit = 0.0_pReal
 constitutive_titanmod_stateInit(1:ns)                             = rho_edge0
 constitutive_titanmod_stateInit(1_pInt*ns+1_pInt:2_pInt*ns)       = rho_screw0
 constitutive_titanmod_stateInit(2_pInt*ns+1_pInt:3_pInt*ns)       = shear_system0
 constitutive_titanmod_stateInit(3_pInt*ns+1_pInt:3_pInt*ns+nt)    = twingamma_dot0
 constitutive_titanmod_stateInit(3_pInt*ns+nt+1_pInt:4_pInt*ns+nt) = segment_edge0
 constitutive_titanmod_stateInit(4_pInt*ns+nt+1_pInt:5_pInt*ns+nt) = segment_screw0        
 constitutive_titanmod_stateInit(5_pInt*ns+nt+1_pInt:6_pInt*ns+nt) = resistance_edge0
 constitutive_titanmod_stateInit(6_pInt*ns+nt+1_pInt:7_pInt*ns+nt) = resistance_screw0
 constitutive_titanmod_stateInit(7_pInt*ns+nt+1_pInt:7_pInt*ns+2_pInt*nt)=resistance_twin0
 
end function constitutive_titanmod_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_aTolState(myInstance)

 implicit none
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(constitutive_titanmod_sizeState(myInstance)) :: constitutive_titanmod_aTolState
 
 constitutive_titanmod_aTolState = constitutive_titanmod_aTolRho(myInstance)
 
endfunction constitutive_titanmod_aTolState


pure function constitutive_titanmod_homogenizedC(state,ipc,ip,el)
use prec,     only: p_vec
use mesh,     only: mesh_NcpElems,mesh_maxNips
use material, only: homogenization_maxNgrains,material_phase,phase_plasticityInstance

implicit none
!* Input-Output variables
integer(pInt), intent(in) :: ipc,ip,el
type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
real(pReal), dimension(6,6) :: constitutive_titanmod_homogenizedC
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   volumefraction_PerTwinSys
!* Local variables 
 integer(pInt) myInstance,ns,nt,i
 real(pReal) sumf
  
 !* Shortened notation
 myInstance = phase_plasticityInstance(material_phase(ipc,ip,el))
 ns = constitutive_titanmod_totalNslip(myInstance)
 nt = constitutive_titanmod_totalNtwin(myInstance)
 
 !* Total twin volume fraction
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state(ipc,ip,el)%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,myInstance)
 enddo
 !sumf = sum(state(ipc,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 !* Homogenized elasticity matrix
 constitutive_titanmod_homogenizedC = (1.0_pReal-sumf)*constitutive_titanmod_Cslip_66(:,:,myInstance)
 do i=1_pInt,nt
    constitutive_titanmod_homogenizedC = &
 !   constitutive_titanmod_homogenizedC + state(ipc,ip,el)%p(6*ns+7*nt+i)*constitutive_titanmod_Ctwin_66(:,:,i,myInstance)
    constitutive_titanmod_homogenizedC + volumefraction_PerTwinSys(i)*constitutive_titanmod_Ctwin_66(:,:,i,myInstance)
 
 enddo 
 
end function constitutive_titanmod_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_titanmod_microstructure(temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase,&
   phase_plasticityInstance

 implicit none
 !* Input-Output variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: Temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: state
 !* Local variables
 integer(pInt) myInstance,myStructure,ns,nt,s,t,i
 real(pReal) sumf,sfe
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
              volumefraction_PerTwinSys
  
 !* Shortened notation
 myInstance = phase_plasticityInstance(material_phase(ipc,ip,el))
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
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state(ipc,ip,el)%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,myInstance)
 
 enddo
 
 !sumf = sum(state(ipc,ip,el)%p((6*ns+7*nt+1):(6*ns+8*nt))) ! safe for nt == 0
 
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 !* Stacking fault energy
 sfe = 0.0002_pReal*Temperature-0.0396_pReal

!--------------------------------------------------------------------------------------------------    
! average segment length for edge dislocations in matrix
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(3_pInt*ns+nt+s) = constitutive_titanmod_CeLambdaSlipPerSlipSys(s,myInstance)/ &
     sqrt(dot_product(state(ipc,ip,el)%p(1:ns), &
     constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
     dot_product(state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns), &
     constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
!--------------------------------------------------------------------------------------------------    
! average segment length for screw dislocations in matrix
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(4_pInt*ns+nt+s) = constitutive_titanmod_CsLambdaSlipPerSlipSys(s,myInstance)/ &
     sqrt(dot_product(state(ipc,ip,el)%p(1:ns), &
     constitutive_titanmod_forestProjectionEdge(1:ns,s,myInstance))+ &
     dot_product(state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns), &
     constitutive_titanmod_forestProjectionScrew(1:ns,s,myInstance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for edge dislocation motion
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(5_pInt*ns+nt+s) = &
     constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSys(s,myInstance)*&
     sqrt(dot_product((state(ipc,ip,el)%p(1:ns)),&
     constitutive_titanmod_interactionMatrix_ee(1:ns,s,myInstance))+ &
     dot_product((state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)),&
     constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for screw dislocation motion
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(6_pInt*ns+nt+s) = &
     constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerSlipSys(s,myInstance)*&
     sqrt(dot_product((state(ipc,ip,el)%p(1:ns)),&
     constitutive_titanmod_interactionMatrix_es(1:ns,s,myInstance))+ &
     dot_product((state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)),&
      constitutive_titanmod_interactionMatrix_ss(1:ns,s,myInstance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for dislocation motion in twin
 forall (t = 1_pInt:nt) &
   state(ipc,ip,el)%p(7_pInt*ns+nt+t) = &
     constitutive_titanmod_Gmod(myInstance)*constitutive_titanmod_burgersPerTwinSys(t,myInstance)*&
     (dot_product((abs(state(ipc,ip,el)%p(2_pInt*ns+1_pInt:2_pInt*ns+nt))),&
     constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,myInstance)))
 
end subroutine constitutive_titanmod_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine constitutive_titanmod_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,&
                                                                     temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   NnonSchmid
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3),                                                  intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9),                                                  intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                                                                  intent(in) :: &
   temperature                                                                                      !< temperature at IP 
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state
 integer(pInt) myInstance,myStructure,ns,nt,f,i,j,k,l,m,n,index_myFamily
 real(pReal) sumf,StressRatio_edge_p,minusStressRatio_edge_p,StressRatio_edge_pminus1,StressRatio_screw_p, &
         StressRatio_screw_pminus1, BoltzmannRatioedge, minusStressRatio_screw_p, &
         screwvelocity_prefactor,twinStressRatio_p,twinminusStressRatio_p,twinStressRatio_pminus1, &
         twinDotGamma0,BoltzmannRatioscrew,BoltzmannRatiotwin,bottomstress_edge,bottomstress_screw
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(constitutive_titanmod_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_slip,dgdot_dtauslip,tau_slip, edge_velocity, screw_velocity,gdot_slip_edge,gdot_slip_screw
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_twin,dgdot_dtautwin,tau_twin, volumefraction_PerTwinSys
 
 !* Shortened notation
 myInstance  = phase_plasticityInstance(material_phase(ipc,ip,el))
 myStructure = constitutive_titanmod_structure(myInstance) 
 ns = constitutive_titanmod_totalNslip(myInstance)
 nt = constitutive_titanmod_totalNtwin(myInstance)
 
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state(ipc,ip,el)%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,myInstance)
 
 enddo
 
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 
 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal
 
 !* Dislocation glide part
 gdot_slip = 0.0_pReal
 gdot_slip_edge = 0.0_pReal
 gdot_slip_screw = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,myStructure)) ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Nslip(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,myStructure)) 
      !*************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!      if(myStructure>=3.and.j>3) then ! for all non-basal slip systems
      if(myStructure==3_pInt) then ! only for prismatic and pyr <a> systems in hex
      screwvelocity_prefactor=constitutive_titanmod_debyefrequency(myInstance)* &
        state(ipc,ip,el)%p(4_pInt*ns+nt+j)*(constitutive_titanmod_burgersPerSlipSys(j,myInstance)/ &
        constitutive_titanmod_kinkcriticallength_PerSlipSys(j,myInstance))**2

     !* Stress ratio for screw ! No slip resistance for screw dislocations, only Peierls stress
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,myInstance)
         StressRatio_screw_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**constitutive_titanmod_ps_PerSlipSys(j,myInstance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif
        
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,myInstance)
      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**(constitutive_titanmod_ps_PerSlipSys(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_kinkf0(myInstance)/(kB*Temperature)

        else  ! if the structure is not hex or the slip family is basal
        screwvelocity_prefactor=constitutive_titanmod_v0s_PerSlipSys(j,myInstance)
        bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,myInstance)+state(ipc,ip,el)%p(6*ns+nt+j)
        StressRatio_screw_p = ((abs(tau_slip(j)))/( bottomstress_screw ))**constitutive_titanmod_ps_PerSlipSys(j,myInstance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif

      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/( bottomstress_screw))** &
      (constitutive_titanmod_ps_PerSlipSys(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_f0_PerSlipSys(j,myInstance)/(kB*Temperature)

        endif

        !* Stress ratio for edge
         bottomstress_edge=constitutive_titanmod_tau0e_PerSlipSys(j,myInstance)+state(ipc,ip,el)%p(5*ns+nt+j)
         StressRatio_edge_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_edge) &
        )**constitutive_titanmod_pe_PerSlipSys(j,myInstance)
                                
        if((1.0_pReal-StressRatio_edge_p)>0.001_pReal) then
        minusStressRatio_edge_p=1.0_pReal-StressRatio_edge_p
        else
        minusStressRatio_edge_p=0.001_pReal
        endif
        
      StressRatio_edge_pminus1 = ((abs(tau_slip(j)))/( bottomstress_edge))** &
      (constitutive_titanmod_pe_PerSlipSys(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio for edge. For screws it is defined above
      BoltzmannRatioedge = constitutive_titanmod_f0_PerSlipSys(j,myInstance)/(kB*Temperature)
            
        screw_velocity(j) =screwvelocity_prefactor * & ! there is no v0 for screw now because it is included in the prefactor
                exp(-BoltzmannRatioscrew*(minusStressRatio_screw_p)** &
        constitutive_titanmod_qs_PerSlipSys(j,myInstance))

        edge_velocity(j) =constitutive_titanmod_v0e_PerSlipSys(j,myInstance)*exp(-BoltzmannRatioedge* &
        (minusStressRatio_edge_p)** &
        constitutive_titanmod_qe_PerSlipSys(j,myInstance))

                !* Shear rates due to edge slip
       gdot_slip_edge(j) = constitutive_titanmod_burgersPerSlipSys(j,myInstance)*(state(ipc,ip,el)%p(j)* &
                edge_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !* Shear rates due to screw slip
       gdot_slip_screw(j) = constitutive_titanmod_burgersPerSlipSys(j,myInstance)*(state(ipc,ip,el)%p(ns+j) * &
                screw_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !Total shear rate
            
       gdot_slip(j) = gdot_slip_edge(j) + gdot_slip_screw(j)
                
       state(ipc,ip,el)%p(7*ns+2*nt+j)=edge_velocity(j)
       state(ipc,ip,el)%p(8*ns+2*nt+j)=screw_velocity(j)
       state(ipc,ip,el)%p(9*ns+2*nt+j)=tau_slip(j)
       state(ipc,ip,el)%p(10*ns+2*nt+j)=gdot_slip_edge(j)
       state(ipc,ip,el)%p(11*ns+2*nt+j)=gdot_slip_screw(j)
       state(ipc,ip,el)%p(12*ns+2*nt+j)=StressRatio_edge_p
       state(ipc,ip,el)%p(13*ns+2*nt+j)=StressRatio_screw_p
                
      !* Derivatives of shear rates
      dgdot_dtauslip(j) = constitutive_titanmod_burgersPerSlipSys(j,myInstance)*(( &
                ( &
                ( &
                ( &
                (edge_velocity(j)*state(ipc,ip,el)%p(j))) * &
                BoltzmannRatioedge*&
        constitutive_titanmod_pe_PerSlipSys(j,myInstance)* &
                constitutive_titanmod_qe_PerSlipSys(j,myInstance) &
                )/ &
                bottomstress_edge &
                )*&
        StressRatio_edge_pminus1*(minusStressRatio_edge_p)** &
                (constitutive_titanmod_qe_PerSlipSys(j,myInstance)-1.0_pReal) &
                ) + &
                ( &
                ( &
                ( &
                (state(ipc,ip,el)%p(ns+j) * screw_velocity(j)) * &
                BoltzmannRatioscrew* &
        constitutive_titanmod_ps_PerSlipSys(j,myInstance)* &
                constitutive_titanmod_qs_PerSlipSys(j,myInstance) &
                )/ &
                bottomstress_screw &
                )*&
        StressRatio_screw_pminus1*(minusStressRatio_screw_p)**(constitutive_titanmod_qs_PerSlipSys(j,myInstance)-1.0_pReal) &
                ) &
                ) !* sign(1.0_pReal,tau_slip(j))
                

                
!*************************************************                
!sumf=0.0_pReal
      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + (1.0_pReal - sumf)*gdot_slip(j)*lattice_Sslip(:,:,index_myFamily+i,myStructure)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
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
do f = 1_pInt,lattice_maxNtwinFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,myStructure)) ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Ntwin(f,myInstance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))        
     
!**************************************************************************************
      !* Stress ratios
!      StressRatio_r = (state(ipc,ip,el)%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_titanmod_r(myInstance)      
      
          !* Shear rates and their derivatives due to twin
!      if ( tau_twin(j) > 0.0_pReal ) !then          
!        gdot_twin(j) =  0.0_pReal!&
!          (constitutive_titanmod_MaxTwinFraction(myInstance)-sumf)*lattice_shearTwin(index_myFamily+i,myStructure)*&
!          state(ipc,ip,el)%p(6*ns+4*nt+j)*constitutive_titanmod_Ndot0PerTwinSys(f,myInstance)*exp(-StressRatio_r) 
!        dgdot_dtautwin(j) = ((gdot_twin(j)*constitutive_titanmod_r(myInstance))/tau_twin(j))*StressRatio_r
!      endif
!**************************************************************************************
   
     !* Stress ratio for edge
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,myInstance)+state(ipc,ip,el)%p(7*ns+nt+j)) &
        )**constitutive_titanmod_twinp_PerTwinSys(j,myInstance)
               
        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
              
      twinStressRatio_pminus1 = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,myInstance)+state(ipc,ip,el)%p(7*ns+nt+j)) &
        )**(constitutive_titanmod_twinp_PerTwinSys(j,myInstance)-1.0_pReal)

      !* Boltzmann ratio
      BoltzmannRatiotwin = constitutive_titanmod_twinf0_PerTwinSys(j,myInstance)/(kB*Temperature)

      !* Initial twin shear rates
      TwinDotGamma0 = &
        constitutive_titanmod_twingamma0_PerTwinSys(j,myInstance)

      !* Shear rates due to twin
         gdot_twin(j) =sign(1.0_pReal,tau_twin(j))*constitutive_titanmod_twingamma0_PerTwinSys(j,myInstance)* &
         exp(-BoltzmannRatiotwin*(twinminusStressRatio_p)**constitutive_titanmod_twinq_PerTwinSys(j,myInstance))
         
                                      
      !* Derivatives of shear rates in twin
      dgdot_dtautwin(j) = ( &
                ( &
                ( &
                (abs(gdot_twin(j))) * &
                BoltzmannRatiotwin*&
        constitutive_titanmod_twinp_PerTwinSys(j,myInstance)* &
                constitutive_titanmod_twinq_PerTwinSys(j,myInstance) &
                )/ &
                constitutive_titanmod_twintau0_PerTwinSys(j,myInstance) &
                )*&
        twinStressRatio_pminus1*(twinminusStressRatio_p)** &
                (constitutive_titanmod_twinq_PerTwinSys(j,myInstance)-1.0_pReal) &
                ) !* sign(1.0_pReal,tau_slip(j))
                
      !* Plastic velocity gradient for mechanical twinning                                                      
!      Lp = Lp + sumf*gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,myStructure)
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,myStructure)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,myStructure)*&
                                  lattice_Stwin(m,n,index_myFamily+i,myStructure)
   enddo
enddo

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

end subroutine constitutive_titanmod_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_titanmod_dotState(Tstar_v,Temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
use lattice,  only: lattice_maxNslipFamily,lattice_maxNtwinFamily, &
                    lattice_NslipSystem,lattice_NtwinSystem, lattice_Stwin_v

implicit none
 real(pReal), dimension(6),                                                    intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                                                                  intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state
real(pReal), dimension(constitutive_titanmod_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
  constitutive_titanmod_dotState

integer(pInt) MyInstance,MyStructure,ns,nt,f,i,j,index_myFamily
real(pReal) sumf,BoltzmannRatio,&
            twinStressRatio_p,twinminusStressRatio_p
real(pReal), dimension(constitutive_titanmod_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
DotRhoEdgeGeneration,DotRhoEdgeAnnihilation,DotRhoScrewAnnihilation,&
DotRhoScrewGeneration
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: gdot_twin, &
tau_twin, &
volumefraction_PerTwinSys
   
!* Shortened notation
myInstance  = phase_plasticityInstance(material_phase(ipc,ip,el))
MyStructure = constitutive_titanmod_structure(myInstance) 
ns = constitutive_titanmod_totalNslip(myInstance)
nt = constitutive_titanmod_totalNtwin(myInstance)

do i=1_pInt,nt
volumefraction_PerTwinSys(i)=state(ipc,ip,el)%p(3_pInt*ns+i)/ &
        constitutive_titanmod_twinshearconstant_PerTwinSys(i,myInstance)

enddo

sumf = sum(abs(volumefraction_PerTwinSys(1_pInt:nt))) ! safe for nt == 0

constitutive_titanmod_dotState = 0.0_pReal

    j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                             ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,myStructure))                 ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Nslip(f,myInstance)                        ! process each (active) slip system in family
     j = j+1_pInt

      !* Multiplication of edge dislocations
      DotRhoEdgeGeneration(j) = (state(ipc,ip,el)%p(ns+j)*state(ipc,ip,el)%p(8*ns+2*nt+j)/state(ipc,ip,el)%p(4*ns+nt+j))
      !* Multiplication of screw dislocations
      DotRhoScrewGeneration(j) = (state(ipc,ip,el)%p(j)*state(ipc,ip,el)%p(7*ns+2*nt+j)/state(ipc,ip,el)%p(3*ns+nt+j))

      !* Annihilation of edge dislocations
      DotRhoEdgeAnnihilation(j) = -((state(ipc,ip,el)%p(j))**2)* &
                constitutive_titanmod_capre_PerSlipSys(j,myInstance)*state(ipc,ip,el)%p(7*ns+2*nt+j)/2.0_pReal

      !* Annihilation of screw dislocations
      DotRhoScrewAnnihilation(j) = -((state(ipc,ip,el)%p(ns+j))**2)* &
                constitutive_titanmod_caprs_PerSlipSys(j,myInstance)*state(ipc,ip,el)%p(8*ns+2*nt+j)/2.0_pReal
       
      !* Edge dislocation density rate of change
      constitutive_titanmod_dotState(j) = &
        DotRhoEdgeGeneration(j)+DotRhoEdgeAnnihilation(j)

      !* Screw dislocation density rate of change
      constitutive_titanmod_dotState(ns+j) = &
        DotRhoScrewGeneration(j)+DotRhoScrewAnnihilation(j)

      constitutive_titanmod_dotState(2*ns+j) = &
                                  state(ipc,ip,el)%p(10*ns+2*nt+j)+state(ipc,ip,el)%p(11*ns+2*nt+j) ! sum of shear due to edge and screw 
    enddo
  enddo
  
!* Twin fraction evolution
j = 0_pInt
do f = 1_pInt,lattice_maxNtwinFamily                                 ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,MyStructure)) ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Ntwin(f,myInstance)          ! process each (active) twin system in family
      j = j+1_pInt

      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,myStructure))

     !* Stress ratio for edge
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,myInstance)+state(ipc,ip,el)%p(7*ns+nt+j)) &
        )**(constitutive_titanmod_twinp_PerTwinSys(j,myInstance))
        

        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
        
        !* Boltzmann ratio
      BoltzmannRatio = constitutive_titanmod_twinf0_PerTwinSys(j,myInstance)/(kB*Temperature)
       
            gdot_twin(j) =constitutive_titanmod_twingamma0_PerTwinSys(j,myInstance)*exp(-BoltzmannRatio* &
              (twinminusStressRatio_p)** &
            constitutive_titanmod_twinq_PerTwinSys(j,myInstance))*sign(1.0_pReal,tau_twin(j))
             
      constitutive_titanmod_dotState(3*ns+j)=gdot_twin(j)

    enddo
 enddo

end function constitutive_titanmod_dotState


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!> @details dummy function, returns 0.0
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_deltaState(Tstar_v,temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6),                                                    intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                                                                  intent(in) :: &
   Temperature                                                                                      !< temperature at integration point
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state
 
 real(pReal), dimension(constitutive_titanmod_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                             constitutive_titanmod_deltaState

 constitutive_titanmod_deltaState = 0.0_pReal
 
end function constitutive_titanmod_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of temperature
!> @details dummy function, returns 0.0
!--------------------------------------------------------------------------------------------------
real(pReal) pure function constitutive_titanmod_dotTemperature(Tstar_v,temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains
 
 implicit none
 real(pReal), dimension(6),                                                    intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                                                                  intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),                                                                intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
   state                                                                                            !< microstructure state

 constitutive_titanmod_dotTemperature = 0.0_pReal

end function constitutive_titanmod_dotTemperature


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_postResults(Tstar_v,Temperature,dt,state,ipc,ip,el)
 use prec,     only: pReal,pInt,p_vec
 use mesh,     only: mesh_NcpElems,mesh_maxNips
 use material, only: homogenization_maxNgrains,material_phase,phase_plasticityInstance,phase_Noutput
 
 implicit none
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt,Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt) myInstance,myStructure,ns,nt,o,i,c
 real(pReal) sumf
 real(pReal), dimension(constitutive_titanmod_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
 constitutive_titanmod_postResults
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
         volumefraction_PerTwinSys
 
 !* Shortened notation
 myInstance  = phase_plasticityInstance(material_phase(ipc,ip,el))
 myStructure = constitutive_titanmod_structure(myInstance) 
 ns = constitutive_titanmod_totalNslip(myInstance)
 nt = constitutive_titanmod_totalNtwin(myInstance)
 
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state(ipc,ip,el)%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,myInstance)
 enddo
 
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 
 !* Required output 
 c = 0_pInt
 constitutive_titanmod_postResults = 0.0_pReal

 do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_titanmod_output(o,myInstance))
     case ('rhoedge')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(1_pInt:ns)
       c = c + ns
     case ('rhoscrew')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)
       c = c + ns
     case ('segment_edge')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(3_pInt*ns+nt+1_pInt:4_pInt*ns+nt)
       c = c + ns
     case ('segment_screw')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(4_pInt*ns+nt+1_pInt:5_pInt*ns+nt)
       c = c + ns
     case ('resistance_edge')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(5_pInt*ns+nt+1_pInt:6_pInt*ns+nt)
       c = c + ns
     case ('resistance_screw')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(6_pInt*ns+nt+1_pInt:7_pInt*ns+nt)
       c = c + ns
     case ('velocity_edge')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(7*ns+2*nt+1:8*ns+2*nt)
       c = c + ns
     case ('velocity_screw')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(8*ns+2*nt+1:9*ns+2*nt)
       c = c + ns
     case ('tau_slip')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(9*ns+2*nt+1:10*ns+2*nt))
       c = c + ns
     case ('gdot_slip_edge')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(10*ns+2*nt+1:11*ns+2*nt))
       c = c + ns
     case ('gdot_slip_screw')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(11*ns+2*nt+1:12*ns+2*nt))
       c = c + ns
     case ('gdot_slip')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(10*ns+2*nt+1:11*ns+2*nt)) + &
                                                  abs(state(ipc,ip,el)%p(11*ns+2*nt+1:12*ns+2*nt))
       c = c + ns
     case ('stressratio_edge_p')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(12*ns+2*nt+1:13*ns+2*nt))
       c = c + ns
     case ('stressratio_screw_p')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(13*ns+2*nt+1:14*ns+2*nt))
       c = c + ns
     case ('shear_system')
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state(ipc,ip,el)%p(2*ns+1:3*ns))
       c = c + ns
     case ('shear_basal')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state(ipc,ip,el)%p(2*ns+1:2*ns+3)))
       c = c + 1_pInt
     case ('shear_prism')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state(ipc,ip,el)%p(2*ns+4:2*ns+6)))
       c = c + 1_pInt
     case ('shear_pyra')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state(ipc,ip,el)%p(2*ns+7:2*ns+12)))
       c = c + 1_pInt
     case ('shear_pyrca')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state(ipc,ip,el)%p(2*ns+13:2*ns+24)))
       c = c + 1_pInt

     case ('rhoedge_basal')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(1:3))
       c = c + 1_pInt
     case ('rhoedge_prism')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(4:6))
       c = c + 1_pInt
     case ('rhoedge_pyra')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(7:12))
       c = c + 1_pInt
     case ('rhoedge_pyrca')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(13:24))
       c = c + 1_pInt

     case ('rhoscrew_basal')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(ns+1:ns+3))
       c = c + 1_pInt
     case ('rhoscrew_prism')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(ns+4:ns+6))
       c = c + 1_pInt
     case ('rhoscrew_pyra')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(ns+7:ns+12))
       c = c + 1_pInt
     case ('rhoscrew_pyrca')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state(ipc,ip,el)%p(ns+13:ns+24))
       c = c + 1_pInt
     case ('shear_total')
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state(ipc,ip,el)%p(2*ns+1:3*ns)))
       c = c + 1_pInt
     case ('twin_fraction')
       constitutive_titanmod_postResults(c+1_pInt:c+nt) = abs(volumefraction_PerTwinSys(1:nt))
       c = c + nt
   end select
 enddo

end function constitutive_titanmod_postResults

end module constitutive_titanmod
