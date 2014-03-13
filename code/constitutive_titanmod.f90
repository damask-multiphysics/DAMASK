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
 character(len=18), dimension(3),          parameter,           private :: &
   CONSTITUTIVE_TITANMOD_listBasicSlipStates = & 
   ['rho_edge    ', 'rho_screw   ', 'shear_system']
 character(len=18), dimension(1),          parameter,           private :: &
   CONSTITUTIVE_TITANMOD_listBasicTwinStates = ['gdot_twin']
 character(len=19), dimension(11),         parameter,           private :: &
   CONSTITUTIVE_TITANMOD_listDependentSlipStates = &
   ['segment_edge       ', 'segment_screw      ', &
    'resistance_edge    ', 'resistance_screw   ', &
    'tau_slip           ', &
    'velocity_edge      ', 'velocity_screw     ', &
    'gdot_slip_edge     ', 'gdot_slip_screw    ', &
    'stressratio_edge_p ', 'stressratio_screw_p'  ]
 character(len=18), dimension(2),          parameter,           private :: &
   constitutive_titanmod_listDependentTwinStates = &
   ['twin_fraction', 'tau_twin     ']
 real(pReal),                              parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin


 integer(pInt), dimension(:), allocatable, public, protected :: &
   constitutive_titanmod_sizeState, &                                                               !< total number of microstructural state variables
   constitutive_titanmod_sizeDotState, &                                                            !< number of dotStates
   constitutive_titanmod_sizePostResults                                                            !<  cumulative size of post results

 integer(pInt), dimension(:,:), allocatable, target, public :: &
   constitutive_titanmod_sizePostResult                                                             !<  size of each post result output

 character(len=64), dimension(:,:), allocatable, target, public :: &
   constitutive_titanmod_output                                                                     !<  name of each post result output

 integer(pInt),     dimension(:),          allocatable,         private :: & 
   constitutive_titanmod_Noutput                                                                    !<  number of outputs per instance of this plasticity                                                            !< ID of the lattice structure

 integer(pInt),     dimension(:),          allocatable,         private :: & 
   constitutive_titanmod_totalNslip, &                                                              !<  total number of active slip systems for each instance
   constitutive_titanmod_totalNtwin                                                                 !<  total number of active twin systems for each instance

 integer(pInt),     dimension(:,:),        allocatable,         private :: &
   constitutive_titanmod_Nslip, &                                                                   !<  number of active slip systems for each family and instance
   constitutive_titanmod_Ntwin, &                                                                   !< number of active twin systems for each family and instance
   constitutive_titanmod_slipFamily, &                                                              !< lookup table relating active slip system to slip family for each instance
   constitutive_titanmod_twinFamily, &                                                              !< lookup table relating active twin system to twin family for each instance
   constitutive_titanmod_slipSystemLattice, &                                                       !< lookup table relating active slip system index to lattice slip system index for each instance
   constitutive_titanmod_twinSystemLattice                                                          !< lookup table relating active twin system index to lattice twin system index for each instance

 real(pReal),       dimension(:),          allocatable,         private :: &
   constitutive_titanmod_debyefrequency, &                                                          !< Debye frequency
   constitutive_titanmod_kinkf0, &                                                                  !<
   constitutive_titanmod_CAtomicVolume, &                                                           !< atomic volume in Bugers vector unit
   constitutive_titanmod_dc, &                                                                      !< prefactor for self-diffusion coefficient
   constitutive_titanmod_twinhpconstant, &                                                          !< activation energy for dislocation climb
   constitutive_titanmod_GrainSize, &                                                               !< grain size - Not being used
   constitutive_titanmod_MaxTwinFraction, &                                                         !< maximum allowed total twin volume fraction
   constitutive_titanmod_r, &                                                                       !< r-exponent in twin nucleation rate
   constitutive_titanmod_CEdgeDipMinDistance, &                                                     !< Not being used
   constitutive_titanmod_Cmfptwin, &                                                                !< Not being used
   constitutive_titanmod_Cthresholdtwin, &                                                          !< Not being used
   constitutive_titanmod_aTolRho                                                                    !< absolute tolerance for integration of dislocation density

 real(pReal),       dimension(:,:),        allocatable,         private :: &
   constitutive_titanmod_rho_edge0, &                                                               !< initial edge dislocation density per slip system for each family and instance
   constitutive_titanmod_rho_screw0, &                                                              !< initial screw dislocation density per slip system for each family and instance
   constitutive_titanmod_shear_system0, &                                                           !< accumulated shear on each system
   constitutive_titanmod_burgersPerSlipFam, &                                                       !< absolute length of burgers vector [m] for each slip family and instance
   constitutive_titanmod_burgersPerSlipSys, &                                                       !< absolute length of burgers vector [m] for each slip system and instance
   constitutive_titanmod_burgersPerTwinFam, &                                                       !< absolute length of burgers vector [m] for each twin family and instance
   constitutive_titanmod_burgersPerTwinSys, &                                                       !< absolute length of burgers vector [m] for each twin system and instance
   constitutive_titanmod_f0_PerSlipFam, &                                                           !< activation energy for glide [J] for each slip family and instance
   constitutive_titanmod_f0_PerSlipSys, &                                                           !< activation energy for glide [J] for each slip system and instance
   constitutive_titanmod_twinf0_PerTwinFam, &                                                       !< activation energy for glide [J] for each twin family and instance
   constitutive_titanmod_twinf0_PerTwinSys, &                                                       !< activation energy for glide [J] for each twin system and instance
   constitutive_titanmod_twinshearconstant_PerTwinFam, &                                            !< activation energy for glide [J] for each twin family and instance
   constitutive_titanmod_twinshearconstant_PerTwinSys, &                                            !< activation energy for glide [J] for each twin system and instance
   constitutive_titanmod_tau0e_PerSlipFam, &                                                        !< Initial yield stress for edge dislocations per slip family
   constitutive_titanmod_tau0e_PerSlipSys, &                                                        !< Initial yield stress for edge dislocations per slip system
   constitutive_titanmod_tau0s_PerSlipFam, &                                                        !< Initial yield stress for screw dislocations per slip family
   constitutive_titanmod_tau0s_PerSlipSys, &                                                        !< Initial yield stress for screw dislocations per slip system
   constitutive_titanmod_twintau0_PerTwinFam, &                                                     !< Initial yield stress for edge dislocations per twin family
   constitutive_titanmod_twintau0_PerTwinSys, &                                                     !< Initial yield stress for edge dislocations per twin system
   constitutive_titanmod_capre_PerSlipFam, &                                                        !< Capture radii for edge dislocations per slip family
   constitutive_titanmod_capre_PerSlipSys, &                                                        !< Capture radii for edge dislocations per slip system
   constitutive_titanmod_caprs_PerSlipFam, &                                                        !< Capture radii for screw dislocations per slip family
   constitutive_titanmod_caprs_PerSlipSys, &                                                        !< Capture radii for screw dislocations per slip system
   constitutive_titanmod_pe_PerSlipFam, &                                                           !< p-exponent in glide velocity
   constitutive_titanmod_ps_PerSlipFam, &                                                           !< p-exponent in glide velocity
   constitutive_titanmod_qe_PerSlipFam, &                                                           !< q-exponent in glide velocity
   constitutive_titanmod_qs_PerSlipFam, &                                                           !< q-exponent in glide velocity
   constitutive_titanmod_pe_PerSlipSys, &                                                           !< p-exponent in glide velocity
   constitutive_titanmod_ps_PerSlipSys, &                                                           !< p-exponent in glide velocity
   constitutive_titanmod_qe_PerSlipSys, &                                                           !< q-exponent in glide velocity
   constitutive_titanmod_qs_PerSlipSys, &                                                           !< q-exponent in glide velocity
   constitutive_titanmod_twinp_PerTwinFam, &                                                        !< p-exponent in glide velocity
   constitutive_titanmod_twinq_PerTwinFam, &                                                        !< q-exponent in glide velocity
   constitutive_titanmod_twinp_PerTwinSys, &                                                        !< p-exponent in glide velocity
   constitutive_titanmod_twinq_PerTwinSys, &                                                        !< p-exponent in glide velocity
   constitutive_titanmod_v0e_PerSlipFam, &                                                          !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_v0e_PerSlipSys, &                                                          !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_v0s_PerSlipFam, &                                                          !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_v0s_PerSlipSys, &                                                          !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_twingamma0_PerTwinFam, &                                                   !< edge dislocation velocity prefactor [m/s] for each family and instance
   constitutive_titanmod_twingamma0_PerTwinSys, &                                                   !< screw dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_titanmod_kinkcriticallength_PerSlipFam, &                                           !< screw dislocation mobility prefactor for kink-pairs per slip family
   constitutive_titanmod_kinkcriticallength_PerSlipSys, &                                           !< screw dislocation mobility prefactor for kink-pairs per slip system
   constitutive_titanmod_twinsizePerTwinFam, &                                                      !< twin thickness [m] for each twin family and instance
   constitutive_titanmod_twinsizePerTwinSys, &                                                      !< twin thickness [m] for each twin system and instance
   constitutive_titanmod_CeLambdaSlipPerSlipFam, &                                                  !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_CeLambdaSlipPerSlipSys, &                                                  !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_CsLambdaSlipPerSlipFam, &                                                  !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_CsLambdaSlipPerSlipSys, &                                                  !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_twinLambdaSlipPerTwinFam, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_titanmod_twinLambdaSlipPerTwinSys, &                                                !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_titanmod_interactionSlipSlip, &                                                     !< coefficients for slip-slip interaction for each interaction type and instance
   constitutive_titanmod_interaction_ee, &                                                          !< coefficients for e-e interaction for each interaction type and instance
   constitutive_titanmod_interaction_ss, &                                                          !< coefficients for s-s interaction for each interaction type and instance
   constitutive_titanmod_interaction_es, &                                                          !< coefficients for e-s-twin interaction for each interaction type and instance
   constitutive_titanmod_interactionSlipTwin, &                                                     !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_titanmod_interactionTwinSlip, &                                                     !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_titanmod_interactionTwinTwin                                                        !< coefficients for twin-twin interaction for each interaction type and instance

 real(pReal),       dimension(:,:,:),      allocatable,         private :: &
   constitutive_titanmod_interactionMatrixSlipSlip, &                                               !< interaction matrix of the different slip systems for each instance
   constitutive_titanmod_interactionMatrix_ee, &                                                    !< interaction matrix of e-e for each instance
   constitutive_titanmod_interactionMatrix_ss, &                                                    !< interaction matrix of s-s for each instance
   constitutive_titanmod_interactionMatrix_es, &                                                    !< interaction matrix of e-s for each instance
   constitutive_titanmod_interactionMatrixSlipTwin, &                                               !< interaction matrix of slip systems with twin systems for each instance
   constitutive_titanmod_interactionMatrixTwinSlip, &                                               !< interaction matrix of twin systems with slip systems for each instance
   constitutive_titanmod_interactionMatrixTwinTwin, &                                               !< interaction matrix of the different twin systems for each instance                                                          
   constitutive_titanmod_forestProjectionEdge, &                                                    !< matrix of forest projections of edge dislocations for each instance  
   constitutive_titanmod_forestProjectionScrew, &                                                   !< matrix of forest projections of screw dislocations for each instance  
   constitutive_titanmod_TwinforestProjectionEdge, &                                                !< matrix of forest projections of edge dislocations in twin system for each instance  
   constitutive_titanmod_TwinforestProjectionScrew                                                  !< matrix of forest projections of screw dislocations in twin system for each instance  

 real(pReal),      dimension(:,:,:,:),     allocatable,         private :: &
  constitutive_titanmod_Ctwin66                                                                    !< twin elasticity matrix in Mandel notation for each instance

 real(pReal),      dimension(:,:,:,:,:,:), allocatable,         private :: &
   constitutive_titanmod_Ctwin3333                                                                 !< twin elasticity matrix for each instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 rhoedge_ID,            rhoscrew_ID, &
                 segment_edge_ID,       segment_screw_ID, &
                 resistance_edge_ID,    resistance_screw_ID, &
                 velocity_edge_ID,      velocity_screw_ID, &
                 tau_slip_ID, &
                 gdot_slip_edge_ID,     gdot_slip_screw_ID, &
                 gdot_slip_ID, &
                 stressratio_edge_p_ID, stressratio_screw_p_ID, &
                 shear_system_ID, &
                 twin_fraction_ID, &
                 shear_basal_ID,    shear_prism_ID,    shear_pyra_ID,    shear_pyrca_ID, &
                 rhoedge_basal_ID,  rhoedge_prism_ID,  rhoedge_pyra_ID,  rhoedge_pyrca_ID, &
                 rhoscrew_basal_ID, rhoscrew_prism_ID, rhoscrew_pyra_ID, rhoscrew_pyrca_ID, &
                 shear_total_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   constitutive_titanmod_outputID                                                                     !< ID of each post result output

 public :: &
   constitutive_titanmod_microstructure, &
   constitutive_titanmod_stateInit, &
   constitutive_titanmod_init, &
   constitutive_titanmod_LpAndItsTangent, &
   constitutive_titanmod_dotState, &
   constitutive_titanmod_postResults, &
   constitutive_titanmod_homogenizedC, &
   constitutive_titanmod_aTolState


 contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_titanmod_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                        ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66,&
   math_Voigt66to3333,&
   math_mul3x3
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   homogenization_maxNgrains, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_TITANMOD_label, &
   PLASTICITY_TITANMOD_ID, &
   MATERIAL_partPhase
 use lattice
 
 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = LATTICE_maxNinteraction + 1_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt) :: &
   phase, &
   instance, j, k, l, m, n, p, q, r, &
   f, o, &
   s, s1, s2, &
   t, t1, t2, &
   ns, nt, &
   Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
   Nchunks_SlipFamilies, Nchunks_TwinFamilies, &
   mySize, &
   maxTotalNslip,maxTotalNtwin, maxNinstance
 character(len=65536) :: &
   tag  = '', &
   line = ''  
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_TITANMOD_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == PLASTICITY_TITANMOD_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(constitutive_titanmod_sizeDotState(maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_sizeState(maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_sizePostResults(maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_sizePostResult(maxval(phase_Noutput),maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_output(maxval(phase_Noutput),maxNinstance))
          constitutive_titanmod_output = ''
 allocate(constitutive_titanmod_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(constitutive_titanmod_Noutput(maxNinstance), source=0_pInt) 
 
 allocate(constitutive_titanmod_Nslip(lattice_maxNslipFamily,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_Ntwin(lattice_maxNtwinFamily,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_slipFamily(lattice_maxNslip,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_twinFamily(lattice_maxNtwin,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_slipSystemLattice(lattice_maxNslip,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_twinSystemLattice(lattice_maxNtwin,maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_totalNslip(maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_totalNtwin(maxNinstance), source=0_pInt) 
 allocate(constitutive_titanmod_debyefrequency(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_kinkf0(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CAtomicVolume(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_dc(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinhpconstant(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_GrainSize(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_MaxTwinFraction(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_r(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CEdgeDipMinDistance(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_Cmfptwin(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_Cthresholdtwin(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_aTolRho(maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_rho_edge0(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_rho_screw0(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_shear_system0(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_burgersPerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_burgersPerTwinFam(lattice_maxNtwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_f0_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_tau0e_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_tau0s_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_capre_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_caprs_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_pe_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_ps_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_qe_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_qs_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_v0e_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_v0s_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_kinkcriticallength_PerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinsizePerTwinFam(lattice_maxNtwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CeLambdaSlipPerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CsLambdaSlipPerSlipFam(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 
 allocate(constitutive_titanmod_twinf0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinshearconstant_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twintau0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinp_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinq_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twingamma0_PerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinLambdaSlipPerTwinFam(lattice_maxNTwinFamily,maxNinstance), source=0.0_pReal) 
 
 allocate(constitutive_titanmod_interactionSlipSlip(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interaction_ee(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interaction_ss(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interaction_es(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionSlipTwin(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionTwinSlip(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionTwinTwin(lattice_maxNinteraction,maxNinstance), source=0.0_pReal) 
 
 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                                    ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     phase = phase + 1_pInt                                                                         ! advance section counter
     if (phase_plasticity(phase) == PLASTICITY_TITANMOD_ID) then
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_TwinFamilies = count(lattice_NtwinSystem(:,phase) > 0_pInt)
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_SlipTwin =     maxval(lattice_interactionSlipTwin(:,:,phase))
       Nchunks_TwinSlip =     maxval(lattice_interactionTwinSlip(:,:,phase))
       Nchunks_TwinTwin =     maxval(lattice_interactionTwinTwin(:,:,phase))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_TITANMOD_ID) then           ! one of my sections. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('plasticity','elasticity','lattice_structure', &
           'covera_ratio','c/a_ratio','c/a', &
           'c11','c12','c13','c22','c23','c33','c44','c55','c66')
       case ('(output)')
         constitutive_titanmod_Noutput(instance) = constitutive_titanmod_Noutput(instance) + 1_pInt
         constitutive_titanmod_output(constitutive_titanmod_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('rhoedge')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoedge_ID
           case ('rhoscrew')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoscrew_ID
           case ('segment_edge')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = segment_edge_ID
           case ('segment_screw')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = segment_screw_ID
           case ('resistance_edge')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = resistance_edge_ID
           case ('resistance_screw')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = resistance_screw_ID
           case ('velocity_edge')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = velocity_edge_ID
           case ('velocity_screw')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = velocity_screw_ID
           case ('tau_slip')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = tau_slip_ID
           case ('gdot_slip_edge')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = gdot_slip_edge_ID
           case ('gdot_slip_screw')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = gdot_slip_screw_ID
           case ('gdot_slip')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = gdot_slip_ID
           case ('stressratio_edge_p')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = stressratio_edge_p_ID
           case ('stressratio_screw_p')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = stressratio_screw_p_ID
           case ('shear_system')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_system_ID
           case ('twin_fraction')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = twin_fraction_ID
           case ('shear_basal')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_basal_ID
           case ('shear_prism')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_prism_ID
           case ('shear_pyra')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_pyra_ID
           case ('shear_pyrca')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_pyrca_ID
           case ('rhoedge_basal')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoedge_basal_ID
           case ('rhoedge_prism')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoedge_prism_ID
           case ('rhoedge_pyra')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoedge_pyra_ID
           case ('rhoedge_pyrca')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoedge_pyrca_ID
           case ('rhoscrew_basal')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoscrew_basal_ID
           case ('rhoscrew_prism')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoscrew_prism_ID
           case ('rhoscrew_pyra')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoscrew_pyra_ID
           case ('rhoscrew_pyrca')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = rhoscrew_pyrca_ID
           case ('shear_total')
             constitutive_titanmod_outputID(constitutive_titanmod_Noutput(instance),instance) = shear_total_ID
           case default
             call IO_error(105_pInt,ext_msg=IO_stringValue(line,positions,2_pInt)//' ('//PLASTICITY_TITANMOD_label//')')
         end select
       case ('debyefrequency')
         constitutive_titanmod_debyefrequency(instance) = IO_floatValue(line,positions,2_pInt)
       case ('kinkf0')
         constitutive_titanmod_kinkf0(instance) = IO_floatValue(line,positions,2_pInt)
       case ('nslip')
         if (positions(1) < 1_pInt + Nchunks_SlipFamilies) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('ntwin')
         if (positions(1) < 1_pInt + Nchunks_TwinFamilies) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_Ntwin(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('rho_edge0')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_rho_edge0(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('rho_screw0')
         do j = 1_pInt, Nchunks_SlipFamilies 
           constitutive_titanmod_rho_screw0(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('slipburgers')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_burgersPerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinburgers')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_burgersPerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('f0')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_f0_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinf0')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinf0_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('tau0e')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_tau0e_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twintau0')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twintau0_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('tau0s')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_tau0s_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('capre')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_capre_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('caprs')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_caprs_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('v0e')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_v0e_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twingamma0')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twingamma0_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('v0s')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_v0s_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('kinkcriticallength')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_kinkcriticallength_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinsize')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinsizePerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('celambdaslip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_CeLambdaSlipPerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinlambdaslip')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinlambdaslipPerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('cslambdaslip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_CsLambdaSlipPerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('grainsize')
         constitutive_titanmod_GrainSize(instance) = IO_floatValue(line,positions,2_pInt)
       case ('maxtwinfraction')
         constitutive_titanmod_MaxTwinFraction(instance) = IO_floatValue(line,positions,2_pInt)
       case ('pe')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_pe_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinp')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinp_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('ps')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_ps_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('qe')
         do j = 1_pInt, Nchunks_SlipFamilies 
           constitutive_titanmod_qe_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinq')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinq_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('qs')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_titanmod_qs_PerSlipFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinshearconstant')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_titanmod_twinshearconstant_PerTwinFam(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('dc')
         constitutive_titanmod_dc(instance) = IO_floatValue(line,positions,2_pInt)
       case ('twinhpconstant')
         constitutive_titanmod_twinhpconstant(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_rho')
         constitutive_titanmod_aTolRho(instance) = IO_floatValue(line,positions,2_pInt)
       case ('interactionee')
         do j = 1_pInt, lattice_maxNinteraction
           constitutive_titanmod_interaction_ee(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interactionss')
         do j = 1_pInt, lattice_maxNinteraction
           constitutive_titanmod_interaction_ss(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interactiones')
         do j = 1_pInt, lattice_maxNinteraction
           constitutive_titanmod_interaction_es(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_slipslip','interactionslipslip')
         if (positions(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           constitutive_titanmod_interactionSlipSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_sliptwin','interactionsliptwin')
         if (positions(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           constitutive_titanmod_interactionSlipTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twinslip','interactiontwinslip')
         if (positions(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           constitutive_titanmod_interactionTwinSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twintwin','interactiontwintwin')
         if (positions(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           constitutive_titanmod_interactionTwinTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case default
         call IO_error(210_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_TITANMOD_label//')')
     end select
   endif; endif
 enddo parsingFile
 
 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
   myPhase: if (phase_plasticity(phase) == PLASTICITY_TITANMOD_ID) then
     instance = phase_plasticityInstance(phase) 
     if (sum(constitutive_titanmod_Nslip(:,instance)) <= 0_pInt) &
       call IO_error(211_pInt,el=instance,ext_msg='nslip ('//PLASTICITY_TITANMOD_label//')')
     if (sum(constitutive_titanmod_Ntwin(:,instance)) < 0_pInt) &
       call IO_error(211_pInt,el=instance,ext_msg='ntwin ('//PLASTICITY_TITANMOD_label//')')
     do f = 1_pInt,lattice_maxNslipFamily
       if (constitutive_titanmod_Nslip(f,instance) > 0_pInt) then   
         if (constitutive_titanmod_rho_edge0(f,instance) < 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='rho_edge0 ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_rho_screw0(f,instance) < 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='rho_screw0 ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_burgersPerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='slipburgers ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_f0_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='f0 ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_tau0e_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='tau0e ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_tau0s_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='tau0s ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_capre_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='capre ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_caprs_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='caprs ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_v0e_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='v0e ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_v0s_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='v0s ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_kinkcriticallength_PerSlipFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='kinkCriticalLength ('//PLASTICITY_TITANMOD_label//')')
       endif
     enddo
     do f = 1_pInt,lattice_maxNtwinFamily
       if (constitutive_titanmod_Ntwin(f,instance) > 0_pInt) then   
         if (constitutive_titanmod_burgersPerTwinFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='twinburgers ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_twinf0_PerTwinFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='twinf0 ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_twinshearconstant_PerTwinFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='twinshearconstant ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_twintau0_PerTwinFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='twintau0 ('//PLASTICITY_TITANMOD_label//')')
         if (constitutive_titanmod_twingamma0_PerTwinFam(f,instance) <= 0.0_pReal) &
           call IO_error(211_pInt,el=instance,ext_msg='twingamma0 ('//PLASTICITY_TITANMOD_label//')')
       endif
     enddo
     if (constitutive_titanmod_dc(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='dc ('//PLASTICITY_TITANMOD_label//')')
     if (constitutive_titanmod_twinhpconstant(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='twinhpconstant ('//PLASTICITY_TITANMOD_label//')')
     if (constitutive_titanmod_aTolRho(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_TITANMOD_label//')')
     
!--------------------------------------------------------------------------------------------------
! determine total number of active slip or twin systems
     constitutive_titanmod_Nslip(:,instance) = min(lattice_NslipSystem(:,phase),constitutive_titanmod_Nslip(:,instance))
     constitutive_titanmod_Ntwin(:,instance) = min(lattice_NtwinSystem(:,phase),constitutive_titanmod_Ntwin(:,instance))
     constitutive_titanmod_totalNslip(instance) = sum(constitutive_titanmod_Nslip(:,instance))
     constitutive_titanmod_totalNtwin(instance) = sum(constitutive_titanmod_Ntwin(:,instance))
   endif myPhase
 enddo sanityChecks

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip = maxval(constitutive_titanmod_totalNslip)
 maxTotalNtwin = maxval(constitutive_titanmod_totalNtwin)
 
 allocate(constitutive_titanmod_burgersPerSlipSys(maxTotalNslip, maxNinstance), source=0.0_pReal) 
 
 allocate(constitutive_titanmod_f0_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_tau0e_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_tau0s_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_capre_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_caprs_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_pe_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_ps_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_qe_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_qs_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_v0e_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_v0s_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_kinkcriticallength_PerSlipSys(maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CeLambdaSlipPerSlipSys(maxTotalNslip, maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_CsLambdaSlipPerSlipSys(maxTotalNslip, maxNinstance), source=0.0_pReal) 
 
 allocate(constitutive_titanmod_burgersPerTwinSys(maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinf0_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinshearconstant_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twintau0_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinp_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinq_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twingamma0_PerTwinSys(maxTotalNTwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinsizePerTwinSys(maxTotalNtwin, maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_twinLambdaSlipPerTwinSys(maxTotalNtwin, maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_Ctwin66                 (6,6,maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_Ctwin3333           (3,3,3,3,maxTotalNtwin,maxNinstance), source=0.0_pReal) 

 allocate(constitutive_titanmod_interactionMatrixSlipSlip(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrix_ee(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrix_ss(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrix_es(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrixSlipTwin(maxTotalNslip,maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrixTwinSlip(maxTotalNtwin,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_interactionMatrixTwinTwin(maxTotalNtwin,maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_forestProjectionScrew(maxTotalNslip,maxTotalNslip,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_TwinforestProjectionEdge(maxTotalNtwin,maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 allocate(constitutive_titanmod_TwinforestProjectionScrew(maxTotalNtwin,maxTotalNtwin,maxNinstance), source=0.0_pReal) 
 
 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(phase) == PLASTICITY_TITANMOD_ID) then
     instance = phase_plasticityInstance(phase) 

!--------------------------------------------------------------------------------------------------
! inverse lookup of slip system family
     l = 0_pInt
     do f = 1_pInt,lattice_maxNslipFamily
       do s = 1_pInt,constitutive_titanmod_Nslip(f,instance)
         l = l + 1_pInt
         constitutive_titanmod_slipFamily(l,instance) = f
         constitutive_titanmod_slipSystemLattice(l,instance) = sum(lattice_NslipSystem(1:f-1_pInt,phase)) + s
     enddo; enddo

!--------------------------------------------------------------------------------------------------
! inverse lookup of twin system family
     l = 0_pInt
     do f = 1_pInt,lattice_maxNtwinFamily
       do t = 1_pInt,constitutive_titanmod_Ntwin(f,instance)
         l = l + 1_pInt
         constitutive_titanmod_twinFamily(l,instance) = f
         constitutive_titanmod_twinSystemLattice(l,instance) = sum(lattice_NtwinSystem(1:f-1_pInt,phase)) + t
     enddo; enddo
   
!--------------------------------------------------------------------------------------------------
! determine size of state array  
     ns = constitutive_titanmod_totalNslip(instance)
     nt = constitutive_titanmod_totalNtwin(instance)
     constitutive_titanmod_sizeDotState(instance) = &
       size(constitutive_titanmod_listBasicSlipStates)*ns + &
       size(constitutive_titanmod_listBasicTwinStates)*nt
     constitutive_titanmod_sizeState(instance) = &
     constitutive_titanmod_sizeDotState(instance)+ &
       size(constitutive_titanmod_listDependentSlipStates)*ns + &
       size(constitutive_titanmod_listDependentTwinStates)*nt

!--------------------------------------------------------------------------------------------------
! determine size of postResults array  
     outputsLoop: do o = 1_pInt,constitutive_titanmod_Noutput(instance)
       mySize = 0_pInt
       select case(constitutive_titanmod_outputID(o,instance))
         case(rhoedge_ID,            rhoscrew_ID, &
              segment_edge_ID,       segment_screw_ID, &
              resistance_edge_ID,    resistance_screw_ID, &
              velocity_edge_ID,      velocity_screw_ID, &
              tau_slip_ID, &
              gdot_slip_edge_ID,     gdot_slip_screw_ID, &
              gdot_slip_ID, &
              stressratio_edge_p_ID, stressratio_screw_p_ID, &
              shear_system_ID)
            mySize = constitutive_titanmod_totalNslip(instance)
          case(twin_fraction_ID)
            mySize = constitutive_titanmod_totalNtwin(instance)
          case(shear_basal_ID,    shear_prism_ID,    shear_pyra_ID,    shear_pyrca_ID, &                  ! use only if all 4 slip families in hex are considered
               rhoedge_basal_ID,  rhoedge_prism_ID,  rhoedge_pyra_ID,  rhoedge_pyrca_ID, &
               rhoscrew_basal_ID, rhoscrew_prism_ID, rhoscrew_pyra_ID, rhoscrew_pyrca_ID, &
               shear_total_ID)
            mySize = 1_pInt
          case default
            call IO_error(105_pInt,ext_msg=constitutive_titanmod_output(o,instance)// &
                                                              ' ('//PLASTICITY_TITANMOD_label//')')
        end select

        outputFound: if (mySize > 0_pInt) then 
          constitutive_titanmod_sizePostResult(o,instance) = mySize
          constitutive_titanmod_sizePostResults(instance)  = constitutive_titanmod_sizePostResults(instance) + mySize
        endif outputFound
      enddo outputsLoop
   
!--------------------------------------------------------------------------------------------------
! construction of the twin elasticity matrices
     do j=1_pInt,lattice_maxNtwinFamily
       do k=1_pInt,constitutive_titanmod_Ntwin(j,instance)           
         do l=1_pInt,3_pInt ; do m=1_pInt,3_pInt ; do n=1_pInt,3_pInt ; do o=1_pInt,3_pInt
           do p=1_pInt,3_pInt ; do q=1_pInt,3_pInt ; do r=1_pInt,3_pInt ; do s=1_pInt,3_pInt
             constitutive_titanmod_Ctwin3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1_pInt,instance))+k,instance) = &
               constitutive_titanmod_Ctwin3333(l,m,n,o,sum(constitutive_titanmod_Nslip(1:j-1_pInt,instance))+k,instance) + &
               lattice_C3333(p,q,r,s,phase)*&
               lattice_Qtwin(l,p,sum(lattice_NslipSystem(1:j-1_pInt,phase))+k,phase)* &
               lattice_Qtwin(m,q,sum(lattice_NslipSystem(1:j-1_pInt,phase))+k,phase)* &
               lattice_Qtwin(n,r,sum(lattice_NslipSystem(1:j-1_pInt,phase))+k,phase)* &
               lattice_Qtwin(o,s,sum(lattice_NslipSystem(1:j-1_pInt,phase))+k,phase)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo ; enddo
         constitutive_titanmod_Ctwin66(1:6,1:6,k,instance) =  &
           math_Mandel3333to66(constitutive_titanmod_Ctwin3333(1:3,1:3,1:3,1:3,k,instance))
     enddo; enddo

!--------------------------------------------------------------------------------------------------
! Burgers vector, dislocation velocity prefactor for each slip system 
     do s = 1_pInt,constitutive_titanmod_totalNslip(instance)   
       f = constitutive_titanmod_slipFamily(s,instance)    

       constitutive_titanmod_burgersPerSlipSys(s,instance) = &
       constitutive_titanmod_burgersPerSlipFam(f,instance)

       constitutive_titanmod_f0_PerSlipSys(s,instance) = &
       constitutive_titanmod_f0_PerSlipFam(f,instance)

       constitutive_titanmod_tau0e_PerSlipSys(s,instance) = &
       constitutive_titanmod_tau0e_PerSlipFam(f,instance)

       constitutive_titanmod_tau0s_PerSlipSys(s,instance) = &
       constitutive_titanmod_tau0s_PerSlipFam(f,instance)

       constitutive_titanmod_capre_PerSlipSys(s,instance) = &
       constitutive_titanmod_capre_PerSlipFam(f,instance)

       constitutive_titanmod_caprs_PerSlipSys(s,instance) = &
       constitutive_titanmod_caprs_PerSlipFam(f,instance)

       constitutive_titanmod_v0e_PerSlipSys(s,instance) = &
       constitutive_titanmod_v0e_PerSlipFam(f,instance)

       constitutive_titanmod_v0s_PerSlipSys(s,instance) = &
       constitutive_titanmod_v0s_PerSlipFam(f,instance)

       constitutive_titanmod_kinkcriticallength_PerSlipSys(s,instance) = &
       constitutive_titanmod_kinkcriticallength_PerSlipFam(f,instance)

       constitutive_titanmod_pe_PerSlipSys(s,instance) = &
       constitutive_titanmod_pe_PerSlipFam(f,instance)

       constitutive_titanmod_ps_PerSlipSys(s,instance) = &
       constitutive_titanmod_ps_PerSlipFam(f,instance)

       constitutive_titanmod_qe_PerSlipSys(s,instance) = &
       constitutive_titanmod_qe_PerSlipFam(f,instance)

       constitutive_titanmod_qs_PerSlipSys(s,instance) = &
       constitutive_titanmod_qs_PerSlipFam(f,instance)

       constitutive_titanmod_CeLambdaSlipPerSlipSys(s,instance) = &
       constitutive_titanmod_CeLambdaSlipPerSlipFam(f,instance)

       constitutive_titanmod_CsLambdaSlipPerSlipSys(s,instance) = &
       constitutive_titanmod_CsLambdaSlipPerSlipFam(f,instance)
     enddo   
   
!--------------------------------------------------------------------------------------------------
! Burgers vector, nucleation rate prefactor and twin size for each twin system 
     do t = 1_pInt,constitutive_titanmod_totalNtwin(instance)   
       f = constitutive_titanmod_twinFamily(t,instance)
     
       constitutive_titanmod_burgersPerTwinSys(t,instance) = &
       constitutive_titanmod_burgersPerTwinFam(f,instance)

       constitutive_titanmod_twinsizePerTwinSys(t,instance) = &
       constitutive_titanmod_twinsizePerTwinFam(f,instance)

       constitutive_titanmod_twinf0_PerTwinSys(t,instance) = &
       constitutive_titanmod_twinf0_PerTwinFam(f,instance)

       constitutive_titanmod_twinshearconstant_PerTwinSys(t,instance) = &
       constitutive_titanmod_twinshearconstant_PerTwinFam(f,instance)

       constitutive_titanmod_twintau0_PerTwinSys(t,instance) = &
       constitutive_titanmod_twintau0_PerTwinFam(f,instance)

       constitutive_titanmod_twingamma0_PerTwinSys(t,instance) = &
       constitutive_titanmod_twingamma0_PerTwinFam(f,instance)

       constitutive_titanmod_twinp_PerTwinSys(t,instance) = &
       constitutive_titanmod_twinp_PerTwinFam(f,instance)

       constitutive_titanmod_twinq_PerTwinSys(t,instance) = &
       constitutive_titanmod_twinq_PerTwinFam(f,instance)

       constitutive_titanmod_twinLambdaSlipPerTwinSys(t,instance) = &
       constitutive_titanmod_twinLambdaSlipPerTwinFam(f,instance)
     enddo   
     
!--------------------------------------------------------------------------------------------------
! Construction of interaction matrices
     do s1 = 1_pInt,constitutive_titanmod_totalNslip(instance)
        do s2 = 1_pInt,constitutive_titanmod_totalNslip(instance)     
           constitutive_titanmod_interactionMatrixSlipSlip(s1,s2,instance) = &
             constitutive_titanmod_interactionSlipSlip(lattice_interactionSlipSlip( &
               constitutive_titanmod_slipSystemLattice(s1,instance),&
               constitutive_titanmod_slipSystemLattice(s2,instance),phase),instance)

           constitutive_titanmod_interactionMatrix_ee(s1,s2,instance) = &
             constitutive_titanmod_interaction_ee(lattice_interactionSlipSlip ( &
               constitutive_titanmod_slipSystemLattice(s1,instance), &
               constitutive_titanmod_slipSystemLattice(s2,instance), phase),instance)

           constitutive_titanmod_interactionMatrix_ss(s1,s2,instance) = &
             constitutive_titanmod_interaction_ss(lattice_interactionSlipSlip( &
               constitutive_titanmod_slipSystemLattice(s1,instance), &
               constitutive_titanmod_slipSystemLattice(s2,instance), phase),instance)

           constitutive_titanmod_interactionMatrix_es(s1,s2,instance) = &
             constitutive_titanmod_interaction_es(lattice_interactionSlipSlip( &
               constitutive_titanmod_slipSystemLattice(s1,instance), &
               constitutive_titanmod_slipSystemLattice(s2,instance), phase),instance)
     enddo; enddo
     
     do s1 = 1_pInt,constitutive_titanmod_totalNslip(instance)
        do t2 = 1_pInt,constitutive_titanmod_totalNtwin(instance)     
           constitutive_titanmod_interactionMatrixSlipTwin(s1,t2,instance) = &
             constitutive_titanmod_interactionSlipTwin(lattice_interactionSlipTwin( &
               constitutive_titanmod_slipSystemLattice(s1,instance), &
               constitutive_titanmod_twinSystemLattice(t2,instance), phase),instance)         
     enddo; enddo
     
     do t1 = 1_pInt,constitutive_titanmod_totalNtwin(instance)
        do s2 = 1_pInt,constitutive_titanmod_totalNslip(instance)     
           constitutive_titanmod_interactionMatrixTwinSlip(t1,s2,instance) = &
             constitutive_titanmod_interactionTwinSlip(lattice_interactionTwinSlip( &
               constitutive_titanmod_twinSystemLattice(t1,instance), &
               constitutive_titanmod_slipSystemLattice(s2,instance), phase),instance)         
     enddo; enddo

     do t1 = 1_pInt,constitutive_titanmod_totalNtwin(instance)
        do t2 = 1_pInt,constitutive_titanmod_totalNtwin(instance)     
           constitutive_titanmod_interactionMatrixTwinTwin(t1,t2,instance) = &
             constitutive_titanmod_interactionTwinTwin(lattice_interactionTwinTwin( &
             constitutive_titanmod_twinSystemLattice(t1,instance), &
             constitutive_titanmod_twinSystemLattice(t2,instance), phase),instance)         
     enddo; enddo
      
     do s1 = 1_pInt,constitutive_titanmod_totalNslip(instance)
        do s2 = 1_pInt,constitutive_titanmod_totalNslip(instance)      
!--------------------------------------------------------------------------------------------------
! calculation of forest projections for edge dislocations      
           constitutive_titanmod_forestProjectionEdge(s1,s2,instance) = &
           abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,instance),phase), &
                           lattice_st(:,constitutive_titanmod_slipSystemLattice(s2,instance),phase))) 

!--------------------------------------------------------------------------------------------------
! calculation of forest projections for screw dislocations 
           constitutive_titanmod_forestProjectionScrew(s1,s2,instance) = &
           abs(math_mul3x3(lattice_sn(:,constitutive_titanmod_slipSystemLattice(s1,instance),phase), &
                           lattice_sd(:,constitutive_titanmod_slipSystemLattice(s2,instance),phase))) 
     enddo; enddo
  
!--------------------------------------------------------------------------------------------------
! calculation of forest projections for edge dislocations in twin system
     do t1 = 1_pInt,constitutive_titanmod_totalNtwin(instance)
        do t2 = 1_pInt,constitutive_titanmod_totalNtwin(instance)      
           constitutive_titanmod_TwinforestProjectionEdge(t1,t2,instance) = &
           abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,instance),phase), &
                           lattice_tt(:,constitutive_titanmod_twinSystemLattice(t2,instance),phase))) 

!--------------------------------------------------------------------------------------------------
! calculation of forest projections for screw dislocations in twin system
           constitutive_titanmod_TwinforestProjectionScrew(t1,t2,instance) = &
           abs(math_mul3x3(lattice_tn(:,constitutive_titanmod_twinSystemLattice(t1,instance),phase), &
                           lattice_td(:,constitutive_titanmod_twinSystemLattice(t2,instance),phase))) 
     enddo; enddo

   endif
 enddo initializeInstances

end subroutine constitutive_titanmod_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_stateInit(instance,phase)
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_mu

 implicit none
 integer(pInt), intent(in) :: instance                                                              !< number specifying the instance of the plasticity
 integer(pInt), intent(in) :: phase                                                                 !< number specifying the phase of the plasticity
 real(pReal), dimension(constitutive_titanmod_sizeState(instance)) :: &
   constitutive_titanmod_stateInit

 integer(pInt) :: &
   s,s0,s1, &
   t,t0,t1, & 
   ns,nt,f
 real(pReal), dimension(constitutive_titanmod_totalNslip(instance)) ::  &
   rho_edge0, &
   rho_screw0, &
   shear_system0, &
   segment_edge0, &
   segment_screw0, &
   resistance_edge0, &
   resistance_screw0
 real(pReal), dimension(constitutive_titanmod_totalNtwin(instance)) ::  &
  twingamma_dot0, &
  resistance_twin0
 
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)

 
!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables for slip
 s1 = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily
    s0 = s1 + 1_pInt
    s1 = s0 + constitutive_titanmod_Nslip(f,instance) - 1_pInt 
    do s = s0,s1
       rho_edge0(s)    = constitutive_titanmod_rho_edge0(f,instance)
       rho_screw0(s) = constitutive_titanmod_rho_screw0(f,instance)
       shear_system0(s) = 0.0_pReal
    enddo 
 enddo
 
!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables for twin
 t1 = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily
    t0 = t1 + 1_pInt
    t1 = t0 + constitutive_titanmod_Ntwin(f,instance) - 1_pInt 
    do t = t0,t1
       twingamma_dot0(t)=0.0_pReal
    enddo 
 enddo
 
!--------------------------------------------------------------------------------------------------
! initialize dependent slip microstructural variables
 forall (s = 1_pInt:ns)
   segment_edge0(s) = constitutive_titanmod_CeLambdaSlipPerSlipSys(s,instance)/ &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,instance))+ &
     dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,instance)))
   segment_screw0(s) = constitutive_titanmod_CsLambdaSlipPerSlipSys(s,instance)/ &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_forestProjectionEdge(1:ns,s,instance))+ &
     dot_product((rho_screw0),constitutive_titanmod_forestProjectionScrew(1:ns,s,instance)))
   resistance_edge0(s) = &
     lattice_mu(phase)*constitutive_titanmod_burgersPerSlipSys(s,instance)* &
     sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_ee(1:ns,s,instance))+ &
     dot_product((rho_screw0),constitutive_titanmod_interactionMatrix_es(1:ns,s,instance)))
   resistance_screw0(s) = &
    lattice_mu(phase)*constitutive_titanmod_burgersPerSlipSys(s,instance)* &
    sqrt(dot_product((rho_edge0),constitutive_titanmod_interactionMatrix_es(1:ns,s,instance))+ &
    dot_product((rho_screw0), constitutive_titanmod_interactionMatrix_ss(1:ns,s,instance)))
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
pure function constitutive_titanmod_aTolState(instance)

 implicit none
 integer(pInt), intent(in) :: instance                                                              !< number specifying the instance of the plasticity
 
 real(pReal), dimension(constitutive_titanmod_sizeState(instance)) :: &
  constitutive_titanmod_aTolState
 
 constitutive_titanmod_aTolState = constitutive_titanmod_aTolRho(instance)
 
end function constitutive_titanmod_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_homogenizedC(state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_C66

implicit none
 real(pReal), dimension(6,6) :: &
   constitutive_titanmod_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec),  intent(in) :: &
   state                                                                                            !< microstructure state
real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   volumefraction_PerTwinSys
 integer(pInt) :: &
   phase, &
   instance, &
   ns, nt, &
   i
 real(pReal) :: &
   sumf
  
!--------------------------------------------------------------------------------------------------
! shortened notation
 phase = material_phase(ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)
 
!--------------------------------------------------------------------------------------------------
! total twin volume fraction
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,instance)
 enddo
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
!--------------------------------------------------------------------------------------------------
! homogenized elasticity matrix
 constitutive_titanmod_homogenizedC = (1.0_pReal-sumf)*lattice_C66(1:6,1:6,phase)
 do i=1_pInt,nt
    constitutive_titanmod_homogenizedC = constitutive_titanmod_homogenizedC &
                                       + volumefraction_PerTwinSys(i)*&
                                                   constitutive_titanmod_Ctwin66(1:6,1:6,i,instance)
 enddo 
 
end function constitutive_titanmod_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine constitutive_titanmod_microstructure(temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase,&
   phase_plasticityInstance
 use lattice, only: &
   lattice_mu

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 
 type(p_vec),   intent(inout) :: &
   state                                                                                            !< microstructure state

 integer(pInt) :: &
   instance, &
   ns, nt, s, t, &
   i, phase
 real(pReal) :: &
   sumf, &
   sfe                                                                                              ! stacking fault energy
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
              volumefraction_PerTwinSys
  
!--------------------------------------------------------------------------------------------------
!Shortened notation
 phase = material_phase(ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)
 
!--------------------------------------------------------------------------------------------------
! total twin volume fraction
 forall (i = 1_pInt:nt) &
   volumefraction_PerTwinSys(i)=state%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,instance) 
 
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 sfe = 0.0002_pReal*Temperature-0.0396_pReal

!--------------------------------------------------------------------------------------------------    
! average segment length for edge dislocations in matrix
 forall (s = 1_pInt:ns) &
   state%p(3_pInt*ns+nt+s) = constitutive_titanmod_CeLambdaSlipPerSlipSys(s,instance)/ &
     sqrt(dot_product(state%p(1:ns), &
     constitutive_titanmod_forestProjectionEdge(1:ns,s,instance))+ &
     dot_product(state%p(ns+1_pInt:2_pInt*ns), &
     constitutive_titanmod_forestProjectionScrew(1:ns,s,instance)))
!--------------------------------------------------------------------------------------------------    
! average segment length for screw dislocations in matrix
 forall (s = 1_pInt:ns) &
   state%p(4_pInt*ns+nt+s) = constitutive_titanmod_CsLambdaSlipPerSlipSys(s,instance)/ &
     sqrt(dot_product(state%p(1:ns), &
     constitutive_titanmod_forestProjectionEdge(1:ns,s,instance))+ &
     dot_product(state%p(ns+1_pInt:2_pInt*ns), &
     constitutive_titanmod_forestProjectionScrew(1:ns,s,instance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for edge dislocation motion
 forall (s = 1_pInt:ns) &
   state%p(5_pInt*ns+nt+s) = &
     lattice_mu(phase)*constitutive_titanmod_burgersPerSlipSys(s,instance)*&
     sqrt(dot_product((state%p(1:ns)),&
     constitutive_titanmod_interactionMatrix_ee(1:ns,s,instance))+ &
     dot_product((state%p(ns+1_pInt:2_pInt*ns)),&
     constitutive_titanmod_interactionMatrix_es(1:ns,s,instance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for screw dislocation motion
 forall (s = 1_pInt:ns) &
   state%p(6_pInt*ns+nt+s) = &
     lattice_mu(phase)*constitutive_titanmod_burgersPerSlipSys(s,instance)*&
     sqrt(dot_product((state%p(1:ns)),&
     constitutive_titanmod_interactionMatrix_es(1:ns,s,instance))+ &
     dot_product((state%p(ns+1_pInt:2_pInt*ns)),&
      constitutive_titanmod_interactionMatrix_ss(1:ns,s,instance)))
!--------------------------------------------------------------------------------------------------    
! threshold stress or slip resistance for dislocation motion in twin
 forall (t = 1_pInt:nt) &
   state%p(7_pInt*ns+nt+t) = &
     lattice_mu(phase)*constitutive_titanmod_burgersPerTwinSys(t,instance)*&
     (dot_product((abs(state%p(2_pInt*ns+1_pInt:2_pInt*ns+nt))),&
     constitutive_titanmod_interactionMatrixTwinTwin(1:nt,t,instance)))
 
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
   lattice_structure, &
   LATTICE_hex_ID
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                 intent(in) :: &
   temperature                                                                                      !< temperature at IP 
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec),                 intent(inout) :: &
   state                                                                                            !< microstructure state
 integer(pInt) :: &
   index_myFamily, instance,phase, &
   ns,nt, &
   f,i,j,k,l,m,n
 real(pReal) :: sumf, &
   StressRatio_edge_p,  minusStressRatio_edge_p,  StressRatio_edge_pminus1,  BoltzmannRatioedge, &
   StressRatio_screw_p, minusStressRatio_screw_p, StressRatio_screw_pminus1, BoltzmannRatioscrew, &
   twinStressRatio_p,   twinminusStressRatio_p,   twinStressRatio_pminus1,   BoltzmannRatiotwin, &
   twinDotGamma0, bottomstress_edge, bottomstress_screw, screwvelocity_prefactor
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(constitutive_titanmod_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_slip,dgdot_dtauslip,tau_slip, &
    edge_velocity, screw_velocity, &
    gdot_slip_edge, gdot_slip_screw
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_twin,dgdot_dtautwin,tau_twin, volumefraction_PerTwinSys
 
!--------------------------------------------------------------------------------------------------
! shortened notation
 phase = material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase) 
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)
 
 do i=1_pInt,nt
 volumefraction_PerTwinSys(i)=state%p(3_pInt*ns+i)/ &
                       constitutive_titanmod_twinshearconstant_PerTwinSys(i,instance)
 
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
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Nslip(f,instance)          ! process each (active) slip system in family
      j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase)) 
      if(lattice_structure(phase)==LATTICE_hex_ID) then ! only for prismatic and pyr <a> systems in hex
      screwvelocity_prefactor=constitutive_titanmod_debyefrequency(instance)* &
        state%p(4_pInt*ns+nt+j)*(constitutive_titanmod_burgersPerSlipSys(j,instance)/ &
        constitutive_titanmod_kinkcriticallength_PerSlipSys(j,instance))**2

     !* Stress ratio for screw ! No slip resistance for screw dislocations, only Peierls stress
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,instance)
         StressRatio_screw_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**constitutive_titanmod_ps_PerSlipSys(j,instance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif
        
         bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,instance)
      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/ &
         ( bottomstress_screw) &
        )**(constitutive_titanmod_ps_PerSlipSys(j,instance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_kinkf0(instance)/(kB*Temperature)

        else  ! if the structure is not hex or the slip family is basal
        screwvelocity_prefactor=constitutive_titanmod_v0s_PerSlipSys(j,instance)
        bottomstress_screw=constitutive_titanmod_tau0s_PerSlipSys(j,instance)+state%p(6*ns+nt+j)
        StressRatio_screw_p = ((abs(tau_slip(j)))/( bottomstress_screw ))**constitutive_titanmod_ps_PerSlipSys(j,instance)

        if((1.0_pReal-StressRatio_screw_p)>0.001_pReal) then
        minusStressRatio_screw_p=1.0_pReal-StressRatio_screw_p
        else
        minusStressRatio_screw_p=0.001_pReal
        endif

      StressRatio_screw_pminus1 = ((abs(tau_slip(j)))/( bottomstress_screw))** &
      (constitutive_titanmod_ps_PerSlipSys(j,instance)-1.0_pReal)

      !* Boltzmann ratio for screw
      BoltzmannRatioscrew = constitutive_titanmod_f0_PerSlipSys(j,instance)/(kB*Temperature)

        endif

        !* Stress ratio for edge
         bottomstress_edge=constitutive_titanmod_tau0e_PerSlipSys(j,instance)+state%p(5*ns+nt+j)
         StressRatio_edge_p = ((abs(tau_slip(j)))/ &
         ( bottomstress_edge) &
        )**constitutive_titanmod_pe_PerSlipSys(j,instance)
                                
        if((1.0_pReal-StressRatio_edge_p)>0.001_pReal) then
        minusStressRatio_edge_p=1.0_pReal-StressRatio_edge_p
        else
        minusStressRatio_edge_p=0.001_pReal
        endif
        
      StressRatio_edge_pminus1 = ((abs(tau_slip(j)))/( bottomstress_edge))** &
      (constitutive_titanmod_pe_PerSlipSys(j,instance)-1.0_pReal)

      !* Boltzmann ratio for edge. For screws it is defined above
      BoltzmannRatioedge = constitutive_titanmod_f0_PerSlipSys(j,instance)/(kB*Temperature)
            
        screw_velocity(j) =screwvelocity_prefactor * & ! there is no v0 for screw now because it is included in the prefactor
                exp(-BoltzmannRatioscrew*(minusStressRatio_screw_p)** &
        constitutive_titanmod_qs_PerSlipSys(j,instance))

        edge_velocity(j) =constitutive_titanmod_v0e_PerSlipSys(j,instance)*exp(-BoltzmannRatioedge* &
        (minusStressRatio_edge_p)** &
        constitutive_titanmod_qe_PerSlipSys(j,instance))

                !* Shear rates due to edge slip
       gdot_slip_edge(j) = constitutive_titanmod_burgersPerSlipSys(j,instance)*(state%p(j)* &
                edge_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !* Shear rates due to screw slip
       gdot_slip_screw(j) = constitutive_titanmod_burgersPerSlipSys(j,instance)*(state%p(ns+j) * &
                screw_velocity(j))* sign(1.0_pReal,tau_slip(j))
                !Total shear rate
            
       gdot_slip(j) = gdot_slip_edge(j) + gdot_slip_screw(j)
                
       state%p(7*ns+2*nt+j)=edge_velocity(j)
       state%p(8*ns+2*nt+j)=screw_velocity(j)
       state%p(9*ns+2*nt+j)=tau_slip(j)
       state%p(10*ns+2*nt+j)=gdot_slip_edge(j)
       state%p(11*ns+2*nt+j)=gdot_slip_screw(j)
       state%p(12*ns+2*nt+j)=StressRatio_edge_p
       state%p(13*ns+2*nt+j)=StressRatio_screw_p
                
      !* Derivatives of shear rates
      dgdot_dtauslip(j) = constitutive_titanmod_burgersPerSlipSys(j,instance)*(( &
                ( &
                ( &
                ( &
                (edge_velocity(j)*state%p(j))) * &
                BoltzmannRatioedge*&
        constitutive_titanmod_pe_PerSlipSys(j,instance)* &
                constitutive_titanmod_qe_PerSlipSys(j,instance) &
                )/ &
                bottomstress_edge &
                )*&
        StressRatio_edge_pminus1*(minusStressRatio_edge_p)** &
                (constitutive_titanmod_qe_PerSlipSys(j,instance)-1.0_pReal) &
                ) + &
                ( &
                ( &
                ( &
                (state%p(ns+j) * screw_velocity(j)) * &
                BoltzmannRatioscrew* &
        constitutive_titanmod_ps_PerSlipSys(j,instance)* &
                constitutive_titanmod_qs_PerSlipSys(j,instance) &
                )/ &
                bottomstress_screw &
                )*&
        StressRatio_screw_pminus1*(minusStressRatio_screw_p)**(constitutive_titanmod_qs_PerSlipSys(j,instance)-1.0_pReal) &
                ) &
                ) !* sign(1.0_pReal,tau_slip(j))
                

                
!*************************************************                
!sumf=0.0_pReal
      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + (1.0_pReal - sumf)*gdot_slip(j)*lattice_Sslip(1:3,1:3,1,index_myFamily+i,phase)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*&
                                  lattice_Sslip(k,l,1,index_myFamily+i,phase)*&
                                  lattice_Sslip(m,n,1,index_myFamily+i,phase) 
   enddo
 enddo slipFamiliesLoop

!* Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
  twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
    index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase)) ! at which index starts my family
    do i = 1_pInt,constitutive_titanmod_Ntwin(f,instance)          ! process each (active) slip system in family
       j = j+1_pInt

      !* Calculation of Lp
      !* Resolved shear stress on twin system
       tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))        
     
!**************************************************************************************
      !* Stress ratios
!      StressRatio_r = (state%p(6*ns+3*nt+j)/tau_twin(j))**constitutive_titanmod_r(instance)      
      
          !* Shear rates and their derivatives due to twin
!      if ( tau_twin(j) > 0.0_pReal ) !then          
!        gdot_twin(j) =  0.0_pReal!&
!          (constitutive_titanmod_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,phase)*&
!          state%p(6*ns+4*nt+j)*constitutive_titanmod_Ndot0PerTwinSys(f,instance)*exp(-StressRatio_r) 
!        dgdot_dtautwin(j) = ((gdot_twin(j)*constitutive_titanmod_r(instance))/tau_twin(j))*StressRatio_r
!      endif
!**************************************************************************************
   
     !* Stress ratio for edge
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,instance)+state%p(7*ns+nt+j)) &
        )**constitutive_titanmod_twinp_PerTwinSys(j,instance)
               
        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
              
      twinStressRatio_pminus1 = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,instance)+state%p(7*ns+nt+j)) &
        )**(constitutive_titanmod_twinp_PerTwinSys(j,instance)-1.0_pReal)

      !* Boltzmann ratio
      BoltzmannRatiotwin = constitutive_titanmod_twinf0_PerTwinSys(j,instance)/(kB*Temperature)

      !* Initial twin shear rates
      TwinDotGamma0 = &
        constitutive_titanmod_twingamma0_PerTwinSys(j,instance)

      !* Shear rates due to twin
         gdot_twin(j) =sign(1.0_pReal,tau_twin(j))*constitutive_titanmod_twingamma0_PerTwinSys(j,instance)* &
         exp(-BoltzmannRatiotwin*(twinminusStressRatio_p)**constitutive_titanmod_twinq_PerTwinSys(j,instance))
         
                                      
      !* Derivatives of shear rates in twin
      dgdot_dtautwin(j) = ( &
                ( &
                ( &
                (abs(gdot_twin(j))) * &
                BoltzmannRatiotwin*&
        constitutive_titanmod_twinp_PerTwinSys(j,instance)* &
                constitutive_titanmod_twinq_PerTwinSys(j,instance) &
                )/ &
                constitutive_titanmod_twintau0_PerTwinSys(j,instance) &
                )*&
        twinStressRatio_pminus1*(twinminusStressRatio_p)** &
                (constitutive_titanmod_twinq_PerTwinSys(j,instance)-1.0_pReal) &
                ) !* sign(1.0_pReal,tau_slip(j))
                
      !* Plastic velocity gradient for mechanical twinning                                                      
!      Lp = Lp + sumf*gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,phase)
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,phase)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,phase)*&
                                  lattice_Stwin(m,n,index_myFamily+i,phase)
   enddo
 enddo twinFamiliesLoop

dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

end subroutine constitutive_titanmod_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_titanmod_dotState(Tstar_v,temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
  use lattice,  only: &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

implicit none
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec),               intent(in) :: &
   state                                                                                            !< microstructure state
 real(pReal), dimension(constitutive_titanmod_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_titanmod_dotState

 integer(pInt) :: &
   index_myFamily, instance,phase, &
   ns,nt,&
   f,i,j
 real(pReal) :: &
    sumf,BoltzmannRatio, &
    twinStressRatio_p,twinminusStressRatio_p
 real(pReal), dimension(constitutive_titanmod_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   DotRhoEdgeGeneration, &
   DotRhoEdgeAnnihilation, &
   DotRhoScrewGeneration, &
   DotRhoScrewAnnihilation
 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_twin, &
   tau_twin, &
   volumefraction_PerTwinSys
   
!--------------------------------------------------------------------------------------------------
! shortened notation
 phase =  material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase)
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)

 do i=1_pInt,nt
   volumefraction_PerTwinSys(i)=state%p(3_pInt*ns+i)/ &
                                   constitutive_titanmod_twinshearconstant_PerTwinSys(i,instance)

 enddo

 sumf = sum(abs(volumefraction_PerTwinSys(1_pInt:nt))) ! safe for nt == 0

 constitutive_titanmod_dotState = 0.0_pReal

 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                 ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Nslip(f,instance)                        ! process each (active) slip system in family
     j = j+1_pInt

      DotRhoEdgeGeneration(j) = &                                                                   ! multiplication of edge dislocations
        state%p(ns+j)*state%p(8*ns+2*nt+j)/state%p(4*ns+nt+j)
      DotRhoScrewGeneration(j) = &                                                                  ! multiplication of screw dislocations
        state%p(j)*state%p(7*ns+2*nt+j)/state%p(3*ns+nt+j)
      DotRhoEdgeAnnihilation(j) = -((state%p(j))**2)* &                                  ! annihilation of edge dislocations
           constitutive_titanmod_capre_PerSlipSys(j,instance)*state%p(7*ns+2*nt+j)*0.5_pReal
      DotRhoScrewAnnihilation(j) = -((state%p(ns+j))**2)* &                              ! annihilation of screw dislocations
           constitutive_titanmod_caprs_PerSlipSys(j,instance)*state%p(8*ns+2*nt+j)*0.5_pReal
      constitutive_titanmod_dotState(j) = &                                                         ! edge dislocation density rate of change
        DotRhoEdgeGeneration(j)+DotRhoEdgeAnnihilation(j)

      constitutive_titanmod_dotState(ns+j) = &                                                      ! screw dislocation density rate of change
        DotRhoScrewGeneration(j)+DotRhoScrewAnnihilation(j)

      constitutive_titanmod_dotState(2*ns+j) = &                                                    ! sum of shear due to edge and screw
        state%p(10*ns+2*nt+j)+state%p(11*ns+2*nt+j) 
    enddo
 enddo slipFamiliesLoop
  
!* Twin fraction evolution
 j = 0_pInt
 twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase)) ! at which index starts my family
   do i = 1_pInt,constitutive_titanmod_Ntwin(f,instance)          ! process each (active) twin system in family
      j = j+1_pInt

      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))

     !* Stress ratio for edge
         twinStressRatio_p = ((abs(tau_twin(j)))/ &
         ( constitutive_titanmod_twintau0_PerTwinSys(j,instance)+state%p(7*ns+nt+j)) &
           )**(constitutive_titanmod_twinp_PerTwinSys(j,instance))
        

        if((1.0_pReal-twinStressRatio_p)>0.001_pReal) then
        twinminusStressRatio_p=1.0_pReal-twinStressRatio_p
        else
        twinminusStressRatio_p=0.001_pReal
        endif
        
      BoltzmannRatio = constitutive_titanmod_twinf0_PerTwinSys(j,instance)/(kB*Temperature)
       
            gdot_twin(j) =constitutive_titanmod_twingamma0_PerTwinSys(j,instance)*exp(-BoltzmannRatio* &
              (twinminusStressRatio_p)** &
            constitutive_titanmod_twinq_PerTwinSys(j,instance))*sign(1.0_pReal,tau_twin(j))
             
      constitutive_titanmod_dotState(3*ns+j)=gdot_twin(j)

    enddo
  enddo twinFamiliesLoop

end function constitutive_titanmod_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_titanmod_postResults(state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance, &
   phase_Noutput
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(p_vec),  intent(in) :: &
   state                                                                                            !< microstructure state
 real(pReal), dimension(constitutive_titanmod_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_titanmod_postResults

 integer(pInt) :: &
   instance, phase,&
   ns,nt,&
   o,i,c
 real(pReal) :: sumf

 real(pReal), dimension(constitutive_titanmod_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
         volumefraction_PerTwinSys
 
!--------------------------------------------------------------------------------------------------
! shortened notation
 phase = material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase)
 ns = constitutive_titanmod_totalNslip(instance)
 nt = constitutive_titanmod_totalNtwin(instance)
 
 do i=1_pInt,nt
   volumefraction_PerTwinSys(i)=state%p(3_pInt*ns+i)/ &
         constitutive_titanmod_twinshearconstant_PerTwinSys(i,instance)
 enddo
 
 sumf = sum(abs(volumefraction_PerTwinSys(1:nt))) ! safe for nt == 0
 
 
!--------------------------------------------------------------------------------------------------
! required output 
 c = 0_pInt
 constitutive_titanmod_postResults = 0.0_pReal

 do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
   select case(constitutive_titanmod_outputID(o,instance))
     case (rhoedge_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(1_pInt:ns)
       c = c + ns
     case (rhoscrew_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(ns+1_pInt:2_pInt*ns)
       c = c + ns
     case (segment_edge_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(3_pInt*ns+nt+1_pInt:4_pInt*ns+nt)
       c = c + ns
     case (segment_screw_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(4_pInt*ns+nt+1_pInt:5_pInt*ns+nt)
       c = c + ns
     case (resistance_edge_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(5_pInt*ns+nt+1_pInt:6_pInt*ns+nt)
       c = c + ns
     case (resistance_screw_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(6_pInt*ns+nt+1_pInt:7_pInt*ns+nt)
       c = c + ns
     case (velocity_edge_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(7*ns+2*nt+1:8*ns+2*nt)
       c = c + ns
     case (velocity_screw_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = state%p(8*ns+2*nt+1:9*ns+2*nt)
       c = c + ns
     case (tau_slip_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(9*ns+2*nt+1:10*ns+2*nt))
       c = c + ns
     case (gdot_slip_edge_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(10*ns+2*nt+1:11*ns+2*nt))
       c = c + ns
     case (gdot_slip_screw_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(11*ns+2*nt+1:12*ns+2*nt))
       c = c + ns
     case (gdot_slip_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(10*ns+2*nt+1:11*ns+2*nt)) + &
                                                  abs(state%p(11*ns+2*nt+1:12*ns+2*nt))
       c = c + ns
     case (stressratio_edge_p_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(12*ns+2*nt+1:13*ns+2*nt))
       c = c + ns
     case (stressratio_screw_p_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(13*ns+2*nt+1:14*ns+2*nt))
       c = c + ns
     case (shear_system_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+ns) = abs(state%p(2*ns+1:3*ns))
       c = c + ns
     case (shear_basal_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state%p(2*ns+1:2*ns+3)))
       c = c + 1_pInt
     case (shear_prism_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state%p(2*ns+4:2*ns+6)))
       c = c + 1_pInt
     case (shear_pyra_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state%p(2*ns+7:2*ns+12)))
       c = c + 1_pInt
     case (shear_pyrca_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state%p(2*ns+13:2*ns+24)))
       c = c + 1_pInt

     case (rhoedge_basal_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(1:3))
       c = c + 1_pInt
     case (rhoedge_prism_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(4:6))
       c = c + 1_pInt
     case (rhoedge_pyra_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(7:12))
       c = c + 1_pInt
     case (rhoedge_pyrca_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(13:24))
       c = c + 1_pInt

     case (rhoscrew_basal_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(ns+1:ns+3))
       c = c + 1_pInt
     case (rhoscrew_prism_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(ns+4:ns+6))
       c = c + 1_pInt
     case (rhoscrew_pyra_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(ns+7:ns+12))
       c = c + 1_pInt
     case (rhoscrew_pyrca_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(state%p(ns+13:ns+24))
       c = c + 1_pInt
     case (shear_total_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+1_pInt) = sum(abs(state%p(2*ns+1:3*ns)))
       c = c + 1_pInt
     case (twin_fraction_ID)
       constitutive_titanmod_postResults(c+1_pInt:c+nt) = abs(volumefraction_PerTwinSys(1:nt))
       c = c + nt
   end select
 enddo

end function constitutive_titanmod_postResults

end module constitutive_titanmod
