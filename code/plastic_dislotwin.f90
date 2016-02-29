!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module plastic_dislotwin
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
    plastic_dislotwin_sizePostResults                                                          !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   plastic_dislotwin_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   plastic_dislotwin_output                                                                    !< name of each post result output
   
 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                         !< Boltzmann constant in J/Kelvin

 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   plastic_dislotwin_Noutput                                                                   !< number of outputs per instance of this plasticity 

 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   plastic_dislotwin_totalNslip, &                                                             !< total number of active slip systems for each instance
   plastic_dislotwin_totalNtwin, &                                                             !< total number of active twin systems for each instance
   plastic_dislotwin_totalNtrans                                                               !< number of active transformation systems 

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   plastic_dislotwin_Nslip, &                                                                  !< number of active slip systems for each family and instance
   plastic_dislotwin_Ntwin, &                                                                  !< number of active twin systems for each family and instance
   plastic_dislotwin_Ntrans                                                                    !< number of active transformation systems for each family and instance

 real(pReal),                         dimension(:),           allocatable,         private :: &
   plastic_dislotwin_CAtomicVolume, &                                                          !< atomic volume in Bugers vector unit
   plastic_dislotwin_D0, &                                                                     !< prefactor for self-diffusion coefficient
   plastic_dislotwin_Qsd, &                                                                    !< activation energy for dislocation climb
   plastic_dislotwin_GrainSize, &                                                              !< grain size
   plastic_dislotwin_pShearBand, &                                                             !< p-exponent in shearband velocity
   plastic_dislotwin_qShearBand, &                                                             !< q-exponent in shearband velocity
   plastic_dislotwin_MaxTwinFraction, &                                                        !< maximum allowed total twin volume fraction
   plastic_dislotwin_CEdgeDipMinDistance, &                                                    !<
   plastic_dislotwin_Cmfptwin, &                                                               !<
   plastic_dislotwin_Cthresholdtwin, &                                                         !<
   plastic_dislotwin_SolidSolutionStrength, &                                                  !< Strength due to elements in solid solution
   plastic_dislotwin_L0_twin, &                                                                !< Length of twin nuclei in Burgers vectors
   plastic_dislotwin_L0_trans, &                                                               !< Length of trans nuclei in Burgers vectors
   plastic_dislotwin_xc_twin, &                                                                !< critical distance for formation of twin nucleus
   plastic_dislotwin_xc_trans, &                                                               !< critical distance for formation of trans nucleus
   plastic_dislotwin_VcrossSlip, &                                                             !< cross slip volume
   plastic_dislotwin_sbResistance, &                                                           !< value for shearband resistance (might become an internal state variable at some point)
   plastic_dislotwin_sbVelocity, &                                                             !< value for shearband velocity_0
   plastic_dislotwin_sbQedge, &                                                                !< value for shearband systems Qedge
   plastic_dislotwin_SFE_0K, &                                                                 !< stacking fault energy at zero K
   plastic_dislotwin_dSFE_dT, &                                                                !< temperature dependance of stacking fault energy
   plastic_dislotwin_dipoleFormationFactor, &                                                  !< scaling factor for dipole formation: 0: off, 1: on. other values not useful
   plastic_dislotwin_aTolRho, &                                                                !< absolute tolerance for integration of dislocation density
   plastic_dislotwin_aTolTwinFrac, &                                                           !< absolute tolerance for integration of twin volume fraction
   plastic_dislotwin_aTolTransFrac, &                                                          !< absolute tolerance for integration of trans volume fraction
   plastic_dislotwin_deltaG, &                                                                 !< Free energy difference between austensite and martensite
   plastic_dislotwin_Cmfptrans, &                                                              !<
   plastic_dislotwin_Cthresholdtrans, &                                                        !<
   plastic_dislotwin_transStackHeight                                                          !< Stack height of hex nucleus

 real(pReal),                         dimension(:,:,:,:),     allocatable,         private :: &
   plastic_dislotwin_Ctwin66                                                                   !< twin elasticity matrix in Mandel notation for each instance
 real(pReal),                         dimension(:,:,:,:,:,:), allocatable,         private :: &
   plastic_dislotwin_Ctwin3333                                                                 !< twin elasticity matrix for each instance
 real(pReal),                         dimension(:,:,:,:),     allocatable,         private :: &
   plastic_dislotwin_Ctrans66                                                                   !< trans elasticity matrix in Mandel notation for each instance
 real(pReal),                         dimension(:,:,:,:,:,:), allocatable,         private :: &
   plastic_dislotwin_Ctrans3333                                                                 !< trans elasticity matrix for each instance
 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   plastic_dislotwin_rhoEdge0, &                                                               !< initial edge dislocation density per slip system for each family and instance
   plastic_dislotwin_rhoEdgeDip0, &                                                            !< initial edge dipole density per slip system for each family and instance
   plastic_dislotwin_burgersPerSlipFamily, &                                                   !< absolute length of burgers vector [m] for each slip family and instance
   plastic_dislotwin_burgersPerSlipSystem, &                                                   !< absolute length of burgers vector [m] for each slip system and instance
   plastic_dislotwin_burgersPerTwinFamily, &                                                   !< absolute length of burgers vector [m] for each twin family and instance
   plastic_dislotwin_burgersPerTwinSystem, &                                                   !< absolute length of burgers vector [m] for each twin system and instance
   plastic_dislotwin_burgersPerTransFamily, &                                                  !< absolute length of burgers vector [m] for each trans family and instance
   plastic_dislotwin_burgersPerTransSystem, &                                                  !< absolute length of burgers vector [m] for each trans system and instance
   plastic_dislotwin_QedgePerSlipFamily, &                                                     !< activation energy for glide [J] for each slip family and instance
   plastic_dislotwin_QedgePerSlipSystem, &                                                     !< activation energy for glide [J] for each slip system and instance
   plastic_dislotwin_v0PerSlipFamily, &                                                        !< dislocation velocity prefactor [m/s] for each family and instance
   plastic_dislotwin_v0PerSlipSystem, &                                                        !< dislocation velocity prefactor [m/s] for each slip system and instance
   plastic_dislotwin_tau_peierlsPerSlipFamily, &                                               !< Peierls stress [Pa] for each family and instance
   plastic_dislotwin_Ndot0PerTwinFamily, &                                                     !< twin nucleation rate [1/m³s] for each twin family and instance
   plastic_dislotwin_Ndot0PerTwinSystem, &                                                     !< twin nucleation rate [1/m³s] for each twin system and instance
   plastic_dislotwin_Ndot0PerTransFamily, &                                                    !< trans nucleation rate [1/m³s] for each trans family and instance
   plastic_dislotwin_Ndot0PerTransSystem, &                                                    !< trans nucleation rate [1/m³s] for each trans system and instance
   plastic_dislotwin_tau_r_twin, &                                                             !< stress to bring partial close together for each twin system and instance
   plastic_dislotwin_tau_r_trans, &                                                            !< stress to bring partial close together for each trans system and instance
   plastic_dislotwin_twinsizePerTwinFamily, &                                                  !< twin thickness [m] for each twin family and instance
   plastic_dislotwin_twinsizePerTwinSystem, &                                                  !< twin thickness [m] for each twin system and instance
   plastic_dislotwin_CLambdaSlipPerSlipFamily, &                                               !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   plastic_dislotwin_CLambdaSlipPerSlipSystem, &                                               !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   plastic_dislotwin_lamellarsizePerTransFamily, &                                             !< martensite lamellar thickness [m] for each trans family and instance
   plastic_dislotwin_lamellarsizePerTransSystem, &                                             !< martensite lamellar thickness [m] for each trans system and instance
   plastic_dislotwin_interaction_SlipSlip, &                                                   !< coefficients for slip-slip interaction for each interaction type and instance
   plastic_dislotwin_interaction_SlipTwin, &                                                   !< coefficients for slip-twin interaction for each interaction type and instance
   plastic_dislotwin_interaction_TwinSlip, &                                                   !< coefficients for twin-slip interaction for each interaction type and instance
   plastic_dislotwin_interaction_TwinTwin, &                                                   !< coefficients for twin-twin interaction for each interaction type and instance
   plastic_dislotwin_interaction_SlipTrans, &                                                  !< coefficients for slip-trans interaction for each interaction type and instance
   plastic_dislotwin_interaction_TransSlip, &                                                  !< coefficients for trans-slip interaction for each interaction type and instance
   plastic_dislotwin_interaction_TransTrans, &                                                 !< coefficients for trans-trans interaction for each interaction type and instance
   plastic_dislotwin_pPerSlipFamily, &                                                         !< p-exponent in glide velocity
   plastic_dislotwin_qPerSlipFamily, &                                                         !< q-exponent in glide velocity
   plastic_dislotwin_rPerTwinFamily, &                                                         !< r-exponent in twin nucleation rate
   plastic_dislotwin_sPerTransFamily                                                           !< s-exponent in trans nucleation rate
 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   plastic_dislotwin_interactionMatrix_SlipSlip, &                                             !< interaction matrix of the different slip systems for each instance
   plastic_dislotwin_interactionMatrix_SlipTwin, &                                             !< interaction matrix of slip systems with twin systems for each instance
   plastic_dislotwin_interactionMatrix_TwinSlip, &                                             !< interaction matrix of twin systems with slip systems for each instance
   plastic_dislotwin_interactionMatrix_TwinTwin, &                                             !< interaction matrix of the different twin systems for each instance
   plastic_dislotwin_interactionMatrix_SlipTrans, &                                            !< interaction matrix of slip systems with trans systems for each instance
   plastic_dislotwin_interactionMatrix_TransSlip, &                                            !< interaction matrix of trans systems with slip systems for each instance
   plastic_dislotwin_interactionMatrix_TransTrans, &                                           !< interaction matrix of the different trans systems for each instance
   plastic_dislotwin_forestProjectionEdge, &                                                   !< matrix of forest projections of edge dislocations for each instance
   plastic_dislotwin_projectionMatrix_Trans                                                    !< matrix for projection of slip system shear on fault band (twin) systems for each instance

 real(pReal),                         dimension(:,:,:,:,:),   allocatable,         private :: &
   plastic_dislotwin_sbSv

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 edge_density_ID, &
                 dipole_density_ID, &
                 shear_rate_slip_ID, &
                 accumulated_shear_slip_ID, &
                 mfp_slip_ID, &
                 resolved_stress_slip_ID, &
                 threshold_stress_slip_ID, &
                 edge_dipole_distance_ID, &
                 stress_exponent_ID, &
                 twin_fraction_ID, &
                 shear_rate_twin_ID, &
                 accumulated_shear_twin_ID, &
                 mfp_twin_ID, &
                 resolved_stress_twin_ID, &
                 threshold_stress_twin_ID, &
                 resolved_stress_shearband_ID, &
                 shear_rate_shearband_ID, &
                 sb_eigenvalues_ID, &
                 sb_eigenvectors_ID, &
                 stress_trans_fraction_ID, &
                 strain_trans_fraction_ID, &
                 trans_fraction_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   plastic_dislotwin_outputID                                                                  !< ID of each post result output
 type, private :: tDislotwinState
   real(pReal), pointer,     dimension(:,:) :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear_slip, &
     twinFraction, &
     accshear_twin, &
     stressTransFraction, &
     strainTransFraction , &
     invLambdaSlip, &
     invLambdaSlipTwin, &
     invLambdaTwin, &
     invLambdaSlipTrans, &
     invLambdaTrans, &
     mfp_slip, &
     mfp_twin, &
     mfp_trans, &
     threshold_stress_slip, &
     threshold_stress_twin, &
     threshold_stress_trans, &
     twinVolume, &
     martensiteVolume
 end type
 type(tDislotwinState), allocatable, dimension(:), private :: &
   state, &
   state0, &
   dotState

 public :: &
   plastic_dislotwin_init, &
   plastic_dislotwin_homogenizedC, &
   plastic_dislotwin_microstructure, &
   plastic_dislotwin_LpAndItsTangent, &
   plastic_dislotwin_dotState, &
   plastic_dislotwin_postResults
 private :: &
   plastic_dislotwin_stateInit, &
   plastic_dislotwin_aTolState

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
   math_mul3x3
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems
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
   PLASTICITY_DISLOTWIN_label, &
   PLASTICITY_DISLOTWIN_ID, &
   material_phase, &  
   plasticState, & 
   MATERIAL_partPhase
 use lattice
 use numerics,only: &
   analyticJaco, &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,maxTotalNslip,maxTotalNtwin,maxTotalNtrans,&
                  f,instance,j,k,l,m,n,o,p,q,r,s,ns,nt,nr, &
                  Nchunks_SlipSlip = 0_pInt, Nchunks_SlipTwin = 0_pInt, &
                  Nchunks_TwinSlip = 0_pInt, Nchunks_TwinTwin = 0_pInt, &
                  Nchunks_SlipTrans = 0_pInt, Nchunks_TransSlip = 0_pInt, Nchunks_TransTrans = 0_pInt, &
                  Nchunks_SlipFamilies = 0_pInt, Nchunks_TwinFamilies = 0_pInt, Nchunks_TransFamilies = 0_pInt, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip, tempPerTwin, tempPerTrans
  
 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOTWIN_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOTWIN_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(plastic_dislotwin_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(plastic_dislotwin_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(plastic_dislotwin_output(maxval(phase_Noutput),maxNinstance))
          plastic_dislotwin_output = ''
 allocate(plastic_dislotwin_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(plastic_dislotwin_Noutput(maxNinstance),                             source=0_pInt)
 allocate(plastic_dislotwin_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(plastic_dislotwin_Ntwin(lattice_maxNtwinFamily,maxNinstance),        source=0_pInt)
 allocate(plastic_dislotwin_Ntrans(lattice_maxNtransFamily,maxNinstance),      source=0_pInt)
 allocate(plastic_dislotwin_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(plastic_dislotwin_totalNtwin(maxNinstance),                          source=0_pInt)
 allocate(plastic_dislotwin_totalNtrans(maxNinstance),                         source=0_pInt)
 allocate(plastic_dislotwin_CAtomicVolume(maxNinstance),                       source=0.0_pReal)
 allocate(plastic_dislotwin_D0(maxNinstance),                                  source=0.0_pReal)
 allocate(plastic_dislotwin_Qsd(maxNinstance),                                 source=0.0_pReal)
 allocate(plastic_dislotwin_GrainSize(maxNinstance),                           source=0.0_pReal)
 allocate(plastic_dislotwin_pShearBand(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_dislotwin_qShearBand(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_dislotwin_MaxTwinFraction(maxNinstance),                     source=0.0_pReal)
 allocate(plastic_dislotwin_CEdgeDipMinDistance(maxNinstance),                 source=0.0_pReal)
 allocate(plastic_dislotwin_Cmfptwin(maxNinstance),                            source=0.0_pReal)
 allocate(plastic_dislotwin_Cthresholdtwin(maxNinstance),                      source=0.0_pReal)
 allocate(plastic_dislotwin_SolidSolutionStrength(maxNinstance),               source=0.0_pReal)
 allocate(plastic_dislotwin_L0_twin(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_dislotwin_L0_trans(maxNinstance),                            source=0.0_pReal)
 allocate(plastic_dislotwin_xc_twin(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_dislotwin_xc_trans(maxNinstance),                            source=0.0_pReal)
 allocate(plastic_dislotwin_VcrossSlip(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_dislotwin_aTolRho(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_dislotwin_aTolTwinFrac(maxNinstance),                        source=0.0_pReal)
 allocate(plastic_dislotwin_aTolTransFrac(maxNinstance),                       source=0.0_pReal)
 allocate(plastic_dislotwin_sbResistance(maxNinstance),                        source=0.0_pReal)
 allocate(plastic_dislotwin_sbVelocity(maxNinstance),                          source=0.0_pReal)
 allocate(plastic_dislotwin_sbQedge(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_dislotwin_SFE_0K(maxNinstance),                              source=0.0_pReal)
 allocate(plastic_dislotwin_dSFE_dT(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_dislotwin_dipoleFormationFactor(maxNinstance),               source=1.0_pReal) !should be on by default
 allocate(plastic_dislotwin_deltaG(maxNinstance),                              source=0.0_pReal)
 allocate(plastic_dislotwin_Cmfptrans(maxNinstance),                           source=0.0_pReal)
 allocate(plastic_dislotwin_Cthresholdtrans(maxNinstance),                     source=0.0_pReal)
 allocate(plastic_dislotwin_transStackHeight(maxNinstance),                    source=0.0_pReal)
 allocate(plastic_dislotwin_rhoEdge0(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(plastic_dislotwin_rhoEdgeDip0(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(plastic_dislotwin_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_burgersPerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_burgersPerTransFamily(lattice_maxNtransFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_QedgePerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_v0PerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_tau_peierlsPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_pPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(plastic_dislotwin_qPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(plastic_dislotwin_Ndot0PerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_Ndot0PerTransFamily(lattice_maxNtransFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_rPerTwinFamily(lattice_maxNtwinFamily,maxNinstance),source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_SlipTrans(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_TransSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_interaction_TransTrans(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_sbSv(6,6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_lamellarsizePerTransFamily(lattice_maxNtransFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(plastic_dislotwin_sPerTransFamily(lattice_maxNtransFamily,maxNinstance),source=0.0_pReal) 


 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     if (phase_plasticity(phase) == PLASTICITY_DISLOTWIN_ID) then
       Nchunks_SlipFamilies  = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_TwinFamilies  = count(lattice_NtwinSystem(:,phase) > 0_pInt)
       Nchunks_TransFamilies = count(lattice_NtransSystem(:,phase)> 0_pInt)
       Nchunks_SlipSlip   = maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_SlipTwin   = maxval(lattice_interactionSlipTwin(:,:,phase))
       Nchunks_TwinSlip   = maxval(lattice_interactionTwinSlip(:,:,phase))
       Nchunks_TwinTwin   = maxval(lattice_interactionTwinTwin(:,:,phase))
       Nchunks_SlipTrans  = maxval(lattice_interactionSlipTrans(:,:,phase))
       Nchunks_TransSlip  = maxval(lattice_interactionTransSlip(:,:,phase))
       Nchunks_TransTrans = maxval(lattice_interactionTransTrans(:,:,phase))
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       if(allocated(tempPerTwin)) deallocate(tempPerTwin)
       if(allocated(tempPerTrans)) deallocate(tempPerTrans)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
       allocate(tempPerTwin(Nchunks_TwinFamilies))
       allocate(tempPerTrans(Nchunks_TransFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif

   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_DISLOTWIN_ID) then          ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('edge_density')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = edge_density_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('dipole_density')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = dipole_density_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shear_rate_slip','shearrate_slip')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = shear_rate_slip_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('accumulated_shear_slip')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = accumulated_shear_slip_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('mfp_slip')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = mfp_slip_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolved_stress_slip')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = resolved_stress_slip_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('threshold_stress_slip')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = threshold_stress_slip_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('edge_dipole_distance')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = edge_dipole_distance_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('stress_exponent')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = stress_exponent_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('twin_fraction')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = twin_fraction_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shear_rate_twin','shearrate_twin')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = shear_rate_twin_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('accumulated_shear_twin')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = accumulated_shear_twin_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('mfp_twin')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = mfp_twin_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolved_stress_twin')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = resolved_stress_twin_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('threshold_stress_twin')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = threshold_stress_twin_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolved_stress_shearband')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = resolved_stress_shearband_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shear_rate_shearband','shearrate_shearband')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = shear_rate_shearband_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('sb_eigenvalues')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = sb_eigenvalues_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('sb_eigenvectors')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = sb_eigenvectors_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('stress_trans_fraction')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = stress_trans_fraction_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('strain_trans_fraction')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = strain_trans_fraction_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('trans_fraction','total_trans_fraction')
             plastic_dislotwin_Noutput(instance) = plastic_dislotwin_Noutput(instance) + 1_pInt
             plastic_dislotwin_outputID(plastic_dislotwin_Noutput(instance),instance) = trans_fraction_ID
             plastic_dislotwin_output(plastic_dislotwin_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip system families
       case ('nslip')
         if (chunkPos(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         if (chunkPos(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           plastic_dislotwin_Nslip(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('rhoedge0','rhoedgedip0','slipburgers','qedge','v0','clambdaslip','tau_peierls','p_slip','q_slip')
         do j = 1_pInt, Nchunks_SlipFamilies
           tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('rhoedge0')
             plastic_dislotwin_rhoEdge0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('rhoedgedip0')
             plastic_dislotwin_rhoEdgeDip0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('slipburgers')
             plastic_dislotwin_burgersPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('qedge')
             plastic_dislotwin_QedgePerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('v0')
             plastic_dislotwin_v0PerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('clambdaslip')
             plastic_dislotwin_CLambdaSlipPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('tau_peierls')
             if (lattice_structure(phase) /= LATTICE_bcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for non-bcc ('//PLASTICITY_DISLOTWIN_label//')')
             plastic_dislotwin_tau_peierlsPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('p_slip')
             plastic_dislotwin_pPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('q_slip')
             plastic_dislotwin_qPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on slip number of twin families
       case ('ntwin')
         if (chunkPos(1) < Nchunks_TwinFamilies + 1_pInt) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         if (chunkPos(1) > Nchunks_TwinFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         Nchunks_TwinFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_TwinFamilies
             plastic_dislotwin_Ntwin(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('ndot0_twin','twinsize','twinburgers','r_twin')
         do j = 1_pInt, Nchunks_TwinFamilies
           tempPerTwin(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('ndot0_twin')
             if (lattice_structure(phase) == LATTICE_fcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for fcc ('//PLASTICITY_DISLOTWIN_label//')')
             plastic_dislotwin_Ndot0PerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('twinsize')
             plastic_dislotwin_twinsizePerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('twinburgers')
             plastic_dislotwin_burgersPerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
           case ('r_twin')
             plastic_dislotwin_rPerTwinFamily(1:Nchunks_TwinFamilies,instance) = tempPerTwin(1:Nchunks_TwinFamilies)
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of transformation system families
       case ('ntrans')
         if (chunkPos(1) < Nchunks_TransFamilies + 1_pInt) &
           call IO_warning(53_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         if (chunkPos(1) > Nchunks_TransFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         Nchunks_TransFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_TransFamilies
           plastic_dislotwin_Ntrans(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('ndot0_trans','lamellarsize','transburgers','s_trans')
         do j = 1_pInt, Nchunks_TransFamilies
           tempPerTrans(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('ndot0_trans')
             if (lattice_structure(phase) == LATTICE_fcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for fcc ('//PLASTICITY_DISLOTWIN_label//')')
             plastic_dislotwin_Ndot0PerTransFamily(1:Nchunks_TransFamilies,instance) = tempPerTrans(1:Nchunks_TransFamilies)
           case ('lamellarsize')
             plastic_dislotwin_lamellarsizePerTransFamily(1:Nchunks_TransFamilies,instance) = tempPerTrans(1:Nchunks_TransFamilies)
           case ('transburgers')
             plastic_dislotwin_burgersPerTransFamily(1:Nchunks_TransFamilies,instance) = tempPerTrans(1:Nchunks_TransFamilies)
           case ('s_trans')
             plastic_dislotwin_sPerTransFamily(1:Nchunks_TransFamilies,instance) = tempPerTrans(1:Nchunks_TransFamilies)
         end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of interactions
       case ('interaction_slipslip','interactionslipslip')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           plastic_dislotwin_interaction_SlipSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_sliptwin','interactionsliptwin')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           plastic_dislotwin_interaction_SlipTwin(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_twinslip','interactiontwinslip')
         if (chunkPos(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           plastic_dislotwin_interaction_TwinSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_twintwin','interactiontwintwin')
         if (chunkPos(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           plastic_dislotwin_interaction_TwinTwin(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_sliptrans','interactionsliptrans')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipTrans) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_SlipTrans
           plastic_dislotwin_interaction_SlipTrans(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_transslip','interactiontransslip')
         if (chunkPos(1) < 1_pInt + Nchunks_TransSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TransSlip
           plastic_dislotwin_interaction_TransSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('interaction_transtrans','interactiontranstrans')
         if (chunkPos(1) < 1_pInt + Nchunks_TransTrans) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TransTrans
           plastic_dislotwin_interaction_TransTrans(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters independent of number of slip/twin/trans systems
       case ('grainsize')
         plastic_dislotwin_GrainSize(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('maxtwinfraction')
         plastic_dislotwin_MaxTwinFraction(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('p_shearband')
         plastic_dislotwin_pShearBand(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('q_shearband')
         plastic_dislotwin_qShearBand(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('d0')
         plastic_dislotwin_D0(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('qsd')
         plastic_dislotwin_Qsd(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_rho')
         plastic_dislotwin_aTolRho(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_twinfrac')
         plastic_dislotwin_aTolTwinFrac(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_transfrac')
         plastic_dislotwin_aTolTransFrac(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cmfptwin')
         plastic_dislotwin_Cmfptwin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cthresholdtwin')
         plastic_dislotwin_Cthresholdtwin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('solidsolutionstrength')
         plastic_dislotwin_SolidSolutionStrength(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('l0_twin')
         plastic_dislotwin_L0_twin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('l0_trans')
         plastic_dislotwin_L0_trans(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('xc_twin')
              plastic_dislotwin_xc_twin(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('xc_trans')
              plastic_dislotwin_xc_trans(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('vcrossslip')
              plastic_dislotwin_VcrossSlip(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cedgedipmindistance')
         plastic_dislotwin_CEdgeDipMinDistance(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('catomicvolume')
         plastic_dislotwin_CAtomicVolume(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('sfe_0k')
         plastic_dislotwin_SFE_0K(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('dsfe_dt')
         plastic_dislotwin_dSFE_dT(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('dipoleformationfactor')
         plastic_dislotwin_dipoleFormationFactor(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('shearbandresistance')
         plastic_dislotwin_sbResistance(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('shearbandvelocity')
         plastic_dislotwin_sbVelocity(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('qedgepersbsystem')
         plastic_dislotwin_sbQedge(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('deltag')
         plastic_dislotwin_deltaG(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cmfptrans')
         plastic_dislotwin_Cmfptrans(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cthresholdtrans')
         plastic_dislotwin_Cthresholdtrans(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('transstackheight')
         plastic_dislotwin_transStackHeight(instance) = IO_floatValue(line,chunkPos,2_pInt)
     end select
   endif; endif
 enddo parsingFile
 
 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
    myPhase: if (phase_plasticity(phase) == PLASTICITY_dislotwin_ID) then
      instance = phase_plasticityInstance(phase)

      if (sum(plastic_dislotwin_Nslip(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Nslip ('//PLASTICITY_DISLOTWIN_label//')')
      if (sum(plastic_dislotwin_Ntwin(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Ntwin ('//PLASTICITY_DISLOTWIN_label//')')
      if (sum(plastic_dislotwin_Ntrans(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Ntrans ('//PLASTICITY_DISLOTWIN_label//')')
      do f = 1_pInt,lattice_maxNslipFamily
        if (plastic_dislotwin_Nslip(f,instance) > 0_pInt) then
          if (plastic_dislotwin_rhoEdge0(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOTWIN_label//')')
          if (plastic_dislotwin_rhoEdgeDip0(f,instance) < 0.0_pReal) & 
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOTWIN_label//')')
          if (plastic_dislotwin_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOTWIN_label//')')
          if (plastic_dislotwin_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOTWIN_label//')')
          if (plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOTWIN_label//')')
        endif
      enddo
      do f = 1_pInt,lattice_maxNtwinFamily
        if (plastic_dislotwin_Ntwin(f,instance) > 0_pInt) then
          if (plastic_dislotwin_burgersPerTwinFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='twinburgers ('//PLASTICITY_DISLOTWIN_label//')')
          if (plastic_dislotwin_Ndot0PerTwinFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='ndot0_twin ('//PLASTICITY_DISLOTWIN_label//')')
        endif
      enddo
      if (plastic_dislotwin_CAtomicVolume(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
      if (plastic_dislotwin_D0(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOTWIN_label//')')
      if (plastic_dislotwin_Qsd(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOTWIN_label//')')
      if (sum(plastic_dislotwin_Ntwin(:,instance)) > 0_pInt) then
        if (abs(plastic_dislotwin_SFE_0K(instance))  <= tiny(0.0_pReal) .and. &
            abs(plastic_dislotwin_dSFE_dT(instance)) <= tiny(0.0_pReal) .and. &
                            lattice_structure(phase) == LATTICE_fcc_ID) &
          call IO_error(211_pInt,el=instance,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
        if (plastic_dislotwin_aTolRho(instance) <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOTWIN_label//')')   
        if (plastic_dislotwin_aTolTwinFrac(instance) <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolTwinFrac ('//PLASTICITY_DISLOTWIN_label//')')
      endif
      if (sum(plastic_dislotwin_Ntrans(:,instance)) > 0_pInt) then
        if (abs(plastic_dislotwin_SFE_0K(instance))  <= tiny(0.0_pReal) .and. &
            abs(plastic_dislotwin_dSFE_dT(instance)) <= tiny(0.0_pReal) .and. &
                            lattice_structure(phase) == LATTICE_fcc_ID) &
          call IO_error(211_pInt,el=instance,ext_msg='SFE0K ('//PLASTICITY_DISLOTWIN_label//')')
        if (plastic_dislotwin_aTolTransFrac(instance) <= 0.0_pReal) &
          call IO_error(211_pInt,el=instance,ext_msg='aTolTransFrac ('//PLASTICITY_DISLOTWIN_label//')')
      endif
      if (plastic_dislotwin_sbResistance(instance) < 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='sbResistance ('//PLASTICITY_DISLOTWIN_label//')')
      if (plastic_dislotwin_sbVelocity(instance) < 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='sbVelocity ('//PLASTICITY_DISLOTWIN_label//')')
      if (plastic_dislotwin_sbVelocity(instance) > 0.0_pReal .and. &
          plastic_dislotwin_pShearBand(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='pShearBand ('//PLASTICITY_DISLOTWIN_label//')')
      if (abs(plastic_dislotwin_dipoleFormationFactor(instance)) >  tiny(0.0_pReal) .and. &
               plastic_dislotwin_dipoleFormationFactor(instance) /= 1.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='dipoleFormationFactor ('//PLASTICITY_DISLOTWIN_label//')')
      if (plastic_dislotwin_sbVelocity(instance) > 0.0_pReal .and. &
          plastic_dislotwin_qShearBand(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='qShearBand ('//PLASTICITY_DISLOTWIN_label//')')

!--------------------------------------------------------------------------------------------------
! Determine total number of active slip or twin systems
      plastic_dislotwin_Nslip(:,instance) = min(lattice_NslipSystem(:,phase),plastic_dislotwin_Nslip(:,instance))
      plastic_dislotwin_Ntwin(:,instance) = min(lattice_NtwinSystem(:,phase),plastic_dislotwin_Ntwin(:,instance))
      plastic_dislotwin_Ntrans(:,instance)= min(lattice_NtransSystem(:,phase),plastic_dislotwin_Ntrans(:,instance))
      plastic_dislotwin_totalNslip(instance) = sum(plastic_dislotwin_Nslip(:,instance))
      plastic_dislotwin_totalNtwin(instance) = sum(plastic_dislotwin_Ntwin(:,instance))
      plastic_dislotwin_totalNtrans(instance) = sum(plastic_dislotwin_Ntrans(:,instance))
   endif myPhase
 enddo sanityChecks
 
!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip  = maxval(plastic_dislotwin_totalNslip)
 maxTotalNtwin  = maxval(plastic_dislotwin_totalNtwin)
 maxTotalNtrans = maxval(plastic_dislotwin_totalNtrans)

 allocate(plastic_dislotwin_burgersPerSlipSystem(maxTotalNslip, maxNinstance),    source=0.0_pReal)
 allocate(plastic_dislotwin_burgersPerTwinSystem(maxTotalNtwin, maxNinstance),    source=0.0_pReal)
 allocate(plastic_dislotwin_burgersPerTransSystem(maxTotalNtrans, maxNinstance),  source=0.0_pReal)
 allocate(plastic_dislotwin_QedgePerSlipSystem(maxTotalNslip, maxNinstance),      source=0.0_pReal)
 allocate(plastic_dislotwin_v0PerSlipSystem(maxTotalNslip, maxNinstance),         source=0.0_pReal)
 allocate(plastic_dislotwin_Ndot0PerTwinSystem(maxTotalNtwin, maxNinstance),      source=0.0_pReal)
 allocate(plastic_dislotwin_Ndot0PerTransSystem(maxTotalNtrans, maxNinstance),    source=0.0_pReal)
 allocate(plastic_dislotwin_tau_r_twin(maxTotalNtwin, maxNinstance),              source=0.0_pReal)
 allocate(plastic_dislotwin_tau_r_trans(maxTotalNtrans, maxNinstance),            source=0.0_pReal)
 allocate(plastic_dislotwin_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance),   source=0.0_pReal)
 allocate(plastic_dislotwin_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance),source=0.0_pReal)
 allocate(plastic_dislotwin_lamellarsizePerTransSystem(maxTotalNtrans, maxNinstance),source=0.0_pReal)

 allocate(plastic_dislotwin_interactionMatrix_SlipSlip(maxval(plastic_dislotwin_totalNslip),&  ! slip resistance from slip activity
                                                            maxval(plastic_dislotwin_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_SlipTwin(maxval(plastic_dislotwin_totalNslip),&  ! slip resistance from twin activity
                                                            maxval(plastic_dislotwin_totalNtwin),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_TwinSlip(maxval(plastic_dislotwin_totalNtwin),&  ! twin resistance from slip activity
                                                            maxval(plastic_dislotwin_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_TwinTwin(maxval(plastic_dislotwin_totalNtwin),&  ! twin resistance from twin activity
                                                            maxval(plastic_dislotwin_totalNtwin),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_SlipTrans(maxval(plastic_dislotwin_totalNslip),&  ! slip resistance from trans activity
                                                            maxval(plastic_dislotwin_totalNtrans),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_TransSlip(maxval(plastic_dislotwin_totalNtrans),&  ! trans resistance from slip activity
                                                            maxval(plastic_dislotwin_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_interactionMatrix_TransTrans(maxval(plastic_dislotwin_totalNtrans),&  ! trans resistance from trans activity
                                                            maxval(plastic_dislotwin_totalNtrans),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_dislotwin_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(plastic_dislotwin_projectionMatrix_Trans(maxTotalNtrans,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(plastic_dislotwin_Ctwin66(6,6,maxTotalNtwin,maxNinstance),              source=0.0_pReal)
 allocate(plastic_dislotwin_Ctwin3333(3,3,3,3,maxTotalNtwin,maxNinstance),        source=0.0_pReal)
 allocate(plastic_dislotwin_Ctrans66(6,6,maxTotalNtrans,maxNinstance),            source=0.0_pReal)
 allocate(plastic_dislotwin_Ctrans3333(3,3,3,3,maxTotalNtrans,maxNinstance),      source=0.0_pReal)
 
 allocate(state(maxNinstance))
 allocate(state0(maxNinstance))
 allocate(dotState(maxNinstance))

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
    myPhase2: if (phase_plasticity(phase) == PLASTICITY_dislotwin_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_plasticityInstance(phase)
 
     ns = plastic_dislotwin_totalNslip(instance)
     nt = plastic_dislotwin_totalNtwin(instance)
     nr = plastic_dislotwin_totalNtrans(instance)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,plastic_dislotwin_Noutput(instance)
       select case(plastic_dislotwin_outputID(o,instance))
         case(edge_density_ID, &
              dipole_density_ID, &
              shear_rate_slip_ID, &
              accumulated_shear_slip_ID, &
              mfp_slip_ID, &
              resolved_stress_slip_ID, &
              threshold_stress_slip_ID, &
              edge_dipole_distance_ID, &
              stress_exponent_ID &
              )
           mySize = ns
         case(twin_fraction_ID, &
              shear_rate_twin_ID, &
              accumulated_shear_twin_ID, &
              mfp_twin_ID, &
              resolved_stress_twin_ID, &
              threshold_stress_twin_ID &
              )
           mySize = nt
         case(resolved_stress_shearband_ID, &
              shear_rate_shearband_ID &
              )
           mySize = 6_pInt
         case(sb_eigenvalues_ID)
           mySize = 3_pInt  
         case(sb_eigenvectors_ID)
           mySize = 9_pInt  
         case(stress_trans_fraction_ID, &
              strain_trans_fraction_ID, &
              trans_fraction_ID &
              )
           mySize = nr
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          plastic_dislotwin_sizePostResult(o,instance) = mySize
          plastic_dislotwin_sizePostResults(instance)  = plastic_dislotwin_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

!--------------------------------------------------------------------------------------------------
! allocate state arrays

     sizeDotState     = int(size(['rhoEdge     ','rhoEdgeDip  ','accshearslip']),pInt) * ns &
                      + int(size(['twinFraction','accsheartwin']),pInt) * nt &
                      + int(size(['stressTransFraction','strainTransFraction']),pInt) * nr
     sizeDeltaState   =  0_pInt
     sizeState        = sizeDotState &
                      + int(size(['invLambdaSlip     ','invLambdaSlipTwin ','invLambdaSlipTrans',&
                                  'meanFreePathSlip  ','tauSlipThreshold  ']),pInt) * ns &
                      + int(size(['invLambdaTwin   ','meanFreePathTwin','tauTwinThreshold',&
                                  'twinVolume      ']),pInt) * nt &
                      + int(size(['invLambdaTrans   ','meanFreePathTrans','tauTransThreshold', &
                                  'martensiteVolume ']),pInt) * nr
                
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_dislotwin_sizePostResults(instance)
     plasticState(phase)%nSlip = plastic_dislotwin_totalNslip(instance)
     plasticState(phase)%nTwin = plastic_dislotwin_totalNtwin(instance)
     plasticState(phase)%nTrans= plastic_dislotwin_totalNtrans(instance)
     allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     if (.not. analyticJaco) then
       allocate(plasticState(phase)%state_backup      (sizeState,NofMyPhase),     source=0.0_pReal)
       allocate(plasticState(phase)%dotState_backup   (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
     offset_slip = 2_pInt*plasticState(phase)%nslip
     plasticState(phase)%slipRate        => &
         plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nslip,1:NofMyPhase)
     plasticState(phase)%accumulatedSlip => &
         plasticState(phase)%state   (offset_slip+1:offset_slip+plasticState(phase)%nslip,1:NofMyPhase)

    !* Process slip related parameters ------------------------------------------------ 
     slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
       index_myFamily = sum(plastic_dislotwin_Nslip(1:f-1_pInt,instance))                      ! index in truncated slip system list
       slipSystemsLoop: do j = 1_pInt,plastic_dislotwin_Nslip(f,instance)

       !* Burgers vector, 
       !  dislocation velocity prefactor,
       !  mean free path prefactor,
       !  and minimum dipole distance
 
         plastic_dislotwin_burgersPerSlipSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_burgersPerSlipFamily(f,instance)
 
         plastic_dislotwin_QedgePerSlipSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_QedgePerSlipFamily(f,instance)
 
         plastic_dislotwin_v0PerSlipSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_v0PerSlipFamily(f,instance)
 
         plastic_dislotwin_CLambdaSlipPerSlipSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_CLambdaSlipPerSlipFamily(f,instance)
  
       !* Calculation of forest projections for edge dislocations
       !* Interaction matrices
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_dislotwin_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (slip)
             plastic_dislotwin_forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,phase))+j,phase), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,phase))+k,phase)))
             plastic_dislotwin_interactionMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
  
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(plastic_dislotwin_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Ntwin(o,instance)                                   ! loop over (active) systems in other family (twin)
             plastic_dislotwin_interactionMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo

         do o = 1_pInt,lattice_maxNtransFamily
           index_otherFamily = sum(plastic_dislotwin_Ntrans(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Ntrans(o,instance)                                  ! loop over (active) systems in other family (trans)
             plastic_dislotwin_interactionMatrix_SlipTrans(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_SlipTrans(lattice_interactionSlipTrans( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtransSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
  
       enddo slipSystemsLoop
     enddo slipFamiliesLoop
  
    !* Process twin related parameters ------------------------------------------------
     twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
       index_myFamily = sum(plastic_dislotwin_Ntwin(1:f-1_pInt,instance))                      ! index in truncated twin system list
       twinSystemsLoop: do j = 1_pInt,plastic_dislotwin_Ntwin(f,instance)
 
       !* Burgers vector,
       !  nucleation rate prefactor,
       !  and twin size
 
         plastic_dislotwin_burgersPerTwinSystem(index_myFamily+j,instance)  = &
         plastic_dislotwin_burgersPerTwinFamily(f,instance)

         plastic_dislotwin_Ndot0PerTwinSystem(index_myFamily+j,instance)  = &
         plastic_dislotwin_Ndot0PerTwinFamily(f,instance)

         plastic_dislotwin_twinsizePerTwinSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_twinsizePerTwinFamily(f,instance)
 
       !* Rotate twin elasticity matrices
         index_otherFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase))                             ! index in full lattice twin list
         do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt
           do p = 1_pInt,3_pInt; do q = 1_pInt,3_pInt; do r = 1_pInt,3_pInt; do s = 1_pInt,3_pInt
             plastic_dislotwin_Ctwin3333(l,m,n,o,index_myFamily+j,instance) = &
             plastic_dislotwin_Ctwin3333(l,m,n,o,index_myFamily+j,instance) + &
               lattice_C3333(p,q,r,s,instance) * &
               lattice_Qtwin(l,p,index_otherFamily+j,phase) * &
               lattice_Qtwin(m,q,index_otherFamily+j,phase) * &
               lattice_Qtwin(n,r,index_otherFamily+j,phase) * &
               lattice_Qtwin(o,s,index_otherFamily+j,phase)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo; enddo
         plastic_dislotwin_Ctwin66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(plastic_dislotwin_Ctwin3333(1:3,1:3,1:3,1:3,index_myFamily+j,instance))
 
      !* Interaction matrices
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_dislotwin_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (slip)
             plastic_dislotwin_interactionMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
 
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(plastic_dislotwin_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Ntwin(o,instance)                                   ! loop over (active) systems in other family (twin)
             plastic_dislotwin_interactionMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
 
       enddo twinSystemsLoop
     enddo twinFamiliesLoop

    !* Process transformation related parameters ------------------------------------------------
     transFamiliesLoop: do f = 1_pInt,lattice_maxNtransFamily
       index_myFamily = sum(plastic_dislotwin_Ntrans(1:f-1_pInt,instance))                      ! index in truncated trans system list
       transSystemsLoop: do j = 1_pInt,plastic_dislotwin_Ntrans(f,instance)

       !* Burgers vector,
       !  nucleation rate prefactor,
       !  and martensite size

         plastic_dislotwin_burgersPerTransSystem(index_myFamily+j,instance)  = &
         plastic_dislotwin_burgersPerTransFamily(f,instance)

         plastic_dislotwin_Ndot0PerTransSystem(index_myFamily+j,instance)  = &
         plastic_dislotwin_Ndot0PerTransFamily(f,instance)

         plastic_dislotwin_lamellarsizePerTransSystem(index_myFamily+j,instance) = &
         plastic_dislotwin_lamellarsizePerTransFamily(f,instance)

       !* Rotate trans elasticity matrices
       index_otherFamily = sum(lattice_NtransSystem(1:f-1_pInt,phase))                               ! index in full lattice trans list
         do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt
           do p = 1_pInt,3_pInt; do q = 1_pInt,3_pInt; do r = 1_pInt,3_pInt; do s = 1_pInt,3_pInt
             plastic_dislotwin_Ctrans3333(l,m,n,o,index_myFamily+j,instance) = &
             plastic_dislotwin_Ctrans3333(l,m,n,o,index_myFamily+j,instance) + &
               lattice_trans_C3333(p,q,r,s,instance) * &
               lattice_Qtrans(l,p,index_otherFamily+j,phase) * &
               lattice_Qtrans(m,q,index_otherFamily+j,phase) * &
               lattice_Qtrans(n,r,index_otherFamily+j,phase) * &
               lattice_Qtrans(o,s,index_otherFamily+j,phase)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo; enddo
         plastic_dislotwin_Ctrans66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(plastic_dislotwin_Ctrans3333(1:3,1:3,1:3,1:3,index_myFamily+j,instance))

       !* Interaction matrices
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_dislotwin_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (slip)
             plastic_dislotwin_interactionMatrix_TransSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_TransSlip(lattice_interactionTransSlip( &
                                                                 sum(lattice_NtransSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo

         do o = 1_pInt,lattice_maxNtransFamily
           index_otherFamily = sum(plastic_dislotwin_Ntrans(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_dislotwin_Ntrans(o,instance)                                  ! loop over (active) systems in other family (trans)
             plastic_dislotwin_interactionMatrix_TransTrans(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_dislotwin_interaction_TransTrans(lattice_interactionTransTrans( &
                                                                 sum(lattice_NtransSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtransSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo

       !* Projection matrices for shear from slip systems to fault-band (twin) systems for strain-induced martensite nucleation
         select case(trans_lattice_structure(phase))
           case (LATTICE_bcc_ID)
             do o = 1_pInt,lattice_maxNtransFamily
               index_otherFamily = sum(plastic_dislotwin_Nslip(1:o-1_pInt,instance))
               do k = 1_pInt,plastic_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (trans)
                 plastic_dislotwin_projectionMatrix_Trans(index_myFamily+j,index_otherFamily+k,instance) = &
                       lattice_projectionTrans( sum(lattice_NtransSystem(1:f-1,phase))+j, &
                                                sum(lattice_NslipSystem(1:o-1,phase))+k, phase)
             enddo; enddo
         end select 

       enddo transSystemsLoop
     enddo transFamiliesLoop

     startIndex=1_pInt
     endIndex=ns
     state(instance)%rhoEdge=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%rhoEdge=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%rhoEdge=>plasticState(phase)%dotState(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%rhoEdgeDip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%rhoEdgeDip=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%rhoEdgeDip=>plasticState(phase)%dotState(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%accshear_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%accshear_slip=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%accshear_slip=>plasticState(phase)%dotState(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%twinFraction=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%twinFraction=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%twinFraction=>plasticState(phase)%dotState(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%accshear_twin=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%accshear_twin=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%accshear_twin=>plasticState(phase)%dotState(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%stressTransFraction=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%stressTransFraction=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%stressTransFraction=>plasticState(phase)%dotState(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%strainTransFraction=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%strainTransFraction=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%strainTransFraction=>plasticState(phase)%dotState(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%invLambdaSlip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaSlip=>plasticState(phase)%state0(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%invLambdaSlipTwin=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaSlipTwin=>plasticState(phase)%state0(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%invLambdaTwin=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaTwin=>plasticState(phase)%state0(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%invLambdaSlipTrans=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaSlipTrans=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%invLambdaTrans=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaTrans=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%mfp_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%mfp_slip=>plasticState(phase)%state0(startIndex:endIndex,:)
     
     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%mfp_twin=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%mfp_twin=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%mfp_trans=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%mfp_trans=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+ns
     state(instance)%threshold_stress_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%threshold_stress_slip=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%threshold_stress_twin=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%threshold_stress_twin=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%threshold_stress_trans=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%threshold_stress_trans=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nt
     state(instance)%twinVolume=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%twinVolume=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1
     endIndex=endIndex+nr
     state(instance)%martensiteVolume=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%martensiteVolume=>plasticState(phase)%state0(startIndex:endIndex,:)
     
     call plastic_dislotwin_stateInit(phase,instance)
     call plastic_dislotwin_aTolState(phase,instance)
   endif myPhase2
 
 enddo initializeInstances
end subroutine plastic_dislotwin_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_stateInit(ph,instance)
 use math, only: &
   pi
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_mu
 use material, only: &
   plasticState

 implicit none
 integer(pInt), intent(in) :: &
   instance, &                                                                                      !< number specifying the instance of the plasticity
   ph 

  real(pReal), dimension(plasticState(ph)%sizeState) :: tempState

 integer(pInt) :: i,j,f,ns,nt,nr, index_myFamily
 real(pReal), dimension(plastic_dislotwin_totalNslip(instance)) :: &
   rhoEdge0, &
   rhoEdgeDip0, &
   invLambdaSlip0, &
   MeanFreePathSlip0, &
   tauSlipThreshold0
 real(pReal), dimension(plastic_dislotwin_totalNtwin(instance)) :: &
   MeanFreePathTwin0,TwinVolume0
 real(pReal), dimension(plastic_dislotwin_totalNtrans(instance)) :: &
   MeanFreePathTrans0,MartensiteVolume0
 tempState = 0.0_pReal
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily   = sum(plastic_dislotwin_Nslip(1:f-1_pInt,instance))                        ! index in truncated slip system list
   rhoEdge0(index_myFamily+1_pInt: &
            index_myFamily+plastic_dislotwin_Nslip(f,instance)) = &
     plastic_dislotwin_rhoEdge0(f,instance)
   rhoEdgeDip0(index_myFamily+1_pInt: &
               index_myFamily+plastic_dislotwin_Nslip(f,instance)) = &
     plastic_dislotwin_rhoEdgeDip0(f,instance)
 enddo
 
 tempState(1_pInt:ns)           = rhoEdge0
 tempState(ns+1_pInt:2_pInt*ns) = rhoEdgeDip0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent slip microstructural variables
 forall (i = 1_pInt:ns) &
   invLambdaSlip0(i) = sqrt(dot_product((rhoEdge0+rhoEdgeDip0),plastic_dislotwin_forestProjectionEdge(1:ns,i,instance)))/ &
                       plastic_dislotwin_CLambdaSlipPerSlipSystem(i,instance)
 tempState(3_pInt*ns+2_pInt*nt+2_pInt*nr+1:4_pInt*ns+2_pInt*nt+2_pInt*nr) = invLambdaSlip0
 
 forall (i = 1_pInt:ns) &
   MeanFreePathSlip0(i) = &
     plastic_dislotwin_GrainSize(instance)/(1.0_pReal+invLambdaSlip0(i)*plastic_dislotwin_GrainSize(instance))
 tempState(6_pInt*ns+3_pInt*nt+3_pInt*nr+1:7_pInt*ns+3_pInt*nt+3_pInt*nr) = MeanFreePathSlip0
 
 forall (i = 1_pInt:ns) &
   tauSlipThreshold0(i) = &
     lattice_mu(ph)*plastic_dislotwin_burgersPerSlipSystem(i,instance) * &
     sqrt(dot_product((rhoEdge0+rhoEdgeDip0),plastic_dislotwin_interactionMatrix_SlipSlip(i,1:ns,instance)))

 tempState(7_pInt*ns+4_pInt*nt+4_pInt*nr+1:8_pInt*ns+4_pInt*nt+4_pInt*nr) = tauSlipThreshold0

!--------------------------------------------------------------------------------------------------
! initialize dependent twin microstructural variables
 forall (j = 1_pInt:nt) &
   MeanFreePathTwin0(j) = plastic_dislotwin_GrainSize(instance)
 tempState(7_pInt*ns+3_pInt*nt+3_pInt*nr+1_pInt:7_pInt*ns+4_pInt*nt+3_pInt*nr) = MeanFreePathTwin0
 
 forall (j = 1_pInt:nt) &
   TwinVolume0(j) = &
     (pi/4.0_pReal)*plastic_dislotwin_twinsizePerTwinSystem(j,instance)*MeanFreePathTwin0(j)**(2.0_pReal)
 tempState(8_pInt*ns+5_pInt*nt+5_pInt*nr+1_pInt:8_pInt*ns+6_pInt*nt+5_pInt*nr) = TwinVolume0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent trans microstructural variables
 forall (j = 1_pInt:nr) &
   MeanFreePathTrans0(j) = plastic_dislotwin_GrainSize(instance)
 tempState(7_pInt*ns+4_pInt*nt+3_pInt*nr+1_pInt:7_pInt*ns+4_pInt*nt+4_pInt*nr) = MeanFreePathTrans0
 
 forall (j = 1_pInt:nr) &
   MartensiteVolume0(j) = &
     (pi/4.0_pReal)*plastic_dislotwin_lamellarsizePerTransSystem(j,instance)*MeanFreePathTrans0(j)**(2.0_pReal)
 tempState(8_pInt*ns+6_pInt*nt+5_pInt*nr+1_pInt:8_pInt*ns+6_pInt*nt+6_pInt*nr) = MartensiteVolume0

plasticState(ph)%state0 = spread(tempState,2,size(plasticState(ph)%state(1,:)))

end subroutine plastic_dislotwin_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_aTolState(ph,instance)
 use material, only: &
  plasticState

 implicit none
 integer(pInt), intent(in) ::  &
   ph, &
   instance                                                                                         ! number specifying the current instance of the plasticity
 
 integer(pInt) :: ns, nt, nr
 
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance) 

 ! Tolerance state for dislocation densities
 plasticState(ph)%aTolState(1_pInt: &
                            2_pInt*ns) = plastic_dislotwin_aTolRho(instance)

 ! Tolerance state for accumulated shear due to slip 
 plasticState(ph)%aTolState(2_pInt*ns+1_pInt: &
                            3_pInt*ns)=1.0e6_pReal
 
 ! Tolerance state for twin volume fraction
 plasticState(ph)%aTolState(3_pInt*ns+1_pInt: &
                            3_pInt*ns+nt) = plastic_dislotwin_aTolTwinFrac(instance)

 ! Tolerance state for accumulated shear due to twin
 plasticState(ph)%aTolState(3_pInt*ns+nt+1_pInt: &
                            3_pInt*ns+2_pInt*nt) = 1.0e6_pReal

! Tolerance state for stress-assisted martensite volume fraction
 plasticState(ph)%aTolState(3_pInt*ns+2_pInt*nt+1_pInt: &
                            3_pInt*ns+2_pInt*nt+nr) = plastic_dislotwin_aTolTransFrac(instance)

! Tolerance state for strain-induced martensite volume fraction
 plasticState(ph)%aTolState(3_pInt*ns+2_pInt*nt+nr+1_pInt: &
                            3_pInt*ns+2_pInt*nt+2_pInt*nr) = plastic_dislotwin_aTolTransFrac(instance)

end subroutine plastic_dislotwin_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_homogenizedC(ipc,ip,el)
 use material, only: &
  phase_plasticityInstance, &
  phaseAt, phasememberAt
 use lattice, only: &
  lattice_C66
 
  implicit none
  real(pReal), dimension(6,6) :: &
    plastic_dislotwin_homogenizedC
  integer(pInt), intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

 integer(pInt) :: instance,ns,nt,nr,i, &
                  ph, &
                  of
 real(pReal) :: sumf, sumftr

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:nt,of))           ! safe for nt == 0

 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:nr,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:nr,of))

 !* Homogenized elasticity matrix
 plastic_dislotwin_homogenizedC = (1.0_pReal-sumf-sumftr)*lattice_C66(1:6,1:6,ph)
 do i=1_pInt,nt
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                   + state(instance)%twinFraction(i,of)*plastic_dislotwin_Ctwin66(1:6,1:6,i,instance)
 enddo
 do i=1_pInt,nr
    plastic_dislotwin_homogenizedC = plastic_dislotwin_homogenizedC &
                   + (state(instance)%stressTransFraction(i,of) + state(instance)%strainTransFraction(i,of))*&
                     plastic_dislotwin_Ctrans66(1:6,1:6,i,instance)
 enddo

 end function plastic_dislotwin_homogenizedC
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_microstructure(temperature,ipc,ip,el)
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   !plasticState, &   !!!!delete
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_mu, &
   lattice_nu

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 

 integer(pInt) :: &
   instance, &
   ns,nt,nr,s,t,r, &
   ph, &
   of
 real(pReal) :: &
   sumf,sfe,x0,sumftr
 real(pReal), dimension(plastic_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: fOverStacksize
 real(pReal), dimension(plastic_dislotwin_totalNtrans(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   ftransOverLamellarSize

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:nt,of)) ! safe for nt == 0
 
 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:nr,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:nr,of))

 !* Stacking fault energy
 sfe = plastic_dislotwin_SFE_0K(instance) + & 
       plastic_dislotwin_dSFE_dT(instance) * Temperature
 
 !* rescaled twin volume fraction for topology
 forall (t = 1_pInt:nt) &
   fOverStacksize(t) = &
     state(instance)%twinFraction(t,of)/plastic_dislotwin_twinsizePerTwinSystem(t,instance)

 !* rescaled trans volume fraction for topology
 forall (r = 1_pInt:nr) &
   ftransOverLamellarSize(r) = &
     (state(instance)%stressTransFraction(r,of)+state(instance)%strainTransFraction(r,of))/&
      plastic_dislotwin_lamellarsizePerTransSystem(r,instance)
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:ns) &
   state(instance)%invLambdaSlip(s,of) = &
     sqrt(dot_product((state(instance)%rhoEdge(1_pInt:ns,of)+state(instance)%rhoEdgeDip(1_pInt:ns,of)),&
                      plastic_dislotwin_forestProjectionEdge(1:ns,s,instance)))/ &
     plastic_dislotwin_CLambdaSlipPerSlipSystem(s,instance)

 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 state(instance)%invLambdaSlipTwin(1_pInt:ns,of) = 0.0_pReal
 if (nt > 0_pInt .and. ns > 0_pInt) &
   state(instance)%invLambdaSlipTwin(1_pInt:ns,of) = &
     matmul(plastic_dislotwin_interactionMatrix_SlipTwin(1:ns,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)
 
 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
 !$OMP CRITICAL (evilmatmul)
 if (nt > 0_pInt) &
   state(instance)%invLambdaTwin(1_pInt:nt,of) = &
     matmul(plastic_dislotwin_interactionMatrix_TwinTwin(1:nt,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)

 !* 1/mean free distance between 2 martensite lamellar from different systems seen by a moving dislocation
 state(instance)%invLambdaSlipTrans(1_pInt:ns,of) = 0.0_pReal
 if (nr > 0_pInt .and. ns > 0_pInt) &
   state(instance)%invLambdaSlipTrans(1_pInt:ns,of) = &
      matmul(plastic_dislotwin_interactionMatrix_SlipTrans(1:ns,1:nr,instance),ftransOverLamellarSize(1:nr))/(1.0_pReal-sumftr)

 !* 1/mean free distance between 2 martensite stacks from different systems seen by a growing martensite (1/lambda_trans)
 if (nr > 0_pInt) &
   state(instance)%invLambdaTrans(1_pInt:nr,of) = &
     matmul(plastic_dislotwin_interactionMatrix_TransTrans(1:nr,1:nr,instance),ftransOverLamellarSize(1:nr))/(1.0_pReal-sumftr)

 !* mean free path between 2 obstacles seen by a moving dislocation
 do s = 1_pInt,ns
    if ((nt > 0_pInt) .or. (nr > 0_pInt)) then
       state(instance)%mfp_slip(s,of) = &
         plastic_dislotwin_GrainSize(instance)/(1.0_pReal+plastic_dislotwin_GrainSize(instance)*&
         (state(instance)%invLambdaSlip(s,of) + &
          state(instance)%invLambdaSlipTwin(s,of) + &
          state(instance)%invLambdaSlipTrans(s,of)))
    else
       state(instance)%mfp_slip(s,of) = &  
         plastic_dislotwin_GrainSize(instance)/&
         (1.0_pReal+plastic_dislotwin_GrainSize(instance)*(state(instance)%invLambdaSlip(s,of))) !!!!!! correct?
    endif
 enddo

 !* mean free path between 2 obstacles seen by a growing twin
 forall (t = 1_pInt:nt) &
   state(instance)%mfp_twin(t,of) = &
     plastic_dislotwin_Cmfptwin(instance)*plastic_dislotwin_GrainSize(instance)/&
       (1.0_pReal+plastic_dislotwin_GrainSize(instance)*state(ph)%invLambdaTwin(t,of))
 
 !* mean free path between 2 obstacles seen by a growing martensite
 forall (r = 1_pInt:nr) &
   state(instance)%mfp_trans(r,of) = &
     plastic_dislotwin_Cmfptrans(instance)*plastic_dislotwin_GrainSize(instance)/&
      (1.0_pReal+plastic_dislotwin_GrainSize(instance)*state(instance)%invLambdaTrans(r,of))

 !* threshold stress for dislocation motion
 forall (s = 1_pInt:ns) &
   state(instance)%threshold_stress_slip(s,of) = &
     lattice_mu(ph)*plastic_dislotwin_burgersPerSlipSystem(s,instance)*&
     sqrt(dot_product((state(instance)%rhoEdge(1_pInt:ns,of)+state(instance)%rhoEdgeDip(1_pInt:ns,of)),&
                      plastic_dislotwin_interactionMatrix_SlipSlip(s,1:ns,instance)))

 !* threshold stress for growing twin
 forall (t = 1_pInt:nt) &
   state(instance)%threshold_stress_twin(t,of) = &
     plastic_dislotwin_Cthresholdtwin(instance)* &
       (sfe/(3.0_pReal*plastic_dislotwin_burgersPerTwinSystem(t,instance)) &
        + 3.0_pReal*plastic_dislotwin_burgersPerTwinSystem(t,instance)*lattice_mu(ph)/&
            (plastic_dislotwin_L0_twin(instance)*plastic_dislotwin_burgersPerSlipSystem(t,instance)) &
        )

 !* threshold stress for growing martensite
 forall (r = 1_pInt:nr) &
   state(instance)%threshold_stress_trans(r,of) = &
     plastic_dislotwin_Cthresholdtrans(instance)* &
       (sfe/(3.0_pReal*plastic_dislotwin_burgersPerTransSystem(r,instance)) &
        + 3.0_pReal*plastic_dislotwin_burgersPerTransSystem(r,instance)*lattice_mu(ph)/&
            (plastic_dislotwin_L0_trans(instance)*plastic_dislotwin_burgersPerSlipSystem(r,instance))&
        + plastic_dislotwin_transStackHeight(instance)*plastic_dislotwin_deltaG(instance)/ &
            (3.0_pReal*plastic_dislotwin_burgersPerTransSystem(r,instance)) &
        )    

 !* final twin volume after growth
 forall (t = 1_pInt:nt) &
   state(instance)%twinVolume(t,of) = &
     (pi/4.0_pReal)*plastic_dislotwin_twinsizePerTwinSystem(t,instance)*&
     state(instance)%mfp_twin(t,of)**(2.0_pReal)

 !* final martensite volume after growth
 forall (r = 1_pInt:nr) &
   state(instance)%martensiteVolume(r,of) = &
     (pi/4.0_pReal)*plastic_dislotwin_lamellarsizePerTransSystem(r,instance)*&
     state(instance)%mfp_trans(r,of)**(2.0_pReal)

 !* equilibrium separation of partial dislocations (twin)
 do t = 1_pInt,nt
   x0 = lattice_mu(ph)*plastic_dislotwin_burgersPerTwinSystem(t,instance)**(2.0_pReal)/&
     (sfe*8.0_pReal*pi)*(2.0_pReal+lattice_nu(ph))/(1.0_pReal-lattice_nu(ph))
   plastic_dislotwin_tau_r_twin(t,instance)= &
        lattice_mu(ph)*plastic_dislotwin_burgersPerTwinSystem(t,instance)/(2.0_pReal*pi)*&     
        (1/(x0+plastic_dislotwin_xc_twin(instance))+cos(pi/3.0_pReal)/x0)
 enddo

 !* equilibrium separation of partial dislocations (trans)
 do r = 1_pInt,nr
   x0 = lattice_mu(ph)*plastic_dislotwin_burgersPerTransSystem(r,instance)**(2.0_pReal)/&
     (sfe*8.0_pReal*pi)*(2.0_pReal+lattice_nu(ph))/(1.0_pReal-lattice_nu(ph))
   plastic_dislotwin_tau_r_trans(r,instance)= &
        lattice_mu(ph)*plastic_dislotwin_burgersPerTransSystem(r,instance)/(2.0_pReal*pi)*&     
        (1/(x0+plastic_dislotwin_xc_trans(instance))+cos(pi/3.0_pReal)/x0)
 enddo

end subroutine plastic_dislotwin_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_eigenValuesVectorsSym, &
   math_tensorproduct33, &
   math_symmetric33, &
   math_mul33x3
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_Strans, &
   lattice_Strans_v, &
   lattice_maxNslipFamily,&
   lattice_maxNtwinFamily, &
   lattice_maxNtransFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NtransSystem, &
   lattice_shearTwin, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID
 
 implicit none
 integer(pInt), intent(in)                  :: ipc,ip,el
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(6),   intent(in)    :: Tstar_v
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(9,9), intent(out)   :: dLp_dTstar99

 integer(pInt) :: instance,ph,of,ns,nt,nr,f,i,j,k,l,m,n,index_myFamily,s1,s2
 real(pReal) :: sumf,sumftr,StressRatio_p,StressRatio_pminus1,StressRatio_r,BoltzmannRatio,DotGamma0,Ndot0_twin,stressRatio, &
    Ndot0_trans,StressRatio_s
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(plastic_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(plastic_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_twin,dgdot_dtautwin,tau_twin
 real(pReal), dimension(plastic_dislotwin_totalNtrans(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_trans,dgdot_dtautrans,tau_trans
 real(pReal), dimension(6) :: gdot_sb,dgdot_dtausb,tau_sb
 real(pReal), dimension(3,3) :: eigVectors, sb_Smatrix
 real(pReal), dimension(3)   :: eigValues, sb_s, sb_m
 logical :: error
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
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! Dislocation glide part
 gdot_slip = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystemsLoop: do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))

      if((abs(tau_slip(j))-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
      !* Stress ratios
        stressRatio =((abs(tau_slip(j))- state(instance)%threshold_stress_slip(j,of))/&
          (plastic_dislotwin_SolidSolutionStrength(instance)+plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))
        StressRatio_p       = stressRatio** plastic_dislotwin_pPerSlipFamily(f,instance)
        StressRatio_pminus1 = stressRatio**(plastic_dislotwin_pPerSlipFamily(f,instance)-1.0_pReal)
      !* Boltzmann ratio
        BoltzmannRatio = plastic_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
      !* Initial shear rates
        DotGamma0 = &
          state(instance)%rhoEdge(j,of)*plastic_dislotwin_burgersPerSlipSystem(j,instance)*&
          plastic_dislotwin_v0PerSlipSystem(j,instance)
 
        !* Shear rates due to slip
        gdot_slip(j) = DotGamma0 &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_dislotwin_qPerSlipFamily(f,instance)) &
                     * sign(1.0_pReal,tau_slip(j))
 
        !* Derivatives of shear rates
        dgdot_dtauslip(j) = &
          abs(gdot_slip(j))*BoltzmannRatio*plastic_dislotwin_pPerSlipFamily(f,instance)&
          *plastic_dislotwin_qPerSlipFamily(f,instance)/&
          (plastic_dislotwin_SolidSolutionStrength(instance)+plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance))*&
          StressRatio_pminus1*(1-StressRatio_p)**(plastic_dislotwin_qPerSlipFamily(f,instance)-1.0_pReal)
      endif
 
      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + gdot_slip(j)*lattice_Sslip(:,:,1,index_myFamily+i,ph)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*&
                                  lattice_Sslip(k,l,1,index_myFamily+i,ph)*&
                                  lattice_Sslip(m,n,1,index_myFamily+i,ph)
   enddo slipSystemsLoop
 enddo slipFamiliesLoop
 
!--------------------------------------------------------------------------------------------------
! correct Lp and dLp_dTstar3333 for twinned and transformed fraction 
 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:nt,of)) ! safe for nt == 0

 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:nr,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:nr,of))
 Lp = Lp * (1.0_pReal - sumf - sumftr)
 dLp_dTstar3333 = dLp_dTstar3333 * (1.0_pReal - sumf - sumftr)
 
!--------------------------------------------------------------------------------------------------
! Shear banding (shearband) part
 if(abs(plastic_dislotwin_sbVelocity(instance))   > tiny(0.0_pReal) .and. &
    abs(plastic_dislotwin_sbResistance(instance)) > tiny(0.0_pReal)) then
   gdot_sb = 0.0_pReal
   dgdot_dtausb = 0.0_pReal
   call math_eigenValuesVectorsSym(math_Mandel6to33(Tstar_v),eigValues,eigVectors,error)
   do j = 1_pInt,6_pInt
     sb_s = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_sComposition(1:3,j))
     sb_m = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_mComposition(1:3,j))
     sb_Smatrix = math_tensorproduct33(sb_s,sb_m)
     plastic_dislotwin_sbSv(1:6,j,ipc,ip,el) = math_Mandel33to6(math_symmetric33(sb_Smatrix))
   
     !* Calculation of Lp
     !* Resolved shear stress on shear banding system
     tau_sb(j) = dot_product(Tstar_v,plastic_dislotwin_sbSv(1:6,j,ipc,ip,el))
   
     !* Stress ratios
      if (abs(tau_sb(j)) < tol_math_check) then
        StressRatio_p = 0.0_pReal
        StressRatio_pminus1 = 0.0_pReal
      else
       StressRatio_p = (abs(tau_sb(j))/plastic_dislotwin_sbResistance(instance))&
                 **plastic_dislotwin_pShearBand(instance)
       StressRatio_pminus1 = (abs(tau_sb(j))/plastic_dislotwin_sbResistance(instance))&
                 **(plastic_dislotwin_pShearBand(instance)-1.0_pReal)
      endif

     !* Boltzmann ratio
     BoltzmannRatio = plastic_dislotwin_sbQedge(instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = plastic_dislotwin_sbVelocity(instance)
 
     !* Shear rates due to shearband
     gdot_sb(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                  plastic_dislotwin_qShearBand(instance))*sign(1.0_pReal,tau_sb(j))
                  
     !* Derivatives of shear rates
     dgdot_dtausb(j) = &
       ((abs(gdot_sb(j))*BoltzmannRatio*&
       plastic_dislotwin_pShearBand(instance)*plastic_dislotwin_qShearBand(instance))/&
       plastic_dislotwin_sbResistance(instance))*&
       StressRatio_pminus1*(1_pInt-StressRatio_p)**(plastic_dislotwin_qShearBand(instance)-1.0_pReal)
 
     !* Plastic velocity gradient for shear banding
     Lp = Lp + gdot_sb(j)*sb_Smatrix
 
     !* Calculation of the tangent of Lp
     forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
       dLp_dTstar3333(k,l,m,n) = &
       dLp_dTstar3333(k,l,m,n) + dgdot_dtausb(j)*&
                                 sb_Smatrix(k,l)*&
                                 sb_Smatrix(m,n)
   enddo
 end if
 
!--------------------------------------------------------------------------------------------------
! Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
 twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph)) ! at which index starts my family
   twinSystemsLoop: do i = 1_pInt,plastic_dislotwin_Ntwin(f,instance)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_twin(j) > tol_math_check) then
        StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/tau_twin(j))**plastic_dislotwin_rPerTwinFamily(f,instance)
        !* Shear rates and their derivatives due to twin
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < plastic_dislotwin_tau_r_twin(j,instance)) then
              Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(ph)%rhoEdgeDip(s2,of))+&  !!!!! correct?
                     abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                    (plastic_dislotwin_L0_twin(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-plastic_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                    (plastic_dislotwin_tau_r_twin(j,instance)-tau_twin(j))))
            else
              Ndot0_twin=0.0_pReal
            end if
          case default
            Ndot0_twin=plastic_dislotwin_Ndot0PerTwinSystem(j,instance)
        end select
        gdot_twin(j) = &
          (1.0_pReal-sumf-sumftr)*lattice_shearTwin(index_myFamily+i,ph)*&
          state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
        dgdot_dtautwin(j) = ((gdot_twin(j)*plastic_dislotwin_rPerTwinFamily(f,instance))/tau_twin(j))*StressRatio_r
      endif
 
      !* Plastic velocity gradient for mechanical twinning
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,ph)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,ph)*&
                                  lattice_Stwin(m,n,index_myFamily+i,ph)
   enddo twinSystemsLoop
 enddo twinFamiliesLoop

 !* Phase transformation part
 gdot_trans = 0.0_pReal
 dgdot_dtautrans = 0.0_pReal
 j = 0_pInt
 transFamiliesLoop: do f = 1_pInt,lattice_maxNtransFamily
   index_myFamily = sum(lattice_NtransSystem(1:f-1_pInt,ph)) ! at which index starts my family
   transSystemsLoop: do i = 1_pInt,plastic_dislotwin_Ntrans(f,instance)
      j = j+1_pInt

      !* Resolved shear stress on transformation system
      tau_trans(j) = dot_product(Tstar_v,lattice_Strans_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_trans(j) > tol_math_check) then
        StressRatio_s = (state(instance)%threshold_stress_trans(j,of)/tau_trans(j))**plastic_dislotwin_sPerTransFamily(f,instance)
        !* Shear rates and their derivatives due to transformation
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_trans(j) < plastic_dislotwin_tau_r_trans(j,instance)) then
              Ndot0_trans=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&  !!!!! correct?
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (plastic_dislotwin_L0_trans(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
                          (1.0_pReal-exp(-plastic_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                          (plastic_dislotwin_tau_r_trans(j,instance)-tau_trans(j))))
            else
              Ndot0_trans=0.0_pReal
            end if
          case default
            Ndot0_trans=plastic_dislotwin_Ndot0PerTransSystem(j,instance)
        end select
        gdot_trans(j) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%martensiteVolume(j,of)*Ndot0_trans*exp(-StressRatio_s)
        dgdot_dtautrans(j) = ((gdot_trans(j)*plastic_dislotwin_sPerTransFamily(f,instance))/tau_trans(j))*StressRatio_s
      endif

      !* Plastic velocity gradient for phase transformation
      Lp = Lp + gdot_trans(j)*lattice_Strans(:,:,index_myFamily+i,ph)

      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautrans(j)*&
                                  lattice_Strans(k,l,index_myFamily+i,ph)*&
                                  lattice_Strans(m,n,index_myFamily+i,ph)

   enddo transSystemsLoop
 enddo transFamiliesLoop

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)
 
end subroutine plastic_dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_dislotwin_dotState(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   phaseAt, phasememberAt
 use lattice,  only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_Strans_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_maxNtransFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NtransSystem, &   
   lattice_sheartwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   lattice_fccTobcc_transNucleationTwinPair, &
   lattice_fccTobcc_shearCritTrans, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: instance,ns,nt,nr,f,i,j,index_myFamily,s1,s2, &
                  ph, &
                  of
 real(pReal) :: sumf,sumftr,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,&
             EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r,Ndot0_twin,stressRatio,&
             Ndot0_trans,StressRatio_s,EdgeDipDistance, ClimbVelocity,DotRhoEdgeDipClimb,DotRhoEdgeDipAnnihilation, &
             DotRhoDipFormation,DotRhoMultiplication,DotRhoEdgeEdgeAnnihilation
 real(pReal), dimension(plastic_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
 gdot_slip,tau_slip

 real(pReal), dimension(plastic_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
              tau_twin
 real(pReal), dimension(plastic_dislotwin_totalNtrans(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
              tau_trans

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:nt,of)) ! safe for nt == 0
 plasticState(ph)%dotState(:,of) = 0.0_pReal
 
 !* Total transformed volume fraction
 sumftr = sum(state(instance)%stressTransFraction(1_pInt:nr,of)) + &
          sum(state(instance)%strainTransFraction(1_pInt:nr,of))
 
 !* Dislocation density evolution
 gdot_slip = 0.0_pReal
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                         ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))   ! at which index starts my family
   do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
      j = j+1_pInt
 
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))

      if((abs(tau_slip(j))-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
      !* Stress ratios
        stressRatio =((abs(tau_slip(j))- state(instance)%threshold_stress_slip(j,of))/&
          (plastic_dislotwin_SolidSolutionStrength(instance)+plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))
        StressRatio_p       = stressRatio** plastic_dislotwin_pPerSlipFamily(f,instance)
        StressRatio_pminus1 = stressRatio**(plastic_dislotwin_pPerSlipFamily(f,instance)-1.0_pReal)
      !* Boltzmann ratio
        BoltzmannRatio = plastic_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
      !* Initial shear rates
        DotGamma0 = &
          plasticState(ph)%state(j, of)*plastic_dislotwin_burgersPerSlipSystem(j,instance)*&
          plastic_dislotwin_v0PerSlipSystem(j,instance)
 
      !* Shear rates due to slip
        gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)** &
                     plastic_dislotwin_qPerSlipFamily(f,instance))*sign(1.0_pReal,tau_slip(j))
      endif
      !* Multiplication
      DotRhoMultiplication = abs(gdot_slip(j))/&
                                (plastic_dislotwin_burgersPerSlipSystem(j,instance)*state(instance)%mfp_slip(j,of))
      !* Dipole formation
      EdgeDipMinDistance = &
        plastic_dislotwin_CEdgeDipMinDistance(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance)
      if (abs(tau_slip(j)) <= tiny(0.0_pReal)) then
        DotRhoDipFormation = 0.0_pReal
      else
        EdgeDipDistance = &
          (3.0_pReal*lattice_mu(ph)*plastic_dislotwin_burgersPerSlipSystem(j,instance))/&
          (16.0_pReal*pi*abs(tau_slip(j)))
        if (EdgeDipDistance>state(instance)%mfp_slip(j,of)) EdgeDipDistance=state(instance)%mfp_slip(j,of)
        if (EdgeDipDistance<EdgeDipMinDistance)       EdgeDipDistance=EdgeDipMinDistance
        DotRhoDipFormation = &
          ((2.0_pReal*(EdgeDipDistance-EdgeDipMinDistance))/plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
          state(instance)%rhoEdge(j,of)*abs(gdot_slip(j))*plastic_dislotwin_dipoleFormationFactor(instance)
      endif
 
      !* Spontaneous annihilation of 2 single edge dislocations
      DotRhoEdgeEdgeAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
        state(instance)%rhoEdge(j,of)*abs(gdot_slip(j))
 
      !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
      DotRhoEdgeDipAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
        state(instance)%rhoEdgeDip(j,of)*abs(gdot_slip(j))
 
      !* Dislocation dipole climb
      AtomicVolume = &
        plastic_dislotwin_CAtomicVolume(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance)**(3.0_pReal)
      VacancyDiffusion = &
        plastic_dislotwin_D0(instance)*exp(-plastic_dislotwin_Qsd(instance)/(kB*Temperature))
      if (abs(tau_slip(j)) <= tiny(0.0_pReal)) then
        DotRhoEdgeDipClimb = 0.0_pReal
      else
        if (EdgeDipDistance-EdgeDipMinDistance <= tiny(0.0_pReal)) then
          DotRhoEdgeDipClimb = 0.0_pReal
        else
          ClimbVelocity = 3.0_pReal*lattice_mu(ph)*VacancyDiffusion*AtomicVolume/ &
                          (2.0_pReal*pi*kB*Temperature*(EdgeDipDistance+EdgeDipMinDistance))
          DotRhoEdgeDipClimb = 4.0_pReal*ClimbVelocity*state(instance)%rhoEdgeDip(j,of)/ &
                          (EdgeDipDistance-EdgeDipMinDistance)
        endif
      endif
 
      !* Edge dislocation density rate of change
      dotState(instance)%rhoEdge(j,of) = &
        DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
 
      !* Edge dislocation dipole density rate of change
      dotState(instance)%rhoEdgeDip(j,of) = &
        DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
 
      !* Dotstate for accumulated shear due to slip
      dotState(instance)%accshear_slip(j,of) = abs(gdot_slip(j))
 
   enddo
 enddo
 
 !* Twin volume fraction evolution
 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                         ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))   ! at which index starts my family
   do i = 1_pInt,plastic_dislotwin_Ntwin(f,instance)          ! process each (active) twin system in family
      j = j+1_pInt
 
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
      !* Stress ratios
      if (tau_twin(j) > tol_math_check) then
        StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/&
                         tau_twin(j))**plastic_dislotwin_rPerTwinFamily(f,instance)
        !* Shear rates and their derivatives due to twin
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < plastic_dislotwin_tau_r_twin(j,instance)) then
              Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                     abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                    (plastic_dislotwin_L0_twin(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-plastic_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                    (plastic_dislotwin_tau_r_twin(j,instance)-tau_twin(j))))
            else
              Ndot0_twin=0.0_pReal
            end if
          case default
            Ndot0_twin=plastic_dislotwin_Ndot0PerTwinSystem(j,instance)
        end select
        dotState(instance)%twinFraction(j,of) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
        !* Dotstate for accumulated shear due to twin
        dotState(instance)%accshear_twin(j,of) = dotState(instance)%twinFraction(j,of) * &
                                                          lattice_sheartwin(index_myfamily+i,ph)
      endif
   enddo
 enddo

 !* Transformation volume fraction evolution
 j = 0_pInt
 do f = 1_pInt,lattice_maxNtransFamily                        ! loop over all trans families
   index_myFamily = sum(lattice_NtransSystem(1:f-1_pInt,ph))  ! at which index starts my family
   do i = 1_pInt,plastic_dislotwin_Ntrans(f,instance)         ! process each (active) trans system in family
      j = j+1_pInt

      !* Resolved shear stress on transformation system
      tau_trans(j) = dot_product(Tstar_v,lattice_Strans_v(:,index_myFamily+i,ph))

      !* Stress ratios
      if (tau_trans(j) > tol_math_check) then
        StressRatio_s = (state(instance)%threshold_stress_trans(j,of)/&
                         tau_trans(j))**plastic_dislotwin_sPerTransFamily(f,instance)
        !* Shear rates and their derivatives due to transformation
        select case(lattice_structure(ph))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_trans(j) < plastic_dislotwin_tau_r_trans(j,instance)) then
              Ndot0_trans=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (plastic_dislotwin_L0_trans(instance)*plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
                          (1.0_pReal-exp(-plastic_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                          (plastic_dislotwin_tau_r_trans(j,instance)-tau_trans(j))))
            else
              Ndot0_trans=0.0_pReal
            end if
          case default
            Ndot0_trans=plastic_dislotwin_Ndot0PerTransSystem(j,instance)
        end select
        dotState(instance)%strainTransFraction(j,of) = &
          (1.0_pReal-sumf-sumftr)*&
          state(instance)%martensiteVolume(j,of)*Ndot0_trans*exp(-StressRatio_s)
        !* Dotstate for accumulated shear due to transformation
        !dotState(instance)%accshear_trans(j,of) = dotState(instance)%strainTransFraction(j,of) * &
        !                                                  lattice_sheartrans(index_myfamily+i,ph)
      endif

   enddo
 enddo

end subroutine plastic_dislotwin_dotState

 
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_dislotwin_postResults(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   pi, &
   math_Mandel6to33, &
   math_eigenValuesSym33, &
   math_eigenValuesVectorsSym33
 use material, only: &
   material_phase, &
   phase_plasticityInstance,& 
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(plastic_dislotwin_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           plastic_dislotwin_postResults
 integer(pInt) :: &
   instance,&
   ns,nt,nr,&
   f,o,i,c,j,index_myFamily,&
   s1,s2, &
   ph, &
   of
 real(pReal) :: sumf,tau,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r,Ndot0_twin,dgdot_dtauslip, &
   stressRatio
 real(preal), dimension(plastic_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip
 real(pReal), dimension(3,3) :: eigVectors
 real(pReal), dimension (3) :: eigValues
 
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_dislotwin_totalNslip(instance)
 nt = plastic_dislotwin_totalNtwin(instance)
 nr = plastic_dislotwin_totalNtrans(instance)

 !* Total twin volume fraction
 sumf = sum(state(instance)%twinFraction(1_pInt:nt,of))                          ! safe for nt == 0
 
 !* Required output
 c = 0_pInt
 plastic_dislotwin_postResults = 0.0_pReal
 do o = 1_pInt,plastic_dislotwin_Noutput(instance)
    select case(plastic_dislotwin_outputID(o,instance))
 
      case (edge_density_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+ns) = state(instance)%rhoEdge(1_pInt:ns,of)
        c = c + ns
      case (dipole_density_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+ns) = state(instance)%rhoEdgeDip(1_pInt:ns,of)
        c = c + ns
      case (shear_rate_slip_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)                                        ! process each (active) slip system in family
              j = j + 1_pInt                                                                        ! could be taken from state by now!
 
              !* Resolved shear stress on slip system
              tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
              !* Stress ratios
              if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
              !* Stress ratios
                stressRatio = ((abs(tau)-state(ph)%threshold_stress_slip(j,of))/&
                                (plastic_dislotwin_SolidSolutionStrength(instance)+&
                                 plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))
                StressRatio_p       = stressRatio** plastic_dislotwin_pPerSlipFamily(f,instance)
                StressRatio_pminus1 = stressRatio**(plastic_dislotwin_pPerSlipFamily(f,instance)-1.0_pReal)
              !* Boltzmann ratio
                BoltzmannRatio = plastic_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
              !* Initial shear rates
                DotGamma0 = &
                  state(instance)%rhoEdge(j,of)*plastic_dislotwin_burgersPerSlipSystem(j,instance)* &
                  plastic_dislotwin_v0PerSlipSystem(j,instance)
 
              !* Shear rates due to slip
                plastic_dislotwin_postResults(c+j) = &
                  DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                               plastic_dislotwin_qPerSlipFamily(f,instance))*sign(1.0_pReal,tau)
              else
                plastic_dislotwin_postResults(c+j) = 0.0_pReal
              endif
              
        enddo ; enddo
        c = c + ns
      case (accumulated_shear_slip_ID)
       plastic_dislotwin_postResults(c+1_pInt:c+ns) = &
                      state(instance)%accshear_slip(1_pInt:ns,of)
        c = c + ns
      case (mfp_slip_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+ns) =&
                      state(instance)%mfp_slip(1_pInt:ns,of)
        c = c + ns
      case (resolved_stress_slip_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)                                        ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) =&
                                dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
        enddo; enddo
        c = c + ns
      case (threshold_stress_slip_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+ns) = &
                                state(instance)%threshold_stress_slip(1_pInt:ns,of)
        c = c + ns
      case (edge_dipole_distance_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                                        ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)                                        ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) = &
                (3.0_pReal*lattice_mu(ph)*plastic_dislotwin_burgersPerSlipSystem(j,instance))/&
                (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))))
              plastic_dislotwin_postResults(c+j)=min(plastic_dislotwin_postResults(c+j),&
                                                            state(instance)%mfp_slip(j,of))
 !            plastic_dislotwin_postResults(c+j)=max(plastic_dislotwin_postResults(c+j),&
 !                                                            plasticState(ph)%state(4*ns+2*nt+2*nr+j, of))
        enddo; enddo
        c = c + ns
       case (resolved_stress_shearband_ID)
         do j = 1_pInt,6_pInt                                                                       ! loop over all shearband families
            plastic_dislotwin_postResults(c+j) = dot_product(Tstar_v, &
                                                       plastic_dislotwin_sbSv(1:6,j,ipc,ip,el))
         enddo
         c = c + 6_pInt
       case (shear_rate_shearband_ID)
         do j = 1_pInt,6_pInt                                                                       ! loop over all shearbands
              !* Resolved shear stress on shearband system
              tau = dot_product(Tstar_v,plastic_dislotwin_sbSv(1:6,j,ipc,ip,el))
              !* Stress ratios
              if (abs(tau) < tol_math_check) then
                StressRatio_p = 0.0_pReal
                StressRatio_pminus1 = 0.0_pReal
              else
                StressRatio_p = (abs(tau)/plastic_dislotwin_sbResistance(instance))**&
                                 plastic_dislotwin_pShearBand(instance)
                StressRatio_pminus1 = (abs(tau)/plastic_dislotwin_sbResistance(instance))**&
                                (plastic_dislotwin_pShearBand(instance)-1.0_pReal)
              endif
              !* Boltzmann ratio
              BoltzmannRatio = plastic_dislotwin_sbQedge(instance)/(kB*Temperature)
              !* Initial shear rates
              DotGamma0 = plastic_dislotwin_sbVelocity(instance)
              ! Shear rate due to shear band
              plastic_dislotwin_postResults(c+j) = &
          DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**plastic_dislotwin_qShearBand(instance))*&
                sign(1.0_pReal,tau)
         enddo 
        c = c + 6_pInt
      case (twin_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nt) = state(instance)%twinFraction(1_pInt:nt,of)
        c = c + nt
      case (shear_rate_twin_ID)
        if (nt > 0_pInt) then
        
          j = 0_pInt
          do f = 1_pInt,lattice_maxNslipFamily                                                      ! loop over all slip families
             index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                               ! at which index starts my family
             do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)                                      ! process each (active) slip system in family
                j = j + 1_pInt
 
               !* Resolved shear stress on slip system
               tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
               !* Stress ratios
               if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
               !* Stress ratios
                 StressRatio_p = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                 (plastic_dislotwin_SolidSolutionStrength(instance)+&
                                  plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))&
                                                    **plastic_dislotwin_pPerSlipFamily(f,instance)
                 StressRatio_pminus1 = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                       (plastic_dislotwin_SolidSolutionStrength(instance)+&
                                        plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))&
                                                    **(plastic_dislotwin_pPerSlipFamily(f,instance)-1.0_pReal)
               !* Boltzmann ratio
                 BoltzmannRatio = plastic_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
               !* Initial shear rates
                 DotGamma0 = &
                   state(instance)%rhoEdge(j,of)*plastic_dislotwin_burgersPerSlipSystem(j,instance)* &
                   plastic_dislotwin_v0PerSlipSystem(j,instance)
 
               !* Shear rates due to slip
                 gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                                plastic_dislotwin_qPerSlipFamily(f,instance))*sign(1.0_pReal,tau)
               else
                 gdot_slip(j) = 0.0_pReal
               endif
          enddo;enddo

          j = 0_pInt
          do f = 1_pInt,lattice_maxNtwinFamily                                                      ! loop over all twin families
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
            do i = 1,plastic_dislotwin_Ntwin(f,instance)                                            ! process each (active) twin system in family
              j = j + 1_pInt

              tau = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))


              !* Shear rates due to twin
              if ( tau > 0.0_pReal ) then
                select case(lattice_structure(ph))
                  case (LATTICE_fcc_ID)
                  s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
                  s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
                  if (tau < plastic_dislotwin_tau_r_twin(j,instance)) then
                    Ndot0_twin=(abs(gdot_slip(s1))*(state(instance)%rhoEdge(s2,of)+state(instance)%rhoEdgeDip(s2,of))+&
                           abs(gdot_slip(s2))*(state(instance)%rhoEdge(s1,of)+state(instance)%rhoEdgeDip(s1,of)))/&
                          (plastic_dislotwin_L0_twin(instance)*&
                           plastic_dislotwin_burgersPerSlipSystem(j,instance))*&
                          (1.0_pReal-exp(-plastic_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                          (plastic_dislotwin_tau_r_twin(j,instance)-tau)))
                  else
                    Ndot0_twin=0.0_pReal
                  end if
                  case default
                    Ndot0_twin=plastic_dislotwin_Ndot0PerTwinSystem(j,instance)
                end select
                StressRatio_r = (state(instance)%threshold_stress_twin(j,of)/tau) &
                                                   **plastic_dislotwin_rPerTwinFamily(f,instance)
                plastic_dislotwin_postResults(c+j) = &
                  (plastic_dislotwin_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,ph)*&
                  state(instance)%twinVolume(j,of)*Ndot0_twin*exp(-StressRatio_r)
              endif
 
          enddo ; enddo
        endif
        c = c + nt
      case (accumulated_shear_twin_ID)
       plastic_dislotwin_postResults(c+1_pInt:c+nt) = state(instance)%accshear_twin(1_pInt:nt,of)
        c = c + nt     
      case (mfp_twin_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nt) = state(instance)%mfp_twin(1_pInt:nt,of)
        c = c + nt
      case (resolved_stress_twin_ID)
        if (nt > 0_pInt) then
          j = 0_pInt
          do f = 1_pInt,lattice_maxNtwinFamily                                                      ! loop over all slip families
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
            do i = 1_pInt,plastic_dislotwin_Ntwin(f,instance)                                  ! process each (active) slip system in family
              j = j + 1_pInt
              plastic_dislotwin_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,ph))
          enddo; enddo
        endif
        c = c + nt
      case (threshold_stress_twin_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nt) = state(instance)%threshold_stress_twin(1_pInt:nt,of)
        c = c + nt
      case (stress_exponent_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                                        ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                  ! at which index starts my family
          do i = 1_pInt,plastic_dislotwin_Nslip(f,instance)                                    ! process each (active) slip system in family
             j = j + 1_pInt

             !* Resolved shear stress on slip system
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
             if((abs(tau)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
             !* Stress ratios
               StressRatio_p = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                               (plastic_dislotwin_SolidSolutionStrength(instance)+&
                                plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))&
                                                  **plastic_dislotwin_pPerSlipFamily(f,instance)
               StressRatio_pminus1 = ((abs(tau)-state(instance)%threshold_stress_slip(j,of))/&
                                     (plastic_dislotwin_SolidSolutionStrength(instance)+&
                                      plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance)))&
                                                  **(plastic_dislotwin_pPerSlipFamily(f,instance)-1.0_pReal)
             !* Boltzmann ratio
               BoltzmannRatio = plastic_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
             !* Initial shear rates
               DotGamma0 = &
                 state(instance)%rhoEdge(j,of)*plastic_dislotwin_burgersPerSlipSystem(j,instance)* &
                 plastic_dislotwin_v0PerSlipSystem(j,instance)

             !* Shear rates due to slip
               gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                            plastic_dislotwin_qPerSlipFamily(f,instance))*sign(1.0_pReal,tau)

             !* Derivatives of shear rates
               dgdot_dtauslip = &
                 abs(gdot_slip(j))*BoltzmannRatio*plastic_dislotwin_pPerSlipFamily(f,instance)&
                 *plastic_dislotwin_qPerSlipFamily(f,instance)/&
                 (plastic_dislotwin_SolidSolutionStrength(instance)+&
                  plastic_dislotwin_tau_peierlsPerSlipFamily(f,instance))*&
                 StressRatio_pminus1*(1-StressRatio_p)**(plastic_dislotwin_qPerSlipFamily(f,instance)-1.0_pReal)
             
             else
               gdot_slip(j) = 0.0_pReal
               dgdot_dtauslip = 0.0_pReal
             endif

             !* Stress exponent
             if (abs(gdot_slip(j))<=tiny(0.0_pReal)) then
               plastic_dislotwin_postResults(c+j) = 0.0_pReal
             else
               plastic_dislotwin_postResults(c+j) = (tau/gdot_slip(j))*dgdot_dtauslip
             endif
         enddo ; enddo
         c = c + ns
      case (sb_eigenvalues_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+3_pInt) = math_eigenvaluesSym33(math_Mandel6to33(Tstar_v))
        c = c + 3_pInt
      case (sb_eigenvectors_ID)
        call math_eigenValuesVectorsSym33(math_Mandel6to33(Tstar_v),eigValues,eigVectors)
        plastic_dislotwin_postResults(c+1_pInt:c+9_pInt) = reshape(eigVectors,[9])
        c = c + 9_pInt
      case (stress_trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nr) = &
          state(instance)%stressTransFraction(1_pInt:nr,of)
        c = c + nr
      case (strain_trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nr) = &
          state(instance)%strainTransFraction(1_pInt:nr,of)
        c = c + nr
      case (trans_fraction_ID)
        plastic_dislotwin_postResults(c+1_pInt:c+nr) = &
          state(instance)%stressTransFraction(1_pInt:nr,of) + &
          state(instance)%strainTransFraction(1_pInt:nr,of)
        c = c + nr
    end select
 enddo
end function plastic_dislotwin_postResults

end module plastic_dislotwin
