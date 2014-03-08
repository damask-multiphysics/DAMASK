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
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module constitutive_dislotwin
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   constitutive_dislotwin_sizeDotState, &                                                           !< number of dotStates
   constitutive_dislotwin_sizeState, &                                                              !< total number of microstructural state variables
   constitutive_dislotwin_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   constitutive_dislotwin_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   constitutive_dislotwin_output                                                                    !< name of each post result output
   
 character(len=12),                   dimension(3),           parameter,           private :: &
   CONSTITUTIVE_DISLOTWIN_listBasicSlipStates = &
   ['rhoEdge     ',      'rhoEdgeDip  ',      'accshearslip']

 character(len=12),                   dimension(2),           parameter,           private :: &
   CONSTITUTIVE_DISLOTWIN_listBasicTwinStates = & 
   ['twinFraction',      'accsheartwin']

 character(len=17),                   dimension(4),           parameter,           private :: &
   CONSTITUTIVE_DISLOTWIN_listDependentSlipStates = &
   ['invLambdaSlip    ', 'invLambdaSlipTwin', 'meanFreePathSlip ', 'tauSlipThreshold ']

 character(len=16),                   dimension(4),           parameter,           private :: &
   CONSTITUTIVE_DISLOTWIN_listDependentTwinStates = & 
   ['invLambdaTwin   ',  'meanFreePathTwin',  'tauTwinThreshold',  'twinVolume      ']

 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 integer(pInt),                       dimension(:),           allocatable,         private :: &
   constitutive_dislotwin_Noutput                                                                   !< number of outputs per instance of this plasticity 

 integer(pInt),                       dimension(:),           allocatable,         private :: &
   constitutive_dislotwin_totalNslip, &                                                             !< total number of active slip systems for each instance
   constitutive_dislotwin_totalNtwin                                                                !< total number of active twin systems for each instance

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   constitutive_dislotwin_Nslip, &                                                                  !< number of active slip systems for each family and instance
   constitutive_dislotwin_Ntwin                                                                     !< number of active twin systems for each family and instance

 real(pReal),                         dimension(:),           allocatable,         private :: &
   constitutive_dislotwin_CAtomicVolume, &                                                          !< atomic volume in Bugers vector unit
   constitutive_dislotwin_D0, &                                                                     !< prefactor for self-diffusion coefficient
   constitutive_dislotwin_Qsd, &                                                                    !< activation energy for dislocation climb
   constitutive_dislotwin_GrainSize, &                                                              !< grain size
   constitutive_dislotwin_p, &                                                                      !< p-exponent in glide velocity
   constitutive_dislotwin_q, &                                                                      !< q-exponent in glide velocity
   constitutive_dislotwin_MaxTwinFraction, &                                                        !< maximum allowed total twin volume fraction
   constitutive_dislotwin_r, &                                                                      !< r-exponent in twin nucleation rate
   constitutive_dislotwin_CEdgeDipMinDistance, &                                                    !<
   constitutive_dislotwin_Cmfptwin, &                                                               !<
   constitutive_dislotwin_Cthresholdtwin, &                                                         !<
   constitutive_dislotwin_SolidSolutionStrength, &                                                  !< Strength due to elements in solid solution
   constitutive_dislotwin_L0, &                                                                     !< Length of twin nuclei in Burgers vectors
   constitutive_dislotwin_xc, &                                                                     !< critical distance for formation of twin nucleus
   constitutive_dislotwin_VcrossSlip, &                                                             !< cross slip volume
   constitutive_dislotwin_sbResistance, &                                                           !< value for shearband resistance (might become an internal state variable at some point)
   constitutive_dislotwin_sbVelocity, &                                                             !< value for shearband velocity_0
   constitutive_dislotwin_sbQedge, &                                                                !< value for shearband systems Qedge
   constitutive_dislotwin_SFE_0K, &                                                                 !< stacking fault energy at zero K
   constitutive_dislotwin_dSFE_dT, &                                                                !< temperature dependance of stacking fault energy
   constitutive_dislotwin_aTolRho, &                                                                !< absolute tolerance for integration of dislocation density
   constitutive_dislotwin_aTolTwinFrac                                                              !< absolute tolerance for integration of twin volume fraction

 real(pReal),                         dimension(:,:,:,:),     allocatable,         private :: &
   constitutive_dislotwin_Ctwin_66                                                                  !< twin elasticity matrix in Mandel notation for each instance

 real(pReal),                         dimension(:,:,:,:,:,:), allocatable,         private :: &
   constitutive_dislotwin_Ctwin_3333                                                                !< twin elasticity matrix for each instance

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   constitutive_dislotwin_rhoEdge0, &                                                               !< initial edge dislocation density per slip system for each family and instance
   constitutive_dislotwin_rhoEdgeDip0, &                                                            !< initial edge dipole density per slip system for each family and instance
   constitutive_dislotwin_burgersPerSlipFamily, &                                                   !< absolute length of burgers vector [m] for each slip family and instance
   constitutive_dislotwin_burgersPerSlipSystem, &                                                   !< absolute length of burgers vector [m] for each slip system and instance
   constitutive_dislotwin_burgersPerTwinFamily, &                                                   !< absolute length of burgers vector [m] for each twin family and instance
   constitutive_dislotwin_burgersPerTwinSystem, &                                                   !< absolute length of burgers vector [m] for each twin system and instance
   constitutive_dislotwin_QedgePerSlipFamily, &                                                     !< activation energy for glide [J] for each slip family and instance
   constitutive_dislotwin_QedgePerSlipSystem, &                                                     !< activation energy for glide [J] for each slip system and instance
   constitutive_dislotwin_v0PerSlipFamily, &                                                        !< dislocation velocity prefactor [m/s] for each family and instance
   constitutive_dislotwin_v0PerSlipSystem, &                                                        !< dislocation velocity prefactor [m/s] for each slip system and instance
   constitutive_dislotwin_Ndot0PerTwinFamily, &                                                     !< twin nucleation rate [1/m³s] for each twin family and instance
   constitutive_dislotwin_Ndot0PerTwinSystem, &                                                     !< twin nucleation rate [1/m³s] for each twin system and instance
   constitutive_dislotwin_tau_r, &                                                                  !< stress to bring partial close together for each twin system and instance
   constitutive_dislotwin_twinsizePerTwinFamily, &                                                  !< twin thickness [m] for each twin family and instance
   constitutive_dislotwin_twinsizePerTwinSystem, &                                                  !< twin thickness [m] for each twin system and instance
   constitutive_dislotwin_CLambdaSlipPerSlipFamily, &                                               !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   constitutive_dislotwin_CLambdaSlipPerSlipSystem, &                                               !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   constitutive_dislotwin_interaction_SlipSlip, &                                                   !< coefficients for slip-slip interaction for each interaction type and instance
   constitutive_dislotwin_interaction_SlipTwin, &                                                   !< coefficients for slip-twin interaction for each interaction type and instance
   constitutive_dislotwin_interaction_TwinSlip, &                                                   !< coefficients for twin-slip interaction for each interaction type and instance
   constitutive_dislotwin_interaction_TwinTwin                                                      !< coefficients for twin-twin interaction for each interaction type and instance
 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   constitutive_dislotwin_interactionMatrix_SlipSlip, &                                             !< interaction matrix of the different slip systems for each instance
   constitutive_dislotwin_interactionMatrix_SlipTwin, &                                             !< interaction matrix of slip systems with twin systems for each instance
   constitutive_dislotwin_interactionMatrix_TwinSlip, &                                             !< interaction matrix of twin systems with slip systems for each instance
   constitutive_dislotwin_interactionMatrix_TwinTwin, &                                             !< interaction matrix of the different twin systems for each instance
   constitutive_dislotwin_forestProjectionEdge                                                      !< matrix of forest projections of edge dislocations for each instance
 real(pReal),                         dimension(:,:,:,:,:),   allocatable,         private :: &
   constitutive_dislotwin_sbSv

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
                 sb_eigenvectors_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   constitutive_dislotwin_outputID                                                                  !< ID of each post result output


 public :: &
   constitutive_dislotwin_init, &
   constitutive_dislotwin_stateInit, &
   constitutive_dislotwin_aTolState, &
   constitutive_dislotwin_homogenizedC, &
   constitutive_dislotwin_microstructure, &
   constitutive_dislotwin_LpAndItsTangent, &
   constitutive_dislotwin_dotState, &
   constitutive_dislotwin_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislotwin_init(fileUnit)
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
   MATERIAL_partPhase
 use lattice
 
 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = LATTICE_maxNinteraction + 1_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,maxTotalNslip,maxTotalNtwin,&
                  f,instance,j,k,l,m,n,o,p,q,r,s,ns,nt, &
                  Nchunks_SlipSlip, Nchunks_SlipTwin, Nchunks_TwinSlip, Nchunks_TwinTwin, &
                  Nchunks_SlipFamilies, Nchunks_TwinFamilies, &
                  index_myFamily, index_otherFamily
 character(len=65536) :: &
   tag  = '', &
   line = ''
  
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOTWIN_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOTWIN_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(constitutive_dislotwin_sizeDotState(maxNinstance),                        source=0_pInt)
 allocate(constitutive_dislotwin_sizeState(maxNinstance),                           source=0_pInt)
 allocate(constitutive_dislotwin_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(constitutive_dislotwin_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(constitutive_dislotwin_output(maxval(phase_Noutput),maxNinstance))
          constitutive_dislotwin_output = ''
 allocate(constitutive_dislotwin_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(constitutive_dislotwin_Noutput(maxNinstance),                             source=0_pInt)
 allocate(constitutive_dislotwin_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(constitutive_dislotwin_Ntwin(lattice_maxNtwinFamily,maxNinstance),        source=0_pInt)
 allocate(constitutive_dislotwin_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(constitutive_dislotwin_totalNtwin(maxNinstance),                          source=0_pInt)
 allocate(constitutive_dislotwin_CAtomicVolume(maxNinstance),                       source=0.0_pReal)
 allocate(constitutive_dislotwin_D0(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislotwin_Qsd(maxNinstance),                                 source=0.0_pReal)
 allocate(constitutive_dislotwin_GrainSize(maxNinstance),                           source=0.0_pReal)
 allocate(constitutive_dislotwin_p(maxNinstance),                                   source=0.0_pReal)
 allocate(constitutive_dislotwin_q(maxNinstance),                                   source=0.0_pReal)
 allocate(constitutive_dislotwin_MaxTwinFraction(maxNinstance),                     source=0.0_pReal)
 allocate(constitutive_dislotwin_r(maxNinstance),                                   source=0.0_pReal)
 allocate(constitutive_dislotwin_CEdgeDipMinDistance(maxNinstance),                 source=0.0_pReal)
 allocate(constitutive_dislotwin_Cmfptwin(maxNinstance),                            source=0.0_pReal)
 allocate(constitutive_dislotwin_Cthresholdtwin(maxNinstance),                      source=0.0_pReal)
 allocate(constitutive_dislotwin_SolidSolutionStrength(maxNinstance),               source=0.0_pReal)
 allocate(constitutive_dislotwin_L0(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislotwin_xc(maxNinstance),                                  source=0.0_pReal)
 allocate(constitutive_dislotwin_VcrossSlip(maxNinstance),                          source=0.0_pReal)
 allocate(constitutive_dislotwin_aTolRho(maxNinstance),                             source=0.0_pReal)
 allocate(constitutive_dislotwin_aTolTwinFrac(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_dislotwin_sbResistance(maxNinstance),                        source=0.0_pReal)
 allocate(constitutive_dislotwin_sbVelocity(maxNinstance),                          source=0.0_pReal)
 allocate(constitutive_dislotwin_sbQedge(maxNinstance),                             source=0.0_pReal)
 allocate(constitutive_dislotwin_SFE_0K(maxNinstance),                              source=0.0_pReal)
 allocate(constitutive_dislotwin_dSFE_dT(maxNinstance),                             source=0.0_pReal)
 allocate(constitutive_dislotwin_rhoEdge0(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(constitutive_dislotwin_rhoEdgeDip0(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(constitutive_dislotwin_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_burgersPerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_QedgePerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_v0PerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_Ndot0PerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_twinsizePerTwinFamily(lattice_maxNtwinFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_interaction_SlipTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_interaction_TwinSlip(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_interaction_TwinTwin(lattice_maxNinteraction,maxNinstance), &
                                                                                    source=0.0_pReal)
 allocate(constitutive_dislotwin_sbSv(6,6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), &
                                                                                    source=0.0_pReal)
 

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
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_TwinFamilies = count(lattice_NtwinSystem(:,phase) > 0_pInt)
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_SlipTwin =     maxval(lattice_interactionSlipTwin(:,:,phase))
       Nchunks_TwinSlip =     maxval(lattice_interactionTwinSlip(:,:,phase))
       Nchunks_TwinTwin =     maxval(lattice_interactionTwinTwin(:,:,phase))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_DISLOTWIN_ID) then          ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('plasticity','elasticity','lattice_structure', &                                       ! already known
           'covera_ratio','c/a_ratio','c/a', &
           'c11','c12','c13','c22','c23','c33','c44','c55','c66')
         cycle
       case ('(output)')
         constitutive_dislotwin_Noutput(instance) = constitutive_dislotwin_Noutput(instance) + 1_pInt
         constitutive_dislotwin_output(constitutive_dislotwin_Noutput(instance),instance) = &
                                                   IO_lc(IO_stringValue(line,positions,2_pInt))
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('edge_density')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = edge_density_ID
           case ('dipole_density')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = dipole_density_ID
           case ('shear_rate_slip')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = shear_rate_slip_ID
           case ('accumulated_shear_slip')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = accumulated_shear_slip_ID
           case ('mfp_slip')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = mfp_slip_ID
           case ('resolved_stress_slip')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = resolved_stress_slip_ID
           case ('edge_dipole_distance')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = edge_dipole_distance_ID
           case ('stress_exponent')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = stress_exponent_ID
           case ('twin_fraction')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = twin_fraction_ID
           case ('shear_rate_twin')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = shear_rate_twin_ID
           case ('accumulated_shear_twin')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = accumulated_shear_twin_ID
           case ('mfp_twin')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = mfp_twin_ID
           case ('resolved_stress_twin')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = resolved_stress_twin_ID
           case ('threshold_stress_twin')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = threshold_stress_twin_ID
           case ('resolved_stress_shearband')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = resolved_stress_shearband_ID
           case ('shear_rate_shearband')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = shear_rate_shearband_ID
           case ('sb_eigenvalues')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = sb_eigenvalues_ID
           case ('sb_eigenvectors')
             constitutive_dislotwin_outputID(constitutive_dislotwin_Noutput(instance),instance) = sb_eigenvectors_ID
           case default
             call IO_error(105_pInt,ext_msg=IO_stringValue(line,positions,2_pInt)//' ('//PLASTICITY_DISLOTWIN_label//')')
          end select
       case ('nslip')
         if (positions(1) < 1_pInt + Nchunks_SlipFamilies) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('ntwin')
         if (positions(1) < 1_pInt + Nchunks_TwinFamilies) &
           call IO_warning(51_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         Nchunks_TwinFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_dislotwin_Ntwin(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
       case ('rhoedge0')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_rhoEdge0(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('rhoedgedip0')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_rhoEdgeDip0(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('slipburgers')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_burgersPerSlipFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinburgers')
         do j = 1_pInt, Nchunks_TwinFamilies 
           constitutive_dislotwin_burgersPerTwinFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('qedge')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_QedgePerSlipFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('v0')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_v0PerSlipFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('ndot0')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_dislotwin_Ndot0PerTwinFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('twinsize')
         do j = 1_pInt, Nchunks_TwinFamilies
           constitutive_dislotwin_twinsizePerTwinFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('clambdaslip')
         do j = 1_pInt, Nchunks_SlipFamilies
           constitutive_dislotwin_CLambdaSlipPerSlipFamily(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('grainsize')
         constitutive_dislotwin_GrainSize(instance) = IO_floatValue(line,positions,2_pInt)
       case ('maxtwinfraction')
         constitutive_dislotwin_MaxTwinFraction(instance) = IO_floatValue(line,positions,2_pInt)
       case ('pexponent')
         constitutive_dislotwin_p(instance) = IO_floatValue(line,positions,2_pInt)
       case ('qexponent')
         constitutive_dislotwin_q(instance) = IO_floatValue(line,positions,2_pInt)
       case ('rexponent')
         constitutive_dislotwin_r(instance) = IO_floatValue(line,positions,2_pInt)
       case ('d0')
         constitutive_dislotwin_D0(instance) = IO_floatValue(line,positions,2_pInt)
       case ('qsd')
         constitutive_dislotwin_Qsd(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_rho')
         constitutive_dislotwin_aTolRho(instance) = IO_floatValue(line,positions,2_pInt)
       case ('atol_twinfrac')
         constitutive_dislotwin_aTolTwinFrac(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cmfptwin')
         constitutive_dislotwin_Cmfptwin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cthresholdtwin')
         constitutive_dislotwin_Cthresholdtwin(instance) = IO_floatValue(line,positions,2_pInt)
       case ('solidsolutionstrength')
         constitutive_dislotwin_SolidSolutionStrength(instance) = IO_floatValue(line,positions,2_pInt)
       case ('l0')
         constitutive_dislotwin_L0(instance) = IO_floatValue(line,positions,2_pInt)
       case ('xc')
              constitutive_dislotwin_xc(instance) = IO_floatValue(line,positions,2_pInt)
       case ('vcrossslip')
              constitutive_dislotwin_VcrossSlip(instance) = IO_floatValue(line,positions,2_pInt)
       case ('cedgedipmindistance')
         constitutive_dislotwin_CEdgeDipMinDistance(instance) = IO_floatValue(line,positions,2_pInt)
       case ('catomicvolume')
         constitutive_dislotwin_CAtomicVolume(instance) = IO_floatValue(line,positions,2_pInt)
       case ('interaction_slipslip','interactionslipslip')
         if (positions(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           constitutive_dislotwin_interaction_SlipSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_sliptwin','interactionsliptwin')
         if (positions(1) < 1_pInt + Nchunks_SlipTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_SlipTwin
           constitutive_dislotwin_interaction_SlipTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twinslip','interactiontwinslip')
         if (positions(1) < 1_pInt + Nchunks_TwinSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TwinSlip
           constitutive_dislotwin_interaction_TwinSlip(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('interaction_twintwin','interactiontwintwin')
         if (positions(1) < 1_pInt + Nchunks_TwinTwin) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
         do j = 1_pInt, Nchunks_TwinTwin
           constitutive_dislotwin_interaction_TwinTwin(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
       case ('sfe_0k')
         constitutive_dislotwin_SFE_0K(instance) = IO_floatValue(line,positions,2_pInt)
       case ('dsfe_dt')
         constitutive_dislotwin_dSFE_dT(instance) = IO_floatValue(line,positions,2_pInt)
       case ('shearbandresistance')
         constitutive_dislotwin_sbResistance(instance) = IO_floatValue(line,positions,2_pInt)
       case ('shearbandvelocity')
         constitutive_dislotwin_sbVelocity(instance) = IO_floatValue(line,positions,2_pInt)
       case ('qedgepersbsystem')
         constitutive_dislotwin_sbQedge(instance) = IO_floatValue(line,positions,2_pInt)
       case default
         call IO_error(210_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOTWIN_label//')')
     end select
   endif; endif
 enddo parsingFile
 
 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
    myPhase: if (phase_plasticity(phase) == PLASTICITY_dislotwin_ID) then
      instance = phase_plasticityInstance(phase)
      if (sum(constitutive_dislotwin_Nslip(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Nslip ('//PLASTICITY_DISLOTWIN_label//')')
      if (sum(constitutive_dislotwin_Ntwin(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Ntwin ('//PLASTICITY_DISLOTWIN_label//')')
      do f = 1_pInt,lattice_maxNslipFamily
        if (constitutive_dislotwin_Nslip(f,instance) > 0_pInt) then
          if (constitutive_dislotwin_rhoEdge0(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOTWIN_label//')')
          if (constitutive_dislotwin_rhoEdgeDip0(f,instance) < 0.0_pReal) & 
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOTWIN_label//')')
          if (constitutive_dislotwin_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOTWIN_label//')')
          if (constitutive_dislotwin_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOTWIN_label//')')
        endif
      enddo
      do f = 1_pInt,lattice_maxNtwinFamily
        if (constitutive_dislotwin_Ntwin(f,instance) > 0_pInt) then
          if (constitutive_dislotwin_burgersPerTwinFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='twinburgers ('//PLASTICITY_DISLOTWIN_label//')')
          if (constitutive_dislotwin_Ndot0PerTwinFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='ndot0 ('//PLASTICITY_DISLOTWIN_label//')')
        endif
      enddo
      if (constitutive_dislotwin_CAtomicVolume(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_D0(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_Qsd(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_SFE_0K(instance) == 0.0_pReal .and. constitutive_dislotwin_dSFE_dT(instance) == 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='SFE ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_aTolRho(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOTWIN_label//')')   
      if (constitutive_dislotwin_aTolTwinFrac(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='aTolTwinFrac ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_sbResistance(instance) < 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='sbResistance ('//PLASTICITY_DISLOTWIN_label//')')
      if (constitutive_dislotwin_sbVelocity(instance) < 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='sbVelocity ('//PLASTICITY_DISLOTWIN_label//')')
!--------------------------------------------------------------------------------------------------
! Determine total number of active slip or twin systems
      constitutive_dislotwin_Nslip(:,instance) = min(lattice_NslipSystem(:,phase),constitutive_dislotwin_Nslip(:,instance))
      constitutive_dislotwin_Ntwin(:,instance) = min(lattice_NtwinSystem(:,phase),constitutive_dislotwin_Ntwin(:,instance))
      constitutive_dislotwin_totalNslip(instance) = sum(constitutive_dislotwin_Nslip(:,instance))
      constitutive_dislotwin_totalNtwin(instance) = sum(constitutive_dislotwin_Ntwin(:,instance))
   endif myPhase
 enddo sanityChecks
 
!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip = maxval(constitutive_dislotwin_totalNslip)
 maxTotalNtwin = maxval(constitutive_dislotwin_totalNtwin)
 
 allocate(constitutive_dislotwin_burgersPerSlipSystem(maxTotalNslip, maxNinstance),    source=0.0_pReal)
 allocate(constitutive_dislotwin_burgersPerTwinSystem(maxTotalNtwin, maxNinstance),    source=0.0_pReal)
 allocate(constitutive_dislotwin_QedgePerSlipSystem(maxTotalNslip, maxNinstance),      source=0.0_pReal)
 allocate(constitutive_dislotwin_v0PerSlipSystem(maxTotalNslip, maxNinstance),         source=0.0_pReal)
 allocate(constitutive_dislotwin_Ndot0PerTwinSystem(maxTotalNtwin, maxNinstance),      source=0.0_pReal)
 allocate(constitutive_dislotwin_tau_r(maxTotalNtwin, maxNinstance),                   source=0.0_pReal)
 allocate(constitutive_dislotwin_twinsizePerTwinSystem(maxTotalNtwin, maxNinstance),   source=0.0_pReal)
 allocate(constitutive_dislotwin_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance),source=0.0_pReal)
 
 allocate(constitutive_dislotwin_interactionMatrix_SlipSlip(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(constitutive_dislotwin_interactionMatrix_SlipTwin(maxTotalNslip,maxTotalNtwin,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(constitutive_dislotwin_interactionMatrix_TwinSlip(maxTotalNtwin,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(constitutive_dislotwin_interactionMatrix_TwinTwin(maxTotalNtwin,maxTotalNtwin,maxNinstance), &
                                                                                       source=0.0_pReal)
 allocate(constitutive_dislotwin_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(phase) == PLASTICITY_dislotwin_ID) then
     instance = phase_plasticityInstance(phase)
 
     ns = constitutive_dislotwin_totalNslip(instance)
     nt = constitutive_dislotwin_totalNtwin(instance)

!--------------------------------------------------------------------------------------------------
! Determine size of state array
     constitutive_dislotwin_sizeDotState(instance) = int(size(CONSTITUTIVE_DISLOTWIN_listBasicSlipStates),pInt) * ns &
                                                   + int(size(CONSTITUTIVE_DISLOTWIN_listBasicTwinStates),pInt) * nt
     constitutive_dislotwin_sizeState(instance) = constitutive_dislotwin_sizeDotState(instance) &
                                        + int(size(CONSTITUTIVE_DISLOTWIN_listDependentSlipStates),pInt) * ns &
                                        + int(size(CONSTITUTIVE_DISLOTWIN_listDependentTwinStates),pInt) * nt

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,constitutive_dislotwin_Noutput(instance)
       select case(constitutive_dislotwin_outputID(o,instance))
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
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          constitutive_dislotwin_sizePostResult(o,instance) = mySize
          constitutive_dislotwin_sizePostResults(instance)  = constitutive_dislotwin_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
  
    !* Process slip related parameters ------------------------------------------------
 
     slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
       index_myFamily = sum(constitutive_dislotwin_Nslip(1:f-1_pInt,instance))                      ! index in truncated slip system list
       slipSystemsLoop: do j = 1_pInt,constitutive_dislotwin_Nslip(f,instance)

      !* Burgers vector, 
      !  dislocation velocity prefactor,
      !  mean free path prefactor,
      !  and minimum dipole distance
 
         constitutive_dislotwin_burgersPerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislotwin_burgersPerSlipFamily(f,instance)
 
         constitutive_dislotwin_QedgePerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislotwin_QedgePerSlipFamily(f,instance)
 
         constitutive_dislotwin_v0PerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislotwin_v0PerSlipFamily(f,instance)
 
         constitutive_dislotwin_CLambdaSlipPerSlipSystem(index_myFamily+j,instance) = &
         constitutive_dislotwin_CLambdaSlipPerSlipFamily(f,instance)
  
       !* Calculation of forest projections for edge dislocations
       !* Interaction matrices
  
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(constitutive_dislotwin_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,constitutive_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (slip)
             constitutive_dislotwin_forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,phase))+j,phase), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,phase))+k,phase)))
             constitutive_dislotwin_interactionMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislotwin_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
  
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(constitutive_dislotwin_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,constitutive_dislotwin_Ntwin(o,instance)                                   ! loop over (active) systems in other family (twin)
             constitutive_dislotwin_interactionMatrix_SlipTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislotwin_interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
  
       enddo slipSystemsLoop
     enddo slipFamiliesLoop
  
    !* Process twin related parameters ------------------------------------------------
    
     twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
       index_myFamily = sum(constitutive_dislotwin_Ntwin(1:f-1_pInt,instance))                      ! index in truncated twin system list
       twinSystemsLoop: do j = 1_pInt,constitutive_dislotwin_Ntwin(f,instance)
 
       !* Burgers vector,
       !  nucleation rate prefactor,
       !  and twin size
 
         constitutive_dislotwin_burgersPerTwinSystem(index_myFamily+j,instance)  = &
         constitutive_dislotwin_burgersPerTwinFamily(f,instance)

         constitutive_dislotwin_Ndot0PerTwinSystem(index_myFamily+j,instance)  = &
         constitutive_dislotwin_Ndot0PerTwinFamily(f,instance)

         constitutive_dislotwin_twinsizePerTwinSystem(index_myFamily+j,instance) = &
         constitutive_dislotwin_twinsizePerTwinFamily(f,instance)
 
       !* Rotate twin elasticity matrices
 
         index_otherFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase))                             ! index in full lattice twin list
         do l = 1_pInt,3_pInt ; do m = 1_pInt,3_pInt ; do n = 1_pInt,3_pInt ; do o = 1_pInt,3_pInt
           do p = 1_pInt,3_pInt ; do q = 1_pInt,3_pInt ; do r = 1_pInt,3_pInt ; do s = 1_pInt,3_pInt
             constitutive_dislotwin_Ctwin_3333(l,m,n,o,index_myFamily+j,instance) = &
             constitutive_dislotwin_Ctwin_3333(l,m,n,o,index_myFamily+j,instance) + &
               lattice_C3333(p,q,r,s,instance) * &
               lattice_Qtwin(l,p,index_otherFamily+j,phase) * &
               lattice_Qtwin(m,q,index_otherFamily+j,phase) * &
               lattice_Qtwin(n,r,index_otherFamily+j,phase) * &
               lattice_Qtwin(o,s,index_otherFamily+j,phase)
           enddo ; enddo ; enddo ; enddo
         enddo ; enddo ; enddo ; enddo
         constitutive_dislotwin_Ctwin_66(1:6,1:6,index_myFamily+j,instance) = &
           math_Mandel3333to66(constitutive_dislotwin_Ctwin_3333(1:3,1:3,1:3,1:3,index_myFamily+j,instance))
 
      !* Interaction matrices
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(constitutive_dislotwin_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,constitutive_dislotwin_Nslip(o,instance)                                   ! loop over (active) systems in other family (slip)
             constitutive_dislotwin_interactionMatrix_TwinSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislotwin_interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
 
         do o = 1_pInt,lattice_maxNtwinFamily
           index_otherFamily = sum(constitutive_dislotwin_Ntwin(1:o-1_pInt,instance))
           do k = 1_pInt,constitutive_dislotwin_Ntwin(o,instance)                                   ! loop over (active) systems in other family (twin)
             constitutive_dislotwin_interactionMatrix_TwinTwin(index_myFamily+j,index_otherFamily+k,instance) = &
                   constitutive_dislotwin_interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                 phase), instance )
         enddo; enddo
 
       enddo twinSystemsLoop
     enddo twinFamiliesLoop
   endif
 
 enddo initializeInstances
 
end subroutine constitutive_dislotwin_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
function constitutive_dislotwin_stateInit(instance,phase)
 use math, only: &
   pi
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_mu
 
 implicit none
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the plasticity
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the plasticity
 real(pReal), dimension(constitutive_dislotwin_sizeState(instance))  :: &
   constitutive_dislotwin_stateInit

 integer(pInt) :: i,j,f,ns,nt, index_myFamily
 real(pReal), dimension(constitutive_dislotwin_totalNslip(instance)) :: &
   rhoEdge0, &
   rhoEdgeDip0, &
   invLambdaSlip0, &
   MeanFreePathSlip0, &
   tauSlipThreshold0
 real(pReal), dimension(constitutive_dislotwin_totalNtwin(instance)) :: &
   MeanFreePathTwin0,TwinVolume0
 
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 constitutive_dislotwin_stateInit = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily   = sum(constitutive_dislotwin_Nslip(1:f-1_pInt,instance))                        ! index in truncated slip system list
   rhoEdge0(index_myFamily+1_pInt: &
            index_myFamily+constitutive_dislotwin_Nslip(f,instance)) = &
     constitutive_dislotwin_rhoEdge0(f,instance)
   rhoEdgeDip0(index_myFamily+1_pInt: &
               index_myFamily+constitutive_dislotwin_Nslip(f,instance)) = &
     constitutive_dislotwin_rhoEdgeDip0(f,instance)
 enddo
 
 constitutive_dislotwin_stateInit(1_pInt:ns)           = rhoEdge0
 constitutive_dislotwin_stateInit(ns+1_pInt:2_pInt*ns) = rhoEdgeDip0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent slip microstructural variables
 forall (i = 1_pInt:ns) &
   invLambdaSlip0(i) = sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislotwin_forestProjectionEdge(1:ns,i,instance)))/ &
                       constitutive_dislotwin_CLambdaSlipPerSlipSystem(i,instance)
 constitutive_dislotwin_stateInit(3_pInt*ns+2_pInt*nt+1:4_pInt*ns+2_pInt*nt) = invLambdaSlip0
 
 forall (i = 1_pInt:ns) &
   MeanFreePathSlip0(i) = &
     constitutive_dislotwin_GrainSize(instance)/(1.0_pReal+invLambdaSlip0(i)*constitutive_dislotwin_GrainSize(instance))
 constitutive_dislotwin_stateInit(5_pInt*ns+3_pInt*nt+1:6_pInt*ns+3_pInt*nt) = MeanFreePathSlip0
 
 forall (i = 1_pInt:ns) &
   tauSlipThreshold0(i) = constitutive_dislotwin_SolidSolutionStrength(instance) + &
     lattice_mu(phase)*constitutive_dislotwin_burgersPerSlipSystem(i,instance) * &
     sqrt(dot_product((rhoEdge0+rhoEdgeDip0),constitutive_dislotwin_interactionMatrix_SlipSlip(i,1:ns,instance)))
 constitutive_dislotwin_stateInit(6_pInt*ns+4_pInt*nt+1:7_pInt*ns+4_pInt*nt) = tauSlipThreshold0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent twin microstructural variables
 forall (j = 1_pInt:nt) &
   MeanFreePathTwin0(j) = constitutive_dislotwin_GrainSize(instance)
 constitutive_dislotwin_stateInit(6_pInt*ns+3_pInt*nt+1_pInt:6_pInt*ns+4_pInt*nt) = MeanFreePathTwin0
 
 forall (j = 1_pInt:nt) &
   TwinVolume0(j) = &
     (pi/4.0_pReal)*constitutive_dislotwin_twinsizePerTwinSystem(j,instance)*MeanFreePathTwin0(j)**(2.0_pReal)
 constitutive_dislotwin_stateInit(7_pInt*ns+5_pInt*nt+1_pInt:7_pInt*ns+6_pInt*nt) = TwinVolume0
 
end function constitutive_dislotwin_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
pure function constitutive_dislotwin_aTolState(instance)

 implicit none
 integer(pInt), intent(in) ::  &
   instance                                                                                         ! number specifying the current instance of the plasticity
 real(pReal), dimension(constitutive_dislotwin_sizeState(instance)) :: &
   constitutive_dislotwin_aTolState                                                                 ! relevant state values for the current instance of this plasticity

 ! Tolerance state for dislocation densities
 constitutive_dislotwin_aTolState(1_pInt:2_pInt*constitutive_dislotwin_totalNslip(instance)) = &
   constitutive_dislotwin_aTolRho(instance)

 ! Tolerance state for accumulated shear due to slip 
 constitutive_dislotwin_aTolState(2_pInt*constitutive_dislotwin_totalNslip(instance)+1_pInt: &
                                  3_pInt*constitutive_dislotwin_totalNslip(instance))=1e6_pReal
   
 
 ! Tolerance state for twin volume fraction
 constitutive_dislotwin_aTolState(3_pInt*constitutive_dislotwin_totalNslip(instance)+1_pInt: &
                                  3_pInt*constitutive_dislotwin_totalNslip(instance)+&
                                   constitutive_dislotwin_totalNtwin(instance)) = &
   constitutive_dislotwin_aTolTwinFrac(instance)

! Tolerance state for accumulated shear due to twin
 constitutive_dislotwin_aTolState(3_pInt*constitutive_dislotwin_totalNslip(instance)+ &
                                  constitutive_dislotwin_totalNtwin(instance)+1_pInt: &
                                  3_pInt*constitutive_dislotwin_totalNslip(instance)+ &
                                  2_pInt*constitutive_dislotwin_totalNtwin(instance)) = 1e6_pReal
   
end function constitutive_dislotwin_aTolState

!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_dislotwin_homogenizedC(state,ipc,ip,el)
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
    constitutive_dislotwin_homogenizedC
  integer(pInt), intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
    state                                                                                           !< microstructure state

 integer(pInt) :: instance,ns,nt,i,phase
 real(pReal) :: sumf
 
 !* Shortened notation
 phase = material_phase(ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((3_pInt*ns+1_pInt):(3_pInt*ns+nt))) ! safe for nt == 0
 
 !* Homogenized elasticity matrix
 constitutive_dislotwin_homogenizedC = (1.0_pReal-sumf)*lattice_C66(1:6,1:6,phase)
 do i=1_pInt,nt
    constitutive_dislotwin_homogenizedC = &
    constitutive_dislotwin_homogenizedC + state(ipc,ip,el)%p(3_pInt*ns+i)*lattice_C66(1:6,1:6,phase)
 enddo
 
 end function constitutive_dislotwin_homogenizedC
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislotwin_microstructure(temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   pi
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
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
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: &
   state                                                                                            !< microstructure state
 
 integer(pInt) :: &
   instance,phase,&
   ns,nt,s,t
 real(pReal) :: &
   sumf,sfe,x0
 real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: fOverStacksize
 
 !* Shortened notation
 phase = material_phase(ipc,ip,el)
 instance = phase_plasticityInstance(phase)
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 !* State: 1           :  ns         rho_edge
 !* State: ns+1        :  2*ns       rho_dipole
 !* State: 2*ns+1      :  3*ns       accumulated shear due to slip
 !* State: 3*ns+1      :  3*ns+nt    f
 !* State: 3*ns+nt+1   :  3*ns+2*nt  accumulated shear due to twin
 !* State: 3*ns+2*nt+1 :  4*ns+2*nt  1/lambda_slip
 !* State: 4*ns+2*nt+1 :  5*ns+2*nt  1/lambda_sliptwin
 !* State: 5*ns+2*nt+1 :  5*ns+3*nt  1/lambda_twin
 !* State: 5*ns+3*nt+1 :  6*ns+3*nt  mfp_slip
 !* State: 6*ns+3*nt+1 :  6*ns+4*nt  mfp_twin
 !* State: 6*ns+4*nt+1 :  7*ns+4*nt  threshold_stress_slip
 !* State: 7*ns+4*nt+1 :  7*ns+5*nt  threshold_stress_twin
 !* State: 7*ns+5*nt+1 :  7*ns+6*nt  twin volume
 
 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((3*ns+1):(3*ns+nt))) ! safe for nt == 0
 
 !* Stacking fault energy
 sfe = constitutive_dislotwin_SFE_0K(instance) + & 
       constitutive_dislotwin_dSFE_dT(instance) * Temperature
 
 !* rescaled twin volume fraction for topology
 forall (t = 1_pInt:nt) &
   fOverStacksize(t) = &
     state(ipc,ip,el)%p(3_pInt*ns+t)/constitutive_dislotwin_twinsizePerTwinSystem(t,instance)
 
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(3_pInt*ns+2_pInt*nt+s) = &
     sqrt(dot_product((state(ipc,ip,el)%p(1:ns)+state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)),&
                      constitutive_dislotwin_forestProjectionEdge(1:ns,s,instance)))/ &
     constitutive_dislotwin_CLambdaSlipPerSlipSystem(s,instance)
 
 !* 1/mean free distance between 2 twin stacks from different systems seen by a moving dislocation
 !$OMP CRITICAL (evilmatmul)
 state(ipc,ip,el)%p((4_pInt*ns+2_pInt*nt+1_pInt):(5_pInt*ns+2_pInt*nt)) = 0.0_pReal
 if (nt > 0_pInt .and. ns > 0_pInt) &
   state(ipc,ip,el)%p((4_pInt*ns+2_pInt*nt+1):(5_pInt*ns+2_pInt*nt)) = &
     matmul(constitutive_dislotwin_interactionMatrix_SlipTwin(1:ns,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)
 
 !* 1/mean free distance between 2 twin stacks from different systems seen by a growing twin
 !$OMP CRITICAL (evilmatmul)
 if (nt > 0_pInt) &
   state(ipc,ip,el)%p((5_pInt*ns+2_pInt*nt+1_pInt):(5_pInt*ns+3_pInt*nt)) = &
     matmul(constitutive_dislotwin_interactionMatrix_TwinTwin(1:nt,1:nt,instance),fOverStacksize(1:nt))/(1.0_pReal-sumf)
 !$OMP END CRITICAL (evilmatmul)
 
 !* mean free path between 2 obstacles seen by a moving dislocation
 do s = 1_pInt,ns
    if (nt > 0_pInt) then
       state(ipc,ip,el)%p(5_pInt*ns+3_pInt*nt+s) = &
         constitutive_dislotwin_GrainSize(instance)/(1.0_pReal+constitutive_dislotwin_GrainSize(instance)*&
         (state(ipc,ip,el)%p(3_pInt*ns+2_pInt*nt+s)+state(ipc,ip,el)%p(4_pInt*ns+2_pInt*nt+s)))
    else
       state(ipc,ip,el)%p(5_pInt*ns+s) = &
         constitutive_dislotwin_GrainSize(instance)/&
         (1.0_pReal+constitutive_dislotwin_GrainSize(instance)*(state(ipc,ip,el)%p(3_pInt*ns+s)))
    endif
 enddo
 
 !* mean free path between 2 obstacles seen by a growing twin
 forall (t = 1_pInt:nt) &
   state(ipc,ip,el)%p(6_pInt*ns+3_pInt*nt+t) = &
     (constitutive_dislotwin_Cmfptwin(instance)*constitutive_dislotwin_GrainSize(instance))/&
     (1.0_pReal+constitutive_dislotwin_GrainSize(instance)*state(ipc,ip,el)%p(5_pInt*ns+2_pInt*nt+t))
 
 !* threshold stress for dislocation motion
 forall (s = 1_pInt:ns) &
   state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+s) = constitutive_dislotwin_SolidSolutionStrength(instance)+ &
     lattice_mu(phase)*constitutive_dislotwin_burgersPerSlipSystem(s,instance)*&
     sqrt(dot_product((state(ipc,ip,el)%p(1:ns)+state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)),&
                      constitutive_dislotwin_interactionMatrix_SlipSlip(s,1:ns,instance)))
 
 !* threshold stress for growing twin
 forall (t = 1_pInt:nt) &
   state(ipc,ip,el)%p(7_pInt*ns+4_pInt*nt+t) = &
     constitutive_dislotwin_Cthresholdtwin(instance)*&
     (sfe/(3.0_pReal*constitutive_dislotwin_burgersPerTwinSystem(t,instance))+&
     3.0_pReal*constitutive_dislotwin_burgersPerTwinSystem(t,instance)*lattice_mu(phase)/&
     (constitutive_dislotwin_L0(instance)*constitutive_dislotwin_burgersPerSlipSystem(t,instance)))
 
 !* final twin volume after growth
 forall (t = 1_pInt:nt) &
   state(ipc,ip,el)%p(7_pInt*ns+5_pInt*nt+t) = &
     (pi/4.0_pReal)*constitutive_dislotwin_twinsizePerTwinSystem(t,instance)*state(ipc,ip,el)%p(6*ns+3*nt+t)**(2.0_pReal)
     
 !* equilibrium seperation of partial dislocations
 do t = 1_pInt,nt
   x0 = lattice_mu(phase)*constitutive_dislotwin_burgersPerTwinSystem(t,instance)**(2.0_pReal)/&
     (sfe*8.0_pReal*pi)*(2.0_pReal+lattice_nu(phase))/(1.0_pReal-lattice_nu(phase))
   constitutive_dislotwin_tau_r(t,instance)= &
        lattice_mu(phase)*constitutive_dislotwin_burgersPerTwinSystem(t,instance)/(2.0_pReal*pi)*&
        (1/(x0+constitutive_dislotwin_xc(instance))+cos(pi/3.0_pReal)/x0)
 enddo
 
end subroutine constitutive_dislotwin_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_spectralDecompositionSym33, &
   math_tensorproduct, &
   math_symmetric33, &
   math_mul33x3
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily,&
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID
 
 implicit none
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: Temperature
 real(pReal), dimension(6), intent(in) :: Tstar_v
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(inout) :: state
 real(pReal), dimension(3,3), intent(out) :: Lp
 real(pReal), dimension(9,9), intent(out) :: dLp_dTstar

 integer(pInt) :: instance,phase,ns,nt,f,i,j,k,l,m,n,index_myFamily,s1,s2
 real(pReal) :: sumf,StressRatio_p,StressRatio_pminus1,StressRatio_r,BoltzmannRatio,DotGamma0,Ndot0
 real(pReal), dimension(3,3,3,3) :: dLp_dTstar3333
 real(pReal), dimension(constitutive_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_slip,dgdot_dtauslip,tau_slip
 real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
    gdot_twin,dgdot_dtautwin,tau_twin
 real(pReal), dimension(6) :: gdot_sb,dgdot_dtausb,tau_sb
 real(pReal), dimension(3,3) :: eigVectors, sb_Smatrix
 real(pReal), dimension(3)   :: eigValues, sb_s, sb_m
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
 logical error
 
 !* Shortened notation
 phase = material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase)
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((3_pInt*ns+1_pInt):(3_pInt*ns+nt))) ! safe for nt == 0
 
 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar = 0.0_pReal
 
 !* Dislocation glide part
 gdot_slip = 0.0_pReal
 dgdot_dtauslip = 0.0_pReal
 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
   slipSystemsLoop: do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
 
      !* Stress ratios
      StressRatio_p = (abs(tau_slip(j))/state(ipc,ip,el)%p(6*ns+4*nt+j))**constitutive_dislotwin_p(instance)
      StressRatio_pminus1 = (abs(tau_slip(j))/state(ipc,ip,el)%p(6*ns+4*nt+j))**(constitutive_dislotwin_p(instance)-1.0_pReal)
      !* Boltzmann ratio
      BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
      !* Initial shear rates
      DotGamma0 = &
        state(ipc,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)*&
        constitutive_dislotwin_v0PerSlipSystem(j,instance)
 
      !* Shear rates due to slip
      gdot_slip(j) = (1.0_pReal - sumf) * DotGamma0 &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** constitutive_dislotwin_q(instance)) &
                   * sign(1.0_pReal,tau_slip(j))
 
      !* Derivatives of shear rates
      dgdot_dtauslip(j) = &
        ((abs(gdot_slip(j))*BoltzmannRatio*&
        constitutive_dislotwin_p(instance)*constitutive_dislotwin_q(instance))/state(ipc,ip,el)%p(6*ns+4*nt+j))*&
        StressRatio_pminus1*(1-StressRatio_p)**(constitutive_dislotwin_q(instance)-1.0_pReal)
 
      !* Plastic velocity gradient for dislocation glide
      Lp = Lp + gdot_slip(j)*lattice_Sslip(:,:,1,index_myFamily+i,phase)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtauslip(j)*&
                                  lattice_Sslip(k,l,1,index_myFamily+i,phase)*&
                                  lattice_Sslip(m,n,1,index_myFamily+i,phase)
   enddo slipSystemsLoop
 enddo slipFamiliesLoop
 
 !* Shear banding (shearband) part
 if(constitutive_dislotwin_sbVelocity(instance) /= 0.0_pReal .and. &
    constitutive_dislotwin_sbResistance(instance) /= 0.0_pReal) then
   gdot_sb = 0.0_pReal
   dgdot_dtausb = 0.0_pReal
   call math_spectralDecompositionSym33(math_Mandel6to33(Tstar_v),eigValues,eigVectors, error)
   do j = 1_pInt,6_pInt
     sb_s = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_sComposition(1:3,j))
     sb_m = 0.5_pReal*sqrt(2.0_pReal)*math_mul33x3(eigVectors,sb_mComposition(1:3,j))
     sb_Smatrix = math_tensorproduct(sb_s,sb_m)
     constitutive_dislotwin_sbSv(1:6,j,ipc,ip,el) = math_Mandel33to6(math_symmetric33(sb_Smatrix))
   
     !* Calculation of Lp
     !* Resolved shear stress on shear banding system
     tau_sb(j) = dot_product(Tstar_v,constitutive_dislotwin_sbSv(1:6,j,ipc,ip,el))
   
     !* Stress ratios
     StressRatio_p = (abs(tau_sb(j))/constitutive_dislotwin_sbResistance(instance))**constitutive_dislotwin_p(instance)
     StressRatio_pminus1 = (abs(tau_sb(j))/constitutive_dislotwin_sbResistance(instance))&
                                                                          **(constitutive_dislotwin_p(instance)-1.0_pReal)
     !* Boltzmann ratio
     BoltzmannRatio = constitutive_dislotwin_sbQedge(instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = constitutive_dislotwin_sbVelocity(instance)
 
     !* Shear rates due to shearband
     gdot_sb(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**constitutive_dislotwin_q(instance))*&
                  sign(1.0_pReal,tau_sb(j))
                  
     !* Derivatives of shear rates
     dgdot_dtausb(j) = &
       ((abs(gdot_sb(j))*BoltzmannRatio*&
       constitutive_dislotwin_p(instance)*constitutive_dislotwin_q(instance))/constitutive_dislotwin_sbResistance(instance))*&
       StressRatio_pminus1*(1_pInt-StressRatio_p)**(constitutive_dislotwin_q(instance)-1.0_pReal)
 
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
 
 !* Mechanical twinning part
 gdot_twin = 0.0_pReal
 dgdot_dtautwin = 0.0_pReal
 j = 0_pInt
 twinFamiliesLoop: do f = 1_pInt,lattice_maxNtwinFamily
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase)) ! at which index starts my family
   twinSystemsLoop: do i = 1_pInt,constitutive_dislotwin_Ntwin(f,instance)
      j = j+1_pInt
 
      !* Calculation of Lp
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))
 
      !* Stress ratios
      StressRatio_r = (state(ipc,ip,el)%p(7*ns+4*nt+j)/tau_twin(j))**constitutive_dislotwin_r(instance)
 
      !* Shear rates and their derivatives due to twin
      if ( tau_twin(j) > 0.0_pReal ) then
        select case(lattice_structure(phase))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < constitutive_dislotwin_tau_r(j,instance)) then
              Ndot0=(abs(gdot_slip(s1))*(state(ipc,ip,el)%p(s2)+state(ipc,ip,el)%p(ns+s2))+&
                     abs(gdot_slip(s2))*(state(ipc,ip,el)%p(s1)+state(ipc,ip,el)%p(ns+s1)))/&
                    (constitutive_dislotwin_L0(instance)*constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-constitutive_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                    (constitutive_dislotwin_tau_r(j,instance)-tau_twin(j))))
            else
              Ndot0=0.0_pReal
            end if
          case default
            Ndot0=constitutive_dislotwin_Ndot0PerTwinSystem(j,instance)
        end select
        gdot_twin(j) = &
          (constitutive_dislotwin_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,phase)*&
          state(ipc,ip,el)%p(7*ns+5*nt+j)*Ndot0*exp(-StressRatio_r)
        dgdot_dtautwin(j) = ((gdot_twin(j)*constitutive_dislotwin_r(instance))/tau_twin(j))*StressRatio_r
      endif
 
      !* Plastic velocity gradient for mechanical twinning
      Lp = Lp + gdot_twin(j)*lattice_Stwin(:,:,index_myFamily+i,phase)
 
      !* Calculation of the tangent of Lp
      forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + dgdot_dtautwin(j)*&
                                  lattice_Stwin(k,l,index_myFamily+i,phase)*&
                                  lattice_Stwin(m,n,index_myFamily+i,phase)
   enddo twinSystemsLoop
 enddo twinFamiliesLoop
 
 dLp_dTstar = math_Plain3333to99(dLp_dTstar3333)
 
end subroutine constitutive_dislotwin_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
pure function constitutive_dislotwin_dotState(Tstar_v,Temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   pi
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 use lattice,  only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_sheartwin, &
   lattice_mu, &
   lattice_structure, &
   lattice_fcc_twinNucleationSlipPair, &
   LATTICE_fcc_ID

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
 real(pReal), dimension(constitutive_dislotwin_sizeDotState(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   constitutive_dislotwin_dotState

 integer(pInt) :: instance,phase,ns,nt,f,i,j,index_myFamily,s1,s2
 real(pReal) :: sumf,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,&
             EdgeDipMinDistance,AtomicVolume,VacancyDiffusion,StressRatio_r,Ndot0
 real(pReal), dimension(constitutive_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
 gdot_slip,tau_slip,DotRhoMultiplication,EdgeDipDistance,DotRhoEdgeEdgeAnnihilation,DotRhoEdgeDipAnnihilation,&
 ClimbVelocity,DotRhoEdgeDipClimb,DotRhoDipFormation
 real(pReal), dimension(constitutive_dislotwin_totalNtwin(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
              tau_twin
 
 !* Shortened notation
 phase = material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase)
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((3_pInt*ns+1_pInt):(3_pInt*ns+nt))) ! safe for nt == 0
 
 constitutive_dislotwin_dotState = 0.0_pReal
 
 !* Dislocation density evolution
 gdot_slip = 0.0_pReal
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
   do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
      j = j+1_pInt
 
 
      !* Resolved shear stress on slip system
      tau_slip(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
      !* Stress ratios
      StressRatio_p = (abs(tau_slip(j))/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                               constitutive_dislotwin_p(instance)
      StressRatio_pminus1 = (abs(tau_slip(j))/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                   (constitutive_dislotwin_p(instance)-1.0_pReal)
      !* Boltzmann ratio
      BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
      !* Initial shear rates
      DotGamma0 = &
        state(ipc,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)*&
        constitutive_dislotwin_v0PerSlipSystem(j,instance)
 
      !* Shear rates due to slip
      gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**constitutive_dislotwin_q(instance))*&
                     sign(1.0_pReal,tau_slip(j))
 
      !* Multiplication
      DotRhoMultiplication(j) = abs(gdot_slip(j))/&
                                (constitutive_dislotwin_burgersPerSlipSystem(j,instance)*state(ipc,ip,el)%p(5*ns+3*nt+j))
 
      !* Dipole formation
      EdgeDipMinDistance = &
        constitutive_dislotwin_CEdgeDipMinDistance(instance)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)
      if (tau_slip(j) == 0.0_pReal) then
        DotRhoDipFormation(j) = 0.0_pReal
      else
        EdgeDipDistance(j) = &
          (3.0_pReal*lattice_mu(phase)*constitutive_dislotwin_burgersPerSlipSystem(j,instance))/&
          (16.0_pReal*pi*abs(tau_slip(j)))
      if (EdgeDipDistance(j)>state(ipc,ip,el)%p(5*ns+3*nt+j)) EdgeDipDistance(j)=state(ipc,ip,el)%p(5*ns+3*nt+j)
      if (EdgeDipDistance(j)<EdgeDipMinDistance) EdgeDipDistance(j)=EdgeDipMinDistance
        DotRhoDipFormation(j) = &
          ((2.0_pReal*EdgeDipDistance(j))/constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
          state(ipc,ip,el)%p(j)*abs(gdot_slip(j))
      endif
 
      !* Spontaneous annihilation of 2 single edge dislocations
      DotRhoEdgeEdgeAnnihilation(j) = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
        state(ipc,ip,el)%p(j)*abs(gdot_slip(j))
 
      !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
      DotRhoEdgeDipAnnihilation(j) = &
        ((2.0_pReal*EdgeDipMinDistance)/constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
        state(ipc,ip,el)%p(ns+j)*abs(gdot_slip(j))
 
      !* Dislocation dipole climb
      AtomicVolume = &
        constitutive_dislotwin_CAtomicVolume(instance)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)**(3.0_pReal)
      VacancyDiffusion = &
        constitutive_dislotwin_D0(instance)*exp(-constitutive_dislotwin_Qsd(instance)/(kB*Temperature))
      if (tau_slip(j) == 0.0_pReal) then
        DotRhoEdgeDipClimb(j) = 0.0_pReal
      else
        ClimbVelocity(j) = &
          ((3.0_pReal*lattice_mu(phase)*VacancyDiffusion*AtomicVolume)/(2.0_pReal*pi*kB*Temperature))*&
          (1/(EdgeDipDistance(j)+EdgeDipMinDistance))
        DotRhoEdgeDipClimb(j) = &
          (4.0_pReal*ClimbVelocity(j)*state(ipc,ip,el)%p(ns+j))/(EdgeDipDistance(j)-EdgeDipMinDistance)
      endif
 
      !* Edge dislocation density rate of change
      constitutive_dislotwin_dotState(j) = &
        DotRhoMultiplication(j)-DotRhoDipFormation(j)-DotRhoEdgeEdgeAnnihilation(j)
 
      !* Edge dislocation dipole density rate of change
      constitutive_dislotwin_dotState(ns+j) = &
        DotRhoDipFormation(j)-DotRhoEdgeDipAnnihilation(j)-DotRhoEdgeDipClimb(j)
 
      !* Dotstate for accumulated shear due to slip
      constitutive_dislotwin_dotstate(2_pInt*ns+j) = gdot_slip(j)
 
   enddo
 enddo
 
 !* Twin volume fraction evolution
 j = 0_pInt
 do f = 1_pInt,lattice_maxNtwinFamily                                 ! loop over all twin families
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase)) ! at which index starts my family
   do i = 1_pInt,constitutive_dislotwin_Ntwin(f,instance)          ! process each (active) twin system in family
      j = j+1_pInt
 
      !* Resolved shear stress on twin system
      tau_twin(j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))
      !* Stress ratios
      StressRatio_r = (state(ipc,ip,el)%p(7*ns+4*nt+j)/tau_twin(j))**constitutive_dislotwin_r(instance)
 
      !* Shear rates and their derivatives due to twin
      if ( tau_twin(j) > 0.0_pReal ) then
        select case(lattice_structure(phase))
          case (LATTICE_fcc_ID)
            s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
            s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
            if (tau_twin(j) < constitutive_dislotwin_tau_r(j,instance)) then
              Ndot0=(abs(gdot_slip(s1))*(state(ipc,ip,el)%p(s2)+state(ipc,ip,el)%p(ns+s2))+&
                     abs(gdot_slip(s2))*(state(ipc,ip,el)%p(s1)+state(ipc,ip,el)%p(ns+s1)))/&
                    (constitutive_dislotwin_L0(instance)*constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
                    (1.0_pReal-exp(-constitutive_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                    (constitutive_dislotwin_tau_r(j,instance)-tau_twin(j))))
            else
              Ndot0=0.0_pReal
            end if
          case default
            Ndot0=constitutive_dislotwin_Ndot0PerTwinSystem(j,instance)
        end select
        constitutive_dislotwin_dotState(3_pInt*ns+j) = &
          (constitutive_dislotwin_MaxTwinFraction(instance)-sumf)*&
          state(ipc,ip,el)%p(7_pInt*ns+5_pInt*nt+j)*Ndot0*exp(-StressRatio_r)
        !* Dotstate for accumulated shear due to twin
        constitutive_dislotwin_dotstate(3_pInt*ns+nt+j) = constitutive_dislotwin_dotState(3_pInt*ns+j) * &
                                                          lattice_sheartwin(index_myfamily+i,phase)
      endif
   enddo
 enddo
 
end function constitutive_dislotwin_dotState

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_dislotwin_postResults(Tstar_v,Temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use math, only: &
   pi, &
   math_Mandel6to33, &
   math_spectralDecompositionSym33
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains,&
   material_phase, &
   phase_plasticityInstance,& 
   phase_Noutput
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

 real(pReal), dimension(constitutive_dislotwin_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           constitutive_dislotwin_postResults

 integer(pInt) :: &
   instance,phase,&
   ns,nt,&
   f,o,i,c,j,index_myFamily,&
   s1,s2
 real(pReal) :: sumf,tau,StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0,StressRatio_r,Ndot0,dgdot_dtauslip
 real(preal), dimension(constitutive_dislotwin_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip
 real(pReal), dimension(3,3) :: eigVectors
 real(pReal), dimension (3) :: eigValues
 logical :: error
 
 !* Shortened notation
 phase = material_phase(ipc,ip,el)
 instance  = phase_plasticityInstance(phase)
 ns = constitutive_dislotwin_totalNslip(instance)
 nt = constitutive_dislotwin_totalNtwin(instance)
 
 !* Total twin volume fraction
 sumf = sum(state(ipc,ip,el)%p((3_pInt*ns+1_pInt):(3_pInt*ns+nt))) ! safe for nt == 0
 
 !* Required output
 c = 0_pInt
 constitutive_dislotwin_postResults = 0.0_pReal
 
 !* Spectral decomposition of stress
 call math_spectralDecompositionSym33(math_Mandel6to33(Tstar_v),eigValues,eigVectors, error)
 
 do o = 1_pInt,phase_Noutput(material_phase(ipc,ip,el))
    select case(constitutive_dislotwin_outputID(o,instance))
 
      case (edge_density_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(1_pInt:ns)
        c = c + ns
      case (dipole_density_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+ns) = state(ipc,ip,el)%p(ns+1_pInt:2_pInt*ns)
        c = c + ns
      case (shear_rate_slip_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
           do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
              j = j + 1_pInt
 
              !* Resolved shear stress on slip system
              tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
              !* Stress ratios
              StressRatio_p = (abs(tau)/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                                constitutive_dislotwin_p(instance)
              StressRatio_pminus1 = (abs(tau)/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                    (constitutive_dislotwin_p(instance)-1.0_pReal)
              !* Boltzmann ratio
              BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
              !* Initial shear rates
              DotGamma0 = &
                state(ipc,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)* &
                constitutive_dislotwin_v0PerSlipSystem(j,instance)
 
              !* Shear rates due to slip
              constitutive_dislotwin_postResults(c+j) = &
                DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                                           constitutive_dislotwin_q(instance))*sign(1.0_pReal,tau)
        enddo ; enddo
        c = c + ns
      case (accumulated_shear_slip_ID)
       constitutive_dislotwin_postResults(c+1_pInt:c+ns) = &
                                state(ipc,ip,el)%p((2_pInt*ns+1_pInt):(3_pInt*ns))
        c = c + ns
      case (mfp_slip_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+ns) =&
                                state(ipc,ip,el)%p((5_pInt*ns+3_pInt*nt+1_pInt):(6_pInt*ns+3_pInt*nt))
        c = c + ns
      case (resolved_stress_slip_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
           do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
              j = j + 1_pInt
              constitutive_dislotwin_postResults(c+j) =&
                                dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
        enddo; enddo
        c = c + ns
      case (threshold_stress_slip_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+ns) = &
                                state(ipc,ip,el)%p((6_pInt*ns+4_pInt*nt+1_pInt):(7_pInt*ns+4_pInt*nt))
        c = c + ns
      case (edge_dipole_distance_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
           do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
              j = j + 1_pInt
              constitutive_dislotwin_postResults(c+j) = &
                (3.0_pReal*lattice_mu(phase)*constitutive_dislotwin_burgersPerSlipSystem(j,instance))/&
                (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))))
              constitutive_dislotwin_postResults(c+j) = min(constitutive_dislotwin_postResults(c+j),state(ipc,ip,el)%p(5*ns+3*nt+j))
 !            constitutive_dislotwin_postResults(c+j) = max(constitutive_dislotwin_postResults(c+j),state(ipc,ip,el)%p(4*ns+2*nt+j))
        enddo; enddo
        c = c + ns
       case (resolved_stress_shearband_ID)
         do j = 1_pInt,6_pInt                                                ! loop over all shearband families
            constitutive_dislotwin_postResults(c+j) = dot_product(Tstar_v, constitutive_dislotwin_sbSv(1:6,j,ipc,ip,el))
         enddo
         c = c + 6_pInt
       case (shear_rate_shearband_ID)
         do j = 1_pInt,6_pInt                                                ! loop over all shearbands
              !* Resolved shear stress on shearband system
              tau = dot_product(Tstar_v,constitutive_dislotwin_sbSv(1:6,j,ipc,ip,el))
              !* Stress ratios
              StressRatio_p = (abs(tau)/constitutive_dislotwin_sbResistance(instance))**constitutive_dislotwin_p(instance)
              StressRatio_pminus1 = (abs(tau)/constitutive_dislotwin_sbResistance(instance))&
                                                                            **(constitutive_dislotwin_p(instance)-1.0_pReal)
              !* Boltzmann ratio
              BoltzmannRatio = constitutive_dislotwin_sbQedge(instance)/(kB*Temperature)
              !* Initial shear rates
              DotGamma0 = constitutive_dislotwin_sbVelocity(instance)
 
              !* Shear rates due to slip
              constitutive_dislotwin_postResults(c+j) = &
          DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**constitutive_dislotwin_q(instance))*sign(1.0_pReal,tau)
         enddo 
        c = c + 6_pInt
      case (twin_fraction_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+nt) = state(ipc,ip,el)%p((3_pInt*ns+1_pInt):(3_pInt*ns+nt))
        c = c + nt
      case (shear_rate_twin_ID)
        if (nt > 0_pInt) then
        
          j = 0_pInt
          do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
             index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
             do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
                j = j + 1_pInt
 
               !* Resolved shear stress on slip system
               tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
               !* Stress ratios
               StressRatio_p = (abs(tau)/state(ipc,ip,el)%p(5_pInt*ns+3_pInt*nt+j))**&
                                                                constitutive_dislotwin_p(instance)
               StressRatio_pminus1 = (abs(tau)/state(ipc,ip,el)%p(5_pInt*ns+3_pInt*nt+j))**&
                                                    (constitutive_dislotwin_p(instance)-1.0_pReal)
               !* Boltzmann ratio
               BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
               !* Initial shear rates
               DotGamma0 = &
                 state(ipc,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)* &
                 constitutive_dislotwin_v0PerSlipSystem(j,instance)
 
               !* Shear rates due to slip
               gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                                           constitutive_dislotwin_q(instance))*sign(1.0_pReal,tau)
          enddo;enddo
        
          j = 0_pInt
          do f = 1_pInt,lattice_maxNtwinFamily                                 ! loop over all twin families
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase))  ! at which index starts my family
            do i = 1,constitutive_dislotwin_Ntwin(f,instance)                ! process each (active) twin system in family
              j = j + 1_pInt
 
              !* Resolved shear stress on twin system
              tau = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))
              !* Stress ratios
              StressRatio_r = (state(ipc,ip,el)%p(7_pInt*ns+4_pInt*nt+j)/tau)**constitutive_dislotwin_r(instance)
 
              !* Shear rates due to twin
              if ( tau > 0.0_pReal ) then
                select case(lattice_structure(phase))
                  case (LATTICE_fcc_ID)
                  s1=lattice_fcc_twinNucleationSlipPair(1,index_myFamily+i)
                  s2=lattice_fcc_twinNucleationSlipPair(2,index_myFamily+i)
                  if (tau < constitutive_dislotwin_tau_r(j,instance)) then
                    Ndot0=(abs(gdot_slip(s1))*(state(ipc,ip,el)%p(s2)+state(ipc,ip,el)%p(ns+s2))+&
                           abs(gdot_slip(s2))*(state(ipc,ip,el)%p(s1)+state(ipc,ip,el)%p(ns+s1)))/&
                          (constitutive_dislotwin_L0(instance)*&
                           constitutive_dislotwin_burgersPerSlipSystem(j,instance))*&
                          (1.0_pReal-exp(-constitutive_dislotwin_VcrossSlip(instance)/(kB*Temperature)*&
                          (constitutive_dislotwin_tau_r(j,instance)-tau)))
                  else
                    Ndot0=0.0_pReal
                  end if
                  case default
                    Ndot0=constitutive_dislotwin_Ndot0PerTwinSystem(j,instance)
                end select
                constitutive_dislotwin_postResults(c+j) = &
                  (constitutive_dislotwin_MaxTwinFraction(instance)-sumf)*lattice_shearTwin(index_myFamily+i,phase)*&
                  state(ipc,ip,el)%p(7_pInt*ns+5_pInt*nt+j)*Ndot0*exp(-StressRatio_r)
              endif
 
          enddo ; enddo
        endif
        c = c + nt
      case (accumulated_shear_twin_ID)
       constitutive_dislotwin_postResults(c+1_pInt:c+nt) = state(ipc,ip,el)%p((3_pInt*ns+nt+1_pInt):(3_pInt*ns+2_pInt*nt))
        c = c + nt     
      case (mfp_twin_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+nt) = state(ipc,ip,el)%p((6_pInt*ns+3_pInt*nt+1_pInt):(6_pInt*ns+4_pInt*nt))
        c = c + nt
      case (resolved_stress_twin_ID)
        if (nt > 0_pInt) then
          j = 0_pInt
          do f = 1_pInt,lattice_maxNtwinFamily                                 ! loop over all slip families
            index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,phase))  ! at which index starts my family
            do i = 1_pInt,constitutive_dislotwin_Ntwin(f,instance)           ! process each (active) slip system in family
              j = j + 1_pInt
              constitutive_dislotwin_postResults(c+j) = dot_product(Tstar_v,lattice_Stwin_v(:,index_myFamily+i,phase))
          enddo; enddo
        endif
        c = c + nt
      case (threshold_stress_twin_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+nt) = state(ipc,ip,el)%p((7_pInt*ns+4_pInt*nt+1_pInt):(7_pInt*ns+5_pInt*nt))
        c = c + nt
      case (stress_exponent_ID)
        j = 0_pInt
        do f = 1_pInt,lattice_maxNslipFamily                                 ! loop over all slip families
          index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase)) ! at which index starts my family
          do i = 1_pInt,constitutive_dislotwin_Nslip(f,instance)          ! process each (active) slip system in family
             j = j + 1_pInt
 
             !* Resolved shear stress on slip system
             tau = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,phase))
             !* Stress ratios
             StressRatio_p = (abs(tau)/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                               constitutive_dislotwin_p(instance)
             StressRatio_pminus1 = (abs(tau)/state(ipc,ip,el)%p(6_pInt*ns+4_pInt*nt+j))**&
                                                   (constitutive_dislotwin_p(instance)-1.0_pReal)
             !* Boltzmann ratio
             BoltzmannRatio = constitutive_dislotwin_QedgePerSlipSystem(j,instance)/(kB*Temperature)
             !* Initial shear rates
             DotGamma0 = &
               state(ipc,ip,el)%p(j)*constitutive_dislotwin_burgersPerSlipSystem(j,instance)* &
               constitutive_dislotwin_v0PerSlipSystem(j,instance)
 
             !* Shear rates due to slip
             gdot_slip(j) = DotGamma0*exp(-BoltzmannRatio*(1_pInt-StressRatio_p)**&
                                          constitutive_dislotwin_q(instance))*sign(1.0_pReal,tau)
 
             !* Derivatives of shear rates
             dgdot_dtauslip = &
               ((abs(gdot_slip(j))*BoltzmannRatio*&
               constitutive_dislotwin_p(instance)*constitutive_dislotwin_q(instance))/state(ipc,ip,el)%p(6*ns+4*nt+j))*&
               StressRatio_pminus1*(1_pInt-StressRatio_p)**(constitutive_dislotwin_q(instance)-1.0_pReal)
 
             !* Stress exponent
             if (gdot_slip(j)==0.0_pReal) then
               constitutive_dislotwin_postResults(c+j) = 0.0_pReal
             else
               constitutive_dislotwin_postResults(c+j) = (tau/gdot_slip(j))*dgdot_dtauslip
             endif
         enddo ; enddo
         c = c + ns
      case (sb_eigenvalues_ID)
        forall (j = 1_pInt:3_pInt) &
        constitutive_dislotwin_postResults(c+j) = eigValues(j)
        c = c + 3_pInt
      case (sb_eigenvectors_ID)
        constitutive_dislotwin_postResults(c+1_pInt:c+9_pInt) = reshape(eigVectors,[9])
        c = c + 9_pInt
    end select
 enddo
end function constitutive_dislotwin_postResults

end module constitutive_dislotwin
