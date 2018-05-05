!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module plastic_disloUCLA
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   plastic_disloUCLA_sizePostResults                                                                !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   plastic_disloUCLA_sizePostResult                                                                 !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   plastic_disloUCLA_output                                                                         !< name of each post result output

 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   plastic_disloUCLA_Noutput                                                                        !< number of outputs per instance of this plasticity 

 integer(pInt),                       dimension(:),           allocatable,         private :: &
   plastic_disloUCLA_totalNslip                                                                     !< total number of active slip systems for each instance

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   plastic_disloUCLA_Nslip                                                                          !< number of active slip systems for each family and instance


 real(pReal),                         dimension(:),           allocatable,         private :: &
   plastic_disloUCLA_CAtomicVolume, &                                                               !< atomic volume in Bugers vector unit
   plastic_disloUCLA_D0, &                                                                          !< prefactor for self-diffusion coefficient
   plastic_disloUCLA_Qsd, &                                                                         !< activation energy for dislocation climb
   plastic_disloUCLA_GrainSize, &                                                                   !< grain size
   plastic_disloUCLA_CEdgeDipMinDistance, &                                                         !<
   plastic_disloUCLA_SolidSolutionStrength, &                                                       !< Strength due to elements in solid solution
   plastic_disloUCLA_dipoleFormationFactor, &                                                       !< scaling factor for dipole formation: 0: off, 1: on. other values not useful
   plastic_disloUCLA_aTolRho                                                                        !< absolute tolerance for integration of dislocation density

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   plastic_disloUCLA_rhoEdge0, &                                                                    !< initial edge dislocation density per slip system for each family and instance
   plastic_disloUCLA_rhoEdgeDip0, &                                                                 !< initial edge dipole density per slip system for each family and instance
   plastic_disloUCLA_burgersPerSlipFamily, &                                                        !< absolute length of burgers vector [m] for each slip family and instance
   plastic_disloUCLA_burgersPerSlipSystem, &                                                        !< absolute length of burgers vector [m] for each slip system and instance
   plastic_disloUCLA_QedgePerSlipFamily, &                                                          !< activation energy for glide [J] for each slip family and instance
   plastic_disloUCLA_QedgePerSlipSystem, &                                                          !< activation energy for glide [J] for each slip system and instance
   plastic_disloUCLA_v0PerSlipFamily, &                                                             !< dislocation velocity prefactor [m/s] for each family and instance
   plastic_disloUCLA_v0PerSlipSystem, &                                                             !< dislocation velocity prefactor [m/s] for each slip system and instance
   plastic_disloUCLA_tau_peierlsPerSlipFamily, &                                                    !< Peierls stress [Pa] for each family and instance
   plastic_disloUCLA_CLambdaSlipPerSlipFamily, &                                                    !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   plastic_disloUCLA_CLambdaSlipPerSlipSystem, &                                                    !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   plastic_disloUCLA_interaction_SlipSlip, &                                                        !< coefficients for slip-slip interaction for each interaction type and instance
   plastic_disloUCLA_pPerSlipFamily, &                                                              !< p-exponent in glide velocity
   plastic_disloUCLA_qPerSlipFamily, &                                                              !< q-exponent in glide velocity       
   !* mobility law parameters                                                                                                                                                      
   plastic_disloUCLA_kinkheight, &                                                                  !< height of the kink pair                                                      
   plastic_disloUCLA_omega, &                                                                       !< attempt frequency for kink pair nucleation                                   
   plastic_disloUCLA_kinkwidth, &                                                                   !< width of the kink pair                                                       
   plastic_disloUCLA_dislolength, &                                                                 !< dislocation length (lamda)                                                   
   plastic_disloUCLA_friction, &                                                                    !< friction coeff. B (kMC)
   !*    
   plastic_disloUCLA_nonSchmidCoeff                                                                 !< non-Schmid coefficients (bcc)
 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   plastic_disloUCLA_interactionMatrix_SlipSlip, &                                                  !< interaction matrix of the different slip systems for each instance
   plastic_disloUCLA_forestProjectionEdge                                                           !< matrix of forest projections of edge dislocations for each instance

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
                 stress_exponent_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   plastic_disloUCLA_outputID                                                                       !< ID of each post result output
 
 type, private :: tDisloUCLAState 
     real(pReal), pointer, dimension(:,:) :: &
       rhoEdge, &
       rhoEdgeDip, &
       accshear_slip, &
       invLambdaSlip, &
       mfp_slip, &
       threshold_stress_slip
 end type 
 type(tDisloUCLAState ), allocatable, dimension(:), private :: &
   state, &
   state0, &
   dotState
 
 public :: &
   plastic_disloUCLA_init, &
   plastic_disloUCLA_microstructure, &
   plastic_disloUCLA_LpAndItsTangent, &
   plastic_disloUCLA_dotState, &
   plastic_disloUCLA_postResults
 private :: &
   plastic_disloUCLA_stateInit, &
   plastic_disloUCLA_aTolState

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
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
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_DISLOUCLA_ID, &
   material_phase, &  
   plasticState, & 
   MATERIAL_partPhase
 use lattice
 use numerics,only: &
   numerics_integrator
 
 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,maxTotalNslip,&
                  f,instance,j,k,o,ns, &
                  Nchunks_SlipSlip = 0_pInt, &
                  Nchunks_SlipFamilies = 0_pInt,Nchunks_nonSchmid = 0_pInt, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip
  
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOUCLA_label//' init  -+>>>'
 write(6,'(/,a)')   ' Cereceda et al., International Journal of Plasticity 78, 2016, 242-256'
 write(6,'(/,a)')   ' http://dx.doi.org/10.1016/j.ijplas.2015.09.002'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOUCLA_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_disloUCLA_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(plastic_disloUCLA_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(plastic_disloUCLA_output(maxval(phase_Noutput),maxNinstance))
          plastic_disloUCLA_output = ''
 allocate(plastic_disloUCLA_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(plastic_disloUCLA_Noutput(maxNinstance),                             source=0_pInt)
 allocate(plastic_disloUCLA_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(plastic_disloUCLA_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(plastic_disloUCLA_CAtomicVolume(maxNinstance),                       source=0.0_pReal)
 allocate(plastic_disloUCLA_D0(maxNinstance),                                  source=0.0_pReal)
 allocate(plastic_disloUCLA_Qsd(maxNinstance),                                 source=0.0_pReal)
 allocate(plastic_disloUCLA_GrainSize(maxNinstance),                           source=0.0_pReal)
 allocate(plastic_disloUCLA_CEdgeDipMinDistance(maxNinstance),                 source=0.0_pReal)
 allocate(plastic_disloUCLA_SolidSolutionStrength(maxNinstance),               source=0.0_pReal)
 allocate(plastic_disloUCLA_aTolRho(maxNinstance),                             source=0.0_pReal)
 allocate(plastic_disloUCLA_dipoleFormationFactor(maxNinstance),               source=1.0_pReal) !should be on by default
 allocate(plastic_disloUCLA_rhoEdge0(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(plastic_disloUCLA_rhoEdgeDip0(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(plastic_disloUCLA_burgersPerSlipFamily(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal)
 allocate(plastic_disloUCLA_kinkheight(lattice_maxNslipFamily,maxNinstance),   source=0.0_pReal)
 allocate(plastic_disloUCLA_omega(lattice_maxNslipFamily,maxNinstance),        source=0.0_pReal)
 allocate(plastic_disloUCLA_kinkwidth(lattice_maxNslipFamily,maxNinstance),    source=0.0_pReal)
 allocate(plastic_disloUCLA_dislolength(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(plastic_disloUCLA_friction(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(plastic_disloUCLA_QedgePerSlipFamily(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal)
 allocate(plastic_disloUCLA_v0PerSlipFamily(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)
 allocate(plastic_disloUCLA_tau_peierlsPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                      source=0.0_pReal)
 allocate(plastic_disloUCLA_pPerSlipFamily(lattice_maxNslipFamily,maxNinstance),       source=0.0_pReal)
 allocate(plastic_disloUCLA_qPerSlipFamily(lattice_maxNslipFamily,maxNinstance),       source=0.0_pReal)

 allocate(plastic_disloUCLA_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                      source=0.0_pReal)

 allocate(plastic_disloUCLA_interaction_SlipSlip(lattice_maxNinteraction,maxNinstance),source=0.0_pReal)

 allocate(plastic_disloUCLA_nonSchmidCoeff(lattice_maxNnonSchmid,maxNinstance),        source=0.0_pReal)
 

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
     if (phase_plasticity(phase) == PLASTICITY_DISLOUCLA_ID) then
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_nonSchmid =    lattice_NnonSchmid(phase)
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_DISLOUCLA_ID) then           ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
      select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('edge_density')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = edge_density_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('dipole_density')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = dipole_density_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('shear_rate_slip','shearrate_slip')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = shear_rate_slip_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('accumulated_shear_slip')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = accumulated_shear_slip_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('mfp_slip')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = mfp_slip_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('resolved_stress_slip')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = resolved_stress_slip_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('threshold_stress_slip')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = threshold_stress_slip_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('edge_dipole_distance')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = edge_dipole_distance_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
           case ('stress_exponent')
             plastic_disloUCLA_Noutput(instance) = plastic_disloUCLA_Noutput(instance) + 1_pInt
             plastic_disloUCLA_outputID(plastic_disloUCLA_Noutput(instance),instance) = stress_exponent_ID
             plastic_disloUCLA_output(plastic_disloUCLA_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select
!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip system families
       case ('nslip')
         if (chunkPos(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOUCLA_label//')')
         if (chunkPos(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOUCLA_label//')')
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
             plastic_disloUCLA_Nslip(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('rhoedge0','rhoedgedip0','slipburgers','qedge','v0','clambdaslip','tau_peierls','p_slip','q_slip',&
             'kink_height','omega','kink_width','dislolength','friction_coeff')
         do j = 1_pInt, Nchunks_SlipFamilies
           tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('rhoedge0')
             plastic_disloUCLA_rhoEdge0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('rhoedgedip0')
             plastic_disloUCLA_rhoEdgeDip0(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('slipburgers')
             plastic_disloUCLA_burgersPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('qedge')
             plastic_disloUCLA_QedgePerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('v0')
             plastic_disloUCLA_v0PerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('clambdaslip')
             plastic_disloUCLA_CLambdaSlipPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('tau_peierls')
             if (lattice_structure(phase) /= LATTICE_bcc_ID) &
               call IO_warning(42_pInt,ext_msg=trim(tag)//' for non-bcc ('//PLASTICITY_DISLOUCLA_label//')')
             plastic_disloUCLA_tau_peierlsPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('p_slip')
             plastic_disloUCLA_pPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('q_slip')
             plastic_disloUCLA_qPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies) 
           case ('kink_height')
             plastic_disloUCLA_kinkheight(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)
           case ('omega')
             plastic_disloUCLA_omega(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)
           case ('kink_width')
             plastic_disloUCLA_kinkwidth(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)
           case ('dislolength')
             plastic_disloUCLA_dislolength(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('friction_coeff')
             plastic_disloUCLA_friction(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)  
         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of interactions
       case ('interaction_slipslip','interactionslipslip')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOUCLA_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           plastic_disloUCLA_interaction_SlipSlip(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('nonschmid_coefficients')
         if (chunkPos(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_DISLOUCLA_label//')')
         do j = 1_pInt,Nchunks_nonSchmid
           plastic_disloUCLA_nonSchmidCoeff(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
!--------------------------------------------------------------------------------------------------
! parameters independent of number of slip systems
       case ('grainsize')
         plastic_disloUCLA_GrainSize(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('d0')
         plastic_disloUCLA_D0(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('qsd')
         plastic_disloUCLA_Qsd(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('atol_rho')
         plastic_disloUCLA_aTolRho(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('solidsolutionstrength')
         plastic_disloUCLA_SolidSolutionStrength(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('cedgedipmindistance')
         plastic_disloUCLA_CEdgeDipMinDistance(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('catomicvolume')
         plastic_disloUCLA_CAtomicVolume(instance) = IO_floatValue(line,chunkPos,2_pInt)
       case ('dipoleformationfactor')
         plastic_disloUCLA_dipoleFormationFactor(instance) = IO_floatValue(line,chunkPos,2_pInt)
     end select
   endif; endif
 enddo parsingFile
 
 sanityChecks: do phase = 1_pInt, size(phase_plasticity)
    myPhase: if (phase_plasticity(phase) == PLASTICITY_disloUCLA_ID) then
      instance = phase_plasticityInstance(phase)
      if (sum(plastic_disloUCLA_Nslip(:,instance)) < 0_pInt) &
        call IO_error(211_pInt,el=instance,ext_msg='Nslip ('//PLASTICITY_DISLOUCLA_label//')')
      do f = 1_pInt,lattice_maxNslipFamily
        if (plastic_disloUCLA_Nslip(f,instance) > 0_pInt) then
          if (plastic_disloUCLA_rhoEdge0(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOUCLA_label//')')
          if (plastic_disloUCLA_rhoEdgeDip0(f,instance) < 0.0_pReal) & 
            call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOUCLA_label//')')
          if (plastic_disloUCLA_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOUCLA_label//')')
          if (plastic_disloUCLA_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOUCLA_label//')')
          if (plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance) < 0.0_pReal) &
            call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOUCLA_label//')')
        endif
      enddo
      if (plastic_disloUCLA_CAtomicVolume(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOUCLA_label//')')
      if (plastic_disloUCLA_D0(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOUCLA_label//')')
      if (plastic_disloUCLA_Qsd(instance) <= 0.0_pReal) &
        call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOUCLA_label//')')
     ! if (plastic_disloUCLA_aTolRho(instance) <= 0.0_pReal) &
     !   call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOUCLA_label//')')   

!--------------------------------------------------------------------------------------------------
! Determine total number of active slip systems
      plastic_disloUCLA_Nslip(:,instance) = min(lattice_NslipSystem(:,phase),plastic_disloUCLA_Nslip(:,instance))
      plastic_disloUCLA_totalNslip(instance) = sum(plastic_disloUCLA_Nslip(:,instance))
   endif myPhase
 enddo sanityChecks
 
!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 maxTotalNslip = maxval(plastic_disloUCLA_totalNslip)
 
 allocate(plastic_disloUCLA_burgersPerSlipSystem(maxTotalNslip, maxNinstance),    source=0.0_pReal)
 allocate(plastic_disloUCLA_QedgePerSlipSystem(maxTotalNslip, maxNinstance),      source=0.0_pReal)
 allocate(plastic_disloUCLA_v0PerSlipSystem(maxTotalNslip, maxNinstance),         source=0.0_pReal)
 allocate(plastic_disloUCLA_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance),source=0.0_pReal)
 
 allocate(plastic_disloUCLA_interactionMatrix_SlipSlip(maxval(plastic_disloUCLA_totalNslip),&  ! slip resistance from slip activity
                                                            maxval(plastic_disloUCLA_totalNslip),&
                                                            maxNinstance), source=0.0_pReal)
 allocate(plastic_disloUCLA_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)
 
 allocate(state(maxNinstance))
 allocate(state0(maxNinstance))
 allocate(dotState(maxNinstance))

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   myPhase2: if (phase_plasticity(phase) == PLASTICITY_disloUCLA_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_plasticityInstance(phase)
     ns = plastic_disloUCLA_totalNslip(instance)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputs: do o = 1_pInt,plastic_disloUCLA_Noutput(instance)
       select case(plastic_disloUCLA_outputID(o,instance))
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
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          plastic_disloUCLA_sizePostResult(o,instance) = mySize
          plastic_disloUCLA_sizePostResults(instance)  = plastic_disloUCLA_sizePostResults(instance) + mySize
       endif
     enddo outputs

!--------------------------------------------------------------------------------------------------
! allocate state arrays

     sizeDotState     = int(size(['rhoEdge     ','rhoEdgeDip  ','accshearslip']),pInt) * ns
     sizeDeltaState   =  0_pInt
     sizeState        = sizeDotState &
                      + int(size(['invLambdaSlip     ',&
                                  'meanFreePathSlip  ','tauSlipThreshold  ']),pInt) * ns

     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_disloUCLA_sizePostResults(instance)
     plasticState(phase)%nSlip = plastic_disloucla_totalNslip(instance)
     plasticState(phase)%nTwin = 0_pInt
     plasticState(phase)%nTrans= 0_pInt
     allocate(plasticState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(plasticState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(plasticState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(plasticState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(plasticState(phase)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)
     offset_slip = 2_pInt*plasticState(phase)%nSlip
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NofMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state   (offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NofMyPhase)
    !* Process slip related parameters ------------------------------------------------ 
 
     mySlipFamilies: do f = 1_pInt,lattice_maxNslipFamily
       index_myFamily = sum(plastic_disloUCLA_Nslip(1:f-1_pInt,instance))                      ! index in truncated slip system list
       mySlipSystems: do j = 1_pInt,plastic_disloUCLA_Nslip(f,instance)

      !* Burgers vector, 
      !  dislocation velocity prefactor,
      !  mean free path prefactor,
      !  and minimum dipole distance
 
         plastic_disloUCLA_burgersPerSlipSystem(index_myFamily+j,instance) = &
         plastic_disloUCLA_burgersPerSlipFamily(f,instance)
 
         plastic_disloUCLA_QedgePerSlipSystem(index_myFamily+j,instance) = &
         plastic_disloUCLA_QedgePerSlipFamily(f,instance)
 
         plastic_disloUCLA_v0PerSlipSystem(index_myFamily+j,instance) = &
         plastic_disloUCLA_v0PerSlipFamily(f,instance)
 
         plastic_disloUCLA_CLambdaSlipPerSlipSystem(index_myFamily+j,instance) = &
         plastic_disloUCLA_CLambdaSlipPerSlipFamily(f,instance)
  
       !* Calculation of forest projections for edge dislocations
       !* Interaction matrices
         otherSlipFamilies: do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_disloUCLA_Nslip(1:o-1_pInt,instance))
           otherSlipSystems: do k = 1_pInt,plastic_disloUCLA_Nslip(o,instance)
             plastic_disloUCLA_forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,phase))+j,phase), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,phase))+k,phase)))
             plastic_disloUCLA_interactionMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k,instance) = &
                   plastic_disloUCLA_interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                 phase), instance )
         enddo otherSlipSystems; enddo otherSlipFamilies
 
       enddo mySlipSystems
     enddo mySlipFamilies
  
     startIndex=1_pInt
     endIndex=ns
     state(instance)%rhoEdge=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%rhoEdge=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%rhoEdge=>plasticState(phase)%dotState(startIndex:endIndex,:) 

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     state(instance)%rhoEdgeDip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%rhoEdgeDip=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%rhoEdgeDip=>plasticState(phase)%dotState(startIndex:endIndex,:)

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     state(instance)%accshear_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%accshear_slip=>plasticState(phase)%state0(startIndex:endIndex,:)
     dotState(instance)%accshear_slip=>plasticState(phase)%dotState(startIndex:endIndex,:)

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     state(instance)%invLambdaSlip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%invLambdaSlip=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     state(instance)%mfp_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%mfp_slip=>plasticState(phase)%state0(startIndex:endIndex,:)

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     state(instance)%threshold_stress_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     state0(instance)%threshold_stress_slip=>plasticState(phase)%state0(startIndex:endIndex,:)

     call plastic_disloUCLA_stateInit(phase,instance)
     call plastic_disloUCLA_aTolState(phase,instance)
   endif myPhase2
 
 enddo initializeInstances
 
end subroutine plastic_disloUCLA_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_stateInit(ph,instance)
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

 integer(pInt) :: i,f,ns, index_myFamily
 real(pReal), dimension(plastic_disloUCLA_totalNslip(instance)) :: &
   rhoEdge0, &
   rhoEdgeDip0, &
   invLambdaSlip0, &
   MeanFreePathSlip0, &
   tauSlipThreshold0
 tempState = 0.0_pReal
 ns = plastic_disloUCLA_totalNslip(instance)

!--------------------------------------------------------------------------------------------------
! initialize basic slip state variables
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily   = sum(plastic_disloUCLA_Nslip(1:f-1_pInt,instance))                        ! index in truncated slip system list
   rhoEdge0(index_myFamily+1_pInt: &
            index_myFamily+plastic_disloUCLA_Nslip(f,instance)) = &
     plastic_disloUCLA_rhoEdge0(f,instance)
   rhoEdgeDip0(index_myFamily+1_pInt: &
               index_myFamily+plastic_disloUCLA_Nslip(f,instance)) = &
     plastic_disloUCLA_rhoEdgeDip0(f,instance)
 enddo
 
 tempState(1_pInt:ns)           = rhoEdge0
 tempState(ns+1_pInt:2_pInt*ns) = rhoEdgeDip0
 
!--------------------------------------------------------------------------------------------------
! initialize dependent slip microstructural variables
 forall (i = 1_pInt:ns) &
   invLambdaSlip0(i) = sqrt(dot_product((rhoEdge0+rhoEdgeDip0),plastic_disloUCLA_forestProjectionEdge(1:ns,i,instance)))/ &
                       plastic_disloUCLA_CLambdaSlipPerSlipSystem(i,instance)
 tempState(3_pInt*ns+1:4_pInt*ns) = invLambdaSlip0
 
 forall (i = 1_pInt:ns) &
   MeanFreePathSlip0(i) = &
     plastic_disloUCLA_GrainSize(instance)/(1.0_pReal+invLambdaSlip0(i)*plastic_disloUCLA_GrainSize(instance))
 tempState(4_pInt*ns+1:5_pInt*ns) = MeanFreePathSlip0
 
 forall (i = 1_pInt:ns) &
   tauSlipThreshold0(i) = &
     lattice_mu(ph)*plastic_disloUCLA_burgersPerSlipSystem(i,instance) * &
     sqrt(dot_product((rhoEdge0+rhoEdgeDip0),plastic_disloUCLA_interactionMatrix_SlipSlip(i,1:ns,instance)))

 tempState(5_pInt*ns+1:6_pInt*ns) = tauSlipThreshold0
 
plasticState(ph)%state0 = spread(tempState,2,size(plasticState(ph)%state(1,:)))

end subroutine plastic_disloUCLA_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this plasticity
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_aTolState(ph,instance)
 use material, only: &
  plasticState

 implicit none
 integer(pInt), intent(in) ::  &
   ph, &
   instance                                                                                         ! number specifying the current instance of the plasticity

 ! Tolerance state for dislocation densities
 plasticState(ph)%aTolState(1_pInt:2_pInt*plastic_disloUCLA_totalNslip(instance)) = &
   plastic_disloUCLA_aTolRho(instance)

 ! Tolerance state for accumulated shear due to slip 
 plasticState(ph)%aTolState(2_pInt*plastic_disloUCLA_totalNslip(instance)+1_pInt: &
                            3_pInt*plastic_disloUCLA_totalNslip(instance))=1e6_pReal

end subroutine plastic_disloUCLA_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_microstructure(temperature,ipc,ip,el)
 use math, only: &
   pi
 use material, only: &
   phase_plasticityInstance, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_mu

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),   intent(in) :: &
   temperature                                                                                      !< temperature at IP 

 integer(pInt) :: &
   instance, &
   ns,s, &
   ph, &
   of

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance)

 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:ns) &
   state(instance)%invLambdaSlip(s,of) = &
     sqrt(dot_product((state(instance)%rhoEdge(1_pInt:ns,of)+state(instance)%rhoEdgeDip(1_pInt:ns,of)),&
                      plastic_disloUCLA_forestProjectionEdge(1:ns,s,instance)))/ &
     plastic_disloUCLA_CLambdaSlipPerSlipSystem(s,instance)
 
 !* mean free path between 2 obstacles seen by a moving dislocation
 do s = 1_pInt,ns
   state(instance)%mfp_slip(s,of) = &
      plastic_disloUCLA_GrainSize(instance)/&
       (1.0_pReal+plastic_disloUCLA_GrainSize(instance)*(state(instance)%invLambdaSlip(s,of)))
 enddo
 
 !* threshold stress for dislocation motion
 forall (s = 1_pInt:ns) &
   state(instance)%threshold_stress_slip(s,of) = &
     lattice_mu(ph)*plastic_disloUCLA_burgersPerSlipSystem(s,instance)*&
     sqrt(dot_product((state(instance)%rhoEdge(1_pInt:ns,of)+state(instance)%rhoEdgeDip(1_pInt:ns,of)),&
                      plastic_disloUCLA_interactionMatrix_SlipSlip(s,1:ns,instance)))
 
end subroutine plastic_disloUCLA_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_Mandel33to6, &
   math_symmetric33, &
   math_mul33x3
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_maxNslipFamily,&
   lattice_NslipSystem, &
   lattice_NnonSchmid
 
 implicit none
 integer(pInt), intent(in)                  :: ipc,ip,el
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(6),   intent(in)    :: Tstar_v
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(9,9), intent(out)   :: dLp_dTstar99

 integer(pInt) :: instance,ph,of,ns,f,i,j,k,l,m,n,index_myFamily
 real(pReal) :: StressRatio_p,StressRatio_pminus1,BoltzmannRatio,DotGamma0, &
                tau_slip_pos,tau_slip_neg,vel_slip,dvel_slip,&
                dgdot_dtauslip_pos,dgdot_dtauslip_neg,stressRatio
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,gdot_slip_neg
   
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance)
 
 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 
!--------------------------------------------------------------------------------------------------
! Dislocation glide part
 gdot_slip_pos = 0.0_pReal
 gdot_slip_neg = 0.0_pReal
 dgdot_dtauslip_pos = 0.0_pReal
 dgdot_dtauslip_neg = 0.0_pReal

 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
     j = j+1_pInt
     !* Boltzmann ratio
     BoltzmannRatio = plastic_disloUCLA_QedgePerSlipSystem(j,instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = &
        state(instance)%rhoEdge(j,of)*plastic_disloUCLA_burgersPerSlipSystem(j,instance)*&
        plastic_disloUCLA_v0PerSlipSystem(j,instance)
     !* Resolved shear stress on slip system
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + plastic_disloUCLA_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k,  index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + plastic_disloUCLA_nonSchmidCoeff(k,instance)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo nonSchmidSystems

     significantPostitiveStress: if((abs(tau_slip_pos)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
       !* Stress ratio
       stressRatio = ((abs(tau_slip_pos)-state(instance)%threshold_stress_slip(j,of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p       = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)
       stressRatio_pminus1 = stressRatio**(plastic_disloUCLA_pPerSlipFamily(f,instance)-1.0_pReal)
       !* Shear rates due to slip                                                                                                                                              
       vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_pos  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
 
       gdot_slip_pos(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos) 

       !* Derivatives of shear rates                                                                                                              
       dvel_slip = &
          2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
          * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
          * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
          * ( &
             (exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) &
          + tau_slip_pos &
          * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&    !deltaf(i)                                                  
          *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
          *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
          (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
          StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  ) &!deltaf(f)                                                 
             ) &
          *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos &
          +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
          *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
          * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
             ) &
          -  (tau_slip_pos &
          * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) )  &
          *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal) &
          +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
          *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
          * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&     !deltaf(i)                                                 
          *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
          *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
          (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
          StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  )& !deltaf(f)
             ) &
            )  &
          / (  &
             ( &
           2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos &
           + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
           *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
           * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
             )**2.0_pReal &
            )
       dgdot_dtauslip_pos = DotGamma0 * dvel_slip

     endif significantPostitiveStress
     significantNegativeStress: if((abs(tau_slip_neg)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
       !* Stress ratio
       stressRatio = ((abs(tau_slip_neg)-state(instance)%threshold_stress_slip(j,of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p       = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)
       stressRatio_pminus1 = stressRatio**(plastic_disloUCLA_pPerSlipFamily(f,instance)-1.0_pReal)
       !* Shear rates due to slip                                                                                                                                              
       vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_neg  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
                
       gdot_slip_neg(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg) 

       !* Derivatives of shear rates                                                                                                              
       dvel_slip = &
          2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
          * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
          * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
          * ( &
             (exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) &
          + tau_slip_neg &
          * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&    !deltaf(i)                                                  
          *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
          *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
          (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
          StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  ) &!deltaf(f)                                                 
             ) &
          *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg &
          +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
          *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
          * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
             ) &
          -  (tau_slip_neg &
          * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) )  &
          *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal) &
          +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
          *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
          * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&     !deltaf(i)                                                 
          *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
          *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
          (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
          StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  )& !deltaf(f)                                                 
             ) &
            )  &
          / (  &
             ( &
           2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg &
           + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
           *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
           * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
             )**2.0_pReal &
            )

       dgdot_dtauslip_neg = DotGamma0 * dvel_slip

     endif significantNegativeStress
     !* Plastic velocity gradient for dislocation glide
     Lp = Lp + (gdot_slip_pos(j)+gdot_slip_neg(j))*0.5_pReal*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     !* Calculation of the tangent of Lp
     forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
        dLp_dTstar3333(k,l,m,n) = &
        dLp_dTstar3333(k,l,m,n) + (dgdot_dtauslip_pos*nonSchmid_tensor(m,n,1)+&
                                   dgdot_dtauslip_neg*nonSchmid_tensor(m,n,2))*0.5_pReal*&
                                   lattice_Sslip(k,l,1,index_myFamily+i,ph)
   enddo slipSystems
 enddo slipFamilies

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)
 
end subroutine plastic_disloUCLA_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dotState(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check, &
   dEq0
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   plasticState, &
   phaseAt, phasememberAt
 use lattice,  only: &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid, &
   lattice_mu

 implicit none
 real(pReal), dimension(6),  intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: instance,ns,f,i,j,k,index_myFamily, &
                  ph, &
                  of
 real(pReal) :: &
   stressRatio_p,&
   BoltzmannRatio,&
   DotGamma0,&
   stressRatio, &
   EdgeDipMinDistance,&
   AtomicVolume,&
   VacancyDiffusion,&
   tau_slip_pos,&
   tau_slip_neg,&
   DotRhoMultiplication,&
   EdgeDipDistance, &
   DotRhoEdgeDipAnnihilation, &
   DotRhoEdgeEdgeAnnihilation, &
   ClimbVelocity, &
   DotRhoEdgeDipClimb, &
   DotRhoDipFormation, &
   vel_slip, &
   gdot_slip
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos, gdot_slip_neg

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance) 

 plasticState(ph)%dotState(:,of) = 0.0_pReal
 
 !* Dislocation density evolution
 gdot_slip_pos = 0.0_pReal
 gdot_slip_neg = 0.0_pReal
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
     j = j+1_pInt
     !* Boltzmann ratio
     BoltzmannRatio = plastic_disloUCLA_QedgePerSlipSystem(j,instance)/(kB*Temperature)
     !* Initial shear rates
     DotGamma0 = &
        state(instance)%rhoEdge(j,of)*plastic_disloUCLA_burgersPerSlipSystem(j,instance)*&
        plastic_disloUCLA_v0PerSlipSystem(j,instance)
     !* Resolved shear stress on slip system
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos

     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
       tau_slip_pos = tau_slip_pos + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems

     significantPositiveStress: if((abs(tau_slip_pos)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_pos)-state(instance)%threshold_stress_slip(j,of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)
       !* Shear rates due to slip                                                                                                                                                                                                                                                                           
       vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_pos  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
                       
       gdot_slip_pos(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos) 
     endif significantPositiveStress
     significantNegativeStress: if((abs(tau_slip_neg)-state(instance)%threshold_stress_slip(j,of)) > tol_math_check) then
       !* Stress ratios
       stressRatio = ((abs(tau_slip_neg)-state(instance)%threshold_stress_slip(j,of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
       stressRatio_p = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)

       vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_neg  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
       
       gdot_slip_neg(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg) 
     endif significantNegativeStress
     gdot_slip = (gdot_slip_pos(j)+gdot_slip_neg(j))*0.5_pReal
     !* Multiplication
     DotRhoMultiplication = abs(gdot_slip)/&
                               (plastic_disloUCLA_burgersPerSlipSystem(j,instance)* &
                                state(instance)%mfp_slip(j,of))
 
     !* Dipole formation
     EdgeDipMinDistance = &
       plastic_disloUCLA_CEdgeDipMinDistance(instance)*plastic_disloUCLA_burgersPerSlipSystem(j,instance)
     if (dEq0(tau_slip_pos)) then
       DotRhoDipFormation = 0.0_pReal
     else
       EdgeDipDistance = &
         (3.0_pReal*lattice_mu(ph)*plastic_disloUCLA_burgersPerSlipSystem(j,instance))/&
         (16.0_pReal*pi*abs(tau_slip_pos))
       if (EdgeDipDistance>state(instance)%mfp_slip(j,of)) EdgeDipDistance=state(instance)%mfp_slip(j,of)
       if (EdgeDipDistance<EdgeDipMinDistance) EdgeDipDistance=EdgeDipMinDistance
       DotRhoDipFormation = &
         ((2.0_pReal*EdgeDipDistance)/plastic_disloUCLA_burgersPerSlipSystem(j,instance))*&
         state(instance)%rhoEdge(j,of)*abs(gdot_slip)*plastic_disloUCLA_dipoleFormationFactor(instance)
     endif
 
    !* Spontaneous annihilation of 2 single edge dislocations
    DotRhoEdgeEdgeAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/plastic_disloUCLA_burgersPerSlipSystem(j,instance))*&
        state(instance)%rhoEdge(j,of)*abs(gdot_slip)
 
    !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
    DotRhoEdgeDipAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/plastic_disloUCLA_burgersPerSlipSystem(j,instance))*&
        state(instance)%rhoEdgeDip(j,of)*abs(gdot_slip)
 
      !* Dislocation dipole climb
     AtomicVolume = &
        plastic_disloUCLA_CAtomicVolume(instance)*plastic_disloUCLA_burgersPerSlipSystem(j,instance)**(3.0_pReal)
     VacancyDiffusion = &
        plastic_disloUCLA_D0(instance)*exp(-plastic_disloUCLA_Qsd(instance)/(kB*Temperature))
     if (dEq0(tau_slip_pos)) then
       DotRhoEdgeDipClimb = 0.0_pReal
     else
       ClimbVelocity = &
          ((3.0_pReal*lattice_mu(ph)*VacancyDiffusion*AtomicVolume)/(2.0_pReal*pi*kB*Temperature))*&
          (1/(EdgeDipDistance+EdgeDipMinDistance))
       DotRhoEdgeDipClimb = &
          (4.0_pReal*ClimbVelocity*state(instance)%rhoEdgeDip(j,of))/(EdgeDipDistance-EdgeDipMinDistance)
     endif
 
     !* Edge dislocation density rate of change
     dotState(instance)%rhoEdge(j,of) = &
        DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
 
     !* Edge dislocation dipole density rate of change
     dotState(instance)%rhoEdgeDip(j,of) = &
        DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb
 
     !* Dotstate for accumulated shear due to slip
     dotState(instance)%accshear_slip(j,of) = gdot_slip
 
   enddo slipSystems
 enddo slipFamilies

 
end subroutine plastic_disloUCLA_dotState

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_disloUCLA_postResults(Tstar_v,Temperature,ipc,ip,el)
 use prec, only: &
   tol_math_check, &
   dEq, dNeq0
 use math, only: &
   pi
 use material, only: &
   material_phase, &
   phase_plasticityInstance,& 
   !plasticState, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid, &
   lattice_mu

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(plastic_disloUCLA_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           plastic_disloUCLA_postResults

 integer(pInt) :: &
   instance,&
   ns,&
   f,o,i,c,j,k,index_myFamily,&
   ph, &
   of
 real(pReal) :: StressRatio_p,StressRatio_pminus1,&
                BoltzmannRatio,DotGamma0,stressRatio,&
                dvel_slip, vel_slip
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg
 
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance)
 
 !* Required output
 c = 0_pInt
 plastic_disloUCLA_postResults = 0.0_pReal

 do o = 1_pInt,plastic_disloUCLA_Noutput(instance)
    select case(plastic_disloUCLA_outputID(o,instance))
 
      case (edge_density_ID)
        plastic_disloUCLA_postResults(c+1_pInt:c+ns) = state(instance)%rhoEdge(1_pInt:ns,of)
        c = c + ns
      case (dipole_density_ID)
        plastic_disloUCLA_postResults(c+1_pInt:c+ns) = state(instance)%rhoEdgeDip(1_pInt:ns,of)
        c = c + ns
      case (shear_rate_slip_ID,stress_exponent_ID)
        gdot_slip_pos = 0.0_pReal
        gdot_slip_neg = 0.0_pReal
        dgdot_dtauslip_pos = 0.0_pReal
        dgdot_dtauslip_neg = 0.0_pReal
        j = 0_pInt
        slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
          index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
          slipSystems: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
            j = j + 1_pInt
            !* Boltzmann ratio
            BoltzmannRatio = plastic_disloUCLA_QedgePerSlipSystem(j,instance)/(kB*Temperature)
            !* Initial shear rates
            DotGamma0 = &
              state(instance)%rhoEdge(j,of)*plastic_disloUCLA_burgersPerSlipSystem(j,instance)*&
              plastic_disloUCLA_v0PerSlipSystem(j,instance)
            !* Resolved shear stress on slip system
            tau_slip_pos(j) = dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
            tau_slip_neg(j)  = tau_slip_pos(j)

            nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph) 
              tau_slip_pos = tau_slip_pos + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
              tau_slip_neg = tau_slip_neg + plastic_disloUCLA_nonSchmidCoeff(k,instance)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
            enddo nonSchmidSystems

            significantPositiveTau: if((abs(tau_slip_pos(j))-state(instance)%threshold_stress_slip(j, of)) > tol_math_check) then
              !* Stress ratio
              stressRatio = ((abs(tau_slip_pos(j))-state(instance)%threshold_stress_slip(j, of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
              stressRatio_p       = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)
              stressRatio_pminus1 = stressRatio**(plastic_disloUCLA_pPerSlipFamily(f,instance)-1.0_pReal)
              !* Shear rates due to slip
              vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_pos(j)  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos(j) &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
                       
              gdot_slip_pos(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos(j))
              !* Derivatives of shear rates 

              dvel_slip = &
                   2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                   * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                   * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                   * ( &
                   (exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) &
                   + tau_slip_pos(j) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&    !deltaf(i)                                         
                   *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
                   *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
                   (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  ) &!deltaf(f)                                        
                   ) &
                   *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos(j) &
                   +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                   ) &
                   -  (tau_slip_pos(j) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) )  &
                   *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal) &
                   +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&     !deltaf(i)                                        
                   *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
                   *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
                   (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  )& !deltaf(f)                                        
                   ) &
                   )  &
                   / (  &
                   ( &
                   2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_pos(j) &
                   + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                   )**2.0_pReal &
                   )

              dgdot_dtauslip_pos(j) = DotGamma0 * dvel_slip

            endif significantPositiveTau
            significantNegativeTau: if((abs(tau_slip_neg(j))-state(instance)%threshold_stress_slip(j, of)) > tol_math_check) then
              !* Stress ratios
              stressRatio = ((abs(tau_slip_neg(j))-state(instance)%threshold_stress_slip(j, of))/&
                      (plastic_disloUCLA_SolidSolutionStrength(instance)+&
                       plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance)))
              stressRatio_p       = stressRatio** plastic_disloUCLA_pPerSlipFamily(f,instance)
              stressRatio_pminus1 = stressRatio**(plastic_disloUCLA_pPerSlipFamily(f,instance)-1.0_pReal)
              !* Shear rates due to slip                                                                                                                                                                                                                                                                           
              vel_slip = 2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                     * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                     * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                     * (tau_slip_neg(j)  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) ) &
                     / ( &
                     2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg(j) &
                     + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                     *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                     )
              
              gdot_slip_neg(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg(j))
              !* Derivatives of shear rates 
              dvel_slip = &
                   2.0_pReal*plastic_disloUCLA_burgersPerSlipFamily(f,instance) &
                   * plastic_disloUCLA_kinkheight(f,instance) * plastic_disloUCLA_omega(f,instance)  &
                   * ( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) ) &
                   * ( &
                   (exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) &
                   + tau_slip_neg(j) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&    !deltaf(i)                                         
                   *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
                   *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
                   (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  ) &!deltaf(f)                                        
                   ) &
                   *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg(j) &
                   +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                   ) &
                   -  (tau_slip_neg(j) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)) )  &
                   *  (2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal) &
                   +  plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance)))&     !deltaf(i)                                        
                   *BoltzmannRatio*plastic_disloUCLA_pPerSlipFamily(f,instance)&
                   *plastic_disloUCLA_qPerSlipFamily(f,instance)/&
                   (plastic_disloUCLA_SolidSolutionStrength(instance)+plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(plastic_disloUCLA_qPerSlipFamily(f,instance)-1.0_pReal)  )& !deltaf(f)                                        
                   ) &
                   )  &
                   / (  &
                   ( &
                   2.0_pReal*(plastic_disloUCLA_burgersPerSlipFamily(f,instance)**2.0_pReal)*tau_slip_neg(j) &
                   + plastic_disloUCLA_omega(f,instance) * plastic_disloUCLA_friction(f,instance) &
                   *(( state(instance)%mfp_slip(j,of) - plastic_disloUCLA_kinkwidth(f,instance) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** plastic_disloUCLA_qPerSlipFamily(f,instance))  &
                   )**2.0_pReal &
                   )


              dgdot_dtauslip_neg(j) = DotGamma0 * dvel_slip

            endif significantNegativeTau
          enddo slipSystems
        enddo slipFamilies

        if     (plastic_disloUCLA_outputID(o,instance) == shear_rate_slip_ID) then
          plastic_disloUCLA_postResults(c+1:c+ns) = (gdot_slip_pos + gdot_slip_neg)*0.5_pReal
          c = c + ns
        elseif(plastic_disloUCLA_outputID(o,instance) == stress_exponent_ID) then
          do j = 1_pInt, ns
            if (dEq(gdot_slip_pos(j)+gdot_slip_neg(j),0.0_pReal)) then
              plastic_disloUCLA_postResults(c+j) = 0.0_pReal
            else
              plastic_disloUCLA_postResults(c+j) = (tau_slip_pos(j)+tau_slip_neg(j))/&
                                                       (gdot_slip_pos(j)+gdot_slip_neg(j))*&
                                                       (dgdot_dtauslip_pos(j)+dgdot_dtauslip_neg(j))* 0.5_pReal
            endif
          enddo
           c = c + ns
        endif

      case (accumulated_shear_slip_ID)
       plastic_disloUCLA_postResults(c+1_pInt:c+ns) = &
                      state(instance)%accshear_slip(1_pInt:ns, of)
        c = c + ns
      case (mfp_slip_ID)
        plastic_disloUCLA_postResults(c+1_pInt:c+ns) =&
                      state(instance)%mfp_slip(1_pInt:ns, of)
        c = c + ns
      case (resolved_stress_slip_ID)
        j = 0_pInt
        slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems1: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
              j = j + 1_pInt
              plastic_disloUCLA_postResults(c+j) =&
                                dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))
        enddo slipSystems1; enddo slipFamilies1
        c = c + ns
      case (threshold_stress_slip_ID)
        plastic_disloUCLA_postResults(c+1_pInt:c+ns) = &
                                state(instance)%threshold_stress_slip(1_pInt:ns,of)
        c = c + ns
      case (edge_dipole_distance_ID)
        j = 0_pInt
        slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems2: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
              j = j + 1_pInt
              if (dNeq0(abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))))) then
              plastic_disloUCLA_postResults(c+j) = &
                (3.0_pReal*lattice_mu(ph)*plastic_disloUCLA_burgersPerSlipSystem(j,instance))/&
                (16.0_pReal*pi*abs(dot_product(Tstar_v,lattice_Sslip_v(:,1,index_myFamily+i,ph))))
              else
              plastic_disloUCLA_postResults(c+j) = huge(1.0_pReal)
              endif
              plastic_disloUCLA_postResults(c+j)=min(plastic_disloUCLA_postResults(c+j),&
                                                            state(instance)%mfp_slip(j,of))
        enddo slipSystems2; enddo slipFamilies2
        c = c + ns
    end select
 enddo
end function plastic_disloUCLA_postResults

end module plastic_disloUCLA
