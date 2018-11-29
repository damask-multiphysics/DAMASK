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
   plastic_disloUCLA_CEdgeDipMinDistance, &                                                         !<
   plastic_disloUCLA_dipoleFormationFactor                                                       !< scaling factor for dipole formation: 0: off, 1: on. other values not useful

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   plastic_disloUCLA_CLambdaSlipPerSlipFamily, &                                                    !< Adj. parameter for distance between 2 forest dislocations for each slip family and instance
   plastic_disloUCLA_CLambdaSlipPerSlipSystem, &                                                    !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
   !* mobility law parameters                                                                                                                                                                                                                                                          
   plastic_disloUCLA_friction                                                                    !< friction coeff. B (kMC)

 real(pReal),                         dimension(:,:,:),       allocatable,         private :: &
   plastic_disloUCLA_forestProjectionEdge                                                           !< matrix of forest projections of edge dislocations for each instance

 enum, bind(c) 
   enumerator :: undefined_ID, &
     rho_ID, &
     rhoDip_ID, &
     shearrate_ID, &
     accumulatedshear_ID, &
     mfp_ID, &
     resolvedstress_ID, &
     thresholdstress_ID, &
     dipoledistance_ID, &
     stressexponent_ID
 end enum

 type, private :: tParameters
   real(pReal) :: &
     aTolRho, &
     grainSize, &
SolidSolutionStrength  !< Strength due to elements in solid solution
   real(pReal),                 allocatable, dimension(:) :: &
     rho0, &                                                                        !< initial edge dislocation density per slip system for each family and instance
     rhoDip0, &                                                                     !< initial edge dipole density per slip system for each family and instance
     burgers, &                                                                     !< absolute length of burgers vector [m] for each slip system and instance
     H0kp, &                                                                        !< activation energy for glide [J] for each slip system and instance
     v0, &                                                                          !< dislocation velocity prefactor [m/s] for each family and instance
     CLambda, &                                                                     !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
     p, &                                                                           !< p-exponent in glide velocity
     q, &                                                                           !< q-exponent in glide velocity       
     !* mobility law parameters                                                                                                                                                      
     kink_height, &                                                                  !< height of the kink pair                                                                                         
     kink_width, &                                                                   !< width of the kink pair   
     omega, &                                                                   !<       attempt frequency for kink pair nucleation                                        
     viscosity, &                                                                   !< friction coeff. B (kMC)
     !*    
     tau_Peierls, &
     nonSchmidCoeff
   real(pReal),                 allocatable, dimension(:,:) :: &
     interaction_SlipSlip                                                                           !< slip resistance from slip activity
   real(pReal),                 allocatable, dimension(:,:,:) :: &
     Schmid_slip, &
     Schmid_twin, &
     nonSchmid_pos, &
     nonSchmid_neg
   integer(pInt) :: &
     totalNslip                                                                                     !< total number of active slip system
   integer(pInt),               allocatable, dimension(:) :: &
     Nslip                                                                                          !< number of active slip systems for each family
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
 end type                                                                                           !< container type for internal constitutive parameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   plastic_disloUCLA_outputID                                                                       !< ID of each post result output
 
 type, private :: tDisloUCLAState 
     real(pReal), pointer, dimension(:,:) :: &
       rhoEdge, &
       rhoEdgeDip, &
       accshear_slip
 end type 

 type, private :: tDisloUCLAMicrostructure
   real(pReal), allocatable,     dimension(:,:) :: &
     mfp, &
     threshold_stress
 end type tDisloUCLAMicrostructure

 type(tDisloUCLAState ), allocatable, dimension(:), private :: &
   state, &
   dotState

 type(tDisloUCLAMicrostructure), allocatable, dimension(:), private :: &
   microstructure
 
 public :: &
   plastic_disloUCLA_init, &
   plastic_disloUCLA_microstructure, &
   plastic_disloUCLA_LpAndItsTangent, &
   plastic_disloUCLA_dotState, &
   plastic_disloUCLA_postResults
 private :: &
   kinetics


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
   math_mul3x3, &
   math_expand
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
material_allocatePlasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice
 
 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,maxTotalNslip,&
                  f,instance,j,k,o,ns, i, &
                  Nchunks_SlipSlip = 0_pInt, outputSize, &
                  Nchunks_SlipFamilies = 0_pInt,Nchunks_nonSchmid = 0_pInt, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex, p
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase
 character(len=65536) :: &
   structure = '',&
   tag  = '', &
   line = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip
 character(len=65536), dimension(:), allocatable :: outputs
 integer(kind(undefined_ID))  :: outputID
 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
  
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
 allocate(plastic_disloUCLA_CEdgeDipMinDistance(maxNinstance),                 source=0.0_pReal)
 allocate(plastic_disloUCLA_dipoleFormationFactor(maxNinstance),               source=1.0_pReal) !should be on by default
 allocate(plastic_disloUCLA_friction(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal)

 allocate(plastic_disloUCLA_CLambdaSlipPerSlipFamily(lattice_maxNslipFamily,maxNinstance), &
                                                                                      source=0.0_pReal)

 
 allocate(param(maxNinstance))
 allocate(state(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(microstructure(maxNinstance))


do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_DISLOUCLA_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)))

   structure          = config_phase(p)%getString('lattice_structure')

   prm%aTolRho = config_phase(p)%getFloat('atol_rho')
!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%Nslip      = config_phase(p)%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%Nslip)
   slipActive: if (prm%totalNslip > 0_pInt) then
     prm%Schmid_slip          = lattice_SchmidMatrix_slip(prm%Nslip,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))
     if(structure=='bcc') then
       prm%nonSchmidCoeff     = config_phase(p)%getFloats('nonschmid_coefficients',&
                                                          defaultVal = emptyRealArray)
       prm%nonSchmid_pos      = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1_pInt)
       prm%nonSchmid_neg      = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1_pInt)
     else
       prm%nonSchmid_pos      = prm%Schmid_slip
       prm%nonSchmid_neg      = prm%Schmid_slip
     endif
     prm%interaction_SlipSlip = lattice_interaction_SlipSlip(prm%Nslip, &
                                                             config_phase(p)%getFloats('interaction_slipslip'), &
                                                             structure(1:3))
     prm%rho0        = config_phase(p)%getFloats('rhoedge0')
     prm%rhoDip0     = config_phase(p)%getFloats('rhoedgedip0')
     prm%burgers     = config_phase(p)%getFloats('slipburgers')
     prm%H0kp        = config_phase(p)%getFloats('qedge')
     prm%v0          = config_phase(p)%getFloats('v0')
     !prm%clambda     = config_phase(p)%getFloats('clambda')
     prm%tau_Peierls  = config_phase(p)%getFloats('tau_peierls')
     prm%p           = config_phase(p)%getFloats('p_slip',defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%q           = config_phase(p)%getFloats('q_slip',defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%kink_height  = config_phase(p)%getFloats('kink_height')
     prm%kink_width   = config_phase(p)%getFloats('kink_width')
     prm%omega   = config_phase(p)%getFloats('omega')
     !prm%viscosity   = config_phase(p)%getFloats('viscosity')


     prm%SolidSolutionStrength  = config_phase(p)%getFloat('solidsolutionstrength')

     prm%grainSize  = config_phase(p)%getFloat('grainsize')

     plastic_disloUCLA_D0(phase_plasticityInstance(p)) = config_phase(p)%getFloat('qsd')
     plastic_disloUCLA_Qsd(phase_plasticityInstance(p)) = config_phase(p)%getFloat('qsd')
     plastic_disloUCLA_CEdgeDipMinDistance(phase_plasticityInstance(p)) = config_phase(p)%getFloat('cedgedipmindistance')
     plastic_disloUCLA_CAtomicVolume(phase_plasticityInstance(p)) = config_phase(p)%getFloat('catomicvolume')
     plastic_disloUCLA_dipoleFormationFactor(phase_plasticityInstance(p)) = config_phase(p)%getFloat('dipoleformationfactor')

     ! expand: family => system
     prm%rho0   = math_expand(prm%rho0,  prm%Nslip)
     prm%rhoDip0   = math_expand(prm%rhoDip0,  prm%Nslip)
     prm%q   = math_expand(prm%q,  prm%Nslip)
     prm%p   = math_expand(prm%p,  prm%Nslip)
     prm%H0kp   = math_expand(prm%H0kp,  prm%Nslip)
     prm%burgers   = math_expand(prm%burgers,  prm%Nslip)
     prm%kink_height   = math_expand(prm%kink_height,  prm%Nslip)
     prm%kink_width   = math_expand(prm%kink_width,  prm%Nslip)
     prm%omega   = math_expand(prm%omega,  prm%Nslip)
     prm%tau_Peierls   = math_expand(prm%tau_Peierls,  prm%Nslip)
     prm%v0   = math_expand(prm%v0,  prm%Nslip)
   endif slipActive


!--------------------------------------------------------------------------------------------------
!  phase outputs

#if defined(__GFORTRAN__)
  outputs = ['GfortranBug86277']
  outputs = config_phase(p)%getStrings('(output)',defaultVal=outputs)
  if (outputs(1) == 'GfortranBug86277') outputs = emptyStringArray
#else
  outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
#endif
   allocate(prm%outputID(0))
   
   do i = 1_pInt, size(outputs)
     outputID = undefined_ID
     outputSize = prm%totalNslip
     select case(trim(outputs(i)))
         case ('edge_density')
           outputID  = merge(rho_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('dipole_density')
           outputID = merge(rhoDip_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('shear_rate','shearrate','shear_rate_slip','shearrate_slip')
           outputID = merge(shearrate_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('accumulated_shear','accumulatedshear','accumulated_shear_slip')
           outputID = merge(accumulatedshear_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('mfp','mfp_slip')
           outputID = merge(mfp_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('resolved_stress','resolved_stress_slip')
           outputID = merge(resolvedstress_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('threshold_stress','threshold_stress_slip')
           outputID = merge(thresholdstress_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('edge_dipole_distance')
           outputID = merge(dipoleDistance_ID,undefined_ID,prm%totalNslip>0_pInt)
         case ('stress_exponent')
           outputID = merge(stressexponent_ID,undefined_ID,prm%totalNslip>0_pInt)
     end select

     if (outputID /= undefined_ID) then       
       plastic_disloUCLA_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_disloUCLA_sizePostResult(i,phase_plasticityInstance(p)) = outputSize   
       prm%outputID = [prm%outputID, outputID]
       plastic_disloUCLA_outputID(i,phase_plasticityInstance(p)) = outputID
       plastic_disloUCLA_sizePostResults(phase_plasticityInstance(p)) = &
          plastic_disloUCLA_sizePostResults(phase_plasticityInstance(p)) + outputSize
plastic_disloUCLA_Noutput(phase_plasticityInstance(p)) = plastic_disloUCLA_Noutput(phase_plasticityInstance(p)) + 1_pInt
     endif

   enddo 
   end associate
 enddo

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
       case ('clambdaslip','friction_coeff')
         do j = 1_pInt, Nchunks_SlipFamilies
           tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('clambdaslip')
             plastic_disloUCLA_CLambdaSlipPerSlipFamily(1:Nchunks_SlipFamilies,instance) = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('friction_coeff')
             plastic_disloUCLA_friction(1:Nchunks_SlipFamilies,instance) = &
                  tempPerSlip(1:Nchunks_SlipFamilies)  
         end select
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
          !if (plastic_disloUCLA_rhoEdge0(f,instance) < 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_rhoEdgeDip0(f,instance) < 0.0_pReal) & 
          !  call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance) < 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOUCLA_label//')')
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
 
 allocate(plastic_disloUCLA_CLambdaSlipPerSlipSystem(maxTotalNslip, maxNinstance),source=0.0_pReal)
 
 allocate(plastic_disloUCLA_forestProjectionEdge(maxTotalNslip,maxTotalNslip,maxNinstance), &
                                                                                       source=0.0_pReal)


 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   myPhase2: if (phase_plasticity(phase) == PLASTICITY_disloUCLA_ID) then
     p = phase
     NofMyPhase=count(material_phase==phase)
     instance = phase_plasticityInstance(phase)
     ns = plastic_disloUCLA_totalNslip(instance)

    associate(prm => param(instance), stt=>state(instance),mse => microstructure(phase_plasticityInstance(p)))
!--------------------------------------------------------------------------------------------------
! allocate state arrays

     sizeDotState     = int(size(['rhoEdge     ','rhoEdgeDip  ','accshearslip']),pInt) * ns
     sizeState        = sizeDotState

   call material_allocatePlasticState(phase,NofMyPhase,sizeState,sizeDotState,0_pInt, &
                                      ns,0_pInt,0_pInt)

     plasticState(phase)%sizePostResults = plastic_disloUCLA_sizePostResults(instance)

     offset_slip = 2_pInt*plasticState(phase)%nSlip
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NofMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state   (offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NofMyPhase)
    !* Process slip related parameters ------------------------------------------------ 
 
     mySlipFamilies: do f = 1_pInt,lattice_maxNslipFamily
       index_myFamily = sum(plastic_disloUCLA_Nslip(1:f-1_pInt,instance))                      ! index in truncated slip system list
       mySlipSystems: do j = 1_pInt,plastic_disloUCLA_Nslip(f,instance)

 
         plastic_disloUCLA_CLambdaSlipPerSlipSystem(index_myFamily+j,instance) = &
         plastic_disloUCLA_CLambdaSlipPerSlipFamily(f,instance)
  
       !* Calculation of forest projections for edge dislocations
         otherSlipFamilies: do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_disloUCLA_Nslip(1:o-1_pInt,instance))
           otherSlipSystems: do k = 1_pInt,plastic_disloUCLA_Nslip(o,instance)
             plastic_disloUCLA_forestProjectionEdge(index_myFamily+j,index_otherFamily+k,instance) = &
               abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,phase))+j,phase), &
                               lattice_st(:,sum(lattice_NslipSystem(1:o-1,phase))+k,phase)))
         enddo otherSlipSystems; enddo otherSlipFamilies
 
       enddo mySlipSystems
     enddo mySlipFamilies
  
     startIndex=1_pInt
     endIndex=ns
     stt%rhoEdge=>plasticState(phase)%state(startIndex:endIndex,:)
     stt%rhoEdge= spread(prm%rho0,2,NofMyPhase)
     dotState(instance)%rhoEdge=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     stt%rhoEdgeDip=>plasticState(phase)%state(startIndex:endIndex,:)
     stt%rhoEdgeDip= spread(prm%rhoDip0,2,NofMyPhase)
     dotState(instance)%rhoEdgeDip=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

     startIndex=endIndex+1_pInt
     endIndex=endIndex+ns
     stt%accshear_slip=>plasticState(phase)%state(startIndex:endIndex,:)
     dotState(instance)%accshear_slip=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(p)%aTolState(startIndex:endIndex) = 1e6_pReal


   allocate(mse%mfp(prm%totalNslip,NofMyPhase),source=0.0_pReal)
   allocate(mse%threshold_stress(prm%totalNslip,NofMyPhase),source=0.0_pReal)


     plasticState(p)%state0 = plasticState(p)%state                                                 ! ToDo: this could be done centrally
   end associate
   endif myPhase2

 enddo initializeInstances
 
end subroutine plastic_disloUCLA_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_microstructure(temperature,ipc,ip,el)
 use math, only: &
   pi
 use material, only: &
   phase_plasticityInstance, &
   phaseAt, phasememberAt, &
   material_phase
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
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
  invLambdaSlip
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance)


 associate(prm => param(instance), stt => state(instance),mse => microstructure(instance))
 !* 1/mean free distance between 2 forest dislocations seen by a moving dislocation
 forall (s = 1_pInt:ns) &
   invLambdaSlip(s) = &
     sqrt(dot_product((stt%rhoEdge(1_pInt:ns,of)+stt%rhoEdgeDip(1_pInt:ns,of)),&
                      plastic_disloUCLA_forestProjectionEdge(1:ns,s,instance)))/ &
     plastic_disloUCLA_CLambdaSlipPerSlipSystem(s,instance)
 
 !* mean free path between 2 obstacles seen by a moving dislocation

 mse%mfp(:,of) = prm%grainSize/(1.0_pReal+prm%grainSize*invLambdaSlip)

 !* threshold stress for dislocation motion
 forall (s = 1_pInt:ns) &
   mse%threshold_stress(s,of) = &
     lattice_mu(ph)*prm%burgers(s)*&
     sqrt(dot_product(stt%rhoEdge(1_pInt:ns,of)+stt%rhoEdgeDip(1_pInt:ns,of),&
                      prm%interaction_SlipSlip(s,1:ns)))
 end associate


end subroutine plastic_disloUCLA_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp,Mp,Temperature,ipc,ip,el)
 use material, only: &
   material_phase, &
   phase_plasticityInstance, &
   phaseAt, phasememberAt
 
 implicit none
 integer(pInt), intent(in)                  :: ipc,ip,el
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(3,3),   intent(in)    :: Mp
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(3,3,3,3), intent(out)   :: dLp_dMp

 integer(pInt) :: instance,ph,of,i,k,l,m,n

 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,gdot_slip_neg,tau_slip_pos,tau_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg
   
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 associate(prm => param(instance))
 
 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal
 
 call kinetics(Mp,Temperature,ph,instance,of, &                                                  
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)
 slipSystems: do i = 1_pInt, prm%totalNslip
   Lp = Lp + (gdot_slip_pos(i)+gdot_slip_neg(i))*prm%Schmid_slip(1:3,1:3,i)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtauslip_pos(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                      + dgdot_dtauslip_neg(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_neg(m,n,i)
 enddo slipSystems
end associate

 Lp = 0.5_pReal * Lp
 dLp_dMp = 0.5_pReal * dLp_dMp

end subroutine plastic_disloUCLA_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dotState(Mp,Temperature,ipc,ip,el)
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
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_mu

 implicit none
 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: instance,ns,f,i,j,index_myFamily, &
                  ph, &
                  of
 real(pReal) :: &
   EdgeDipMinDistance,&
   AtomicVolume,&
   VacancyDiffusion,&
   DotRhoMultiplication,&
   EdgeDipDistance, &
   DotRhoEdgeDipAnnihilation, &
   DotRhoEdgeEdgeAnnihilation, &
   ClimbVelocity, &
   DotRhoEdgeDipClimb, &
   DotRhoDipFormation
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos, gdot_slip_neg,&
   tau_slip_pos,&
   tau_slip_neg, &
   dgdot_dtauslip_neg,dgdot_dtauslip_pos

 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance) 

 plasticState(ph)%dotState(:,of) = 0.0_pReal
  associate(prm => param(instance), stt => state(instance),mse => microstructure(instance))
 !* Dislocation density evolution
 call kinetics(Mp,Temperature,ph,instance,of, &                                                
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)
 dotState(instance)%accshear_slip(:,of) = (gdot_slip_pos+gdot_slip_neg)*0.5_pReal

 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph)) ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
     j = j+1_pInt

     !* Multiplication
     DotRhoMultiplication = abs(dotState(instance)%accshear_slip(j,of))/&
                               (prm%burgers(j)* &
                                mse%mfp(j,of))
 
     !* Dipole formation
     EdgeDipMinDistance = &
       plastic_disloUCLA_CEdgeDipMinDistance(instance)*prm%burgers(j)
     if (dEq0(tau_slip_pos(j))) then
       DotRhoDipFormation = 0.0_pReal
     else
       EdgeDipDistance = &
         (3.0_pReal*lattice_mu(ph)*prm%burgers(j))/&
         (16.0_pReal*pi*abs(tau_slip_pos(j)))
       if (EdgeDipDistance>mse%mfp(j,of)) EdgeDipDistance=mse%mfp(j,of)
       if (EdgeDipDistance<EdgeDipMinDistance) EdgeDipDistance=EdgeDipMinDistance
       DotRhoDipFormation = &
         ((2.0_pReal*EdgeDipDistance)/prm%burgers(j))*&
         stt%rhoEdge(j,of)*abs(dotState(instance)%accshear_slip(j,of))*plastic_disloUCLA_dipoleFormationFactor(instance)
     endif
 
    !* Spontaneous annihilation of 2 single edge dislocations
    DotRhoEdgeEdgeAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/prm%burgers(j))*&
        stt%rhoEdge(j,of)*abs(dotState(instance)%accshear_slip(j,of))
 
    !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
    DotRhoEdgeDipAnnihilation = &
        ((2.0_pReal*EdgeDipMinDistance)/prm%burgers(j))*&
        stt%rhoEdgeDip(j,of)*abs(dotState(instance)%accshear_slip(j,of))
 
      !* Dislocation dipole climb
     AtomicVolume = &
        plastic_disloUCLA_CAtomicVolume(instance)*prm%burgers(j)**(3.0_pReal)
     VacancyDiffusion = &
        plastic_disloUCLA_D0(instance)*exp(-plastic_disloUCLA_Qsd(instance)/(kB*Temperature))
     if (dEq0(tau_slip_pos(j))) then
       DotRhoEdgeDipClimb = 0.0_pReal
     else
       ClimbVelocity = &
          ((3.0_pReal*lattice_mu(ph)*VacancyDiffusion*AtomicVolume)/(2.0_pReal*pi*kB*Temperature))*&
          (1/(EdgeDipDistance+EdgeDipMinDistance))
       DotRhoEdgeDipClimb = &
          (4.0_pReal*ClimbVelocity*stt%rhoEdgeDip(j,of))/(EdgeDipDistance-EdgeDipMinDistance)
     endif
 
     !* Edge dislocation density rate of change
     dotState(instance)%rhoEdge(j,of) = &
        DotRhoMultiplication-DotRhoDipFormation-DotRhoEdgeEdgeAnnihilation
 
     !* Edge dislocation dipole density rate of change
     dotState(instance)%rhoEdgeDip(j,of) = &
        DotRhoDipFormation-DotRhoEdgeDipAnnihilation-DotRhoEdgeDipClimb

 
   enddo slipSystems
 enddo slipFamilies
end associate
 
end subroutine plastic_disloUCLA_dotState

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_disloUCLA_postResults(Mp,Temperature,ipc,ip,el) result(postResults)
 use prec, only: &
   tol_math_check, &
   dEq, dNeq0
 use math, only: &
   pi, &
math_mul33xx33
 use material, only: &
   material_phase, &
   phase_plasticityInstance,& 
   !plasticState, &
   phaseAt, phasememberAt
 use lattice, only: &
   lattice_Sslip, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_mu

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(plastic_disloUCLA_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           postResults

 integer(pInt) :: &
   instance,&
   ns,&
   f,o,i,c,j,index_myFamily,&
   ph, &
   of
 real(pReal), dimension(plastic_disloUCLA_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg
 
 !* Shortened notation
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance  = phase_plasticityInstance(ph)
 ns = plastic_disloUCLA_totalNslip(instance)
 
 !* Required output
 c = 0_pInt
 postResults = 0.0_pReal
 associate (prm => param(instance),stt =>state(instance),mse => microstructure(instance))
 do o = 1_pInt,plastic_disloUCLA_Noutput(instance)
    select case(plastic_disloUCLA_outputID(o,instance))
 
      case (rho_ID)
        postResults(c+1_pInt:c+ns) = stt%rhoEdge(1_pInt:ns,of)
        c = c + ns
      case (rhoDip_ID)
        postResults(c+1_pInt:c+ns) = stt%rhoEdgeDip(1_pInt:ns,of)
        c = c + ns
      case (shearrate_ID,stressexponent_ID)
 call kinetics(Mp,Temperature,ph,instance,of, &                                                
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)

        if     (plastic_disloUCLA_outputID(o,instance) == shearrate_ID) then
          postResults(c+1:c+ns) = (gdot_slip_pos + gdot_slip_neg)*0.5_pReal
          c = c + ns
        elseif(plastic_disloUCLA_outputID(o,instance) == stressexponent_ID) then
          do j = 1_pInt, ns
            if (dEq(gdot_slip_pos(j)+gdot_slip_neg(j),0.0_pReal)) then
              postResults(c+j) = 0.0_pReal
            else
              postResults(c+j) = (tau_slip_pos(j)+tau_slip_neg(j))/&
                                                       (gdot_slip_pos(j)+gdot_slip_neg(j))*&
                                                       (dgdot_dtauslip_pos(j)+dgdot_dtauslip_neg(j))* 0.5_pReal
            endif
          enddo
           c = c + ns
        endif

      case (accumulatedshear_ID)
       postResults(c+1_pInt:c+ns) = &
                      stt%accshear_slip(1_pInt:ns, of)
        c = c + ns
      case (mfp_ID)
        postResults(c+1_pInt:c+ns) = mse%mfp(1_pInt:ns, of)
        c = c + ns
      case (resolvedstress_ID)
        j = 0_pInt
        slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems1: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
              j = j + 1_pInt
              postResults(c+j) =&
                                math_mul33xx33(Mp,lattice_Sslip(:,:,1,index_myFamily+i,ph))
        enddo slipSystems1; enddo slipFamilies1
        c = c + ns
      case (thresholdstress_ID)
        postResults(c+1_pInt:c+ns) = mse%threshold_stress(1_pInt:ns,of)
        c = c + ns
      case (dipoleDistance_ID)
        j = 0_pInt
        slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
           index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                 ! at which index starts my family
           slipSystems2: do i = 1_pInt,plastic_disloUCLA_Nslip(f,instance)
              j = j + 1_pInt
              if (dNeq0(abs(math_mul33xx33(Mp,lattice_Sslip(:,:,1,index_myFamily+i,ph))))) then
              postResults(c+j) = &
                (3.0_pReal*lattice_mu(ph)*prm%burgers(j))/&
                (16.0_pReal*pi*abs(math_mul33xx33(Mp,lattice_Sslip(:,:,1,index_myFamily+i,ph))))
              else
              postResults(c+j) = huge(1.0_pReal)
              endif
              postResults(c+j)=min(postResults(c+j),&
                                                            mse%mfp(j,of))
        enddo slipSystems2; enddo slipFamilies2
        c = c + ns
    end select
 enddo
end associate
end function plastic_disloUCLA_postResults


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
subroutine kinetics(Mp,Temperature,ph,instance,of, &
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)
 use prec, only: &
   tol_math_check, &
   dEq, dNeq0
 use math, only: &
   pi, &
math_mul33xx33
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
ph, instance,of

 integer(pInt) :: &
   ns,&
   f,i,j,index_myFamily
 real(pReal) :: StressRatio_p,StressRatio_pminus1,&
                BoltzmannRatio,DotGamma0,stressRatio,&
                dvel_slip, vel_slip
 real(pReal), intent(out), dimension(plastic_disloUCLA_totalNslip(instance)) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg
 associate(prm => param(instance), stt => state(instance),mse => microstructure(instance))
 !* Shortened notation
 ns = plastic_disloUCLA_totalNslip(instance)
 

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
            BoltzmannRatio = prm%H0kp(j)/(kB*Temperature)
            !* Initial shear rates
            DotGamma0 = stt%rhoEdge(j,of)*prm%burgers(j)*prm%v0(j)
            !* Resolved shear stress on slip system
            tau_slip_pos(j) = math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,j))
            tau_slip_neg(j) = math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,j))

            significantPositiveTau: if((abs(tau_slip_pos(j))-mse%threshold_stress(j, of)) > tol_math_check) then
              !* Stress ratio
              stressRatio = ((abs(tau_slip_pos(j))-mse%threshold_stress(j, of))/&
                      (prm%solidSolutionStrength+&
                       prm%tau_Peierls(j)))
              stressRatio_p       = stressRatio** prm%p(j)
              stressRatio_pminus1 = stressRatio**(prm%p(j)-1.0_pReal)
              !* Shear rates due to slip
              vel_slip = 2.0_pReal*prm%burgers(j) &
                     * prm%kink_height(j) * prm%omega(j)  &
                     * ( mse%mfp(j,of) - prm%kink_width(j) ) &
                     * (tau_slip_pos(j)  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) ) &
                     / ( &
                     2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
                     + prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                     *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                     )
                       
              gdot_slip_pos(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_pos(j))
              !* Derivatives of shear rates 

              dvel_slip = &
                   2.0_pReal*prm%burgers(j) &
                   * prm%kink_height(j) * prm%omega(j)  &
                   * ( mse%mfp(j,of) - prm%kink_width(j) ) &
                   * ( &
                   (exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) &
                   + tau_slip_pos(j) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)))&    !deltaf(i)                                         
                   *BoltzmannRatio*prm%p(j)&
                   *prm%q(j)/&
                   (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(prm%q(j)-1.0_pReal)  ) &!deltaf(f)                                        
                   ) &
                   *  (2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
                   +  prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                   ) &
                   -  (tau_slip_pos(j) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) )  &
                   *  (2.0_pReal*(prm%burgers(j)**2.0_pReal) &
                   +  prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)))&     !deltaf(i)                                        
                   *BoltzmannRatio*prm%p(j)&
                   *prm%q(j)/&
                   (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(prm%q(j)-1.0_pReal)  )& !deltaf(f)                                        
                   ) &
                   )  &
                   / (  &
                   ( &
                   2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
                   + prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                   )**2.0_pReal &
                   )

              dgdot_dtauslip_pos(j) = DotGamma0 * dvel_slip

            endif significantPositiveTau
            significantNegativeTau: if((abs(tau_slip_neg(j))-mse%threshold_stress(j, of)) > tol_math_check) then
              !* Stress ratios
              stressRatio = ((abs(tau_slip_neg(j))-mse%threshold_stress(j, of))/&
                      (prm%solidSolutionStrength+&
                       prm%tau_Peierls(j)))
              stressRatio_p       = stressRatio** prm%p(j)
              stressRatio_pminus1 = stressRatio**(prm%p(j)-1.0_pReal)
              !* Shear rates due to slip                                                                                                                                                                                                                                                                           
              vel_slip = 2.0_pReal*prm%burgers(j) &
                     * prm%kink_height(j) * prm%omega(j)  &
                     * ( mse%mfp(j,of) - prm%kink_width(j) ) &
                     * (tau_slip_neg(j)  &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) ) &
                     / ( &
                     2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
                     + prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                     *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                     * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                     )
              
              gdot_slip_neg(j) = DotGamma0 &
                       * vel_slip & 
                       * sign(1.0_pReal,tau_slip_neg(j))
              !* Derivatives of shear rates 
              dvel_slip = &
                   2.0_pReal*prm%burgers(j) &
                   * prm%kink_height(j) * prm%omega(j)  &
                   * ( mse%mfp(j,of) - prm%kink_width(j) ) &
                   * ( &
                   (exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) &
                   + tau_slip_neg(j) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)))&    !deltaf(i)                                         
                   *BoltzmannRatio*prm%p(j)&
                   *prm%q(j)/&
                   (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(prm%q(j)-1.0_pReal)  ) &!deltaf(f)                                        
                   ) &
                   *  (2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
                   +  prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                   ) &
                   -  (tau_slip_neg(j) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)) )  &
                   *  (2.0_pReal*(prm%burgers(j)**2.0_pReal) &
                   +  prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * (abs(exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j)))&     !deltaf(i)                                        
                   *BoltzmannRatio*prm%p(j)&
                   *prm%q(j)/&
                   (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
                   StressRatio_pminus1*(1-StressRatio_p)**(prm%q(j)-1.0_pReal)  )& !deltaf(f)                                        
                   ) &
                   )  &
                   / (  &
                   ( &
                   2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
                   + prm%omega(j) * plastic_disloUCLA_friction(f,instance) &
                   *(( mse%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
                   * exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q(j))  &
                   )**2.0_pReal &
                   )


              dgdot_dtauslip_neg(j) = DotGamma0 * dvel_slip

            endif significantNegativeTau
          enddo slipSystems
        enddo slipFamilies
 end associate
end subroutine kinetics

end module plastic_disloUCLA
