!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power law formulation for the shear rates
!!         and a Voce-type kinematic hardening rule
!--------------------------------------------------------------------------------------------------
module plastic_kinehardening
 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_sizePostResult                                                             !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_output                                                                     !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_kinehardening_Noutput                                                                    !< number of outputs per instance
 
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_kinehardening_totalNslip                                                                 !< no. of slip system used in simulation
  

 integer(pInt),                       dimension(:,:),   allocatable,         private :: &
   plastic_kinehardening_Nslip                                                                      !< active number of slip systems per family (input parameter, per family)
   

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     crss_ID, &                                                                         !< critical resolved stress
     crss_back_ID, &                                                                    !< critical resolved back stress
     sense_ID, &                                                                        !< sense of acting shear stress (-1 or +1)
     chi0_ID, &                                                                         !< backstress at last switch of stress sense (positive?)
     gamma0_ID, &                                                                       !< accumulated shear at last switch of stress sense (at current switch?)
     accshear_ID, &
     shearrate_ID, &
     resolvedstress_ID
 end enum


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     gdot0, &                                                                                       !< reference shear strain rate for slip (input parameter)
     n_slip, &                                                                                      !< stress exponent for slip (input parameter)
     aTolResistance, &
     aTolShear
     
   
   real(pReal),                         dimension(:),   allocatable,          private :: &
     crss0, &                                                                                        !< initial critical shear stress for slip (input parameter, per family)
     theta0, &                                                                                      !< initial hardening rate of forward stress for each slip
     theta1, &                                                                                      !< asymptotic hardening rate of forward stress for each slip >
     theta0_b, &                                                                                    !< initial hardening rate of back stress for each slip > 
     theta1_b, &                                                                                    !< asymptotic hardening rate of back stress for each slip >
     tau1, &
     tau1_b, &
     interaction_slipslip, &                                                                        !< latent hardening matrix
     nonSchmidCoeff

   real(pReal),                 allocatable, dimension(:,:,:) :: &
     Schmid_slip, &
     Schmid_twin, &
     nonSchmid_pos, &
     nonSchmid_neg
        
   real(pReal),                       dimension(:,:),   allocatable,          private :: &
     hardeningMatrix_SlipSlip
   integer(pInt) :: &
     totalNslip                                                                                     !< total number of active slip system
   integer(pInt),               allocatable, dimension(:) :: &
     Nslip                                                                                          !< number of active slip systems for each family
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
 end type

 type, private :: tKinehardeningState
   real(pReal), pointer, dimension(:,:) :: &                                                        !< vectors along NipcMyInstance
     crss, &                                                                                        !< critical resolved stress
     crss_back, &                                                                                   !< critical resolved back stress
     sense, &                                                                                       !< sense of acting shear stress (-1 or +1)
     chi0, &                                                                                        !< backstress at last switch of stress sense
     gamma0, &                                                                                      !< accumulated shear at last switch of stress sense
     accshear                                                                                       !< accumulated (absolute) shear

   real(pReal), pointer, dimension(:) :: &                                                          !< scalars along NipcMyInstance
     sumGamma                                                                                       !< accumulated shear across all systems
 end type

 type(tParameters), dimension(:), allocatable, private :: &
   param, &                                                                                            !< containers of constitutive parameters (len Ninstance)
   paramNew ! temp

 type(tKinehardeningState), allocatable, dimension(:), private :: &
   dotState, &
   deltaState, &
   state, &
   state0

  
 public :: &
   plastic_kinehardening_init, &
   plastic_kinehardening_LpAndItsTangent, &
   plastic_kinehardening_dotState, &
   plastic_kinehardening_deltaState, &
   plastic_kinehardening_postResults
 private :: &
   plastic_kinehardening_shearRates


contains



!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   dEq0
 use debug, only: &
   debug_level, &
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
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
   material_allocatePlasticState, &
   PLASTICITY_kinehardening_label, &
   PLASTICITY_kinehardening_ID, &
   material_phase, &
   plasticState
 use config, only: &
   config_phase, &
   MATERIAL_partPhase
 use lattice

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(kind(undefined_ID)) :: &
   output_ID
 integer(pInt) :: &
   o, i,j, k, f, p,  &
   phase, & 
   instance, &
   maxNinstance, &
   NipcMyPhase, &
   outputSize, &
   Nchunks_SlipSlip = 0_pInt, Nchunks_SlipFamilies = 0_pInt, &
   Nchunks_nonSchmid = 0_pInt, &
   offset_slip, index_myFamily, index_otherFamily, &
   startIndex, endIndex, &
   mySize, nSlip, nSlipFamilies, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState

 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 real(pReal), dimension(:), allocatable :: tempPerSlip
 integer(kind(undefined_ID)) :: &
   outputID                                                                                         !< ID of each post result output
   
 character(len=65536), dimension(:), allocatable :: &
  outputs
 character(len=65536) :: &
   tag       = '', &
   line      = '', &
   extmsg    = '', &
   structure    = ''

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_KINEHARDENING_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == PLASTICITY_KINEHARDENING_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a,1x,i5,/)') '# instances:',maxNinstance                                      
   
 allocate(plastic_kinehardening_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                             source=0_pInt)
 allocate(plastic_kinehardening_output(maxval(phase_Noutput),maxNinstance))
          plastic_kinehardening_output                                             = ''
 allocate(plastic_kinehardening_Noutput(maxNinstance),                       source=0_pInt)
 allocate(plastic_kinehardening_Nslip(lattice_maxNslipFamily,maxNinstance),  source=0_pInt)
 allocate(plastic_kinehardening_totalNslip(maxNinstance),                    source=0_pInt)
 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
 allocate(paramNew(maxNinstance))
 allocate(state(maxNinstance))
 allocate(state0(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(deltaState(maxNinstance))
 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_KINEHARDENING_ID) cycle
     instance = phase_plasticityInstance(p)                                                     ! which instance of my phase
   associate(prm => paramNew(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             delta => deltaState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)))

   structure          = config_phase(p)%getString('lattice_structure')

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
   prm%aTolResistance = config_phase(p)%getFloat('atol_resistance',defaultVal=1.0_pReal)
   prm%aTolShear      = config_phase(p)%getFloat('atol_shear',     defaultVal=1.0e-6_pReal)

   ! sanity checks
   if (prm%aTolResistance <= 0.0_pReal) extmsg = trim(extmsg)//'aTolresistance '
   if (prm%aTolShear      <= 0.0_pReal) extmsg = trim(extmsg)//'aTolShear '

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

     prm%crss0                = config_phase(p)%getFloats('crss0',   requiredShape=shape(prm%Nslip))
     prm%tau1                 = config_phase(p)%getFloats('tau1', requiredShape=shape(prm%Nslip))
     prm%tau1_b                 = config_phase(p)%getFloats('tau1_b', requiredShape=shape(prm%Nslip))
     prm%theta0                 = config_phase(p)%getFloats('theta0', requiredShape=shape(prm%Nslip))
     prm%theta1                 = config_phase(p)%getFloats('theta1', requiredShape=shape(prm%Nslip))
     prm%theta0_b                 = config_phase(p)%getFloats('theta0_b', requiredShape=shape(prm%Nslip))
     prm%theta1_b                 = config_phase(p)%getFloats('theta1_b', requiredShape=shape(prm%Nslip))


     prm%gdot0           = config_phase(p)%getFloat('gdot0')
     prm%n_slip               = config_phase(p)%getFloat('n_slip')



     !prm%interaction_SlipSlip = lattice_interaction_SlipSlip(prm%Nslip, &
     !                                                        config_phase(p)%getFloats('interaction_slipslip'), &
     !                                                        structure(1:3))
   endif slipActive

   
!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
        case ('resistance')
          outputID = merge(crss_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('accumulatedshear')
          outputID = merge(accshear_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('shearrate')
          outputID = merge(shearrate_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('resolvedstress')
          outputID = merge(resolvedstress_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('backstress')
          outputID = merge(crss_back_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('sense')
          outputID = merge(sense_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('chi0')
          outputID = merge(chi0_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('gamma0')
          outputID = merge(gamma0_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip

      end select

      if (outputID /= undefined_ID) then
        plastic_kinehardening_Noutput(instance) = plastic_kinehardening_Noutput(instance) + 1_pInt
        plastic_kinehardening_output(i,phase_plasticityInstance(p)) = outputs(i)
        plastic_kinehardening_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
        prm%outputID = [prm%outputID , outputID]
      endif

   end do
param(instance)%outputID = prm%outputID
   nslip = prm%totalNslip
!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeDotState = nSlip &                                        !< crss
                  + nSlip &                                        !< crss_back
                  + nSlip &                                        !< accumulated (absolute) shear
                  + 1_pInt                                         !< sum(gamma)
               
     sizeDeltaState = nSlip &                                      !< sense of acting shear stress (-1 or +1)
                    + nSlip &                                      !< backstress at last switch of stress sense
                    + nSlip                                        !< accumulated shear at last switch of stress sense

     sizeState = sizeDotState + sizeDeltaState
     NipcMyPhase = count(material_phase == p)                                                   ! number of IPCs containing my phase
     call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,sizeDeltaState, &
                                        nSlip,0_pInt,0_pInt)
     plasticState(p)%sizePostResults = sum(plastic_kinehardening_sizePostResult(:,phase_plasticityInstance(p)))
     plasticState(p)%offsetDeltaState = sizeDotState


     endindex = 0_pInt
     o = endIndex                                                                                           ! offset of dotstate index relative to state index
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     stt%crss          => plasticState(p)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     dot%crss          => plasticState(p)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex-o:endIndex-o) = prm%aTolResistance
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip 
     stt%crss_back          => plasticState(p)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     dot%crss_back          => plasticState(p)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex-o:endIndex-o) = prm%aTolResistance
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     stt%accshear          => plasticState(p)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     dot%accshear          => plasticState(p)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex-o:endIndex-o) = prm%aTolShear
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + 1_pInt
     stt%sumGamma        => plasticState(p)%state    (startIndex             ,1:NipcMyPhase)
     dot%sumGamma        => plasticState(p)%dotState (startIndex-o           ,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex-o:endIndex-o) =prm%aTolShear
     
!----------------------------------------------------------------------------------------------
!locally define deltaState alias
     o = endIndex
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     stt%sense          => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     delta%sense        => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     stt%chi0           => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     delta%chi0         => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)

     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     stt%gamma0         => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)         
     delta%gamma0       => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)


   end associate
 end do

 
 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     if (phase_plasticity(phase) == PLASTICITY_KINEHARDENING_ID) then
       instance = phase_plasticityInstance(phase)                                                   ! count instances of my constitutive law
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)                          ! maximum number of slip families according to lattice type of current phase
       Nchunks_SlipSlip     = maxval(lattice_interactionSlipSlip(:,:,phase))
       Nchunks_nonSchmid    = lattice_NnonSchmid(phase)
       allocate(param(instance)%crss0   (Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%tau1    (Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%tau1_b  (Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%theta0  (Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%theta1  (Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%theta0_b(Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%theta1_b(Nchunks_SlipFamilies), source=0.0_pReal)
       allocate(param(instance)%interaction_slipslip(Nchunks_SlipSlip), source=0.0_pReal)
       allocate(param(instance)%nonSchmidCoeff(Nchunks_nonSchmid),      source=0.0_pReal)
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (phase_plasticity(phase) == PLASTICITY_KINEHARDENING_ID) then      ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)


!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip families 
       case ('nslip')
         if (chunkPos(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_KINEHARDENING_label//')')
         if (chunkPos(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_KINEHARDENING_label//')')
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt                                                 ! user specified number of (possibly) active slip families (e.g. 6 0 6 --> 3)
         do j = 1_pInt, Nchunks_SlipFamilies
           plastic_kinehardening_Nslip(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
      
       case ('crss0','tau1','tau1_b','theta0','theta1','theta0_b','theta1_b')
         tempPerSlip = 0.0_pReal
         do j = 1_pInt, Nchunks_SlipFamilies
           if (plastic_kinehardening_Nslip(j,instance) > 0_pInt) &
             tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('crss0')
             param(instance)%crss0(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('tau1')
             param(instance)%tau1(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('tau1_b')
             param(instance)%tau1_b(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('theta0')
             param(instance)%theta0(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('theta1')
             param(instance)%theta1(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('theta0_b')
             param(instance)%theta0_b(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
           case ('theta1_b')
             param(instance)%theta1_b(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)  
         end select
          
!--------------------------------------------------------------------------------------------------
! parameters depending on number of interactions 
       case ('interaction_slipslip')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_KINEHARDENING_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           param(instance)%interaction_slipslip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
       case ('nonschmidcoeff')
         if (chunkPos(1) < 1_pInt + Nchunks_nonSchmid) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_KINEHARDENING_label//')')
         do j = 1_pInt,Nchunks_nonSchmid
           param(instance)%nonSchmidCoeff(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo  
!--------------------------------------------------------------------------------------------------
       case ('gdot0')
         param(instance)%gdot0                    = IO_floatValue(line,chunkPos,2_pInt)
             
       case ('n_slip')
         param(instance)%n_slip                   = IO_floatValue(line,chunkPos,2_pInt)
               
       case default

     end select
   endif; endif
 enddo parsingFile

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems



 initializeInstances: do phase = 1_pInt, size(phase_plasticity)                                     ! loop through all phases in material.config
   myPhase2: if (phase_plasticity(phase) == PLASTICITY_KINEHARDENING_ID) then                       ! only consider my phase
     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my phase
     plastic_kinehardening_Nslip(1:lattice_maxNslipFamily,instance) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase),&                                    ! limit active slip systems per family to min of available and requested
           plastic_kinehardening_Nslip(1:lattice_maxNslipFamily,instance))

     plastic_kinehardening_totalNslip(instance)  = sum(plastic_kinehardening_Nslip(:,instance))     ! how many slip systems altogether
     nSlipFamilies = count(plastic_kinehardening_Nslip(:,instance) > 0_pInt)
     nSlip = plastic_kinehardening_totalNslip(instance)                                             ! total number of active slip systems
     
!--------------------------------------------------------------------------------------------------
!  sanity checks
 
     if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%crss0 (1:nSlipFamilies)          < 0.0_pReal)) extmsg = trim(extmsg)//' crss0'
     if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%tau1  (1:nSlipFamilies)         <= 0.0_pReal)) extmsg = trim(extmsg)//' tau1'
     if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%tau1_b(1:nSlipFamilies)         < 0.0_pReal)) extmsg = trim(extmsg)//' tau1_b'
     if (param(instance)%gdot0            <= 0.0_pReal) extmsg = trim(extmsg)//' gdot0'
     if (param(instance)%n_slip           <= 0.0_pReal) extmsg = trim(extmsg)//' n_slip'               
     if (extmsg /= '') then 
       extmsg = trim(extmsg)//' ('//PLASTICITY_KINEHARDENING_label//')'                                 ! prepare error message identifier
       call IO_error(211_pInt,ip=instance,ext_msg=extmsg)
     endif
   
    
     offset_slip = plasticState(phase)%nSlip+plasticState(phase)%nTwin+2_pInt
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase) 

     allocate(param(instance)%hardeningMatrix_SlipSlip(nSlip,nSlip),  source=0.0_pReal)
     do f = 1_pInt,lattice_maxNslipFamily                                                                    ! >>> interaction slip -- X
       index_myFamily = sum(plastic_kinehardening_Nslip(1:f-1_pInt,instance))
        do j = 1_pInt,plastic_kinehardening_Nslip(f,instance)                                                ! loop over (active) systems in my family (slip)
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(plastic_kinehardening_Nslip(1:o-1_pInt,instance))
           do k = 1_pInt,plastic_kinehardening_Nslip(o,instance)                                        ! loop over (active) systems in other family (slip)
             param(instance)%hardeningMatrix_SlipSlip(index_myFamily+j,index_otherFamily+k) = &
                 param(instance)%interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                   sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                   sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                   phase))
           enddo; enddo
     enddo; enddo
     
     endindex = 0_pInt
     o = endIndex                                                                                           ! offset of dotstate index relative to state index
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state0  (instance)%crss          => plasticState(phase)%state0   (startIndex  :endIndex  ,1:NipcMyPhase)

     state0(instance)%crss                                  = spread(math_expand(param(instance)%crss0,&
                                                                                 plastic_kinehardening_Nslip(:,instance)), &
                                                                     2, NipcMyPhase)
   endif myPhase2
 enddo initializeInstances

end subroutine plastic_kinehardening_init

!--------------------------------------------------------------------------------------------------
!> @brief calculation of shear rates (\dot \gamma)
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                            Mp,ph,instance,of)

 use math
 use lattice, only: &
   lattice_NslipSystem, &
   lattice_Sslip, &
   lattice_maxNslipFamily, &
   lattice_NnonSchmid

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 integer(pInt),               intent(in) :: &
   ph, &                                                                                           !< phase ID
   instance, &                                                                                     !< instance of that phase
   of                                                                                              !< index of phaseMember
 real(pReal), dimension(plastic_kinehardening_totalNslip(instance)), intent(out) :: &
   gdot_pos, &                                                                                     !< shear rates from positive line segments
   gdot_neg, &                                                                                     !< shear rates from negative line segments
   tau_pos, &                                                                                      !< shear stress on positive line segments
   tau_neg                                                                                         !< shear stress on negative line segments

 integer(pInt) :: &
   index_myFamily, &
   f,i,j,k


 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
     j = j + 1_pInt
     tau_pos(j) = math_mul33xx33(Mp,lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph))
     tau_neg(j) = tau_pos(j)
     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph)
       tau_pos(j) = tau_pos(j) + param(instance)%nonSchmidCoeff(k)* &
                                 math_mul33xx33(Mp,lattice_Sslip(1:3,1:3,2*k+0,index_myFamily+i,ph))
       tau_neg(j) = tau_neg(j) + param(instance)%nonSchmidCoeff(k)* &
                                 math_mul33xx33(Mp,lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems
   enddo slipSystems
 enddo slipFamilies

 gdot_pos = 0.5_pReal * param(instance)%gdot0 * &
            (abs(tau_pos-state(instance)%crss_back(:,of))/ &
            state(instance)%crss(:,of))**param(instance)%n_slip &
            *sign(1.0_pReal,tau_pos-state(instance)%crss_back(:,of)) 
 gdot_neg = 0.5_pReal * param(instance)%gdot0 * &
            (abs(tau_neg-state(instance)%crss_back(:,of))/ &
            state(instance)%crss(:,of))**param(instance)%n_slip &
            *sign(1.0_pReal,tau_neg-state(instance)%crss_back(:,of)) 
            

end subroutine plastic_kinehardening_shearRates


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_LpAndItsTangent(Lp,dLp_dMp, &
                                                 Mp,ipc,ip,el)
 use prec, only: &
   dNeq0
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33, &
   math_transpose33
 use lattice, only: &
   lattice_Sslip, &       !< schmid matrix
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid
 use material, only: &
   phaseAt, phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to the Mandel stress

 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(3,3), intent(in) :: &
   Mp

 integer(pInt) :: &
   instance, &
   index_myFamily, &
   f,i,j,k,l,m,n, &
   of, &
   ph
   
 real(pReal), dimension(plastic_kinehardening_totalNslip(phase_plasticityInstance(phaseAt(ipc,ip,el)))) :: &
   gdot_pos,gdot_neg, &
   tau_pos,tau_neg
 real(pReal) :: &
   dgdot_dtau_pos,dgdot_dtau_neg
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor

 ph = phaseAt(ipc,ip,el)                                                                            !< figures phase for each material point 
 of = phasememberAt(ipc,ip,el)                                                                      !< index of the positions of each constituent of material point, phasememberAt is a function in material that helps figure them out
 instance = phase_plasticityInstance(ph)

 Lp = 0.0_pReal 
 dLp_dMp = 0.0_pReal

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Mp,ph,instance,of)


 j = 0_pInt                                                                                          ! reading and marking the starting index for each slip family
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                          ! at which index starts my family
   slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
     j = j + 1_pInt

     ! build nonSchmid tensor
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     do k = 1,lattice_NnonSchmid(ph)
       nonSchmid_tensor(1:3,1:3,1) = &
       nonSchmid_tensor(1:3,1:3,1) + param(instance)%nonSchmidCoeff(k) * &
                                     lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = &
       nonSchmid_tensor(1:3,1:3,2) + param(instance)%nonSchmidCoeff(k) * &
                                     lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo

     Lp = Lp + (gdot_pos(j)+gdot_neg(j))*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)                ! sum of all gdot*SchmidTensor gives Lp

     ! Calculation of the tangent of Lp                                                              ! sensitivity of Lp
     if (dNeq0(gdot_pos(j))) then
       dgdot_dtau_pos = gdot_pos(j)*param(instance)%n_slip/(tau_pos(j)-state(instance)%crss_back(j,of))
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dMp(k,l,m,n) = &
         dLp_dMp(k,l,m,n) + dgdot_dtau_pos*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                  nonSchmid_tensor(m,n,1)
     endif

     if (dNeq0(gdot_neg(j))) then
       dgdot_dtau_neg = gdot_neg(j)*param(instance)%n_slip/(tau_neg(j)-state(instance)%crss_back(j,of))
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dMp(k,l,m,n) = &
         dLp_dMp(k,l,m,n) + dgdot_dtau_neg*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                  nonSchmid_tensor(m,n,2)
     endif
   enddo slipSystems
 enddo slipFamilies


end subroutine plastic_kinehardening_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_deltaState(Mp,ipc,ip,el)
 use prec, only: &
   dNeq, &
   dEq0
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use material, only: &
   phaseAt, &
   phasememberAt, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plastic_kinehardening_totalNslip(phase_plasticityInstance(phaseAt(ipc,ip,el)))) :: &
   gdot_pos,gdot_neg, &
   tau_pos,tau_neg, &
   sense
 integer(pInt) :: &
   ph, &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of, &
   j                                                                                                !< shortcut notation for offset position in state array

 ph = phaseAt(ipc,ip,el)
 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(ph)

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Mp,ph,instance,of)
 sense = merge(state(instance)%sense(:,of), &                                                       ! keep existing...
               sign(1.0_pReal,gdot_pos+gdot_neg), &                                                 ! ...or have a defined 
               dEq0(gdot_pos+gdot_neg,1e-10_pReal))                                                 ! current sense of shear direction

#ifdef DEBUG
         if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
           write(6,'(a)') '======= kinehardening delta state ======='
         endif
#endif

!--------------------------------------------------------------------------------------------------
! switch in sense of shear?
 do j = 1,plastic_kinehardening_totalNslip(instance)
#ifdef DEBUG
         if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
                   .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
           write(6,'(i2,1x,f7.4,1x,f7.4)') j,sense(j),state(instance)%sense(j,of)
         endif
#endif
   if (dNeq(sense(j),state(instance)%sense(j,of),0.1_pReal)) then
     deltaState(instance)%sense (j,of) = sense(j) - state(instance)%sense(j,of)                              ! switch sense
     deltaState(instance)%chi0  (j,of) = abs(state(instance)%crss_back(j,of)) - state(instance)%chi0(j,of)   ! remember current backstress magnitude
     deltaState(instance)%gamma0(j,of) = state(instance)%accshear(j,of) - state(instance)%gamma0(j,of)       ! remember current accumulated shear
   else
     deltaState(instance)%sense (j,of) = 0.0_pReal                                                           ! no change
     deltaState(instance)%chi0  (j,of) = 0.0_pReal
     deltaState(instance)%gamma0(j,of) = 0.0_pReal
   endif
 enddo

end subroutine plastic_kinehardening_deltaState



!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_dotState(Mp,ipc,ip,el)
 use lattice, only: &
   lattice_maxNslipFamily
 use material, only: &
   material_phase, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element !< microstructure state

 integer(pInt) :: &
   instance,ph, &
   f,i,j, &
   nSlip, &
   of
 
 real(pReal), dimension(plastic_kinehardening_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_pos,gdot_neg, &
   tau_pos,tau_neg
 
 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 nSlip = plastic_kinehardening_totalNslip(instance)
 
 dotState(instance)%sumGamma(of) = 0.0_pReal

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Mp,ph,instance,of)
                                       
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
     j = j+1_pInt    
     dotState(instance)%crss(j,of) = &                                                                               ! evolution of slip resistance j
          dot_product(param(instance)%hardeningMatrix_SlipSlip(j,1:nSlip),abs(gdot_pos+gdot_neg)) * &
          ( param(instance)%theta1(f) + &
           (param(instance)%theta0(f) - param(instance)%theta1(f) &
            + param(instance)%theta0(f)*param(instance)%theta1(f)*state(instance)%sumGamma(of)/param(instance)%tau1(f)) &
           *exp(-state(instance)%sumGamma(of)*param(instance)%theta0(f)/param(instance)%tau1(f)) &                   ! V term depending on the harding law
          )
     dotState(instance)%crss_back(j,of) = &                                                                          ! evolution of back stress resistance j
          state(instance)%sense(j,of)*abs(gdot_pos(j)+gdot_neg(j)) * &
          ( param(instance)%theta1_b(f) + &
           (param(instance)%theta0_b(f) - param(instance)%theta1_b(f) &
            + param(instance)%theta0_b(f)*param(instance)%theta1_b(f)/(param(instance)%tau1_b(f)+state(instance)%chi0(j,of)) &
            *(state(instance)%accshear(j,of)-state(instance)%gamma0(j,of))) &
           *exp(-(state(instance)%accshear(j,of)-state(instance)%gamma0(j,of)) &
                 *param(instance)%theta0_b(f)/(param(instance)%tau1_b(f)+state(instance)%chi0(j,of))) &
          )                                                                                                    ! V term depending on the harding law for back stress
    
     dotState(instance)%accshear(j,of) = abs(gdot_pos(j)+gdot_neg(j))
     dotState(instance)%sumGamma(of) = dotState(instance)%sumGamma(of) + dotState(instance)%accshear(j,of)
   enddo slipSystems
 enddo slipFamilies

end subroutine plastic_kinehardening_dotState

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_kinehardening_postResults(Mp,ipc,ip,el) result(postResults)
 use math
 use material, only: &
   material_phase, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_Sslip, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(sum(plastic_kinehardening_sizePostResult(:,phase_plasticityInstance(material_phase(ipc,ip,el))))) :: &
   postResults
 integer(pInt) :: &
   instance,ph, of, &
   nSlip,&
   o,f,i,c,j,&
   index_myFamily
   
 real(pReal), dimension(plastic_kinehardening_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_pos,gdot_neg, &
   tau_pos,tau_neg

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 nSlip = plastic_kinehardening_totalNslip(instance)
 
 postResults = 0.0_pReal
 c = 0_pInt

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Mp,ph,instance,of) 
 associate( prm => paramNew(instance), stt => state(instance))
 outputsLoop: do o = 1_pInt,plastic_kinehardening_Noutput(instance)
   select case(prm%outputID(o))
     case (crss_ID)
       postResults(c+1_pInt:c+nSlip) = stt%crss(:,of)
       c = c + nSlip
       
     case(crss_back_ID)
       postResults(c+1_pInt:c+nSlip) = stt%crss_back(:,of)
       c = c + nSlip
       
     case (sense_ID)
       postResults(c+1_pInt:c+nSlip) = stt%sense(:,of)
       c = c + nSlip
                                                                        
     case (chi0_ID)
       postResults(c+1_pInt:c+nSlip) = stt%chi0(:,of)
       c = c + nSlip
       
     case (gamma0_ID)
       postResults(c+1_pInt:c+nSlip) = stt%gamma0(:,of)
       c = c + nSlip
     
     case (accshear_ID)
       postResults(c+1_pInt:c+nSlip) = stt%accshear(:,of)
       c = c + nSlip
       
     case (shearrate_ID)
       postResults(c+1_pInt:c+nSlip) = gdot_pos+gdot_neg
       c = c + nSlip

     case (resolvedstress_ID)
       j = 0_pInt
       slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
           j = j + 1_pInt
           postResults(c+j) = &
                             math_mul33xx33(Mp,lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph))
         enddo slipSystems
       enddo slipFamilies
       c = c + nSlip
       
   end select
 enddo outputsLoop
 end associate

end function plastic_kinehardening_postResults


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on slip systems and derivatives with respect to resolved stress
!> @details: Shear rates are calculated only optionally. NOTE: Against the common convention, the
!> result (i.e. intent(out)) variables are the last to have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(prm,stt,of,Mp,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)
 use prec, only: &
  dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 type(tParameters), intent(in) :: &
   prm
 type(tKinehardeningState), intent(in) :: &
   stt
 integer(pInt),     intent(in) :: &
   of
 real(pReal), dimension(prm%totalNslip), intent(out) :: &
   gdot_pos, &
   gdot_neg
 real(pReal), dimension(prm%totalNslip), optional, intent(out) :: &
   dgdot_dtau_pos, &
   dgdot_dtau_neg
 real(pReal), dimension(3,3), intent(in) :: &
   Mp

 real(pReal), dimension(prm%totalNslip) :: &
   tau_pos, &
   tau_neg
 integer(pInt) :: i
 logical       :: nonSchmidActive

 nonSchmidActive = size(prm%nonSchmidCoeff) > 0_pInt

 do i = 1_pInt, prm%totalNslip
   tau_pos(i) =       math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))
   tau_neg(i) = merge(math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,i)), &
                           0.0_pReal, nonSchmidActive)
 enddo

 tau_pos = tau_pos - stt%crss_back(:,of)
 tau_neg = tau_neg - stt%crss_back(:,of)

 where(dNeq0(tau_pos))
   gdot_pos = prm%gdot0 * merge(0.5_pReal,1.0_pReal, nonSchmidActive) &                             ! 1/2 if non-Schmid active
            * sign(abs(tau_pos/stt%crss(:,of))**prm%n_slip,  tau_pos)
 else where
   gdot_pos = 0.0_pReal
 end where

 where(dNeq0(tau_neg))
   gdot_neg = prm%gdot0 * 0.5_pReal &                                                               ! only used if non-Schmid active, always 1/2
            * sign(abs(tau_neg/stt%crss(:,of))**prm%n_slip,  tau_neg)
 else where
   gdot_neg = 0.0_pReal
 end where

 if (present(dgdot_dtau_pos)) then
   where(dNeq0(gdot_pos))
     dgdot_dtau_pos = gdot_pos*prm%n_slip/tau_pos
   else where
     dgdot_dtau_pos = 0.0_pReal
   end where
 endif
 if (present(dgdot_dtau_neg)) then
   where(dNeq0(gdot_neg))
     dgdot_dtau_neg = gdot_neg*prm%n_slip/tau_neg
   else where
     dgdot_dtau_neg = 0.0_pReal
   end where
 endif

end subroutine kinetics

end module plastic_kinehardening
