!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Zhuowen Zhao, Michigan State University
!> @brief Introducing Voce-type kinematic hardening rule into crystal plasticity  
!! formulation using a power law fitting
!--------------------------------------------------------------------------------------------------
module plastic_kinehardening
 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_kinehardening_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_sizePostResult                                                                   !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_output                                                                           !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_kinehardening_Noutput                                                                          !< number of outputs per instance
 
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_kinehardening_totalNslip                                                                       !< no. of slip system used in simulation
  

 integer(pInt),                       dimension(:,:),   allocatable,         private :: &
   plastic_kinehardening_Nslip                                                                            !< active number of slip systems per family (input parameter, per family)
   

 enum, bind(c)
   enumerator :: undefined_ID, &
                 crss_ID, &                                                                               !< critical resolved stress
                 crss_back_ID, &                                                                          !< critical resolved back stress
                 sense_ID, &                                                                              !< sense of acting shear stress (-1 or +1)
                 chi0_ID, &                                                                               !< backstress at last switch of stress sense (positive?)
                 gamma0_ID, &                                                                             !< accumulated shear at last switch of stress sense (at current switch?)
                 accshear_ID, &
                 sumGamma_ID, &
                 shearrate_ID, &
                 resolvedstress_ID
                
 end enum
 
 
 type, private :: tParameters                                                                          !< container type for internal constitutive parameters
   integer(kind(undefined_ID)),         dimension(:),   allocatable,          private :: &
     outputID                                                                                          !< ID of each post result output
    
   real(pReal) :: &
 !     F0, &
!      mu, &
!      mu0, &
!      tau_hat0, &
!      p1, &
!      q1, &
     gdot0, &                                                                                          !< reference shear strain rate for slip (input parameter)
     n_slip, &                                                                                         !< stress exponent for slip (input parameter)
     aTolResistance, &
     aTolShear
     
   
   real(pReal),                         dimension(:),   allocatable,          private :: &
     crss0, &                                                                                           !< initial critical shear stress for slip (input parameter, per family)
     theta0, &                                                                                         !< initial hardening rate of forward stress for each slip
     theta1, &                                                                                         !< asymptotic hardening rate of forward stress for each slip >
     theta0_b, &                                                                                       !< initial hardening rate of back stress for each slip > 
     theta1_b, &                                                                                       !< asymptotic hardening rate of back stress for each slip >
     tau1, &
     tau1_b, &
     interaction_slipslip, &                                                                           !< latent hardening matrix
     nonSchmidCoeff
        
   real(pReal),                       dimension(:,:),   allocatable,          private :: &
     hardeningMatrix_SlipSlip
 end type

 type, private :: tKinehardeningState
   real(pReal), pointer, dimension(:,:) :: &                                                           !< vectors along NipcMyInstance
     crss, &                                                                                           !< critical resolved stress
     crss_back, &                                                                                      !< critical resolved back stress
     sense, &                                                                                          !< sense of acting shear stress (-1 or +1)
     chi0, &                                                                                           !< backstress at last switch of stress sense
     gamma0, &                                                                                         !< accumulated shear at last switch of stress sense
     accshear                                                                                          !< accumulated (absolute) shear

   real(pReal), pointer, dimension(:) :: &                                                             !< scalars along NipcMyInstance
     sumGamma                                                                                          !< accumulated shear across all systems
 end type

 type(tParameters), dimension(:), allocatable, private :: &
   param                                                                                               !< containers of constitutive parameters (len Ninstance)
 
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
   PLASTICITY_kinehardening_label, &
   PLASTICITY_kinehardening_ID, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase
 use lattice
 use numerics,only: &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   o, j, k, f, &
   output_ID, &
   phase, & 
   instance, &
   maxNinstance, &
   NipcMyPhase, &
   Nchunks_SlipSlip = 0_pInt, Nchunks_SlipFamilies = 0_pInt, &
   Nchunks_nonSchmid = 0_pInt, &
   offset_slip, index_myFamily, index_otherFamily, &
   startIndex, endIndex, &
   mySize, nSlip, nSlipFamilies, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState

 real(pReal), dimension(:), allocatable :: tempPerSlip
   
 character(len=65536) :: &
   tag       = '', &
   line      = '', &
   extmsg    = ''
 character(len=64) :: &
   outputtag = ''

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_KINEHARDENING_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == PLASTICITY_KINEHARDENING_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a,1x,i5,/)') '# instances:',maxNinstance                                      
   
 allocate(plastic_kinehardening_sizePostResults(maxNinstance),               source=0_pInt)
 allocate(plastic_kinehardening_sizePostResult(maxval(phase_Noutput),maxNinstance), &
                                                                             source=0_pInt)
 allocate(plastic_kinehardening_output(maxval(phase_Noutput),maxNinstance))
          plastic_kinehardening_output                                             = ''
 allocate(plastic_kinehardening_Noutput(maxNinstance),                       source=0_pInt)
 allocate(plastic_kinehardening_Nslip(lattice_maxNslipFamily,maxNinstance),  source=0_pInt)
 allocate(plastic_kinehardening_totalNslip(maxNinstance),                    source=0_pInt)
 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
   
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
       allocate(param(instance)%outputID(phase_Noutput(phase)), source=undefined_ID)                           ! allocate space for IDs of every requested output
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
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)
       case ('(output)')
         outputtag = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
         output_ID = undefined_ID
         select case(outputtag)
           case ('resistance')
             output_ID = crss_ID
           case ('backstress')
             output_ID = crss_back_ID
           case ('sense')
             output_ID = sense_ID
           case ('chi0')
             output_ID = chi0_ID
           case ('gamma0')
             output_ID = gamma0_ID
           case ('accumulatedshear')
             output_ID = accshear_ID
           case ('totalshear')
             output_ID = sumGamma_ID
           case ('shearrate')
             output_ID = shearrate_ID
           case ('resolvedstress')
             output_ID = resolvedstress_ID
         end select

         if (output_ID /= undefined_ID) then
           plastic_kinehardening_Noutput(instance) = plastic_kinehardening_Noutput(instance) + 1_pInt
           plastic_kinehardening_output(plastic_kinehardening_Noutput(instance),instance) = outputtag
           param(instance)%outputID (plastic_kinehardening_Noutput(instance)) = output_ID
         endif
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
! parameters independent of number of slip families 
      !  case ('F0')
!          param(instance)%F0                    = IO_floatValue(line,chunkPos,2_pInt)
!          
!        case ('mu')
!          param(instance)%mu                    = IO_floatValue(line,chunkPos,2_pInt)
!          
!        case ('mu0')
!          param(instance)%mu0                    = IO_floatValue(line,chunkPos,2_pInt)
!          
!        case ('tau_hat0')
!          param(instance)%tau_hat0                    = IO_floatValue(line,chunkPos,2_pInt)
!          
!        case ('p1')
!          param(instance)%p1                    = IO_floatValue(line,chunkPos,2_pInt)
!          
!        case ('q1')
!          param(instance)%q1                    = IO_floatValue(line,chunkPos,2_pInt)
                                                
       case ('gdot0')
         param(instance)%gdot0                    = IO_floatValue(line,chunkPos,2_pInt)
             
       case ('n_slip')
         param(instance)%n_slip                   = IO_floatValue(line,chunkPos,2_pInt)
             
       case ('atol_resistance')
         param(instance)%aTolResistance           = IO_floatValue(line,chunkPos,2_pInt)
       
       case ('atol_shear')
         param(instance)%aTolShear                = IO_floatValue(line,chunkPos,2_pInt)
               
       case default

     end select
   endif; endif
 enddo parsingFile

!--------------------------------------------------------------------------------------------------
! allocation of variables whose size depends on the total number of active slip systems
 allocate(state(maxNinstance))
 allocate(state0(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(deltaState(maxNinstance))


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
 
     if (any(plastic_kinehardening_Nslip(1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%crss0(1:nSlipFamilies)           < 0.0_pReal)) extmsg = trim(extmsg)//' crss0'
     if (any(plastic_kinehardening_Nslip(1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%tau1(1:nSlipFamilies)          <= 0.0_pReal)) extmsg = trim(extmsg)//' tau1'
     if (any(plastic_kinehardening_Nslip(1:nSlipFamilies,instance) > 0_pInt &
             .and. param(instance)%tau1_b(1:nSlipFamilies)         < 0.0_pReal)) extmsg = trim(extmsg)//' tau1_b'
     if (param(instance)%gdot0            <= 0.0_pReal) extmsg = trim(extmsg)//' gdot0'
     if (param(instance)%n_slip           <= 0.0_pReal) extmsg = trim(extmsg)//' n_slip'
     if (param(instance)%aTolResistance   <= 0.0_pReal) param(instance)%aTolResistance = 1.0_pReal      ! default absolute tolerance 1 Pa
     if (param(instance)%aTolShear        <= 0.0_pReal) param(instance)%aTolShear = 1.0e-6_pReal        ! default absolute tolerance 1e-6                  
     if (extmsg /= '') then 
       extmsg = trim(extmsg)//' ('//PLASTICITY_KINEHARDENING_label//')'                                 ! prepare error message identifier
       call IO_error(211_pInt,ip=instance,ext_msg=extmsg)
     endif
     

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     
     outputsLoop: do o = 1_pInt,plastic_kinehardening_Noutput(instance)
       select case(param(instance)%outputID(o))
         case(crss_ID, &                                                                                           !< critical resolved stress
              crss_back_ID, &                                                                                      !< critical resolved back stress
              sense_ID, &                                                                                          !< sense of acting shear stress (-1 or +1)
              chi0_ID, &                                                                                           !< backstress at last switch of stress sense
              gamma0_ID, &                                                                                         !< accumulated shear at last switch of stress sense
              accshear_ID, &
              shearrate_ID, &
              resolvedstress_ID)
           mySize = nSlip
         case(sumGamma_ID)
           mySize = 1_pInt
         case default
       end select

       outputFound: if (mySize > 0_pInt) then
         plastic_kinehardening_sizePostResult(o,instance) = mySize
         plastic_kinehardening_sizePostResults(instance)  = plastic_kinehardening_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop     
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
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%offsetDeltaState = sizeDotState
     plasticState(phase)%sizePostResults = plastic_kinehardening_sizePostResults(instance)
     plasticState(phase)%nSlip = nSlip
  
     allocate(plasticState(phase)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%aTolState          (sizeDotState),             source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)          ! allocate space for deltaState 
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)

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
     
!----------------------------------------------------------------------------------------------
!locally define dotState alias

     endindex = 0_pInt
     o = endIndex                                                                                           ! offset of dotstate index relative to state index
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state   (instance)%crss          => plasticState(phase)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     state0  (instance)%crss          => plasticState(phase)%state0   (startIndex  :endIndex  ,1:NipcMyPhase)
     dotState(instance)%crss          => plasticState(phase)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)

     state0(instance)%crss                                  = spread(math_expand(param(instance)%crss0,&
                                                                                 plastic_kinehardening_Nslip(:,instance)), &
                                                                     2, NipcMyPhase)
     plasticState(phase)%aTolState(startIndex-o:endIndex-o) = param(instance)%aTolResistance
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip 
     state   (instance)%crss_back          => plasticState(phase)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     state0  (instance)%crss_back          => plasticState(phase)%state0   (startIndex  :endIndex  ,1:NipcMyPhase)
     dotState(instance)%crss_back          => plasticState(phase)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)
     
     state0(instance)%crss_back                             = 0.0_pReal
     plasticState(phase)%aTolState(startIndex-o:endIndex-o) = param(instance)%aTolResistance
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state   (instance)%accshear          => plasticState(phase)%state    (startIndex  :endIndex  ,1:NipcMyPhase)
     state0  (instance)%accshear          => plasticState(phase)%state0   (startIndex  :endIndex  ,1:NipcMyPhase)
     dotState(instance)%accshear          => plasticState(phase)%dotState (startIndex-o:endIndex-o,1:NipcMyPhase)

     state0(instance)%accshear                              = 0.0_pReal
     plasticState(phase)%aTolState(startIndex-o:endIndex-o) = param(instance)%aTolShear
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + 1_pInt
     state   (instance)%sumGamma        => plasticState(phase)%state    (startIndex             ,1:NipcMyPhase)
     state0  (instance)%sumGamma        => plasticState(phase)%state0   (startIndex             ,1:NipcMyPhase)
     dotState(instance)%sumGamma        => plasticState(phase)%dotState (startIndex-o           ,1:NipcMyPhase)

     state0(instance)%sumGamma                              = 0.0_pReal
     plasticState(phase)%aTolState(startIndex-o:endIndex-o) = param(instance)%aTolShear
     
!----------------------------------------------------------------------------------------------
!locally define deltaState alias
     o = endIndex
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state   (instance)%sense          => plasticState(phase)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     state0  (instance)%sense          => plasticState(phase)%state0    (startIndex  :endIndex  ,1:NipcMyPhase)
     deltaState(instance)%sense        => plasticState(phase)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)
     
     state0(instance)%sense          = 0.0_pReal
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state   (instance)%chi0           => plasticState(phase)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     state0  (instance)%chi0           => plasticState(phase)%state0    (startIndex  :endIndex  ,1:NipcMyPhase)
     deltaState(instance)%chi0         => plasticState(phase)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)
     
     state0(instance)%chi0           = 0.0_pReal
     
!    .............................................
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + nSlip
     state   (instance)%gamma0         => plasticState(phase)%state     (startIndex  :endIndex  ,1:NipcMyPhase)         
     state0  (instance)%gamma0         => plasticState(phase)%state0    (startIndex  :endIndex  ,1:NipcMyPhase)
     deltaState(instance)%gamma0       => plasticState(phase)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)

     state0(instance)%gamma0         = 0.0_pReal
     
   endif myPhase2
 enddo initializeInstances

end subroutine plastic_kinehardening_init

!--------------------------------------------------------------------------------------------------
!> @brief calculation of shear rates (\dot \gamma)
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                            Tstar_v,ph,instance,of)

 use lattice, only: &
   lattice_NslipSystem, &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NnonSchmid

 implicit none
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                         !< 2nd Piola Kirchhoff stress tensor in Mandel notation
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
     tau_pos(j) = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_neg(j) = tau_pos(j)
     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph)
       tau_pos(j) = tau_pos(j) + param(instance)%nonSchmidCoeff(k)* &
                                 dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+0,index_myFamily+i,ph))
       tau_neg(j) = tau_neg(j) + param(instance)%nonSchmidCoeff(k)* &
                                 dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems
   enddo slipSystems
 enddo slipFamilies

 gdot_pos = 0.5_pReal * param(instance)%gdot0 * &
            (abs(tau_pos-state(instance)%sense(:,of)*state(instance)%crss_back(:,of))/ &
            state(instance)%crss(:,of))**param(instance)%n_slip &
            *sign(1.0_pReal,tau_pos) 
 gdot_neg = 0.5_pReal * param(instance)%gdot0 * &
            (abs(tau_neg-state(instance)%sense(:,of)*state(instance)%crss_back(:,of))/ &
            state(instance)%crss(:,of))**param(instance)%n_slip &
            *sign(1.0_pReal,tau_neg) 
            
!  gdot_pos = 0.5_pReal * param(instance)%gdot0 * &
!             exp(-param(instance)%F0/(1.38e-23*298.15)* &
!                 (1-((abs(tau_pos-state(instance)%crss_back(:,of)) &
!                       -state(instance)%crss(:,of)*param(instance)%mu/param(instance)%mu) / &
!                     !----------------------------------------------------------------------------
!                       param(instance)%tau_hat0*param(instance)%mu/param(instance)%mu &
!                     )**param(instance)%p1 &
!                 )**param(instance)%q1 &
!                )*sign(1.0_pReal,(tau_pos-state(instance)%crss_back(:,of)))
!             
!             
!              
!  gdot_neg = 0.5_pReal * param(instance)%gdot0 * &
!             exp(-param(instance)%F0/(1.38e-23*298.15)* &
!                 (1-((abs(tau_neg-state(instance)%crss_back(:,of)) &
!                       -state(instance)%crss(:,of)*param(instance)%mu/param(instance)%mu) / &
!                      !---------------------------------------------------------------------------- 
!                       param(instance)%tau_hat0*param(instance)%mu/param(instance)%mu &
!                     )**param(instance)%p1 &
!                 )**param(instance)%q1 &
!                )*sign(1.0_pReal,(tau_neg-state(instance)%crss_back(:,of)))      

  

end subroutine plastic_kinehardening_shearRates


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_LpAndItsTangent(Lp,dLp_dTstar99, &
                                                 Tstar_v,ipc,ip,el)
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
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid
 use material, only: &
   phaseAt, phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation

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
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333                                                                                   !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor

 ph = phaseAt(ipc,ip,el)                                                                            !< figures phase for each material point 
 of = phasememberAt(ipc,ip,el)                                                                      !< index of the positions of each constituent of material point, phasememberAt is a function in material that helps figure them out
 instance = phase_plasticityInstance(ph)

 Lp = 0.0_pReal 
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Tstar_v,ph,instance,of)


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
         dLp_dTstar3333(k,l,m,n) = &
         dLp_dTstar3333(k,l,m,n) + dgdot_dtau_pos*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                  nonSchmid_tensor(m,n,1)
     endif

     if (dNeq0(gdot_neg(j))) then
       dgdot_dtau_neg = gdot_neg(j)*param(instance)%n_slip/(tau_neg(j)-state(instance)%crss_back(j,of))
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = &
         dLp_dTstar3333(k,l,m,n) + dgdot_dtau_neg*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                  nonSchmid_tensor(m,n,2)
     endif
   enddo slipSystems
 enddo slipFamilies

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

end subroutine plastic_kinehardening_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_deltaState(Tstar_v,ipc,ip,el)
 use prec, only: &
   dNeq
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
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric 2nd Piola Kirchhoff stress tensor in Mandel notation
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
                                       Tstar_v,ph,instance,of)

 sense = sign(1.0_pReal,gdot_pos+gdot_neg)                                                          ! current sense of shear direction

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
subroutine plastic_kinehardening_dotState(Tstar_v,ipc,ip,el)
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid
 use material, only: &
   material_phase, &
   phaseAt, phasememberAt, &
   plasticState, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation, vector form
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element !< microstructure state

 integer(pInt) :: &
   instance,ph, &
   f,i,j,k, &
   index_myFamily,index_otherFamily, &
   nSlip, &
   offset_accshear, &
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
                                       Tstar_v,ph,instance,of)
                                       
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
     j = j+1_pInt    
     dotState(instance)%crss(j,of) = &                                                                        ! evolution of slip resistance j
          dot_product(param(instance)%hardeningMatrix_SlipSlip(j,1:nSlip),abs(gdot_pos+gdot_neg)) * &
          ( param(instance)%theta1(f) + &
           (param(instance)%theta0(f) - param(instance)%theta1(f) &
            + param(instance)%theta0(f)*param(instance)%theta1(f)*state(instance)%sumGamma(of)/param(instance)%tau1(f)) &
           *exp(-state(instance)%sumGamma(of)*param(instance)%theta0(f)/param(instance)%tau1(f)) &                  ! V term depending on the harding law
          )
     dotState(instance)%crss_back(j,of) = &                                                                   ! evolution of back stress resistance j
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
function plastic_kinehardening_postResults(Tstar_v,ipc,ip,el)
 use material, only: &
   material_phase, &
   plasticState, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_NnonSchmid

 implicit none
 real(pReal), dimension(6), intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(plastic_kinehardening_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   plastic_kinehardening_postResults

 integer(pInt) :: &
   instance,ph, of, &
   nSlip,&
   o,f,i,c,j,k, &
   index_myFamily
   
 real(pReal), dimension(plastic_kinehardening_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_pos,gdot_neg, &
   tau_pos,tau_neg

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 nSlip = plastic_kinehardening_totalNslip(instance)
 
 plastic_kinehardening_postResults = 0.0_pReal
 c = 0_pInt

 call plastic_kinehardening_shearRates(gdot_pos,gdot_neg,tau_pos,tau_neg, &
                                       Tstar_v,ph,instance,of) 

 outputsLoop: do o = 1_pInt,plastic_kinehardening_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (crss_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%crss(:,of)
       c = c + nSlip
       
     case(crss_back_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%crss_back(:,of)
       c = c + nSlip
       
     case (sense_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%sense(:,of)
       c = c + nSlip
                                                                        
     case (chi0_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%chi0(:,of)
       c = c + nSlip
       
     case (gamma0_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%gamma0(:,of)
       c = c + nSlip
     
     case (accshear_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = state(instance)%accshear(:,of)
       c = c + nSlip
       
     case (sumGamma_ID)
       plastic_kinehardening_postResults(c+1_pInt) = state(instance)%sumGamma(of)
       c = c + 1_pInt  
       
     case (shearrate_ID)
       plastic_kinehardening_postResults(c+1_pInt:c+nSlip) = gdot_pos+gdot_neg
       c = c + nSlip

     case (resolvedstress_ID)
       j = 0_pInt
       slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems: do i = 1_pInt,plastic_kinehardening_Nslip(f,instance)
           j = j + 1_pInt
           plastic_kinehardening_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
         enddo slipSystems
       enddo slipFamilies
       c = c + nSlip
       
   end select
 enddo outputsLoop

end function plastic_kinehardening_postResults

end module plastic_kinehardening
