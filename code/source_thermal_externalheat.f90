!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_thermal_externalheat
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_thermal_externalheat_sizePostResults, &                                                        !< cumulative size of post results
   source_thermal_externalheat_offset, &                                                                 !< which source is my current thermal dissipation mechanism?
   source_thermal_externalheat_instance                                                                  !< instance of thermal dissipation source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_thermal_externalheat_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_thermal_externalheat_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_thermal_externalheat_Noutput                                                                   !< number of outputs per instance of this source 

 integer(pInt),                       dimension(:),           allocatable,        private :: &
   source_thermal_externalheat_nIntervals

 real(pReal),                         dimension(:,:),         allocatable,        private :: &
   source_thermal_externalheat_time, &
   source_thermal_externalheat_rate

 public :: &
   source_thermal_externalheat_init, &
   source_thermal_externalheat_dotState, &
   source_thermal_externalheat_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
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
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_thermal_externalheat_label, &
   SOURCE_thermal_externalheat_ID, &
   material_Nphase, &
   material_phase, &  
   sourceState, &
   MATERIAL_partPhase
 use numerics,only: &
   analyticJaco, &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,source,sourceOffset
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase,interval   
 character(len=65536) :: &
   tag  = '', &
   line = ''
 real(pReal), allocatable, dimension(:,:) :: temp_time, temp_rate  

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_thermal_externalheat_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_source == SOURCE_thermal_externalheat_ID),pInt)
 if (maxNinstance == 0_pInt) return
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_thermal_externalheat_offset(material_Nphase), source=0_pInt)
 allocate(source_thermal_externalheat_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_thermal_externalheat_instance(phase) = count(phase_source(:,1:phase) == SOURCE_thermal_externalheat_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == SOURCE_thermal_externalheat_ID) &
       source_thermal_externalheat_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_thermal_externalheat_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_thermal_externalheat_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_thermal_externalheat_output  (maxval(phase_Noutput),maxNinstance))
          source_thermal_externalheat_output = ''
 allocate(source_thermal_externalheat_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(source_thermal_externalheat_nIntervals(maxNinstance),                          source=0_pInt) 

 allocate(temp_time(maxNinstance,1000), source=0.0_pReal) 
 allocate(temp_rate(maxNinstance,1000), source=0.0_pReal) 

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
     cycle                                                                                          ! skip to next line
   endif

   if (phase > 0_pInt ) then; if (any(phase_source(:,phase) == SOURCE_thermal_externalheat_ID)) then ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran

     instance = source_thermal_externalheat_instance(phase)                                          ! which instance of my source is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('externalheat_time')
         if (chunkPos(1) <= 2_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//SOURCE_thermal_externalheat_label//')')
         source_thermal_externalheat_nIntervals(instance) = chunkPos(1) - 2_pInt
         do interval = 1, source_thermal_externalheat_nIntervals(instance) + 1_pInt
           temp_time(instance, interval) = IO_floatValue(line,chunkPos,1_pInt + interval)
         enddo  

       case ('externalheat_rate')
         do interval = 1, source_thermal_externalheat_nIntervals(instance) + 1_pInt
           temp_rate(instance, interval) = IO_floatValue(line,chunkPos,1_pInt + interval)
         enddo  

     end select
   endif; endif
 enddo parsingFile
 
 allocate(source_thermal_externalheat_time(maxNinstance,maxval(source_thermal_externalheat_nIntervals)+1_pInt), source=0.0_pReal) 
 allocate(source_thermal_externalheat_rate(maxNinstance,maxval(source_thermal_externalheat_nIntervals)+1_pInt), source=0.0_pReal) 

 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_thermal_externalheat_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_thermal_externalheat_instance(phase)
     sourceOffset = source_thermal_externalheat_offset(phase)
     source_thermal_externalheat_time(instance,1:source_thermal_externalheat_nIntervals(instance)+1_pInt) = &
       temp_time(instance,1:source_thermal_externalheat_nIntervals(instance)+1_pInt)
     source_thermal_externalheat_rate(instance,1:source_thermal_externalheat_nIntervals(instance)+1_pInt) = &
       temp_rate(instance,1:source_thermal_externalheat_nIntervals(instance)+1_pInt)

     sizeDotState              =   1_pInt
     sizeDeltaState            =   0_pInt
     sizeState                 =   1_pInt
     sourceState(phase)%p(sourceOffset)%sizeState =       sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState =    sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState =  sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_thermal_externalheat_sizePostResults(instance)
     allocate(sourceState(phase)%p(sourceOffset)%aTolState           (sizeState),                source=0.00001_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state               (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(sourceState(phase)%p(sourceOffset)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     if (.not. analyticJaco) then
       allocate(sourceState(phase)%p(sourceOffset)%state_backup      (sizeState,NofMyPhase),     source=0.0_pReal)
       allocate(sourceState(phase)%p(sourceOffset)%dotState_backup   (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 1_pInt)) then
       allocate(sourceState(phase)%p(sourceOffset)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(sourceState(phase)%p(sourceOffset)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(sourceState(phase)%p(sourceOffset)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(sourceState(phase)%p(sourceOffset)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)      

   endif
 
 enddo initializeInstances
end subroutine source_thermal_externalheat_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_dotState(ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, &
   constituent, &
   sourceOffset

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 sourceOffset = source_thermal_externalheat_offset(phase)
 
 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = 1.0_pReal

end subroutine source_thermal_externalheat_dotState

!--------------------------------------------------------------------------------------------------
!> @brief returns local vacancy generation rate 
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_externalheat_getRateAndItsTangent(TDot, dTDot_dT, ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(out) :: &
   TDot, &
   dTDot_dT
 integer(pInt) :: &
   instance, phase, constituent, sourceOffset, interval
 real(pReal) :: &
   norm_time   

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 instance = source_thermal_externalheat_instance(phase)
 sourceOffset = source_thermal_externalheat_offset(phase)

 do interval = 1, source_thermal_externalheat_nIntervals(instance)
   norm_time = (sourceState(phase)%p(sourceOffset)%state(1,constituent) - &
                source_thermal_externalheat_time(instance,interval)) / &
               (source_thermal_externalheat_time(instance,interval+1) - &
                source_thermal_externalheat_time(instance,interval))
   if (norm_time >= 0.0_pReal .and. norm_time < 1.0_pReal) &
     TDot = source_thermal_externalheat_rate(instance,interval  ) * (1.0_pReal - norm_time) + &
            source_thermal_externalheat_rate(instance,interval+1) * norm_time
 enddo
 dTDot_dT = 0.0
 
end subroutine source_thermal_externalheat_getRateAndItsTangent
 
end module source_thermal_externalheat
