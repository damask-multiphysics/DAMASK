!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic ductile damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoDuctile
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_damage_isoDuctile_sizePostResults, &                                                        !< cumulative size of post results
   source_damage_isoDuctile_offset, &                                                                 !< which source is my current damage mechanism?
   source_damage_isoDuctile_instance                                                                  !< instance of damage source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_damage_isoDuctile_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_damage_isoDuctile_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_damage_isoDuctile_Noutput                                                                   !< number of outputs per instance of this damage 


 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_drivingforce_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   source_damage_isoDuctile_outputID                                                                  !< ID of each post result output


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     critPlasticStrain, &
     N, &
     aTol
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID
 end type tParameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)


 public :: &
   source_damage_isoDuctile_init, &
   source_damage_isoDuctile_dotState, &
   source_damage_isoDuctile_getRateAndItsTangent, &
   source_damage_isoDuctile_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_init(fileUnit)
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pStringLen
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
   material_allocateSourceState, &
   phase_source, &
   phase_Nsources, &
   phase_Noutput, &
   SOURCE_damage_isoDuctile_label, &
   SOURCE_damage_isoDuctile_ID, &
   material_phase, &  
   sourceState
 use config, only: &
   config_phase, &
   material_Nphase, &
   MATERIAL_partPhase

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: Ninstance,mySize=0_pInt,phase,instance,source,sourceOffset,o
 integer(pInt) :: NofMyPhase,p,i
 character(len=65536) :: &
   tag  = '', &
   line = ''
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_DAMAGE_ISODUCTILE_LABEL//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 Ninstance = int(count(phase_source == SOURCE_damage_isoDuctile_ID),pInt)
 if (Ninstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
 allocate(source_damage_isoDuctile_offset(material_Nphase), source=0_pInt)
 allocate(source_damage_isoDuctile_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_damage_isoDuctile_instance(phase) = count(phase_source(:,1:phase) == source_damage_isoDuctile_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_damage_isoDuctile_ID) &
       source_damage_isoDuctile_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_damage_isoDuctile_sizePostResults(Ninstance),                     source=0_pInt)
 allocate(source_damage_isoDuctile_sizePostResult(maxval(phase_Noutput),Ninstance),source=0_pInt)
 allocate(source_damage_isoDuctile_output(maxval(phase_Noutput),Ninstance))
          source_damage_isoDuctile_output = ''
 allocate(source_damage_isoDuctile_outputID(maxval(phase_Noutput),Ninstance),      source=undefined_ID)
 allocate(source_damage_isoDuctile_Noutput(Ninstance),                             source=0_pInt)

 allocate(param(Ninstance))
 
 do p=1, size(config_phase)
   if (all(phase_source(:,p) /= SOURCE_DAMAGE_ISODUCTILE_ID)) cycle
   associate(prm => param(source_damage_isoDuctile_instance(p)), &
             config => config_phase(p))
             
   prm%aTol              = config%getFloat('isoductile_atol',defaultVal = 1.0e-3_pReal)

   prm%N                 = config%getFloat('isoductile_ratesensitivity')
   prm%critPlasticStrain = config%getFloat('isoductile_criticalplasticstrain')
   
   ! sanity checks
   if (prm%aTol                 < 0.0_pReal) extmsg = trim(extmsg)//' isoductile_atol'
   
   if (prm%N                   <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_ratesensitivity'
   if (prm%critPlasticStrain   <= 0.0_pReal) extmsg = trim(extmsg)//' isoductile_criticalplasticstrain'
   
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211_pInt,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ISODUCTILE_LABEL//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case ('isoductile_drivingforce')
       
     end select

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
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (any(phase_source(:,phase) == SOURCE_damage_isoDuctile_ID)) then   ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = source_damage_isoDuctile_instance(phase)                                            ! which instance of my damage is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,chunkPos,2_pInt)))
           case ('isoductile_drivingforce')
             source_damage_isoDuctile_Noutput(instance) = source_damage_isoDuctile_Noutput(instance) + 1_pInt
             source_damage_isoDuctile_outputID(source_damage_isoDuctile_Noutput(instance),instance) = damage_drivingforce_ID
             source_damage_isoDuctile_output(source_damage_isoDuctile_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,chunkPos,2_pInt))
          end select

     end select
   endif; endif
 enddo parsingFile

 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_damage_isoDuctile_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_damage_isoDuctile_instance(phase)
     sourceOffset = source_damage_isoDuctile_offset(phase)
!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,source_damage_isoDuctile_Noutput(instance)
       select case(source_damage_isoDuctile_outputID(o,instance))
         case(damage_drivingforce_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          source_damage_isoDuctile_sizePostResult(o,instance) = mySize
          source_damage_isoDuctile_sizePostResults(instance)  = source_damage_isoDuctile_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop

     call material_allocateSourceState(phase,sourceOffset,NofMyPhase,1_pInt)
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_damage_isoDuctile_sizePostResults(instance)
     sourceState(phase)%p(sourceOffset)%aTolState=param(instance)%aTol
              
   endif
 
 enddo initializeInstances
end subroutine source_damage_isoDuctile_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_dotState(ipc, ip, el)
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   sourceState, &
   material_homog, &
   damage, &
   damageMapping

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, constituent, instance, homog, sourceOffset, damageOffset

 phase = phaseAt(ipc,ip,el)
 constituent = phasememberAt(ipc,ip,el)
 instance = source_damage_isoDuctile_instance(phase)
 sourceOffset = source_damage_isoDuctile_offset(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)

 sourceState(phase)%p(sourceOffset)%dotState(1,constituent) = &
   sum(plasticState(phase)%slipRate(:,constituent))/ &
   ((damage(homog)%p(damageOffset))**param(instance)%N)/ & 
   param(instance)%critPlasticStrain 

end subroutine source_damage_isoDuctile_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)
 use material, only: &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   constituent
 real(pReal),  intent(in) :: &
   phi
 real(pReal),  intent(out) :: &
   localphiDot, &
   dLocalphiDot_dPhi
 integer(pInt) :: &
   sourceOffset

 sourceOffset = source_damage_isoDuctile_offset(phase)
 
 localphiDot = 1.0_pReal - &
               sourceState(phase)%p(sourceOffset)%state(1,constituent)* &
               phi
 
 dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_isoDuctile_getRateAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of local damage results
!--------------------------------------------------------------------------------------------------
function source_damage_isoDuctile_postResults(phase, constituent)
 use material, only: &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   phase, &
   constituent
 real(pReal), dimension(source_damage_isoDuctile_sizePostResults( &
                          source_damage_isoDuctile_instance(phase))) :: &
   source_damage_isoDuctile_postResults

 integer(pInt) :: &
   instance, sourceOffset, o, c
   
 instance = source_damage_isoDuctile_instance(phase)
 sourceOffset = source_damage_isoDuctile_offset(phase)

 c = 0_pInt
 source_damage_isoDuctile_postResults = 0.0_pReal

 do o = 1_pInt,source_damage_isoDuctile_Noutput(instance)
    select case(source_damage_isoDuctile_outputID(o,instance))
      case (damage_drivingforce_ID)
        source_damage_isoDuctile_postResults(c+1_pInt) = sourceState(phase)%p(sourceOffset)%state(1,constituent)
        c = c + 1

    end select
 enddo
end function source_damage_isoDuctile_postResults

end module source_damage_isoDuctile
