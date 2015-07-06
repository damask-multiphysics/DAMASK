!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating isotropic brittle damage source mechanism
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_damage_isoBrittle
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   source_damage_isoBrittle_sizePostResults, &                                                        !< cumulative size of post results
   source_damage_isoBrittle_offset, &                                                                 !< which source is my current damage mechanism?
   source_damage_isoBrittle_instance                                                                  !< instance of damage source mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   source_damage_isoBrittle_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   source_damage_isoBrittle_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   source_damage_isoBrittle_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),     allocatable,         private :: &
   source_damage_isoBrittle_aTol, &
   source_damage_isoBrittle_critStrainEnergy

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 damage_drivingforce_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   source_damage_isoBrittle_outputID                                                                  !< ID of each post result output


 public :: &
   source_damage_isoBrittle_init, &
   source_damage_isoBrittle_deltaState, &
   source_damage_isoBrittle_getRateAndItsTangent, &
   source_damage_isoBrittle_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_init(fileUnit)
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
   SOURCE_damage_isoBrittle_label, &
   SOURCE_damage_isoBrittle_ID, &
   material_Nphase, &
   material_phase, &  
   sourceState, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,source,sourceOffset,o
 integer(pInt) :: sizeState, sizeDotState, sizeDeltaState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  source_'//SOURCE_damage_isoBrittle_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_source == SOURCE_damage_isoBrittle_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(source_damage_isoBrittle_offset(material_Nphase), source=0_pInt)
 allocate(source_damage_isoBrittle_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   source_damage_isoBrittle_instance(phase) = count(phase_source(:,1:phase) == source_damage_isoBrittle_ID)
   do source = 1, phase_Nsources(phase)
     if (phase_source(source,phase) == source_damage_isoBrittle_ID) &
       source_damage_isoBrittle_offset(phase) = source
   enddo    
 enddo
   
 allocate(source_damage_isoBrittle_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(source_damage_isoBrittle_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(source_damage_isoBrittle_output(maxval(phase_Noutput),maxNinstance))
          source_damage_isoBrittle_output = ''
 allocate(source_damage_isoBrittle_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(source_damage_isoBrittle_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(source_damage_isoBrittle_critStrainEnergy(maxNinstance),                    source=0.0_pReal) 
 allocate(source_damage_isoBrittle_aTol(maxNinstance),                                source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (any(phase_source(:,phase) == SOURCE_damage_isoBrittle_ID)) then   ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = source_damage_isoBrittle_instance(phase)                                            ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('isobrittle_drivingforce')
             source_damage_isoBrittle_Noutput(instance) = source_damage_isoBrittle_Noutput(instance) + 1_pInt
             source_damage_isoBrittle_outputID(source_damage_isoBrittle_Noutput(instance),instance) = damage_drivingforce_ID
             source_damage_isoBrittle_output(source_damage_isoBrittle_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('isobrittle_criticalstrainenergy')
         source_damage_isoBrittle_critStrainEnergy(instance) = IO_floatValue(line,positions,2_pInt)

       case ('isobrittle_atol')
         source_damage_isoBrittle_aTol(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile


!--------------------------------------------------------------------------------------------------
!  sanity checks
 sanityChecks: do phase = 1_pInt, material_Nphase 
   myPhase: if (any(phase_source(:,phase) == SOURCE_damage_isoBrittle_ID)) then
     instance = source_damage_isoBrittle_instance(phase)
     if (source_damage_isoBrittle_aTol(instance) < 0.0_pReal) &
       source_damage_isoBrittle_aTol(instance) = 1.0e-3_pReal                                              ! default absolute tolerance 1e-3
     if (source_damage_isoBrittle_critStrainEnergy(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='criticalStrainEnergy ('//SOURCE_damage_isoBrittle_LABEL//')')
   endif myPhase
 enddo sanityChecks

 initializeInstances: do phase = 1_pInt, material_Nphase
   if (any(phase_source(:,phase) == SOURCE_damage_isoBrittle_ID)) then
     NofMyPhase=count(material_phase==phase)
     instance = source_damage_isoBrittle_instance(phase)
     sourceOffset = source_damage_isoBrittle_offset(phase)
!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,source_damage_isoBrittle_Noutput(instance)
       select case(source_damage_isoBrittle_outputID(o,instance))
         case(damage_drivingforce_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          source_damage_isoBrittle_sizePostResult(o,instance) = mySize
          source_damage_isoBrittle_sizePostResults(instance)  = source_damage_isoBrittle_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   1_pInt
     sizeDeltaState            =   1_pInt
     sizeState                 =   1_pInt
                
     sourceState(phase)%p(sourceOffset)%sizeState      = sizeState
     sourceState(phase)%p(sourceOffset)%sizeDotState   = sizeDotState
     sourceState(phase)%p(sourceOffset)%sizeDeltaState = sizeDeltaState
     sourceState(phase)%p(sourceOffset)%sizePostResults = source_damage_isoBrittle_sizePostResults(instance)
     allocate(sourceState(phase)%p(sourceOffset)%aTolState           (sizeState),                &
              source=source_damage_isoBrittle_aTol(instance))
     allocate(sourceState(phase)%p(sourceOffset)%state0              (sizeState,NofMyPhase),     source=.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(sourceState(phase)%p(sourceOffset)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%deltaState        (sizeDeltaState,NofMyPhase),  source=0.0_pReal)
     allocate(sourceState(phase)%p(sourceOffset)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
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
end subroutine source_damage_isoBrittle_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_deltaState(C, Fe, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   sourceState, &
   material_homog, &
   phase_NstiffnessDegradations, &
   phase_stiffnessDegradation, &
   porosity, &
   porosityMapping, &  
   STIFFNESS_DEGRADATION_porosity_ID
 use math, only : &
   math_mul33x33, &
   math_mul66x6, &
   math_Mandel33to6, &
   math_transpose33, &
   math_I3

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(3,3) :: &
   Fe
 real(pReal),  intent(in), dimension(6,6) :: &
   C
 integer(pInt) :: &
   phase, constituent, instance, sourceOffset, mech
 real(pReal) :: &
   strain(6), &
   stiffness(6,6), &
   strainenergy

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = source_damage_isoBrittle_instance(phase)
 sourceOffset = source_damage_isoBrittle_offset(phase)

 stiffness = C
 do mech = 1_pInt, phase_NstiffnessDegradations(phase)
   select case(phase_stiffnessDegradation(mech,phase))
     case (STIFFNESS_DEGRADATION_porosity_ID)
       stiffness = porosity(material_homog(ip,el))%p(porosityMapping(material_homog(ip,el))%p(ip,el))* &
                   porosity(material_homog(ip,el))%p(porosityMapping(material_homog(ip,el))%p(ip,el))* &
                   stiffness 
   end select
 enddo                                                 
 strain = 0.5_pReal*math_Mandel33to6(math_mul33x33(math_transpose33(Fe),Fe)-math_I3)

 strainenergy = 2.0_pReal*sum(strain*math_mul66x6(stiffness,strain))/ &
                source_damage_isoBrittle_critStrainEnergy(instance) 
 if (strainenergy > sourceState(phase)%p(sourceOffset)%subState0(1,constituent)) then
   sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
     strainenergy - sourceState(phase)%p(sourceOffset)%state(1,constituent)
 else
   sourceState(phase)%p(sourceOffset)%deltaState(1,constituent) = &
     sourceState(phase)%p(sourceOffset)%subState0(1,constituent) - &
     sourceState(phase)%p(sourceOffset)%state(1,constituent)
 endif
 
end subroutine source_damage_isoBrittle_deltaState
 
!--------------------------------------------------------------------------------------------------
!> @brief returns local part of nonlocal damage driving force
!--------------------------------------------------------------------------------------------------
subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   sourceState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in) :: &
   phi
 real(pReal),  intent(out) :: &
   localphiDot, &
   dLocalphiDot_dPhi
 integer(pInt) :: &
   phase, constituent, sourceOffset

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 sourceOffset = source_damage_isoBrittle_offset(phase)
 
 localphiDot = 1.0_pReal - phi*sourceState(phase)%p(sourceOffset)%state(1,constituent)
 dLocalphiDot_dPhi = -sourceState(phase)%p(sourceOffset)%state(1,constituent)
 
end subroutine source_damage_isoBrittle_getRateAndItsTangent
 
!--------------------------------------------------------------------------------------------------
!> @brief return array of local damage results
!--------------------------------------------------------------------------------------------------
function source_damage_isoBrittle_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   sourceState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(source_damage_isoBrittle_sizePostResults( &
                         source_damage_isoBrittle_instance(mappingConstitutive(2,ipc,ip,el)))) :: &
   source_damage_isoBrittle_postResults

 integer(pInt) :: &
   instance, phase, constituent, sourceOffset, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = source_damage_isoBrittle_instance(phase)
 sourceOffset = source_damage_isoBrittle_offset(phase)

 c = 0_pInt
 source_damage_isoBrittle_postResults = 0.0_pReal

 do o = 1_pInt,source_damage_isoBrittle_Noutput(instance)
    select case(source_damage_isoBrittle_outputID(o,instance))
      case (damage_drivingforce_ID)
        source_damage_isoBrittle_postResults(c+1_pInt) = sourceState(phase)%p(sourceOffset)%state(1,constituent)
        c = c + 1

    end select
 enddo
end function source_damage_isoBrittle_postResults

end module source_damage_isoBrittle
