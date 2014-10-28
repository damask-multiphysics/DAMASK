!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating ductile damage
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_isoDuctile
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_isoDuctile_sizePostResults                                                                   !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_isoDuctile_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_isoDuctile_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   damage_isoDuctile_Noutput                                                                           !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_isoDuctile_aTol, &
   damage_isoDuctile_critpStrain

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_isoDuctile_outputID                                                                  !< ID of each post result output


 public :: &
   damage_isoDuctile_init, &
   damage_isoDuctile_stateInit, &
   damage_isoDuctile_aTolState, &
   damage_isoDuctile_dotState, &
   damage_isoDuctile_microstructure, &
   damage_isoDuctile_getDamage, &
   damage_isoDuctile_getSlipDamage, &
   damage_isoDuctile_putLocalDamage, &
   damage_isoDuctile_getLocalDamage, &
   damage_isoDuctile_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
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
   phase_damage, &
   phase_damageInstance, &
   phase_Noutput, &
   LOCAL_damage_isoDuctile_label, &
   LOCAL_damage_isoDuctile_ID, &
   material_phase, &  
   damageState, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank, &
   numerics_integrator

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,o
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_damage_isoDuctile_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_damage == LOCAL_damage_isoDuctile_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_isoDuctile_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_isoDuctile_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_isoDuctile_output(maxval(phase_Noutput),maxNinstance))
          damage_isoDuctile_output = ''
 allocate(damage_isoDuctile_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_isoDuctile_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_isoDuctile_critpStrain(maxNinstance),                         source=0.0_pReal) 
 allocate(damage_isoDuctile_aTol(maxNinstance),                                source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_damage_isoDuctile_ID) then              ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_isoDuctile_Noutput(instance) = damage_isoDuctile_Noutput(instance) + 1_pInt
             damage_isoDuctile_outputID(damage_isoDuctile_Noutput(instance),instance) = local_damage_ID
             damage_isoDuctile_output(damage_isoDuctile_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('critical_plastic_strain')
         damage_isoDuctile_critpStrain(instance) = IO_floatValue(line,positions,2_pInt)

       case ('atol_damage')
         damage_isoDuctile_aTol(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_damage_isoDuctile_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_isoDuctile_Noutput(instance)
       select case(damage_isoDuctile_outputID(o,instance))
         case(local_damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_isoDuctile_sizePostResult(o,instance) = mySize
          damage_isoDuctile_sizePostResults(instance)  = damage_isoDuctile_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   1_pInt
     sizeState                 =   2_pInt
                
     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_isoDuctile_sizePostResults(instance)
     allocate(damageState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(damageState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(damageState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(damageState(phase)%deltaState          (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(damageState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(damageState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(damageState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(damageState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(damageState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)

     call damage_isoDuctile_stateInit(phase)
     call damage_isoDuctile_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_isoDuctile_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_stateInit(phase)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1) = 1.0_pReal
 tempState(2) = 1.0_pReal


 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_isoDuctile_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_isoDuctile_aTol(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_isoDuctile_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_dotState(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState
 use math, only: &
   math_norm33
 use lattice, only: &
   lattice_DamageMobility

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, constituent

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 
 damageState(phase)%dotState(1,constituent) = &
   (1.0_pReal/lattice_DamageMobility(phase))* &
   (damageState(phase)%state(2,constituent) - &
    damageState(phase)%state(1,constituent))
  
end subroutine damage_isoDuctile_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_microstructure(nSlip,accumulatedSlip,ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_Mandel6to33, &
   math_mul33x33, &
   math_transpose33, &
   math_I3, &
   math_norm33

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(nSlip), intent(in) :: &
   accumulatedSlip
 integer(pInt) :: &
   phase, constituent, instance

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 damageState(phase)%state(2,constituent) = min(damageState(phase)%state(2,constituent), &
                                                     damage_isoDuctile_critpStrain(instance)/ &
                                                     sum(accumulatedSlip))                         !< akin to damage surface
                                                       
end subroutine damage_isoDuctile_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief returns damage 
!--------------------------------------------------------------------------------------------------
function damage_isoDuctile_getDamage(ipc, ip, el)
 use material, only: &
   material_homog, &
   mappingHomogenization, &
   fieldDamage, &
   field_damage_type, &
   FIELD_DAMAGE_LOCAL_ID, &
   FIELD_DAMAGE_NONLOCAL_ID

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: damage_isoDuctile_getDamage
 
 select case(field_damage_type(material_homog(ip,el)))                                                   
   case (FIELD_DAMAGE_LOCAL_ID)
    damage_isoDuctile_getDamage = damage_isoDuctile_getLocalDamage(ipc, ip, el)
    
   case (FIELD_DAMAGE_NONLOCAL_ID)
    damage_isoDuctile_getDamage =    fieldDamage(material_homog(ip,el))% &
      field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

 end select
 
end function damage_isoDuctile_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns slip damage 
!--------------------------------------------------------------------------------------------------
function damage_isoDuctile_getSlipDamage(ipc, ip, el)

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: damage_isoDuctile_getSlipDamage, damage
 
 damage = damage_isoDuctile_getDamage(ipc, ip, el)
 damage_isoDuctile_getSlipDamage = damage*damage
 
end function damage_isoDuctile_getSlipDamage

!--------------------------------------------------------------------------------------------------
!> @brief puts local damage 
!--------------------------------------------------------------------------------------------------
subroutine damage_isoDuctile_putLocalDamage(ipc, ip, el, localDamage)
 use material, only: &
   mappingConstitutive, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: localDamage
 
 damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el)) = &
   localDamage
 
end subroutine damage_isoDuctile_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage
!--------------------------------------------------------------------------------------------------
function damage_isoDuctile_getLocalDamage(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: damage_isoDuctile_getLocalDamage
 
 damage_isoDuctile_getLocalDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function damage_isoDuctile_getLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_isoDuctile_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_isoDuctile_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_isoDuctile_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_isoDuctile_postResults = 0.0_pReal

 do o = 1_pInt,damage_isoDuctile_Noutput(instance)
    select case(damage_isoDuctile_outputID(o,instance))
      case (local_damage_ID)
        damage_isoDuctile_postResults(c+1_pInt) = damageState(phase)%state(1,constituent)
        c = c + 1

    end select
 enddo
end function damage_isoDuctile_postResults

end module damage_isoDuctile
