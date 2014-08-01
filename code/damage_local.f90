!--------------------------------------------------------------------------------------------------
! $Id: damage_local.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_local
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_local_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_local_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_local_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   damage_local_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),     allocatable,         public :: &
   damage_local_crack_mobility, &
   damage_local_aTol

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_local_outputID                                                                  !< ID of each post result output


 public :: &
   damage_local_init, &
   damage_local_stateInit, &
   damage_local_aTolState, &
   damage_local_dotState, &
   damage_local_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_local_init(fileUnit)
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
   DAMAGE_LOCAL_label, &
   DAMAGE_local_ID, &
   material_phase, &  
   damageState, &
   MATERIAL_partPhase
 use numerics,only: &
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
  
 write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_LOCAL_label//' init  -+>>>'
 write(6,'(a)')     ' $Id: damage_local.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_damage == DAMAGE_local_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_local_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_local_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_local_output(maxval(phase_Noutput),maxNinstance))
          damage_local_output = ''
 allocate(damage_local_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_local_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_local_crack_mobility(maxNinstance),                      source=0.0_pReal) 
 allocate(damage_local_aTol(maxNinstance),                                source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == DAMAGE_local_ID) then               ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                     ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_local_Noutput(instance) = damage_local_Noutput(instance) + 1_pInt
             damage_local_outputID(damage_local_Noutput(instance),instance) = local_damage_ID
             damage_local_output(damage_local_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('crack_mobility')
         damage_local_crack_mobility(instance) = IO_floatValue(line,positions,2_pInt)

       case ('atol_damage')
         damage_local_aTol(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == DAMAGE_local_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_local_Noutput(instance)
       select case(damage_local_outputID(o,instance))
         case(local_damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_local_sizePostResult(o,instance) = mySize
          damage_local_sizePostResults(instance)  = damage_local_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   2_pInt
     sizeState                 =   2_pInt
                
     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_local_sizePostResults(instance)
     allocate(damageState(phase)%aTolState           (sizeState),                source=0.0_pReal)
     allocate(damageState(phase)%state0              (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%partionedState0     (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%subState0           (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%state               (sizeState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%state_backup        (sizeState,NofMyPhase),     source=0.0_pReal)

     allocate(damageState(phase)%dotState            (sizeDotState,NofMyPhase),  source=0.0_pReal)
     allocate(damageState(phase)%deltaState          (sizeDotState,NofMyPhase),     source=0.0_pReal)
     allocate(damageState(phase)%dotState_backup     (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(damageState(phase)%previousDotState  (sizeDotState,NofMyPhase),  source=0.0_pReal)
       allocate(damageState(phase)%previousDotState2 (sizeDotState,NofMyPhase),  source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(damageState(phase)%RK4dotState       (sizeDotState,NofMyPhase),  source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(damageState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase),source=0.0_pReal)

     call damage_local_stateInit(phase,instance)
     call damage_local_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_local_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_local_stateInit(phase,instance)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the damage
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1) = 0.0_pReal
 tempState(2) = 1.0_pReal
 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_local_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_local_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_local_aTol(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_local_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_local_dotState(Tstar_v, Fe, Lp, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_Mandel66to3333, &
   math_mul33x33, &
   math_mul3333xx33, &
   math_transpose33, &
   math_trace33, &
   math_I3
 use lattice, only: &
   lattice_surfaceEnergy33, &
   lattice_C66

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp, &
   Fe
 integer(pInt) :: &
   phase, constituent, instance
 real(pReal) :: &
   trialDamage, strain(3,3)

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 strain = 0.5_pReal*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3)
 trialDamage = min(1.0_pReal, &
                   (math_trace33(lattice_surfaceEnergy33(1:3,1:3,phase))/3.0_pReal)/ & 
                   (abs(sum(strain*math_mul3333xx33(math_Mandel66to3333(lattice_C66(1:6,1:6,phase)),strain))) + &
                    damageState(phase)%state(1,constituent)))
 
 damageState(phase)%dotState(1,constituent) = &
   0.0_pReal
 damageState(phase)%dotState(2,constituent) = &
   damage_local_crack_mobility(instance)* &
   (trialDamage - damageState(phase)%state(2,constituent))
  
end subroutine damage_local_dotState

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_local_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_local_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_local_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_local_postResults = 0.0_pReal

 do o = 1_pInt,damage_local_Noutput(instance)
    select case(damage_local_outputID(o,instance))
      case (local_damage_ID)
        damage_local_postResults(c+1_pInt) = damageState(phase)%state(2,constituent)
        c = c + 1

    end select
 enddo
end function damage_local_postResults

end module damage_local
