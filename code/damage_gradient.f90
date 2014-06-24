!--------------------------------------------------------------------------------------------------
! $Id: damage_gradient.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_gradient
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_gradient_sizeDotState, &                                                           !< number of dotStates
   damage_gradient_sizeState, &                                                              !< total number of microstructural state variables
   damage_gradient_sizePostResults                                                           !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_gradient_sizePostResult                                                            !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_gradient_output                                                                    !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   damage_gradient_Noutput                                                                   !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),     allocatable,         private :: &
   damage_gradient_crack_mobility

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID, &
                 gradient_damage_ID
 end enum
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_gradient_outputID                                                                  !< ID of each post result output


 public :: &
   damage_gradient_init, &
   damage_gradient_stateInit, &
   damage_gradient_aTolState, &
   damage_gradient_microstructure, &
   damage_gradient_dotState, &
   damage_gradient_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_gradient_init(fileUnit)
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
   DAMAGE_GRADIENT_label, &
   DAMAGE_gradient_ID, &
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
  
 write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_GRADIENT_label//' init  -+>>>'
 write(6,'(a)')     ' $Id: damage_gradient.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_damage == DAMAGE_gradient_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_gradient_sizeDotState(maxNinstance),                        source=0_pInt)
 allocate(damage_gradient_sizeState(maxNinstance),                           source=0_pInt)
 allocate(damage_gradient_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_gradient_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_gradient_output(maxval(phase_Noutput),maxNinstance))
          damage_gradient_output = ''
 allocate(damage_gradient_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_gradient_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_gradient_crack_mobility(maxNinstance),                      source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == DAMAGE_gradient_ID) then               ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                     ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_gradient_outputID(damage_gradient_Noutput(instance),instance) = local_damage_ID
             damage_gradient_Noutput(instance) = damage_gradient_Noutput(instance) + 1_pInt
             damage_gradient_output(damage_gradient_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
           case ('gradient_damage')
             damage_gradient_outputID(damage_gradient_Noutput(instance),instance) = gradient_damage_ID
             damage_gradient_Noutput(instance) = damage_gradient_Noutput(instance) + 1_pInt
             damage_gradient_output(damage_gradient_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('crack_mobility')
         damage_gradient_crack_mobility(instance) = IO_floatValue(line,positions,2_pInt)
     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == DAMAGE_gradient_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)
     damage_gradient_sizeDotState(instance) = 1_pInt
     damage_gradient_sizeState(instance) = 3_pInt

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_gradient_Noutput(instance)
       select case(damage_gradient_outputID(o,instance))
         case(local_damage_ID, &
              gradient_damage_ID &
              )
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_gradient_sizePostResult(o,instance) = mySize
          damage_gradient_sizePostResults(instance)  = damage_gradient_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   damage_gradient_sizeDotState(instance)
     sizeState                 =   damage_gradient_sizeState   (instance)
                
     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
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

     call damage_gradient_stateInit(phase,instance)
     call damage_gradient_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_gradient_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_gradient_stateInit(phase,instance)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the damage
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1) = 0.0_pReal
 tempState(2:3) = 1.0_pReal
 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_gradient_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_gradient_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = 0.0_pReal
 damageState(phase)%aTolState = tempTol
end subroutine damage_gradient_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_gradient_microstructure(Tstar_v, Fe, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use mesh, only: &
   charLength
 use math, only: &
   math_Mandel6to33, &
   math_mul33x33, &
   math_transpose33, &
   math_I3
 use lattice, only: &
   lattice_surfaceEnergy33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Fe
 integer(pInt) :: &
   instance, phase, constituent 
 real(pReal) :: &
   damage

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 damage = damageState(phase)%state(3,constituent)*damageState(phase)%state(3,constituent)
 
 damageState(phase)%state(2,constituent) = &
   min(1.0_pReal, &   
       2.0_pReal*charLength*maxval(lattice_surfaceEnergy33(1:3,1:3,phase))/ &
       (0.125_pReal*sum(math_Mandel6to33(Tstar_v/damage)*(math_mul33x33(math_transpose33(Fe),Fe)-math_I3)) + &
        0.5_pReal*damageState(phase)%state(1,constituent)) &
       )  
  
end subroutine damage_gradient_microstructure
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_gradient_dotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_Mandel6to33

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp
 integer(pInt) :: &
   instance, phase, constituent 

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 damageState(phase)%dotState(1,constituent) = &
   sum(abs(math_Mandel6to33(Tstar_v)*Lp))
  
end subroutine damage_gradient_dotState

 
!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_gradient_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_gradient_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_gradient_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_gradient_postResults = 0.0_pReal

 do o = 1_pInt,damage_gradient_Noutput(instance)
    select case(damage_gradient_outputID(o,instance))
 
      case (local_damage_ID)
        damage_gradient_postResults(c+1_pInt) = damageState(phase)%state(2,constituent)
        c = c + 1
      case (gradient_damage_ID)
        damage_gradient_postResults(c+1_pInt) = damageState(phase)%state(3,constituent)
        c = c + 1
    end select
 enddo
end function damage_gradient_postResults

end module damage_gradient
