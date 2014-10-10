!--------------------------------------------------------------------------------------------------
! $Id: damage_gurson.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating gurson damage
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_gurson
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_gurson_sizePostResults                                                                   !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   damage_gurson_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   damage_gurson_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public :: &
   damage_gurson_Noutput                                                                           !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_gurson_aTol, &
   damage_gurson_critpStrain

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_gurson_outputID                                                                  !< ID of each post result output


 public :: &
   damage_gurson_init, &
   damage_gurson_stateInit, &
   damage_gurson_aTolState, &
   damage_gurson_dotState, &
   damage_gurson_microstructure, &
   constitutive_gurson_getDamage, &
   constitutive_gurson_putDamage, &
   damage_gurson_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_init(fileUnit)
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
   LOCAL_DAMAGE_gurson_label, &
   LOCAL_DAMAGE_gurson_ID, &
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
 write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_DAMAGE_gurson_LABEL//' init  -+>>>'
 write(6,'(a)')     ' $Id: damage_gurson.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_damage == LOCAL_DAMAGE_gurson_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_gurson_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_gurson_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_gurson_output(maxval(phase_Noutput),maxNinstance))
          damage_gurson_output = ''
 allocate(damage_gurson_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_gurson_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_gurson_critpStrain(maxNinstance),                    source=0.0_pReal) 
 allocate(damage_gurson_aTol(maxNinstance),                                source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_DAMAGE_gurson_ID) then              ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_gurson_Noutput(instance) = damage_gurson_Noutput(instance) + 1_pInt
             damage_gurson_outputID(damage_gurson_Noutput(instance),instance) = local_damage_ID
             damage_gurson_output(damage_gurson_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('critical_plastic_strain')
         damage_gurson_critpStrain(instance) = IO_floatValue(line,positions,2_pInt)

       case ('atol_damage')
         damage_gurson_aTol(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_DAMAGE_gurson_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_gurson_Noutput(instance)
       select case(damage_gurson_outputID(o,instance))
         case(local_damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_gurson_sizePostResult(o,instance) = mySize
          damage_gurson_sizePostResults(instance)  = damage_gurson_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   4_pInt
     sizeState                 =   6_pInt
                
     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_gurson_sizePostResults(instance)
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

     call damage_gurson_stateInit(phase,instance)
     call damage_gurson_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_gurson_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_stateInit(phase,instance)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: instance                                                 !< number specifying the instance of the damage
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1) = 1.0_pReal
 tempState(2) = 0.0_pReal
 tempState(3) = 1.0_pReal

 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_gurson_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_gurson_aTol(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_gurson_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_dotState(Lp, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState
 use math, only: &
   math_equivStrain33, &
   math_trace33
 use lattice, only: &
   lattice_DamageMobility

 implicit none
 real(pReal), intent(in), dimension(3,3) :: &
   Lp
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
   (damageState(phase)%state(6,constituent) - &
    damageState(phase)%state(1,constituent))

 damageState(phase)%dotState(2,constituent) = &
   damageState(phase)%dotState(3,constituent) + &
   damageState(phase)%dotState(4,constituent)                            ! total rate of void fraction evolution
   
 damageState(phase)%dotState(3,constituent) = &
   damageState(phase)%state(6,constituent) + &
   damageState(phase)%dotState(4,constituent)                            ! void nucleation rate
   
 damageState(phase)%dotState(4,constituent) = &
   (1_pReal - damageState(phase)%state(3,constituent)) * &
    math_trace33(Lp)                                                     ! void growth rate( proportional to hydrostatic part of Lp )
  
end subroutine damage_gurson_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_microstructure(ipc, ip, el)
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
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, constituent

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 
 damageState(phase)%state(5,constituent) = 0.0_pReal                                                       !< a statistical measure of aeformation hetrogeneity
 damageState(phase)%state(6,constituent) = min(damageState(phase)%state(6,constituent), &
                                                     damage_gurson_critpStrain(phase)/ &
                                                      damageState(phase)%state(2,constituent))      !< akin to damage surface
                                                       
end subroutine damage_gurson_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on local damage model state layout 
!--------------------------------------------------------------------------------------------------
function constitutive_gurson_getDamage(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_gurson_getDamage
 
 constitutive_gurson_getDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function constitutive_gurson_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage value based on local damage 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_gurson_putDamage(ipc, ip, el, localDamage)
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
 
end subroutine constitutive_gurson_putDamage

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_gurson_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_gurson_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_gurson_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_gurson_postResults = 0.0_pReal

 do o = 1_pInt,damage_gurson_Noutput(instance)
    select case(damage_gurson_outputID(o,instance))
      case (local_damage_ID)
        damage_gurson_postResults(c+1_pInt) = damageState(phase)%state(1,constituent)
        c = c + 1

    end select
 enddo
end function damage_gurson_postResults

end module damage_gurson
