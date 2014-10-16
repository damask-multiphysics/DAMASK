!--------------------------------------------------------------------------------------------------
! $Id: damage_anisotropic.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_anisotropic
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_anisotropic_sizePostResults                                                                   !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   damage_anisotropic_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   damage_anisotropic_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   damage_anisotropic_Noutput                                                                           !< number of outputs per instance of this damage 
   
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   damage_anisotropic_totalNslip                                                                            !< Todo
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   damage_anisotropic_Nslip                                                                            !< Todo
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_anisotropic_aTol

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   damage_anisotropic_critpStrain

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_anisotropic_outputID                                                                  !< ID of each post result output


 public :: &
   damage_anisotropic_init, &
   damage_anisotropic_stateInit, &
   damage_anisotropic_aTolState, &
   damage_anisotropic_dotState, &
   damage_anisotropic_microstructure, &
   constitutive_anisotropic_getDamage, &
   constitutive_anisotropic_putDamage, &
   damage_anisotropic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_anisotropic_init(fileUnit)
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
   LOCAL_DAMAGE_anisotropic_label, &
   LOCAL_DAMAGE_anisotropic_ID, &
   material_phase, &  
   damageState, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank, &
   numerics_integrator
 use lattice, only: &
   lattice_maxNslipFamily

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,o
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase   
 integer(pInt) :: Nchunks_SlipFamilies, j   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_DAMAGE_anisotropic_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id: damage_anisotropic.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_damage == LOCAL_DAMAGE_anisotropic_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_anisotropic_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_anisotropic_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_anisotropic_output(maxval(phase_Noutput),maxNinstance))
          damage_anisotropic_output = ''
 allocate(damage_anisotropic_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_anisotropic_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_anisotropic_critpStrain(lattice_maxNslipFamily,maxNinstance),  source=0.0_pReal) 
 allocate(damage_anisotropic_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(damage_anisotropic_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(damage_anisotropic_aTol(maxNinstance),                                source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_DAMAGE_anisotropic_ID) then              ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_anisotropic_Noutput(instance) = damage_anisotropic_Noutput(instance) + 1_pInt
             damage_anisotropic_outputID(damage_anisotropic_Noutput(instance),instance) = local_damage_ID
             damage_anisotropic_output(damage_anisotropic_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('atol_damage')
         damage_anisotropic_aTol(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('Nslip')  !
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisotropic_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
         damage_anisotropic_totalNslip(instance) = sum(damage_anisotropic_Nslip(:,instance))

       case ('critical_plastic_strain')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisotropic_critpStrain(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_DAMAGE_anisotropic_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_anisotropic_Noutput(instance)
       select case(damage_anisotropic_outputID(o,instance))
         case(local_damage_ID)
           mySize = damage_anisotropic_totalNslip(instance)
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_anisotropic_sizePostResult(o,instance) = mySize
          damage_anisotropic_sizePostResults(instance)  = damage_anisotropic_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              =   damage_anisotropic_totalNslip(instance)
     sizeState                 =   2_pInt * damage_anisotropic_totalNslip(instance)
                
     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_anisotropic_sizePostResults(instance)
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

     call damage_anisotropic_stateInit(phase)
     call damage_anisotropic_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_anisotropic_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisotropic_stateInit(phase)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState = 1.0_pReal
 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_anisotropic_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisotropic_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_anisotropic_aTol(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_anisotropic_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisotropic_dotState(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
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
   phase, &
   constituent, &
   instance, &
   i

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 do i = 1_pInt,damage_anisotropic_totalNslip(instance)
   damageState(phase)%dotState(i,constituent) = &
     (1.0_pReal/lattice_DamageMobility(phase))* &
     (damageState(phase)%state(i+damage_anisotropic_totalNslip(instance),constituent) - &
      damageState(phase)%state(i,constituent))
 enddo     
  
end subroutine damage_anisotropic_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisotropic_microstructure(nSlip,accumulatedSlip,ipc, ip, el)
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
 use lattice, only: &
   lattice_maxNslipFamily

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(nSlip), intent(in) :: &
   accumulatedSlip
 integer(pInt) :: &
   phase, constituent, instance, i, j, f

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 j = 0_pInt
 do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,damage_anisotropic_Nslip(f,instance)                                          ! process each (active) slip system in family
     j = j+1_pInt
     damageState(phase)%state(j+damage_anisotropic_totalNslip(instance),constituent) = &
       min(damageState(phase)%state(j+damage_anisotropic_totalNslip(instance),constituent), &
           damage_anisotropic_critpStrain(f,instance)/accumulatedSlip(j))
   enddo
 enddo     
                                                        
end subroutine damage_anisotropic_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on local damage model state layout 
!--------------------------------------------------------------------------------------------------
function constitutive_anisotropic_getDamage(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: &
   constitutive_anisotropic_getDamage(damage_anisotropic_totalNslip(phase_damageInstance(mappingConstitutive(2,ipc,ip,el))))
 
 constitutive_anisotropic_getDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))% &
   state(1:damage_anisotropic_totalNslip(phase_damageInstance(mappingConstitutive(2,ipc,ip,el))), &
         mappingConstitutive(1,ipc,ip,el))
 
end function constitutive_anisotropic_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage value based on local damage 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_anisotropic_putDamage(ipc, ip, el, localDamage)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in) :: &
   localDamage(damage_anisotropic_totalNslip(phase_damageInstance(mappingConstitutive(2,ipc,ip,el))))
 integer(pInt) :: &
   phase, constituent, instance  
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 damageState(phase)%state(1:damage_anisotropic_totalNslip(instance),constituent) = &
   localDamage
 
end subroutine constitutive_anisotropic_putDamage

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_anisotropic_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_anisotropic_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_anisotropic_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_anisotropic_postResults = 0.0_pReal

 do o = 1_pInt,damage_anisotropic_Noutput(instance)
    select case(damage_anisotropic_outputID(o,instance))
      case (local_damage_ID)
        damage_anisotropic_postResults(c+1_pInt:c+damage_anisotropic_totalNslip(instance)) = &
          damageState(phase)%state(1,constituent)
        c = c + damage_anisotropic_totalNslip(instance)

    end select
 enddo
end function damage_anisotropic_postResults

end module damage_anisotropic
