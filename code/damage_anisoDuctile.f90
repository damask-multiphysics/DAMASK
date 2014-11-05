!--------------------------------------------------------------------------------------------------
! $Id: damage_anisoDuctile.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_anisoDuctile
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_anisoDuctile_sizePostResults                                                                   !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   damage_anisoDuctile_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   damage_anisoDuctile_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   damage_anisoDuctile_Noutput                                                                           !< number of outputs per instance of this damage 
   
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   damage_anisoDuctile_totalNslip                                                                    !< total number of slip systems
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   damage_anisoDuctile_Nslip                                                                         !< number of slip systems per family
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_anisoDuctile_aTol_damage

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   damage_anisoDuctile_critAccShear

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_anisoDuctile_outputID                                                                  !< ID of each post result output


 public :: &
   damage_anisoDuctile_init, &
   damage_anisoDuctile_stateInit, &
   damage_anisoDuctile_aTolState, &
   damage_anisoDuctile_dotState, &
   damage_anisoDuctile_microstructure, &
   damage_anisoDuctile_getDamage, &
   damage_anisoDuctile_putLocalDamage, &
   damage_anisoDuctile_getLocalDamage, &
   damage_anisoDuctile_getSlipDamage, &
   damage_anisoDuctile_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_init(fileUnit)
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
   LOCAL_damage_anisoDuctile_label, &
   LOCAL_damage_anisoDuctile_ID, &
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
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_damage_anisoDuctile_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id: damage_anisoDuctile.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_damage == LOCAL_damage_anisoDuctile_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_anisoDuctile_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_anisoDuctile_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_anisoDuctile_output(maxval(phase_Noutput),maxNinstance))
          damage_anisoDuctile_output = ''
 allocate(damage_anisoDuctile_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_anisoDuctile_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_anisoDuctile_critAccShear(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(damage_anisoDuctile_Nslip(lattice_maxNslipFamily,maxNinstance),source=0_pInt)
 allocate(damage_anisoDuctile_totalNslip(maxNinstance),                      source=0_pInt)
 allocate(damage_anisoDuctile_aTol_damage(maxNinstance),                         source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_damage_anisoDuctile_ID) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_anisoDuctile_Noutput(instance) = damage_anisoDuctile_Noutput(instance) + 1_pInt
             damage_anisoDuctile_outputID(damage_anisoDuctile_Noutput(instance),instance) = local_damage_ID
             damage_anisoDuctile_output(damage_anisoDuctile_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('atol_damage')
         damage_anisoDuctile_aTol_damage(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('nslip')  !
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoDuctile_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
         damage_anisoDuctile_totalNslip(instance) = sum(damage_anisoDuctile_Nslip(:,instance))

       case ('critical_accshear')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoDuctile_critAccShear(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_damage_anisoDuctile_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_anisoDuctile_Noutput(instance)
       select case(damage_anisoDuctile_outputID(o,instance))
         case(local_damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_anisoDuctile_sizePostResult(o,instance) = mySize
          damage_anisoDuctile_sizePostResults(instance)  = damage_anisoDuctile_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              = 1_pInt ! non-local damage
     sizeState                 = sizeDotState + 1_pInt

     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_anisoDuctile_sizePostResults(instance)
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

     call damage_anisoDuctile_stateInit(phase,instance)
     call damage_anisoDuctile_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_anisoDuctile_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_stateInit(phase,instance)
 use material, only: &
   damageState
 use math, only: &
   math_I3  
 
 implicit none
 integer(pInt),              intent(in) :: phase, instance                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState = 1.0_pReal
 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_anisoDuctile_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_anisoDuctile_aTol_damage(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_anisoDuctile_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_dotState(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState
 use lattice, only: &
   lattice_DamageMobility

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, &
   constituent

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)

 damageState(phase)%dotState(1,constituent) = &
   (damageState(phase)%state(2,constituent) - &
    damageState(phase)%state(1,constituent))/lattice_DamageMobility(phase)

end subroutine damage_anisoDuctile_dotState

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_microstructure(nSlip, accumulatedSlip, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
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
   phase, &
   constituent, &
   instance, &
   index, f, i

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 index = 1_pInt
 damageState(phase)%state(2,constituent) = 1.0_pReal
 do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,damage_anisoDuctile_Nslip(f,instance)                                            ! process each (active) slip system in family
     damageState(phase)%state(2,constituent) = damageState(phase)%state(2,constituent) - &
       accumulatedSlip(index)/damage_anisoDuctile_critAccShear(f,instance)
     index = index + 1_pInt
   enddo
 enddo
 damageState(phase)%state(2,constituent) = max(0.0_pReal, damageState(phase)%state(2,constituent))

end subroutine damage_anisoDuctile_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief returns damage
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_getDamage(ipc, ip, el)
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
 real(pReal) :: damage_anisoDuctile_getDamage
 
 select case(field_damage_type(material_homog(ip,el)))                                                   
   case (FIELD_DAMAGE_LOCAL_ID)
    damage_anisoDuctile_getDamage = damage_anisoDuctile_getLocalDamage(ipc, ip, el)
    
   case (FIELD_DAMAGE_NONLOCAL_ID)
    damage_anisoDuctile_getDamage = fieldDamage(material_homog(ip,el))% &
                                      field(1,mappingHomogenization(1,ip,el))                       ! Taylor type 

 end select
 
end function damage_anisoDuctile_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage value based on local damage 
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_putLocalDamage(ipc, ip, el, localDamage)
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
   localDamage
 
 damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el)) = &
   localDamage
 
end subroutine damage_anisoDuctile_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_getLocalDamage(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal)   :: &
   damage_anisoDuctile_getLocalDamage

 damage_anisoDuctile_getLocalDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(2,mappingConstitutive(1,ipc,ip,el))
 
end function damage_anisoDuctile_getLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns slip system damage
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_getSlipDamage(nSlip, accumulatedSlip, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use lattice, only: &
   lattice_maxNslipFamily

 implicit none
 integer(pInt), intent(in) :: &
   nSlip, &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(nSlip), intent(in) :: &
   accumulatedSlip
 real(pReal) :: &
   damage_anisoDuctile_getSlipDamage(nSlip), &
   nonlocalFactor
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   index, f, i

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 nonlocalFactor = damage_anisoDuctile_getDamage     (ipc, ip, el) - &
                  damage_anisoDuctile_getLocalDamage(ipc, ip, el)
 index = 1_pInt
 do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,damage_anisoDuctile_Nslip(f,instance)                                            ! process each (active) cleavage system in family
     damage_anisoDuctile_getSlipDamage(index) = &
       1.0_pReal/(1.0_pReal + &
                  accumulatedSlip(index)/damage_anisoDuctile_critAccShear(f,instance) - &
                  nonlocalFactor)

     index = index + 1_pInt
   enddo
 enddo

end function damage_anisoDuctile_getSlipDamage

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_anisoDuctile_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_anisoDuctile_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_anisoDuctile_postResults = 0.0_pReal

 do o = 1_pInt,damage_anisoDuctile_Noutput(instance)
    select case(damage_anisoDuctile_outputID(o,instance))
      case (local_damage_ID)
        damage_anisoDuctile_postResults(c+1_pInt) = &
          damage_anisoDuctile_getLocalDamage(ipc, ip, el)
        c = c + 1_pInt

    end select
 enddo
end function damage_anisoDuctile_postResults

end module damage_anisoDuctile
