!--------------------------------------------------------------------------------------------------
! $Id: damage_anisoBrittle.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut fŸr Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut fŸr Eisenforschung GmbH
!> @brief material subroutine incorporating anisotropic ductile damage
!> @details to be done
!--------------------------------------------------------------------------------------------------
module damage_anisoBrittle
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   damage_anisoBrittle_sizePostResults                                                                   !< cumulative size of post results

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   damage_anisoBrittle_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   damage_anisoBrittle_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   damage_anisoBrittle_Noutput                                                                           !< number of outputs per instance of this damage 
   
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   damage_anisoBrittle_totalNslip                                                                            !< Todo
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   damage_anisoBrittle_Nslip                                                                            !< Todo
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_anisoBrittle_aTol, &
   damage_anisoBrittle_sdot_0, &
   damage_anisoBrittle_N

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   damage_anisoBrittle_critDisp, &
   damage_anisoBrittle_critLoad

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_anisoBrittle_outputID                                                                  !< ID of each post result output


 public :: &
   damage_anisoBrittle_init, &
   damage_anisoBrittle_stateInit, &
   damage_anisoBrittle_aTolState, &
   damage_anisoBrittle_dotState, &
   damage_anisoBrittle_getDamage, &
   damage_anisoBrittle_putLocalDamage, &
   damage_anisoBrittle_getLocalDamage, &
   damage_anisoBrittle_getDamageStrain, &
   damage_anisoBrittle_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_init(fileUnit)
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
   LOCAL_damage_anisoBrittle_label, &
   LOCAL_damage_anisoBrittle_ID, &
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
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_damage_anisoBrittle_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id: damage_anisoBrittle.f90 3210 2014-06-17 15:24:44Z MPIE\m.diehl $'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_damage == LOCAL_damage_anisoBrittle_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(damage_anisoBrittle_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(damage_anisoBrittle_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(damage_anisoBrittle_output(maxval(phase_Noutput),maxNinstance))
          damage_anisoBrittle_output = ''
 allocate(damage_anisoBrittle_outputID(maxval(phase_Noutput),maxNinstance),      source=undefined_ID)
 allocate(damage_anisoBrittle_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(damage_anisoBrittle_critDisp(lattice_maxNslipFamily,maxNinstance),   source=0.0_pReal) 
 allocate(damage_anisoBrittle_critLoad  (lattice_maxNslipFamily,maxNinstance),   source=0.0_pReal) 
 allocate(damage_anisoBrittle_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(damage_anisoBrittle_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(damage_anisoBrittle_aTol(maxNinstance),                                source=0.0_pReal) 
 allocate(damage_anisoBrittle_sdot_0(maxNinstance),                              source=0.0_pReal) 
 allocate(damage_anisoBrittle_N(maxNinstance),                                   source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_damage_anisoBrittle_ID) then              ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = phase_damageInstance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
           case ('local_damage')
             damage_anisoBrittle_Noutput(instance) = damage_anisoBrittle_Noutput(instance) + 1_pInt
             damage_anisoBrittle_outputID(damage_anisoBrittle_Noutput(instance),instance) = local_damage_ID
             damage_anisoBrittle_output(damage_anisoBrittle_Noutput(instance),instance) = &
                                                       IO_lc(IO_stringValue(line,positions,2_pInt))
          end select

       case ('atol_damage')
         damage_anisoBrittle_aTol(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('sdot_0')
         damage_anisoBrittle_sdot_0(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('n_damage')
         damage_anisoBrittle_N(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('Nslip')  !
         Nchunks_SlipFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoBrittle_Nslip(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo
         damage_anisoBrittle_totalNslip(instance) = sum(damage_anisoBrittle_Nslip(:,instance))

       case ('critical_displacement')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoBrittle_critDisp(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

       case ('critical_load')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoBrittle_critLoad(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

     end select
   endif; endif
 enddo parsingFile
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_damage_anisoBrittle_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_anisoBrittle_Noutput(instance)
       select case(damage_anisoBrittle_outputID(o,instance))
         case(local_damage_ID)
           mySize = damage_anisoBrittle_totalNslip(instance)
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_anisoBrittle_sizePostResult(o,instance) = mySize
          damage_anisoBrittle_sizePostResults(instance)  = damage_anisoBrittle_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              = 2_pInt + & ! viscous and non-viscous damage values
                                 9_pInt + & ! damage deformation gradient  
                                 damage_anisoBrittle_totalNslip(instance) ! opening on each damage system
     sizeState                 = sizeDotState

     damageState(phase)%sizeState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = damage_anisoBrittle_sizePostResults(instance)
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

     call damage_anisoBrittle_stateInit(phase)
     call damage_anisoBrittle_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_anisoBrittle_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_stateInit(phase)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState     = 0.0_pReal
 tempState(1)  = 1.0_pReal
 tempState(2)  = 1.0_pReal
 tempState(3)  = 1.0_pReal
 tempState(4)  = 0.0_pReal
 tempState(5)  = 0.0_pReal
 tempState(6)  = 0.0_pReal
 tempState(7)  = 1.0_pReal
 tempState(8)  = 0.0_pReal
 tempState(9)  = 0.0_pReal
 tempState(10) = 0.0_pReal
 tempState(11) = 1.0_pReal
 damageState(phase)%state = spread(tempState,2,size(damageState(phase)%state(1,:)))
 damageState(phase)%state0 = damageState(phase)%state
 damageState(phase)%partionedState0 = damageState(phase)%state
end subroutine damage_anisoBrittle_stateInit

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_aTolState(phase,instance)
 use material, only: &
  damageState

 implicit none
 integer(pInt), intent(in) ::  &
   phase, &
   instance                                                                                         ! number specifying the current instance of the damage
 real(pReal), dimension(damageState(phase)%sizeState) :: tempTol

 tempTol = damage_anisoBrittle_aTol(instance)
 damageState(phase)%aTolState = tempTol
end subroutine damage_anisoBrittle_aTolState
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_dotState(Tstar_v,ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_mul33x33
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_DamageMobility

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   j, f, i, index_myFamily
 real(pReal), dimension(3,3) :: &
   Ld
 real(pReal) :: &
   tau, &
   tau_critical, &
   nonLocalFactor    

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 damageState(phase)%dotState(1,constituent) = &
   (1.0_pReal/lattice_DamageMobility(phase))* &
   (damageState(phase)%state(2,constituent) - &
    damageState(phase)%state(1,constituent))

 nonLocalFactor = 1.0_pReal + &
                  (damageState(phase)%state(1,constituent) - &
                   damage_anisoBrittle_getDamage(ipc, ip, el)) 
 Ld = 0.0_pReal
 j = 0_pInt
 slipFamiliesLoop: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,damage_anisoBrittle_Nslip(f,instance)                                            ! process each (active) slip system in family
     j = j+1_pInt
     tau = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,phase))
     tau_critical = (1.0_pReal - damageState(phase)%state(11+j,constituent)/&
                                 damage_anisoBrittle_critDisp(f,instance))* &
                    damage_anisoBrittle_critLoad(f,instance)*nonLocalFactor             
     damageState(phase)%dotState(11+j,constituent) = &
       damage_anisoBrittle_sdot_0(instance)*(tau/tau_critical)**damage_anisoBrittle_N(instance)
     damageState(phase)%dotState(2,constituent) = &
       damageState(phase)%dotState(2,constituent) - &
       2.0_pReal*tau*damageState(phase)%dotState(11+j,constituent)/ &
       (damage_anisoBrittle_critDisp(f,instance)*damage_anisoBrittle_critLoad(f,instance))
     Ld = Ld + damageState(phase)%dotState(11+j,constituent)* &
               lattice_Sslip(1:3,1:3,1,index_myFamily+i,phase)
   enddo
 enddo slipFamiliesLoop
 damageState(phase)%dotState(3:11,constituent) = &
   reshape(math_mul33x33(Ld,reshape(damageState(phase)%state(3:11,constituent),shape=[3,3])),shape=[9])

end subroutine damage_anisoBrittle_dotState
 
!--------------------------------------------------------------------------------------------------
!> @brief returns damage
!--------------------------------------------------------------------------------------------------
function damage_anisoBrittle_getDamage(ipc, ip, el)
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
 real(pReal) :: damage_anisoBrittle_getDamage
 
 select case(field_damage_type(material_homog(ip,el)))                                                   
   case (FIELD_DAMAGE_LOCAL_ID)
    damage_anisoBrittle_getDamage = damage_anisoBrittle_getLocalDamage(ipc, ip, el)
    
   case (FIELD_DAMAGE_NONLOCAL_ID)
    damage_anisoBrittle_getDamage =    fieldDamage(material_homog(ip,el))% &
      field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

 end select
 
end function damage_anisoBrittle_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage value based on local damage 
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_putLocalDamage(ipc, ip, el, localDamage)
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
 
end subroutine damage_anisoBrittle_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage
!--------------------------------------------------------------------------------------------------
function damage_anisoBrittle_getLocalDamage(ipc, ip, el)
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
   damage_anisoBrittle_getLocalDamage
 
 damage_anisoBrittle_getLocalDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function damage_anisoBrittle_getLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
function damage_anisoBrittle_getDamageStrain(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   damage_anisoBrittle_getDamageStrain
 
 damage_anisoBrittle_getDamageStrain = &
   reshape(damageState(mappingConstitutive(2,ipc,ip,el))%state(3:11,mappingConstitutive(1,ipc,ip,el)), &
           shape=[3,3])
 
end function damage_anisoBrittle_getDamageStrain

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_anisoBrittle_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance,& 
   damageState

 implicit none
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(damage_anisoBrittle_sizePostResults(phase_damageInstance(mappingConstitutive(2,ipc,ip,el)))) :: &
   damage_anisoBrittle_postResults

 integer(pInt) :: &
   instance, phase, constituent, o, c
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 c = 0_pInt
 damage_anisoBrittle_postResults = 0.0_pReal

 do o = 1_pInt,damage_anisoBrittle_Noutput(instance)
    select case(damage_anisoBrittle_outputID(o,instance))
      case (local_damage_ID)
        damage_anisoBrittle_postResults(c+1_pInt:c+damage_anisoBrittle_totalNslip(instance)) = &
          damageState(phase)%state(1,constituent)
        c = c + damage_anisoBrittle_totalNslip(instance)

    end select
 enddo
end function damage_anisoBrittle_postResults

end module damage_anisoBrittle
