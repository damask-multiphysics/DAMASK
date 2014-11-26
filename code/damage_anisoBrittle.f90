!--------------------------------------------------------------------------------------------------
! $Id$
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
   damage_anisoBrittle_totalNcleavage                                                                    !< total number of cleavage systems
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   damage_anisoBrittle_Ncleavage                                                                         !< number of cleavage systems per family
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_anisoBrittle_aTol_damage, &
   damage_anisoBrittle_aTol_disp, &
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
   damage_anisoBrittle_microstructure, &
   damage_anisoBrittle_LdAndItsTangent, &
   damage_anisoBrittle_getFd, &
   damage_anisoBrittle_putFd, &
   damage_anisoBrittle_getFd0, &
   damage_anisoBrittle_getPartionedFd0, &
   damage_anisoBrittle_getDamage, &
   damage_anisoBrittle_putLocalDamage, &
   damage_anisoBrittle_getLocalDamage, &
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
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,mySize=0_pInt,phase,instance,o
 integer(pInt) :: sizeState, sizeDotState
 integer(pInt) :: NofMyPhase   
 integer(pInt) :: Nchunks_CleavageFamilies, j   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_damage_anisoBrittle_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
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
 allocate(damage_anisoBrittle_critDisp(lattice_maxNcleavageFamily,maxNinstance), source=0.0_pReal) 
 allocate(damage_anisoBrittle_critLoad(lattice_maxNcleavageFamily,maxNinstance), source=0.0_pReal) 
 allocate(damage_anisoBrittle_Ncleavage(lattice_maxNcleavageFamily,maxNinstance),source=0_pInt)
 allocate(damage_anisoBrittle_totalNcleavage(maxNinstance),                      source=0_pInt)
 allocate(damage_anisoBrittle_aTol_damage(maxNinstance),                         source=0.0_pReal) 
 allocate(damage_anisoBrittle_aTol_disp(maxNinstance),                           source=0.0_pReal) 
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
   if (phase > 0_pInt ) then; if (phase_damage(phase) == LOCAL_damage_anisoBrittle_ID) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
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
         damage_anisoBrittle_aTol_damage(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('atol_disp')
         damage_anisoBrittle_aTol_disp(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('sdot0')
         damage_anisoBrittle_sdot_0(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('damageratesensitivity')
         damage_anisoBrittle_N(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('ncleavage')  !
         Nchunks_CleavageFamilies = positions(1) - 1_pInt
         do j = 1_pInt, Nchunks_CleavageFamilies
           damage_anisoBrittle_Ncleavage(j,instance) = IO_intValue(line,positions,1_pInt+j)
         enddo

       case ('criticaldisplacement')
         do j = 1_pInt, Nchunks_CleavageFamilies
           damage_anisoBrittle_critDisp(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

       case ('criticalload')
         do j = 1_pInt, Nchunks_CleavageFamilies
           damage_anisoBrittle_critLoad(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo

     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_damage)   
   myPhase: if (phase_damage(phase) == LOCAL_damage_anisoBrittle_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)
!  sanity checks
     damage_anisoBrittle_Ncleavage(1:lattice_maxNcleavageFamily,instance) = &
       min(lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,phase),&                            ! limit active cleavage systems per family to min of available and requested
           damage_anisoBrittle_Ncleavage(1:lattice_maxNcleavageFamily,instance))
     damage_anisoBrittle_totalNcleavage(instance)  = sum(damage_anisoBrittle_Ncleavage(:,instance)) ! how many cleavage systems altogether
     if (damage_anisoBrittle_aTol_damage(instance) < 0.0_pReal) &
       damage_anisoBrittle_aTol_damage(instance) = 1.0e-3_pReal                                     ! default absolute tolerance 1e-3
     if (damage_anisoBrittle_aTol_disp(instance) >= 1.0e-3_pReal) &
       damage_anisoBrittle_aTol_disp(instance) = 1.0e-3_pReal                                       ! default absolute tolerance 1e-3
     if (damage_anisoBrittle_sdot_0(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//LOCAL_DAMAGE_anisoBrittle_LABEL//')')
     if (any(damage_anisoBrittle_critDisp(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='critical_displacement ('//LOCAL_DAMAGE_anisoBrittle_LABEL//')')
     if (any(damage_anisoBrittle_critLoad(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='critical_load ('//LOCAL_DAMAGE_anisoBrittle_LABEL//')')
     if (damage_anisoBrittle_N(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity_damage ('//LOCAL_DAMAGE_anisoBrittle_LABEL//')')
   endif myPhase
 enddo sanityChecks
 
 initializeInstances: do phase = 1_pInt, size(phase_damage)
   if (phase_damage(phase) == LOCAL_damage_anisoBrittle_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,damage_anisoBrittle_Noutput(instance)
       select case(damage_anisoBrittle_outputID(o,instance))
         case(local_damage_ID)
           mySize = 1_pInt
       end select
 
       if (mySize > 0_pInt) then  ! any meaningful output found
          damage_anisoBrittle_sizePostResult(o,instance) = mySize
          damage_anisoBrittle_sizePostResults(instance)  = damage_anisoBrittle_sizePostResults(instance) + mySize
       endif
     enddo outputsLoop
! Determine size of state array
     sizeDotState              = 1_pInt + & ! non-local damage
                                 damage_anisoBrittle_totalNcleavage(instance)     ! opening on each damage system
     sizeState                 = sizeDotState + &
                                 damage_anisoBrittle_totalNcleavage(instance) + & ! local damage on each damage system
                                 9_pInt ! Fd

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

     call damage_anisoBrittle_stateInit(phase,instance)
     call damage_anisoBrittle_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_anisoBrittle_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_stateInit(phase,instance)
 use material, only: &
   damageState
 use math, only: &
   math_I3  
 
 implicit none
 integer(pInt),              intent(in) :: phase, instance                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1)                                                 = 1.0_pReal
 tempState(2 : &
           1 +  damage_anisoBrittle_totalNcleavage(instance)) = 0.0_pReal 
 tempState(2 +  damage_anisoBrittle_totalNcleavage(instance): &
           1 +2*damage_anisoBrittle_totalNcleavage(instance)) = 1.0_pReal
 tempState(2 +2*damage_anisoBrittle_totalNcleavage(instance): &
           10+2*damage_anisoBrittle_totalNcleavage(instance)) = reshape(math_I3, shape=[9])
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

 tempTol(1)                                                 = damage_anisoBrittle_aTol_damage(instance)
 tempTol(2 : &
         1 +  damage_anisoBrittle_totalNcleavage(instance)) = damage_anisoBrittle_aTol_disp  (instance) 
 tempTol(2 +  damage_anisoBrittle_totalNcleavage(instance): &
         1 +2*damage_anisoBrittle_totalNcleavage(instance)) = damage_anisoBrittle_aTol_damage(instance)
 tempTol(2 +2*damage_anisoBrittle_totalNcleavage(instance): &
         10+2*damage_anisoBrittle_totalNcleavage(instance)) = damage_anisoBrittle_aTol_damage(instance)
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
 use lattice, only: &
   lattice_Scleavage_v, &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem, &
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
   f, i, index_d, index_o, index_myFamily
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, udott, udotn, &
   nonlocalFactor, localDamage

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 localDamage = max(0.0_pReal, &
                   1.0_pReal - sum(1.0_pReal - damageState(phase)% &
                                    state(2+  damage_anisoBrittle_totalNcleavage(instance): &
                                          1+2*damage_anisoBrittle_totalNcleavage(instance),constituent)))
 damageState(phase)%dotState(1,constituent) = &
   (localDamage - damageState(phase)%state(1,constituent))/lattice_DamageMobility(phase)
 nonlocalFactor = damage_anisoBrittle_getDamage(ipc, ip, el) - localDamage

 index_o = 2_pInt
 index_d = 2_pInt + damage_anisoBrittle_totalNcleavage(instance)
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,damage_anisoBrittle_Ncleavage(f,instance)                                            ! process each (active) cleavage system in family
     traction_d    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,1,index_myFamily+i,phase))
     traction_t    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,2,index_myFamily+i,phase))
     traction_n    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,3,index_myFamily+i,phase))
     traction_crit = damage_anisoBrittle_critLoad(f,instance)* &
                     damageState(phase)%state(index_d,constituent)* &
                     damageState(phase)%state(index_d,constituent)
                    
     udotd = &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_d)/traction_crit)**damage_anisoBrittle_N(instance)
     udott = &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_t)/traction_crit)**damage_anisoBrittle_N(instance)
     udotn = &
       damage_anisoBrittle_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit)**damage_anisoBrittle_N(instance)

     damageState(phase)%dotState(index_o,constituent) = &
       (udotd + udott + udotn)/damage_anisoBrittle_critDisp(f,instance)

     index_d = index_d + 1_pInt; index_o = index_o + 1_pInt
   enddo
 enddo

end subroutine damage_anisoBrittle_dotState

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_microstructure(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use lattice, only: &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   f, i, index_d, index_o, index_myFamily
 real(pReal) :: &
   nonlocalFactor, localDamage

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 localDamage = max(0.0_pReal, &
                   1.0_pReal - sum(1.0_pReal - damageState(phase)% &
                                    state(2+  damage_anisoBrittle_totalNcleavage(instance): &
                                          1+2*damage_anisoBrittle_totalNcleavage(instance),constituent)))
 nonlocalFactor = damage_anisoBrittle_getDamage(ipc, ip, el) - localDamage

 index_o = 2_pInt
 index_d = 2_pInt + damage_anisoBrittle_totalNcleavage(instance)
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,damage_anisoBrittle_Ncleavage(f,instance)                                            ! process each (active) cleavage system in family
     damageState(phase)%state(index_d,constituent) = &
       min(damageState(phase)%state0(index_d,constituent), &
           1.0_pReal/max(0.0_pReal,damageState(phase)%state(index_o,constituent) - &
                                   nonlocalFactor))

     index_d = index_d + 1_pInt; index_o = index_o + 1_pInt
   enddo
 enddo

end subroutine damage_anisoBrittle_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_LdAndItsTangent(Ld, dLd_dTstar, Tstar_v, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_Plain3333to99
 use lattice, only: &
   lattice_Scleavage, &
   lattice_Scleavage_v, &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(out), dimension(3,3) :: &
   Ld                                                                                               !< damage velocity gradient
 real(pReal),   intent(out), dimension(9,9) :: &
   dLd_dTstar                                                                                       !< derivative of Ld with respect to Tstar (2nd-order tensor)
 real(pReal),   dimension(3,3,3,3) :: &
   dLd_dTstar3333                                                                                   !< derivative of Ld with respect to Tstar (4th-order tensor)
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   f, i, index, index_myFamily, k, l, m, n
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 Ld = 0.0_pReal
 dLd_dTstar3333 = 0.0_pReal
 index = 2_pInt + damage_anisoBrittle_totalNcleavage(instance)
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,damage_anisoBrittle_Ncleavage(f,instance)                                            ! process each (active) cleavage system in family
     traction_d    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,1,index_myFamily+i,phase))
     traction_t    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,2,index_myFamily+i,phase))
     traction_n    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,3,index_myFamily+i,phase))
     traction_crit = damage_anisoBrittle_critLoad(f,instance)* &
                     damageState(phase)%state(index,constituent)* &
                     damageState(phase)%state(index,constituent)
     udotd = &
       sign(1.0_pReal,traction_d)* &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_d)/traction_crit)**damage_anisoBrittle_N(instance)
     if (udotd /= 0.0_pReal) then
       Ld = Ld + udotd*lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase)
       dudotd_dt = udotd*damage_anisoBrittle_N(instance)/traction_d
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotd_dt*lattice_Scleavage(k,l,1,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,1,index_myFamily+i,phase)
     endif                

     udott = &
       sign(1.0_pReal,traction_t)* &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_t)/traction_crit)**damage_anisoBrittle_N(instance)
     if (udott /= 0.0_pReal) then
       Ld = Ld + udott*lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase)
       dudott_dt = udott*damage_anisoBrittle_N(instance)/traction_t
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudott_dt*lattice_Scleavage(k,l,2,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,2,index_myFamily+i,phase)
     endif                

     udotn = &
       damage_anisoBrittle_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit)**damage_anisoBrittle_N(instance)
     if (udotn /= 0.0_pReal) then
       Ld = Ld + udotn*lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase)
       dudotn_dt = udotn*damage_anisoBrittle_N(instance)/traction_n
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotn_dt*lattice_Scleavage(k,l,3,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,3,index_myFamily+i,phase)
     endif                

     index = index + 1_pInt
   enddo
 enddo
 dLd_dTstar = math_Plain3333to99(dLd_dTstar3333)
 
end subroutine damage_anisoBrittle_LdAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoBrittle_getFd(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   phase, &
   constituent, &
   instance
 real(pReal), dimension(3,3) :: &
   damage_anisoBrittle_getFd
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoBrittle_getFd = &
   reshape(damageState(phase)% &
     state(2*damage_anisoBrittle_totalNcleavage(instance)+2: &
           2*damage_anisoBrittle_totalNcleavage(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoBrittle_getFd

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_putFd(Tstar_v, dt, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_mul33x33, &
   math_inv33, &
   math_I3
 use lattice, only: &
   lattice_Scleavage, &
   lattice_Scleavage_v, &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(in) :: &
   dt
 real(pReal), dimension(3,3) :: &
   Ld                                                                                               !< damage velocity gradient
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   f, i, index, index_myFamily
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, udott, udotn

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 Ld = 0.0_pReal
 index = 2_pInt + damage_anisoBrittle_totalNcleavage(instance)
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,damage_anisoBrittle_Ncleavage(f,instance)                                            ! process each (active) cleavage system in family
     traction_d    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,1,index_myFamily+i,phase))
     traction_t    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,2,index_myFamily+i,phase))
     traction_n    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,3,index_myFamily+i,phase))
     traction_crit = damage_anisoBrittle_critLoad(f,instance)* &
                     damageState(phase)%state(index,constituent)* &
                     damageState(phase)%state(index,constituent)
     udotd = &
       sign(1.0_pReal,traction_d)* &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_d)/traction_crit)**damage_anisoBrittle_N(instance)
     Ld = Ld + udotd*lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase)

     udott = &
       sign(1.0_pReal,traction_t)* &
       damage_anisoBrittle_sdot_0(instance)* &
       (abs(traction_t)/traction_crit)**damage_anisoBrittle_N(instance)
     Ld = Ld + udott*lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase)

     udotn = &
       damage_anisoBrittle_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit)**damage_anisoBrittle_N(instance)
     Ld = Ld + udotn*lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase)

     index = index + 1_pInt
   enddo
 enddo
 damageState(phase)%state(2*damage_anisoBrittle_totalNcleavage(instance)+2: &
                          2*damage_anisoBrittle_totalNcleavage(instance)+10,constituent) = &
   reshape(math_mul33x33(math_inv33(math_I3 - dt*Ld), &
           damage_anisoBrittle_getFd0(ipc, ip, el)), shape=[9])                             

end subroutine damage_anisoBrittle_putFd 

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoBrittle_getFd0(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   phase, &
   constituent, &
   instance
 real(pReal), dimension(3,3) :: &
   damage_anisoBrittle_getFd0
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoBrittle_getFd0 = &
   reshape(damageState(phase)% &
     subState0(2*damage_anisoBrittle_totalNcleavage(instance)+2: &
               2*damage_anisoBrittle_totalNcleavage(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoBrittle_getFd0

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoBrittle_getPartionedFd0(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   phase, &
   constituent, &
   instance
 real(pReal), dimension(3,3) :: &
   damage_anisoBrittle_getPartionedFd0
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoBrittle_getPartionedFd0 = &
   reshape(damageState(phase)% &
     partionedState0(2*damage_anisoBrittle_totalNcleavage(instance)+2: &
                     2*damage_anisoBrittle_totalNcleavage(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoBrittle_getPartionedFd0

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
    damage_anisoBrittle_getDamage = fieldDamage(material_homog(ip,el))% &
                                      field(1,mappingHomogenization(1,ip,el))                       ! Taylor type 

 end select
 
end function damage_anisoBrittle_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns damage value based on local damage 
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoBrittle_putLocalDamage(ipc, ip, el, localDamage)
 use material, only: &
   mappingConstitutive, &
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
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_anisoBrittle_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance

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
        damage_anisoBrittle_postResults(c+1_pInt) = &
          damage_anisoBrittle_getLocalDamage(ipc, ip, el)
        c = c + 1_pInt

    end select
 enddo
end function damage_anisoBrittle_postResults

end module damage_anisoBrittle
