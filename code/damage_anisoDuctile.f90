!--------------------------------------------------------------------------------------------------
! $Id$
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
   damage_anisoDuctile_critPlasticStrain

 real(pReal),                         dimension(:),           allocatable,         private :: &
   damage_anisoDuctile_sdot_0, &
   damage_anisoDuctile_N

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   damage_anisoDuctile_critLoad
   
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum 
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_anisoDuctile_outputID                                                                  !< ID of each post result output


 public :: &
   damage_anisoDuctile_init, &
   damage_anisoDuctile_stateInit, &
   damage_anisoDuctile_aTolState, &
   damage_anisoDuctile_microstructure, &
   damage_anisoDuctile_LdAndItsTangent, &
   damage_anisoDuctile_getFd, &
   damage_anisoDuctile_putFd, &
   damage_anisoDuctile_getFd0, &
   damage_anisoDuctile_getPartionedFd0, &
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
   LOCAL_damage_anisoDuctile_label, &
   LOCAL_damage_anisoDuctile_ID, &
   material_phase, &  
   damageState, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank, &
   numerics_integrator
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem

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
   write(6,'(a)')     ' $Id$'
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
 allocate(damage_anisoDuctile_critLoad(lattice_maxNslipFamily,maxNinstance), source=0.0_pReal) 
 allocate(damage_anisoDuctile_critPlasticStrain(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal) 
 allocate(damage_anisoDuctile_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(damage_anisoDuctile_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(damage_anisoDuctile_N(maxNinstance),                                   source=0.0_pReal) 
 allocate(damage_anisoDuctile_sdot_0(maxNinstance),                              source=0.0_pReal) 
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

       case ('sdot0')
         damage_anisoDuctile_sdot_0(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('criticalplasticstrain')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoDuctile_critPlasticStrain(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
         
       case ('damageratesensitivity')
         damage_anisoDuctile_N(instance) = IO_floatValue(line,positions,2_pInt)

       case ('criticalload')
         do j = 1_pInt, Nchunks_SlipFamilies
           damage_anisoDuctile_critLoad(j,instance) = IO_floatValue(line,positions,1_pInt+j)
         enddo
         
     end select
   endif; endif
 enddo parsingFile

 sanityChecks: do phase = 1_pInt, size(phase_damage)   
   myPhase: if (phase_damage(phase) == LOCAL_damage_anisoDuctile_ID) then
     NofMyPhase=count(material_phase==phase)
     instance = phase_damageInstance(phase)
!  sanity checks
     damage_anisoDuctile_Nslip(1:lattice_maxNslipFamily,instance) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase),&                                    ! limit active cleavage systems per family to min of available and requested
           damage_anisoDuctile_Nslip(1:lattice_maxNslipFamily,instance))
         damage_anisoDuctile_totalNslip(instance) = sum(damage_anisoDuctile_Nslip(:,instance))
     if (damage_anisoDuctile_aTol_damage(instance) < 0.0_pReal) &
       damage_anisoDuctile_aTol_damage(instance) = 1.0e-3_pReal                                     ! default absolute tolerance 1e-3
     if (damage_anisoDuctile_sdot_0(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//LOCAL_DAMAGE_anisoDuctile_LABEL//')')
     if (any(damage_anisoDuctile_critPlasticStrain(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='criticaPlasticStrain ('//LOCAL_DAMAGE_anisoDuctile_LABEL//')')
     if (damage_anisoDuctile_N(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity_damage ('//LOCAL_DAMAGE_anisoDuctile_LABEL//')')
   endif myPhase
 enddo sanityChecks
  

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
     sizeDotState              = 0_pInt
     sizeState                 = sizeDotState + &
                                 1_pInt + & ! time regularised damage
                                 damage_anisoDuctile_totalNslip(instance) + & ! slip system damages
                                 9 ! Fd
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
subroutine damage_anisoDuctile_stateInit(phase, instance)
 use material, only: &
   damageState
 use math, only: &
   math_I3  
 implicit none
 integer(pInt),              intent(in) :: phase , instance                                         !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1)                                                 = 1.0_pReal
 tempState(2 : &
           1 +  damage_anisoDuctile_totalNslip(instance))     =  1.0_pReal 
 tempState(damage_anisoDuctile_totalNslip(instance)+2: &
           damage_anisoDuctile_totalNslip(instance)+10) = reshape(math_I3, shape=[9])
           
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
subroutine damage_anisoDuctile_microstructure(subdt, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   plasticState, &
   damageState
 use lattice, only: &
   lattice_maxNslipFamily
 use lattice, only: &
   lattice_DamageMobility

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal),  intent(in) :: &
   subdt
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   index, f, i
 real(pReal) :: &
   localDamage, &
   drivingForce

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 localDamage    = minval(damageState(phase)%state(2:1+damage_anisoDuctile_totalNslip(instance),constituent))

 index = 1_pInt
 do f = 1_pInt,lattice_maxNslipFamily
   do i = 1_pInt,damage_anisoDuctile_Nslip(f,instance)                                            ! process each (active) slip system in family
     if (localDamage == damageState(phase)%state(index+1,constituent)) then
       drivingForce = plasticState(phase)%accumulatedSlip(index,constituent)/damage_anisoDuctile_critPlasticStrain(f,instance) - &
                      damage_anisoDuctile_getDamage(ipc, ip, el)
       damageState(phase)%state(index+1,constituent) = &
         min(damageState(phase)%state0(index+1,constituent), &
             (sqrt(drivingForce*drivingForce + 4.0_pReal) - drivingForce)/2.0_pReal)
     else
       drivingForce = plasticState(phase)%accumulatedSlip(index,constituent)/damage_anisoDuctile_critPlasticStrain(f,instance)
       damageState(phase)%state(index+1,constituent) = &
         min(damageState(phase)%state0(index+1,constituent), &
             1.0_pReal/drivingForce)                                                                 ! irreversibility
     endif
     index = index + 1_pInt
   enddo
 enddo

 localDamage = minval(damageState(phase)%state(2:1+damage_anisoDuctile_totalNslip(instance),constituent))

 damageState(phase)%state(1,constituent) = &
   localDamage + &
   (damageState(phase)%subState0(1,constituent) - localDamage)* &
   exp(-subdt/lattice_DamageMobility(phase))

end subroutine damage_anisoDuctile_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_LdAndItsTangent(Ld, dLd_dTstar, Tstar_v, ipc, ip, el)
 use numerics, only: &
   residualStiffness
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_sd, &
   lattice_st, &
   lattice_sn
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_Plain3333to99, &
   math_I3, &
   math_identity4th, &
   math_symmetric33, &
   math_Mandel33to6, &
   math_tensorproduct
 
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
 real(pReal),   dimension(3,3) :: &
   projection_d, projection_t, projection_n                                                         !< projection modes 3x3 tensor
 real(pReal),   dimension(6) :: &
   projection_d_v, projection_t_v, projection_n_v                                                   !< projection modes 3x3 vector
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
 

 index = 2_pInt
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                                      ! at which index starts my family
   do i = 1_pInt,damage_anisoDuctile_Nslip(f,instance)                                              ! process each (active) slip system in family
   
     projection_d = math_tensorproduct(lattice_sd(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_t = math_tensorproduct(lattice_st(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_n = math_tensorproduct(lattice_sn(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))

     projection_d_v(1:6) = math_Mandel33to6(math_symmetric33(projection_d(1:3,1:3)))
     projection_t_v(1:6) = math_Mandel33to6(math_symmetric33(projection_t(1:3,1:3)))
     projection_n_v(1:6) = math_Mandel33to6(math_symmetric33(projection_n(1:3,1:3)))
   
     traction_d    = dot_product(Tstar_v,projection_d_v(1:6))
     traction_t    = dot_product(Tstar_v,projection_t_v(1:6))
     traction_n    = dot_product(Tstar_v,projection_n_v(1:6))
     
     traction_crit = damage_anisoDuctile_critLoad(f,instance)* &
                     (damageState(phase)%state0(index,constituent) + residualStiffness)             ! degrading critical load carrying capacity by damage 

     udotd = &
       sign(1.0_pReal,traction_d)* &
       damage_anisoDuctile_sdot_0(instance)* &
       (abs(traction_d)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udotd /= 0.0_pReal) then
       Ld = Ld + udotd*projection_d
       dudotd_dt = udotd*damage_anisoDuctile_N(instance)/traction_d
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotd_dt*projection_d(k,l)* projection_d(m,n)
     endif                
 
     udott = &
       sign(1.0_pReal,traction_t)* &
       damage_anisoDuctile_sdot_0(instance)* &
       (abs(traction_t)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udott /= 0.0_pReal) then
       Ld = Ld + udott*projection_t
       dudott_dt = udott*damage_anisoDuctile_N(instance)/traction_t
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudott_dt*projection_t(k,l)*projection_t(m,n)
     endif
     udotn = &
       damage_anisoDuctile_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udotn /= 0.0_pReal) then
       Ld = Ld + udotn*projection_n(1:3,1:3)
       dudotn_dt = udotn*damage_anisoDuctile_N(instance)/traction_n
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotn_dt*projection_n(k,l)* projection_n(m,n)
     endif                
     index = index + 1_pInt
   enddo
 enddo
 
 dLd_dTstar = math_Plain3333to99(dLd_dTstar3333)
 
end subroutine damage_anisoDuctile_LdAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoDuctile_getFd(ipc, ip, el)
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
   damage_anisoDuctile_getFd
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoDuctile_getFd = &
   reshape(damageState(phase)% &
     state(damage_anisoDuctile_totalNslip(instance)+2: &
           damage_anisoDuctile_totalNslip(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoDuctile_getFd

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine damage_anisoDuctile_putFd(Tstar_v, dt, ipc, ip, el)
 use numerics, only: &
   residualStiffness
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance, &
   damageState
 use math, only: &
   math_mul33x33, &
   math_inv33, &
   math_Plain3333to99, &
   math_I3, &
   math_identity4th, &
   math_symmetric33, &
   math_Mandel33to6, &
   math_tensorproduct
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_sd, &
   lattice_st, &
   lattice_sn

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
 real(pReal),   dimension(3,3) :: &
   projection_d, projection_t, projection_n                                                         !< projection modes 3x3 tensor
 real(pReal),   dimension(6) :: &
   projection_d_v, projection_t_v, projection_n_v                                                   !< projection modes 3x3 vector
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
 index = 2_pInt
 Ld = 0.0_pReal
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                                      ! at which index starts my family
   do i = 1_pInt,damage_anisoDuctile_Nslip(f,instance)                                              ! process each (active) slip system in family
   
     projection_d = math_tensorproduct(lattice_sd(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_t = math_tensorproduct(lattice_st(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_n = math_tensorproduct(lattice_sn(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))

     projection_d_v(1:6) = math_Mandel33to6(math_symmetric33(projection_d(1:3,1:3)))
     projection_t_v(1:6) = math_Mandel33to6(math_symmetric33(projection_t(1:3,1:3)))
     projection_n_v(1:6) = math_Mandel33to6(math_symmetric33(projection_n(1:3,1:3)))
   
     traction_d    = dot_product(Tstar_v,projection_d_v(1:6))
     traction_t    = dot_product(Tstar_v,projection_t_v(1:6))
     traction_n    = dot_product(Tstar_v,projection_n_v(1:6))
     
     traction_crit = damage_anisoDuctile_critLoad(f,instance)* &
                     (damageState(phase)%state0(index,constituent) + residualStiffness)             ! degrading critical load carrying capacity by damage 

     udotd = &
       sign(1.0_pReal,traction_d)* &
       damage_anisoDuctile_sdot_0(instance)* &
       (abs(traction_d)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udotd /= 0.0_pReal) then
       Ld = Ld + udotd*projection_d
     endif                
 
     udott = &
       sign(1.0_pReal,traction_t)* &
       damage_anisoDuctile_sdot_0(instance)* &
       (abs(traction_t)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udott /= 0.0_pReal) then
       Ld = Ld + udott*projection_t
     endif
                     
     udotn = &
       damage_anisoDuctile_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit)**damage_anisoDuctile_N(instance)
     if (udotn /= 0.0_pReal) then
       Ld = Ld + udotn*projection_n(1:3,1:3)
     endif     
           
     index = index + 1_pInt
   enddo
 enddo
 
 damageState(phase)%state(damage_anisoDuctile_totalNslip(instance)+2: &
                          damage_anisoDuctile_totalNslip(instance)+10,constituent) = &
                            reshape(math_mul33x33(math_inv33(math_I3 - dt*Ld), &
                            damage_anisoDuctile_getFd0(ipc, ip, el)), shape=[9])                             

end subroutine damage_anisoDuctile_putFd 

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoDuctile_getFd0(ipc, ip, el)
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
   damage_anisoDuctile_getFd0
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoDuctile_getFd0 = &
   reshape(damageState(phase)% &
     subState0(damage_anisoDuctile_totalNslip(instance)+2: &
               damage_anisoDuctile_totalNslip(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoDuctile_getFd0

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage deformation gradient
!--------------------------------------------------------------------------------------------------
pure function damage_anisoDuctile_getPartionedFd0(ipc, ip, el)
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
   damage_anisoDuctile_getPartionedFd0
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)

 damage_anisoDuctile_getPartionedFd0 = &
   reshape(damageState(phase)% &
     partionedState0(damage_anisoDuctile_totalNslip(instance)+2: &
                     damage_anisoDuctile_totalNslip(instance)+10,constituent), &
           shape=[3,3])
 
end function damage_anisoDuctile_getPartionedFd0

!--------------------------------------------------------------------------------------------------
!> @brief returns damage
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_getDamage(ipc, ip, el)
 use material, only: &
   material_homog, &
   mappingHomogenization, &
   mappingConstitutive, &
   damageState, &
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
    damage_anisoDuctile_getDamage = damageState(mappingConstitutive(2,ipc,ip,el))% &
      state0(1,mappingConstitutive(1,ipc,ip,el))
    
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
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal)   :: &
   damage_anisoDuctile_getLocalDamage

 damage_anisoDuctile_getLocalDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))

end function damage_anisoDuctile_getLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns slip system damage
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_getSlipDamage(ipc, ip, el)
 use numerics, only: &
   residualStiffness
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
   damage_anisoDuctile_getSlipDamage(damage_anisoDuctile_totalNslip( &
       phase_damageInstance(mappingConstitutive(2,ipc,ip,el))))
 integer(pInt) :: &
   phase, &
   constituent, &
   instance

 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = phase_damageInstance(phase)
 
 damage_anisoDuctile_getSlipDamage = &
   damageState(phase)%state0(2:1+damage_anisoDuctile_totalNslip(instance),constituent) + &
   residualStiffness

end function damage_anisoDuctile_getSlipDamage

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function damage_anisoDuctile_postResults(ipc,ip,el)
 use material, only: &
   mappingConstitutive, &
   phase_damageInstance

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
