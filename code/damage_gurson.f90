!--------------------------------------------------------------------------------------------------
! $Id$
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
   damage_gurson_coeff_torsion, &
   damage_gurson_coeff_ten_comp, &
   damage_gurson_coeff_triaxiality, &
   damage_gurson_fracture_tough, &
   damage_gurson_lengthscale, &
   damage_gurson_crit_void_fraction
   

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 local_damage_ID
 end enum                                                                                          !!!!! ToDo
 
 integer(kind(undefined_ID)),         dimension(:,:),         allocatable,          private :: & 
   damage_gurson_outputID                                                                          !< ID of each post result output


 public :: &
   damage_gurson_init, &
   damage_gurson_stateInit, &
   damage_gurson_aTolState, &
   damage_gurson_dotState, &
   damage_gurson_microstructure, &
   damage_gurson_getDamage, &
   damage_gurson_getSlipDamage, &
   damage_gurson_putLocalDamage, &
   damage_gurson_getLocalDamage, &
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
   write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_DAMAGE_gurson_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

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
 allocate(damage_gurson_coeff_torsion(maxNinstance),                       source=0.0_pReal)
 allocate(damage_gurson_coeff_ten_comp(maxNinstance),                      source=0.0_pReal)
 allocate(damage_gurson_coeff_triaxiality(maxNinstance),                   source=0.0_pReal)
 allocate(damage_gurson_fracture_tough(maxNinstance),                      source=0.0_pReal)
 allocate(damage_gurson_lengthscale(maxNinstance),                         source=0.0_pReal)
 allocate(damage_gurson_crit_void_fraction(maxNinstance),                  source=0.0_pReal)
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
! input parameters 
       case ('coeff_torsion')
         damage_gurson_coeff_torsion(instance)  = IO_floatValue(line,positions,2_pInt)               !> coefficent of torsional stress component
         
       case ('coeff_tension_comp')
         damage_gurson_coeff_ten_comp(instance) = IO_floatValue(line,positions,2_pInt)               !> coefficent of tensile or compressive stress component
         
       case ('coeff_triaxiality')
         damage_gurson_coeff_triaxiality(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('fracture_toughness')
         damage_gurson_fracture_tough(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('lengthscale')
         damage_gurson_lengthscale(instance) = IO_floatValue(line,positions,2_pInt)
         
       case ('critical_voidFraction')
         damage_gurson_crit_void_fraction(instance) = IO_floatValue(line,positions,2_pInt)

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
     sizeDotState              =   3_pInt
     sizeState                 =   4_pInt
                
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

     call damage_gurson_stateInit(phase)
     call damage_gurson_aTolState(phase,instance)
   endif
 
 enddo initializeInstances
end subroutine damage_gurson_init

!--------------------------------------------------------------------------------------------------
!> @brief sets the relevant  NEW state values for a given instance of this damage
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_stateInit(phase)
 use material, only: &
   damageState
 
 implicit none
 integer(pInt),              intent(in) :: phase                                                    !< number specifying the phase of the damage

 real(pReal), dimension(damageState(phase)%sizeState) :: tempState

 tempState(1) = 1.0_pReal
 tempState(2) = 1.0_pReal
 tempState(3) = 1.0_pReal
 tempState(4) = 1.0_pReal

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
subroutine damage_gurson_dotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState
 use math, only: &
   math_equivStrain33, &
   math_norm33, &
   math_j3_33, &
   math_trace33, &
   math_I3, &
   math_Mandel6to33
 use lattice, only: &
   lattice_DamageMobility

 implicit none
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal), intent(in), dimension(3,3) :: &
   Lp
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 integer(pInt) :: &
   phase, constituent
 real(pReal) :: &
   i1, j2, j3
 real(pReal) , dimension(3,3) :: &
   Tstar_dev
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 Tstar_dev = math_Mandel6to33(Tstar_v) - math_trace33(math_Mandel6to33(Tstar_v))/3*math_I3
 i1 = sum(Tstar_v(1:3))
 j2 = 0.5_pReal*(math_norm33(Tstar_dev))**2
 j3 = math_j3_33(math_Mandel6to33(Tstar_v))
 
 damageState(phase)%dotState(1,constituent) = &
   (1.0_pReal/lattice_DamageMobility(phase))* &
   (damageState(phase)%state(4,constituent) - &
    damageState(phase)%state(1,constituent))


 damageState(phase)%dotState(2,constituent) = &                                                     !> void nucleation rate
   math_norm33(Lp)*sqrt(damage_gurson_lengthscale(phase))/damage_gurson_fracture_tough(phase)* &
    damageState(phase)%state(2,constituent) * ( &
             damage_gurson_coeff_torsion(phase) * ((4_pReal/27_pReal) - (j3**(2)/j2**(3))) + &
             damage_gurson_coeff_ten_comp(phase) * (j3/j2**(1.5_pReal)) + &
             damage_gurson_coeff_triaxiality(phase) * abs(i1/sqrt(j2)))                             !> to be coupled with vacancy generation
             
 damageState(phase)%dotState(3,constituent) = &
   ( damageState(phase)%state(4,constituent)) * math_trace33(Lp)                          !> void growth rate

  
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
 integer(pReal) :: &
   voidFraction
 integer(pInt) :: &
   phase, constituent
 
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 voidFraction =  damageState(phase)%state(2,constituent) + damageState(phase)%state(3,constituent)
 
 if(voidFraction < damage_gurson_crit_void_fraction(phase)) then
   damageState(phase)%state(4,constituent) =  1_pReal - voidFraction                                ! damage parameter is 1 when no void present
 else 
   damageState(phase)%state(4,constituent) =  1_pReal - damage_gurson_crit_void_fraction(phase) + &
                                  5_pReal * (voidFraction - damage_gurson_crit_void_fraction(phase)) ! this accelerated void increase models the effect of void coalescence
 endif
                                                       
end subroutine damage_gurson_microstructure

!--------------------------------------------------------------------------------------------------
!> @brief returns damage
!--------------------------------------------------------------------------------------------------
function damage_gurson_getDamage(ipc, ip, el)
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
 real(pReal) :: damage_gurson_getDamage
 
 select case(field_damage_type(material_homog(ip,el)))                                                   
   case (FIELD_DAMAGE_LOCAL_ID)
    damage_gurson_getDamage = damage_gurson_getLocalDamage(ipc, ip, el)
    
   case (FIELD_DAMAGE_NONLOCAL_ID)
    damage_gurson_getDamage =    fieldDamage(material_homog(ip,el))% &
      field(1,mappingHomogenization(1,ip,el))                                                     ! Taylor type 

 end select
 
end function damage_gurson_getDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns slip damage
!--------------------------------------------------------------------------------------------------
function damage_gurson_getSlipDamage(Tstar_v, ipc, ip, el)

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: damage_gurson_getSlipDamage, porosity
 
 porosity = damage_gurson_getDamage(ipc, ip, el)
 damage_gurson_getSlipDamage = porosity*porosity      ! Gurson yield function should go here
 
end function damage_gurson_getSlipDamage

!--------------------------------------------------------------------------------------------------
!> @brief puts local damage 
!--------------------------------------------------------------------------------------------------
subroutine damage_gurson_putLocalDamage(ipc, ip, el, localDamage)
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
 
end subroutine damage_gurson_putLocalDamage

!--------------------------------------------------------------------------------------------------
!> @brief returns local damage
!--------------------------------------------------------------------------------------------------
function damage_gurson_getLocalDamage(ipc, ip, el)
 use material, only: &
   mappingConstitutive, &
   damageState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: damage_gurson_getLocalDamage
 
 damage_gurson_getLocalDamage = &
   damageState(mappingConstitutive(2,ipc,ip,el))%state(1,mappingConstitutive(1,ipc,ip,el))
 
end function damage_gurson_getLocalDamage

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
