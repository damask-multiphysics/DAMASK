!--------------------------------------------------------------------------------------------------
! $Id: damage_none.f90 3148 2014-05-27 14:46:03Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely elastic material
!--------------------------------------------------------------------------------------------------
module damage_none
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   damage_none_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   damage_none_sizePostResult                                                                 !< size of each post result output

 public :: &
   damage_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_none_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use IO, only: &
   IO_timeStamp
 use numerics, only: &
   numerics_integrator
 use material, only: &
   phase_damage, &
   phase_Noutput, &
   LOCAL_DAMAGE_NONE_label, &
   LOCAL_DAMAGE_NONE_ID, &
   material_phase, &
   damageState, &
   MATERIAL_partPhase

 implicit none

 integer(pInt), intent(in) :: fileUnit
 integer(pInt) :: &
   maxNinstance, &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState
 write(6,'(/,a)')   ' <<<+-  damage_'//LOCAL_DAMAGE_NONE_label//' init  -+>>>'
 write(6,'(a)')     ' $Id: damage_none.f90 3148 2014-05-27 14:46:03Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 maxNinstance = int(count(phase_damage == LOCAL_DAMAGE_NONE_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 initializeInstances: do phase = 1_pInt, size(phase_damage)
   NofMyPhase=count(material_phase==phase)
   if (phase_damage(phase) == LOCAL_DAMAGE_none_ID) then
     sizeState    = 0_pInt
     damageState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     damageState(phase)%sizeDotState = sizeDotState
     damageState(phase)%sizePostResults = 0_pInt
     allocate(damageState(phase)%state0         (sizeState,NofMyPhase))
     allocate(damageState(phase)%partionedState0(sizeState,NofMyPhase))
     allocate(damageState(phase)%subState0      (sizeState,NofMyPhase))
     allocate(damageState(phase)%state          (sizeState,NofMyPhase))
     allocate(damageState(phase)%state_backup   (sizeState,NofMyPhase))
     allocate(damageState(phase)%aTolState      (NofMyPhase))
     allocate(damageState(phase)%dotState       (sizeDotState,NofMyPhase))
     allocate(damageState(phase)%dotState_backup(sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(damageState(phase)%previousDotState  (sizeDotState,NofMyPhase))
       allocate(damageState(phase)%previousDotState2 (sizeDotState,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(damageState(phase)%RK4dotState       (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(damageState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase))
   endif
 enddo initializeInstances

 allocate(damage_none_sizePostResults(maxNinstance), source=0_pInt)

end subroutine damage_none_init

end module damage_none
