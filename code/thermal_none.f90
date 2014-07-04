!--------------------------------------------------------------------------------------------------
! $Id: thermal_none.f90 3148 2014-05-27 14:46:03Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely elastic material
!--------------------------------------------------------------------------------------------------
module thermal_none
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   thermal_none_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   thermal_none_sizePostResult                                                                 !< size of each post result output

 public :: &
   thermal_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine thermal_none_init(fileUnit)
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
   phase_thermal, &
   phase_Noutput, &
   THERMAL_NONE_label, &
   material_phase, &
   thermalState, &
   THERMAL_NONE_ID, &
   MATERIAL_partPhase

 implicit none

 integer(pInt), intent(in) :: fileUnit
 integer(pInt) :: &
   maxNinstance, &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState
 
 write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_NONE_label//' init  -+>>>'
 write(6,'(a)')     ' $Id: thermal_none.f90 3148 2014-05-27 14:46:03Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_thermal == THERMAL_NONE_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
   
 initializeInstances: do phase = 1_pInt, size(phase_thermal)
   NofMyPhase=count(material_phase==phase)
   if (phase_thermal(phase) == THERMAL_none_ID .and. NofMyPhase/=0) then
     sizeState    = 0_pInt
     thermalState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     thermalState(phase)%sizeDotState = sizeDotState
     thermalState(phase)%sizePostResults = 0_pInt
     allocate(thermalState(phase)%state0         (sizeState,NofMyPhase))
     allocate(thermalState(phase)%partionedState0(sizeState,NofMyPhase))
     allocate(thermalState(phase)%subState0      (sizeState,NofMyPhase))
     allocate(thermalState(phase)%state          (sizeState,NofMyPhase))
     allocate(thermalState(phase)%state_backup   (sizeState,NofMyPhase))
     allocate(thermalState(phase)%aTolState      (NofMyPhase))
     allocate(thermalState(phase)%dotState       (sizeDotState,NofMyPhase))
     allocate(thermalState(phase)%dotState_backup(sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(thermalState(phase)%previousDotState  (sizeDotState,NofMyPhase))
       allocate(thermalState(phase)%previousDotState2 (sizeDotState,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(thermalState(phase)%RK4dotState       (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(thermalState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase))
   endif
 enddo initializeInstances
 allocate(thermal_none_sizePostResults(maxNinstance), source=0_pInt)

end subroutine thermal_none_init

end module thermal_none
