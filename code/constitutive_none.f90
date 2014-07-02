!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely elastic material
!--------------------------------------------------------------------------------------------------
module constitutive_none
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   constitutive_none_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   constitutive_none_sizePostResult                                                                 !< size of each post result output

 public :: &
   constitutive_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine constitutive_none_init(fileUnit)
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
   phase_plasticity, &
   phase_Noutput, &
   PLASTICITY_NONE_label, &
   material_phase, &
   plasticState, &
   PLASTICITY_none_ID, &
   MATERIAL_partPhase

 implicit none

 integer(pInt), intent(in) :: fileUnit
 integer(pInt) :: &
   maxNinstance, &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_NONE_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_none_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(phase) == PLASTICITY_none_ID) then
   NofMyPhase=count(material_phase==phase)

     sizeState    = 0_pInt
     plasticState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizePostResults = 0_pInt
     allocate(plasticState(phase)%aTolState          (sizeState))
     allocate(plasticState(phase)%state0             (sizeState,NofMyPhase))
     allocate(plasticState(phase)%partionedState0    (sizeState,NofMyPhase))
     allocate(plasticState(phase)%subState0          (sizeState,NofMyPhase))
     allocate(plasticState(phase)%state              (sizeState,NofMyPhase))
     allocate(plasticState(phase)%state_backup       (sizeState,NofMyPhase))

     allocate(plasticState(phase)%dotState           (sizeDotState,NofMyPhase))
     allocate(plasticState(phase)%dotState_backup    (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NofMyPhase))
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NofMyPhase))
   endif
 enddo initializeInstances

 allocate(constitutive_none_sizePostResults(maxNinstance), source=0_pInt)

end subroutine constitutive_none_init

end module constitutive_none
