!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely elastic material
!--------------------------------------------------------------------------------------------------
module plastic_none

 implicit none
 private

 public :: &
   plastic_none_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_none_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pInt
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
   PLASTICITY_NONE_label, &
   material_phase, &
   plasticState, &
   PLASTICITY_none_ID

 implicit none

 integer(pInt) :: &
   maxNinstance, &
   phase, &
   NofMyPhase
 
 write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_NONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_none_ID),pInt)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(phase) == PLASTICITY_none_ID) then
   NofMyPhase=count(material_phase==phase)

     allocate(plasticState(phase)%aTolState          (0_pInt))
     allocate(plasticState(phase)%state0             (0_pInt,NofMyPhase))
     allocate(plasticState(phase)%partionedState0    (0_pInt,NofMyPhase))
     allocate(plasticState(phase)%subState0          (0_pInt,NofMyPhase))
     allocate(plasticState(phase)%state              (0_pInt,NofMyPhase))

     allocate(plasticState(phase)%dotState           (0_pInt,NofMyPhase))
     allocate(plasticState(phase)%deltaState         (0_pInt,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (0_pInt,NofMyPhase))
       allocate(plasticState(phase)%previousDotState2(0_pInt,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (0_pInt,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,0_pInt,NofMyPhase))
   endif
 enddo initializeInstances

end subroutine plastic_none_init

end module plastic_none
