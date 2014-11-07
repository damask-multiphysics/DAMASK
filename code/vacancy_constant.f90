!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant vacancy concentration
!--------------------------------------------------------------------------------------------------
module vacancy_constant
 use prec, only: &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,          public, protected :: &
   vacancy_constant_sizePostResults

 integer(pInt),                       dimension(:,:),   allocatable, target,  public :: &
   vacancy_constant_sizePostResult                                                                 !< size of each post result output

 public :: &
   vacancy_constant_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine vacancy_constant_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use IO, only: &
   IO_timeStamp
 use numerics, only: &
   worldrank, &
   numerics_integrator
 use material, only: &
   phase_vacancy, &
   phase_Noutput, &
   LOCAL_VACANCY_CONSTANT_label, &
   LOCAL_VACANCY_CONSTANT_ID, &
   material_phase, &
   vacancyState, &
   MATERIAL_partPhase

 implicit none

 integer(pInt) :: &
   maxNinstance, &
   phase, &
   NofMyPhase, &
   sizeState, &
   sizeDotState

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  vacancy_'//LOCAL_VACANCY_CONSTANT_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_vacancy == LOCAL_VACANCY_CONSTANT_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
   
 initializeInstances: do phase = 1_pInt, size(phase_vacancy)
   NofMyPhase=count(material_phase==phase)
   
   if (phase_vacancy(phase) == LOCAL_VACANCY_CONSTANT_ID) then
     sizeState    = 0_pInt
     vacancyState(phase)%sizeState = sizeState
     sizeDotState = sizeState
     vacancyState(phase)%sizeDotState = sizeDotState
     vacancyState(phase)%sizePostResults = 0_pInt
     allocate(vacancyState(phase)%state0         (sizeState,NofMyPhase))
     allocate(vacancyState(phase)%partionedState0(sizeState,NofMyPhase))
     allocate(vacancyState(phase)%subState0      (sizeState,NofMyPhase))
     allocate(vacancyState(phase)%state          (sizeState,NofMyPhase))
     allocate(vacancyState(phase)%state_backup   (sizeState,NofMyPhase))
     allocate(vacancyState(phase)%aTolState      (NofMyPhase))
     allocate(vacancyState(phase)%dotState       (sizeDotState,NofMyPhase))
     allocate(vacancyState(phase)%dotState_backup(sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 1_pInt)) then
       allocate(vacancyState(phase)%previousDotState  (sizeDotState,NofMyPhase))
       allocate(vacancyState(phase)%previousDotState2 (sizeDotState,NofMyPhase))
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(vacancyState(phase)%RK4dotState       (sizeDotState,NofMyPhase))
     if (any(numerics_integrator == 5_pInt)) &
       allocate(vacancyState(phase)%RKCK45dotState    (6,sizeDotState,NofMyPhase))
   endif
 enddo initializeInstances
 allocate(vacancy_constant_sizePostResults(maxNinstance), source=0_pInt)

end subroutine vacancy_constant_init

end module vacancy_constant
