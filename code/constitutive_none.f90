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
#ifndef NEWSTATE
   constitutive_none_sizeDotState, &
   constitutive_none_sizeState, &
#endif
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

 use material, only: &
   phase_plasticity, &
   phase_Noutput, &
   PLASTICITY_NONE_label, &
#ifdef NEWSTATE
   material_phase, &
   plasticState, &
#endif
   PLASTICITY_NONE_ID, &
   MATERIAL_partPhase

 implicit none

 integer(pInt), intent(in) :: fileUnit
 integer(pInt) :: &
   maxNinstance, &
   phase
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_NONE_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_NONE_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

#ifdef NEWSTATE
 initializeInstances: do phase = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(phase) == PLASTICITY_none_ID .and. count(material_phase==phase)/=0) &
     plasticState(phase)%sizeState = 0_pInt
 enddo initializeInstances
#else
 allocate(constitutive_none_sizeDotState(maxNinstance),    source=1_pInt)
 allocate(constitutive_none_sizeState(maxNinstance),       source=1_pInt)
#endif
 allocate(constitutive_none_sizePostResults(maxNinstance), source=0_pInt)

end subroutine constitutive_none_init

end module constitutive_none
