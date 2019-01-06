!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Dummy plasticity for purely elastic material
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
 use material, only: &
   phase_plasticity, &
   material_allocatePlasticState, &
   PLASTICITY_NONE_label, &
   PLASTICITY_NONE_ID, &
   material_phase, &
   plasticState

 implicit none
 integer(pInt) :: &
   Ninstance, &
   p, &
   NipcMyPhase

 write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_NONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 Ninstance = int(count(phase_plasticity == PLASTICITY_NONE_ID),pInt)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

 do p = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_NONE_ID) cycle

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phase == p)

   call material_allocatePlasticState(p,NipcMyPhase,0_pInt,0_pInt,0_pInt, &
                                      0_pInt,0_pInt,0_pInt)
   plasticState(p)%sizePostResults = 0_pInt

 enddo

end subroutine plastic_none_init

end module plastic_none
