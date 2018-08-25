!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief dummy homogenization homogenization scheme
!--------------------------------------------------------------------------------------------------
module homogenization_none

 implicit none
 private
 
 public :: &
   homogenization_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_none_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pReal, &
   pInt 
 use IO, only: &
   IO_timeStamp
 use material
 use config
 
 implicit none
 integer(pInt) :: &
   h, &
   NofMyHomog

 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_NONE_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 do h = 1_pInt, size(homogenization_type)
   if (homogenization_type(h) /= HOMOGENIZATION_NONE_ID) cycle
   
   NofMyHomog = count(material_homog == h)
   homogState(h)%sizeState = 0_pInt
   homogState(h)%sizePostResults = 0_pInt
   allocate(homogState(h)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
   allocate(homogState(h)%state    (0_pInt,NofMyHomog), source=0.0_pReal)

 enddo

end subroutine homogenization_none_init

end module homogenization_none
