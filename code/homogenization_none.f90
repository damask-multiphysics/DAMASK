!--------------------------------------------------------------------------------------------------
! $Id$
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
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pReal, &
   pInt 
 use IO, only: &
   IO_timeStamp
 use material
 use numerics, only: &
   worldrank
 
 implicit none
 integer(pInt) :: &
   homog, &
   NofMyHomog

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_NONE_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (homogenization_type(homog) == HOMOGENIZATION_none_ID) then
     NofMyHomog = count(material_homog == homog)
     homogState(homog)%sizeState = 0_pInt
     homogState(homog)%sizePostResults = 0_pInt
     allocate(homogState(homog)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (0_pInt,NofMyHomog), source=0.0_pReal)

   endif myhomog
 enddo initializeInstances


end subroutine homogenization_none_init

end module homogenization_none
