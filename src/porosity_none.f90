!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant porosity
!--------------------------------------------------------------------------------------------------
module porosity_none

 implicit none
 private
 
 public :: &
   porosity_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine porosity_none_init()
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
   write(6,'(/,a)')   ' <<<+-  porosity_'//POROSITY_none_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (porosity_type(homog) == POROSITY_none_ID) then
     NofMyHomog = count(material_homog == homog)
     porosityState(homog)%sizeState = 0_pInt
     porosityState(homog)%sizePostResults = 0_pInt
     allocate(porosityState(homog)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(porosityState(homog)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(porosityState(homog)%state    (0_pInt,NofMyHomog), source=0.0_pReal)
     
     deallocate(porosity(homog)%p)
     allocate  (porosity(homog)%p(1), source=porosity_initialPhi(homog))
     
   endif myhomog
 enddo initializeInstances


end subroutine porosity_none_init

end module porosity_none
