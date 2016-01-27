!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant damage field
!--------------------------------------------------------------------------------------------------
module damage_none

 implicit none
 private
 
 public :: &
   damage_none_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine damage_none_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
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
   write(6,'(/,a)')   ' <<<+-  damage_'//DAMAGE_none_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (damage_type(homog) == DAMAGE_none_ID) then
     NofMyHomog = count(material_homog == homog)
     damageState(homog)%sizeState = 0_pInt
     damageState(homog)%sizePostResults = 0_pInt
     allocate(damageState(homog)%state0   (0_pInt,NofMyHomog))
     allocate(damageState(homog)%subState0(0_pInt,NofMyHomog))
     allocate(damageState(homog)%state    (0_pInt,NofMyHomog))
     
     deallocate(damage(homog)%p)
     allocate  (damage(homog)%p(1), source=damage_initialPhi(homog))
     
   endif myhomog
 enddo initializeInstances


end subroutine damage_none_init

end module damage_none
