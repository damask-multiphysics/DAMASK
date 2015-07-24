!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isothermal temperature field
!--------------------------------------------------------------------------------------------------
module thermal_isothermal

 implicit none
 private
 
 public :: &
   thermal_isothermal_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine thermal_isothermal_init()
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
   NofMyHomog, &
   sizeState

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  thermal_'//THERMAL_isothermal_label//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (thermal_type(homog) == THERMAL_isothermal_ID) then
     NofMyHomog = count(material_homog == homog)
     sizeState = 0_pInt
     thermalState(homog)%sizeState = sizeState
     thermalState(homog)%sizePostResults = sizeState
     allocate(thermalState(homog)%state0   (sizeState,NofMyHomog), source=0.0_pReal)
     allocate(thermalState(homog)%subState0(sizeState,NofMyHomog), source=0.0_pReal)
     allocate(thermalState(homog)%state    (sizeState,NofMyHomog), source=0.0_pReal)
     
     deallocate(temperature    (homog)%p)
     allocate  (temperature    (homog)%p(1), source=thermal_initialT(homog))
     deallocate(temperatureRate(homog)%p)
     allocate  (temperatureRate(homog)%p(1), source=0.0_pReal)

   endif myhomog
 enddo initializeInstances


end subroutine thermal_isothermal_init

end module thermal_isothermal
