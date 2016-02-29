!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant vacancy concentration
!--------------------------------------------------------------------------------------------------
module vacancyflux_isoconc

 implicit none
 private
 
 public :: &
   vacancyflux_isoconc_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine vacancyflux_isoconc_init()
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
   write(6,'(/,a)')   ' <<<+-  vacancyflux_'//VACANCYFLUX_isoconc_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (vacancyflux_type(homog) == VACANCYFLUX_isoconc_ID) then
     NofMyHomog = count(material_homog == homog)
     vacancyfluxState(homog)%sizeState = 0_pInt
     vacancyfluxState(homog)%sizePostResults = 0_pInt
     allocate(vacancyfluxState(homog)%state0   (0_pInt,NofMyHomog))
     allocate(vacancyfluxState(homog)%subState0(0_pInt,NofMyHomog))
     allocate(vacancyfluxState(homog)%state    (0_pInt,NofMyHomog))
     
     deallocate(vacancyConc    (homog)%p)
     allocate  (vacancyConc    (homog)%p(1), source=vacancyflux_initialCv(homog))
     deallocate(vacancyConcRate(homog)%p)
     allocate  (vacancyConcRate(homog)%p(1), source=0.0_pReal)

   endif myhomog
 enddo initializeInstances


end subroutine vacancyflux_isoconc_init

end module vacancyflux_isoconc
