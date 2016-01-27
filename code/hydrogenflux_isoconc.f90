!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for constant hydrogen concentration
!--------------------------------------------------------------------------------------------------
module hydrogenflux_isoconc

 implicit none
 private
 
 public :: &
   hydrogenflux_isoconc_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine hydrogenflux_isoconc_init()
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
   write(6,'(/,a)')   ' <<<+-  hydrogenflux_'//HYDROGENFLUX_isoconc_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

  initializeInstances: do homog = 1_pInt, material_Nhomogenization
   
   myhomog: if (hydrogenflux_type(homog) == HYDROGENFLUX_isoconc_ID) then
     NofMyHomog = count(material_homog == homog)
     hydrogenfluxState(homog)%sizeState = 0_pInt
     hydrogenfluxState(homog)%sizePostResults = 0_pInt
     allocate(hydrogenfluxState(homog)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(hydrogenfluxState(homog)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(hydrogenfluxState(homog)%state    (0_pInt,NofMyHomog), source=0.0_pReal)
     
     deallocate(hydrogenConc    (homog)%p)
     deallocate(hydrogenConcRate(homog)%p)
     allocate  (hydrogenConc    (homog)%p(1), source=hydrogenflux_initialCh(homog))
     allocate  (hydrogenConcRate(homog)%p(1), source=0.0_pReal)

   endif myhomog
 enddo initializeInstances


end subroutine hydrogenflux_isoconc_init

end module hydrogenflux_isoconc
