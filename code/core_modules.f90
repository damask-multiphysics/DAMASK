! Copyright 2011,2012 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$
!##############################################################

module prec

 implicit none 
!    *** Precision of real and integer variables for python interfacing***
 integer, parameter :: pReal = 8
 integer, parameter :: pInt  = 4
 real(pReal), parameter :: DAMASK_NaN = real(Z'7FF0000000000001',pReal)
 real(pReal), parameter :: tol_math_check = 1.0e-8_pReal
 
end module prec


module debug

 use prec, only: pInt
 
 implicit none
 integer(pInt), parameter, public :: &
   debug_levelBasic = 2_pInt**1_pInt, &
   debug_math       = 2_pInt
 integer(pInt), dimension(11+2),  public :: &
   debug_what                    = debug_levelBasic
   
end module debug


module numerics

 use prec, only: pInt, pReal
 
 implicit none
 real(pReal), parameter   :: fftw_timelimit    = -1.0_pReal
 integer(pInt), parameter :: fftw_planner_flag = 32_pInt
 integer(pInt), parameter :: fixedSeed         = 1_pInt
 
end module numerics


module IO
 contains
 
 subroutine IO_error(error_ID,e,i,g,ext_msg)

 use prec, only: pInt
 
 implicit none
 integer(pInt), intent(in) :: error_ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=1024) msg

 select case (error_ID)
   case default
     print*, 'Error messages not supported when interfacing to Python'
 end select
 end subroutine IO_error

end module IO
