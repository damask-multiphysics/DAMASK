! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
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
!* $Id: prec.f90 1033 2011-10-20 16:46:11Z MPIE\m.diehl $
!##############################################################

MODULE prec
 implicit none 
!    *** Precision of real and integer variables for python interfacing***
 integer, parameter :: pReal = selected_real_kind(8)
 integer, parameter :: pInt  = selected_int_kind(9)           ! up to +- 1e9
 real(pReal), parameter :: DAMASK_NaN = Z'7FF0000000000001'
 real(pReal), parameter :: tol_math_check = 1.0e-8_pReal
END MODULE prec

MODULE debug
 use prec, only: pInt
 implicit none
 integer(pInt), parameter :: debug_verbosity = 0_pInt
END MODULE debug

MODULE numerics
 use prec, only: pInt
 implicit none
 real*8, parameter :: fftw_timelimit = -1.0
 integer*8, parameter :: fftw_planner_flag = 32
 integer(pInt), parameter :: fixedSeed = 1_pInt
END MODULE numerics

MODULE IO
 CONTAINS
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

END MODULE IO
