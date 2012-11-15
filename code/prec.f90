! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief setting precision for real and int type depending on makros "FLOAT" and "INT"
!--------------------------------------------------------------------------------------------------

module prec

 implicit none
 private 

#if (FLOAT==4)
#ifdef Spectral
 SPECTRAL SOLVER DOES NOT SUPPORT SINGLE PRECISION, STOPPING COMPILATION
#endif
 integer,     parameter, public :: pReal = 4                                                        !< floating point single precition (was selected_real_kind(6,37), number with 6 significant digits, up to 1e+-37)
#ifdef __INTEL_COMPILER
 real(pReal), parameter, public :: DAMASK_NaN = Z'7F800001'                                         !< quiet NaN for single precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#ifdef __GFORTRAN__
 real(pReal), parameter, public :: DAMASK_NaN = real(Z'7F800001', pReal)                            !< quiet NaN for single precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#elif (FLOAT==8)
 integer,     parameter, public :: pReal = 8                                                        !< floating point double precision (was selected_real_kind(15,300), number with 15 significant digits, up to 1e+-300)
#ifdef __INTEL_COMPILER
 real(pReal), parameter, public :: DAMASK_NaN = Z'7FF8000000000000'                                 !< quiet NaN for double precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#ifdef __GFORTRAN__
 real(pReal), parameter, public :: DAMASK_NaN = real(Z'7FF8000000000000', pReal)                    !< quiet NaN for double precision (from http://www.hpc.unimelb.edu.au/doc/f90lrm/dfum_035.html, copy can be found in documentation/Code/Fortran)
#endif
#else
 NO SUITABLE PRECISION SELECTED, STOPPING COMPILATION
#endif

#if (INT==4)
 integer,     parameter, public :: pInt  = 4                                                        !< integer representation 32 bit (was selected_int_kind(9), number with at least up to +- 1e9)
#elif (INT==8)
 integer,     parameter, public :: pInt  = 8                                                        !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
#else
 NO SUITABLE PRECISION SELECTED, STOPPING COMPILATION
#endif

 integer,     parameter, public :: pLongInt  = 8                                                    !< integer representation 64 bit (was selected_int_kind(12), number with at least up to +- 1e12)
 real(pReal), parameter, public :: tol_math_check = 1.0e-8_pReal
 real(pReal), parameter, public :: tol_gravityNodePos = 1.0e-100_pReal
 

 type, public :: p_vec
   real(pReal), dimension(:), pointer :: p
 end type p_vec

 public :: prec_init
 
contains

!--------------------------------------------------------------------------------------------------
!> @brief reporting precision and checking if DAMASK_NaN is set correctly
!--------------------------------------------------------------------------------------------------
subroutine prec_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 
 implicit none

!$OMP CRITICAL (write2out)

 write(6,*)
 write(6,*) '<<<+-  prec init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 write(6,'(a,i3)')    ' Bytes for pReal:    ',pReal
 write(6,'(a,i3)')    ' Bytes for pInt:     ',pInt
 write(6,'(a,i3)')    ' Bytes for pLongInt: ',pLongInt
 write(6,'(a,e10.3)') ' NaN:           ',     DAMASK_NaN
 write(6,'(a,l3)')    ' NaN /= NaN:         ',DAMASK_NaN/=DAMASK_NaN
 write(6,*)
 if (DAMASK_NaN == DAMASK_NaN) call quit(9000)
!$OMP END CRITICAL (write2out)

end subroutine prec_init

end module prec
