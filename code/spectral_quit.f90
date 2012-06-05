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
!##################################################################################################
!* $Id$
!##################################################################################################
! Material subroutine for BVP solution using spectral method
!
! Run 'DAMASK_spectral.exe --help' to get usage hints
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts,
!            D.D. Tjahjanto,
!            C. Kords,
!            M. Diehl,
!            R. Lebensohn
!
! MPI fuer Eisenforschung, Duesseldorf
!********************************************************************
! quit subroutine to satisfy IO_error
!
!********************************************************************
subroutine quit(stop_id)
 use prec, only: &
   pInt
   
 implicit none
 integer(pInt), intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer

 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') 'DAMASK terminated on:'
 write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                      dateAndTime(2),'/',&
                                                      dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                      dateAndTime(6),':',&
                                                      dateAndTime(7)  
 if (stop_id == 0_pInt) stop 0                                                                      ! normal termination
 if (stop_id <  0_pInt) then                                                                        ! trigger regridding
   write(0,'(a,i6)') 'restart a', stop_id*(-1_pInt)
   stop 2
 endif
 stop 1                                                                                             ! error
end subroutine
