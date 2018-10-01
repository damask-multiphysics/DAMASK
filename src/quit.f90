!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quit subroutine
!> @details exits the program and reports current time and duration. Exit code 0 signals
!> everything is fine. Exit code 1 signals an error, message according to IO_error. Exit code
!> 2 signals no severe problems, but some increments did not converge
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
#include <petsc/finclude/petscsys.h>
#ifdef _OPENMP
 use MPI, only: &
   MPI_finalize
#endif
 use prec, only: &
   pInt
 use PetscSys

 implicit none
 integer(pInt), intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer
 integer(pInt) :: error = 0_pInt
 PetscErrorCode :: ierr = 0
 
 call PETScFinalize(ierr)
 CHKERRQ(ierr)
#ifdef _OPENMP
 call MPI_finalize(error)
 if (error /= 0) write(6,'(a)') ' Error in MPI_finalize'
#endif
 
 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') 'DAMASK terminated on:'
 write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                      dateAndTime(2),'/',&
                                                      dateAndTime(1)
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                      dateAndTime(6),':',&
                                                      dateAndTime(7)

 if (stop_id == 0_pInt .and. ierr == 0_pInt .and. error == 0_pInt) stop 0                           ! normal termination
 if (stop_id == 2_pInt .and. ierr == 0_pInt .and. error == 0_pInt) stop 2                           ! not all incs converged
 stop 1                                                                                             ! error (message from IO_error)

end subroutine quit
