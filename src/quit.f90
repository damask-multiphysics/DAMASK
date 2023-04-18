!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quit subroutine
!> @details exits the program and reports current time and duration. Exit code 0 signals
!> everything is fine. Exit code 1 signals an error, message according to IO_error.
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
  use, intrinsic :: ISO_fortran_env, only: ERROR_UNIT
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif
  use HDF5

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif

  integer, intent(in) :: stop_id

  integer, dimension(8) :: dateAndTime
  integer :: err_HDF5
  integer(MPI_INTEGER_KIND) :: err_MPI, worldsize
  PetscErrorCode :: err_PETSc


  call h5open_f(err_HDF5)                                                                           ! prevents error if not opened yet
  if (err_HDF5 < 0) write(ERROR_UNIT,'(a,i5)') ' Error in h5open_f ',err_HDF5
  call h5close_f(err_HDF5)
  if (err_HDF5 < 0) write(ERROR_UNIT,'(a,i5)') ' Error in h5close_f ',err_HDF5

  call PetscFinalize(err_PETSc)

  call date_and_time(values = dateAndTime)
  write(6,'(/,a)') ' DAMASK terminated on:'
  write(6,'(a,2(i2.2,a),i4.4)') ' Date:               ',dateAndTime(3),'/',&
                                                        dateAndTime(2),'/',&
                                                        dateAndTime(1)
  write(6,'(a,2(i2.2,a),i2.2)') ' Time:               ',dateAndTime(5),':',&
                                                        dateAndTime(6),':',&
                                                        dateAndTime(7)

  if (stop_id == 0 .and. err_HDF5 == 0 .and. err_PETSC == 0) then
    call MPI_Finalize(err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI_Finalize error'
    stop 0                                                                                          ! normal termination
  else
    call MPI_Comm_size(MPI_COMM_WORLD,worldsize,err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI_Comm error'
    if (stop_id /= 0 .and. worldsize > 1) call MPI_Abort(MPI_COMM_WORLD,1,err_MPI)
    stop 1                                                                                          ! error (message from IO_error)
  endif

end subroutine quit
