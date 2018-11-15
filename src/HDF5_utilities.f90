!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module HDF5_Utilities
  use prec
  use IO
  use HDF5
#ifdef PETSc
  use PETSC
#endif

 implicit none
 private
 integer(HID_T), public, protected :: tempCoordinates, tempResults
 integer(HID_T), private :: resultsFile, currentIncID, plist_id
 integer(pInt),  private :: currentInc

!--------------------------------------------------------------------------------------------------
!> @brief reads pInt or pReal data of defined shape from file
!--------------------------------------------------------------------------------------------------
 interface HDF5_read
   module procedure HDF5_read_pReal_1
   module procedure HDF5_read_pReal_2
   module procedure HDF5_read_pReal_3
   module procedure HDF5_read_pReal_4
   module procedure HDF5_read_pReal_5
   module procedure HDF5_read_pReal_6
   module procedure HDF5_read_pReal_7

   module procedure HDF5_read_pInt_1
   module procedure HDF5_read_pInt_2
   module procedure HDF5_read_pInt_3
   module procedure HDF5_read_pInt_4
   module procedure HDF5_read_pInt_5
   module procedure HDF5_read_pInt_6
   module procedure HDF5_read_pInt_7  !ABOVE 8 DIMENSIONS IT GIVES ERROR: THE CALL TO H5DREAD_F DOESNT WORK

 end interface HDF5_read

!--------------------------------------------------------------------------------------------------
!> @brief writes pInt or pReal data of defined shape to file
!--------------------------------------------------------------------------------------------------
 interface HDF5_write
   module procedure HDF5_write_pReal1
   module procedure HDF5_write_pReal2
   module procedure HDF5_write_pReal3
   module procedure HDF5_write_pReal4
   module procedure HDF5_write_pReal5
   module procedure HDF5_write_pReal6
   module procedure HDF5_write_pReal7

   module procedure HDF5_write_pInt1
   module procedure HDF5_write_pInt2
   module procedure HDF5_write_pInt3
   module procedure HDF5_write_pInt4
   module procedure HDF5_write_pInt5
   module procedure HDF5_write_pInt6
   module procedure HDF5_write_pInt7  !ABOVE 8 DIMENSIONS IT GIVES ERROR: THE CALL TO H5DREAD_F DOESNT WORK

 end interface HDF5_write

 public :: &
   HDF5_utilities_init, &
   HDF5_closeGroup ,&
   HDF5_openGroup2, &
   HDF5_closeFile, &
   HDF5_addGroup2, &
   HDF5_openFile, &
   HDF5_read, &
   HDF5_write
contains

subroutine HDF5_Utilities_init
 use, intrinsic :: &
   iso_fortran_env  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

 implicit none
 integer              :: hdferr
 integer(SIZE_T)      :: typeSize

 write(6,'(/,a)') ' <<<+-  HDF5_Utilities init  -+>>>'
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
!initialize HDF5 library and check if integer and float type size match
 call h5open_f(hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5open_f')
 call h5tget_size_f(H5T_NATIVE_INTEGER,typeSize, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5tget_size_f (int)')
 if (int(pInt,SIZE_T)/=typeSize) call IO_error(0_pInt,ext_msg='pInt does not match H5T_NATIVE_INTEGER')
 call h5tget_size_f(H5T_NATIVE_DOUBLE,typeSize, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5tget_size_f (double)')
 if (int(pReal,SIZE_T)/=typeSize) call IO_error(0_pInt,ext_msg='pReal does not match H5T_NATIVE_DOUBLE')

end subroutine HDF5_Utilities_init




!--------------------------------------------------------------------------------------------------
!> @brief creates and initializes HDF5 output files
!--------------------------------------------------------------------------------------------------
 integer(HID_T) function HDF5_createFile(path)
 use hdf5
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer                       :: hdferr
 integer(SIZE_T)               :: typeSize
 character(len=*), intent(in)  :: path
#ifdef PETSc
#include <petsc/finclude/petscsys.h>
#endif
  call h5open_f(hdferr) !############################################################ DANGEROUS
#ifdef PETSc
 call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pcreate_f')
 call h5pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_Utilities_init: h5pset_fapl_mpio_f')
#endif
!--------------------------------------------------------------------------------------------------
! create a file
 !call h5fcreate_f(path,H5F_ACC_TRUNC_F,resultsFile,hdferr)
 call h5fcreate_f(path,H5F_ACC_TRUNC_F,HDF5_createFile,hdferr,access_prp = plist_id)
 if (hdferr < 0) call IO_error(100_pInt,ext_msg=path)
 !call HDF5_addStringAttribute(HDF5_createFile,'createdBy',DAMASKVERSION)
 call h5pclose_f(plist_id, hdferr)  !neu

end function HDF5_createFile



!--------------------------------------------------------------------------------------------------
!> @brief open and initializes HDF5 output file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openFile(fileName,mode)

 implicit none
 character(len=*), intent(in)           :: fileName
 character,        intent(in), optional :: mode
 character                              :: m
 integer :: hdferr

 if (present(mode)) then
   m = mode
 else
   m = 'r'
 endif

 if (m == 'w') then
   call h5fcreate_f(fileName,H5F_ACC_TRUNC_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fcreate_f',el=hdferr)
 elseif(m == 'a') then
   call h5fopen_f(fileName,H5F_ACC_RDWR_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (a)',el=hdferr)
 elseif(m == 'r') then
   call h5fopen_f(fileName,H5F_ACC_RDONLY_F,HDF5_openFile,hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (r)',el=hdferr)
 else
   call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f unknown access mode',el=hdferr)
 endif

end function HDF5_openFile


!--------------------------------------------------------------------------------------------------
!> @brief close the opened HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeFile(fileHandle)

 implicit none
 integer :: hdferr
 integer(HID_T), intent(in) :: fileHandle
 call h5fclose_f(fileHandle,hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_closeFile: h5fclose_f',el=hdferr)

end subroutine HDF5_closeFile


!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the fileHandle (additional to addGroup2)
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup2(fileHandle,groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer(HID_T), intent(in)   :: fileHandle
 integer                      :: hdferr

 call h5gcreate_f(fileHandle, trim(groupName), HDF5_addGroup2, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup2: h5gcreate_f ('//trim(groupName)//')')

end function HDF5_addGroup2


!--------------------------------------------------------------------------------------------------
!> @brief open an existing group of a file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup2(FileReadID,groupName)
 use hdf5

 implicit none
 character(len=*), intent(in) :: groupName
 integer                      :: hdferr
 integer(HID_T), intent(in)   :: FileReadID

 call h5gopen_f(FileReadID, trim(groupName), HDF5_openGroup2, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup2: h5gopen_f ('//trim(groupName)//')')

end function HDF5_openGroup2


!--------------------------------------------------------------------------------------------------
!> @brief close a group
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeGroup(ID)
 use hdf5

 implicit none
 integer(HID_T), intent(in) :: ID
 integer                    :: hdferr

 call h5gclose_f(ID, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_closeGroup: h5gclose_f (el is ID)', el = int(ID,pInt))

end subroutine HDF5_closeGroup


!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_1(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape
 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape1: h5dclose_f')

end subroutine HDF5_read_pReal_1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_2(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape2: h5dclose_f')

end subroutine HDF5_read_pReal_2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_3(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape
 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape3: h5dclose_f')

end subroutine HDF5_read_pReal_3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_4(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape4: h5dclose_f')

end subroutine HDF5_read_pReal_4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_5(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape5: h5dclose_f')

end subroutine HDF5_read_pReal_5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_6(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape6: h5dclose_f')

end subroutine HDF5_read_pReal_6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pReal_7(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pReal_shape7: h5dclose_f')

end subroutine HDF5_read_pReal_7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_1(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape1: h5dclose_f')

end subroutine HDF5_read_pInt_1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_2(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape2: h5dclose_f')

end subroutine HDF5_read_pInt_2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_3(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape3: h5dclose_f')

end subroutine HDF5_read_pInt_3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_4(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape4: h5dclose_f')

end subroutine HDF5_read_pInt_4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_5(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape5: h5dclose_f')

end subroutine HDF5_read_pInt_5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_6(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape6: h5dclose_f')

end subroutine HDF5_read_pInt_6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for reading dataset of the type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_pInt_7(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file
 integer(pInt),dimension(:), allocatable  :: myShape

 integer        :: hdferr
 integer(HID_T) :: dset_id
 myShape = shape(dataset)

 call h5dopen_f(loc_id,datasetName,dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dopen_f')
 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T),hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dread_f')
 call h5dclose_f(dset_id,hdferr)
 if (hdferr /= 0) call IO_error(0_pInt,ext_msg='HDF5_read_pInt__shape7: h5dclose_f')

end subroutine HDF5_read_pInt_7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 1 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal1(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal1: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal1: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal1: h5sclose_f')

end subroutine HDF5_write_pReal1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal2(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal2: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal2: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal2: h5sclose_f')

end subroutine HDF5_write_pReal2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal3(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal3: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal3: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal3: h5sclose_f')

end subroutine HDF5_write_pReal3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal4(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal4: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal4: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal4: h5sclose_f')

end subroutine HDF5_write_pReal4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal5(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal5: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal5: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal5: h5sclose_f')

end subroutine HDF5_write_pReal5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal6(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal6: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal6: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal6: h5sclose_f')

end subroutine HDF5_write_pReal6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pReal7(dataset,loc_id,datasetName)

 implicit none
 real(pReal),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(myShape,HSIZE_T), hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pReal7: h5dcreate_f')

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal7: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pReal7: h5sclose_f')

end subroutine HDF5_write_pReal7

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 1 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt1(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt1: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt1: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt1: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt1: h5sclose_f')

end subroutine HDF5_write_pInt1

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt2(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt2: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt2: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt2: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt2: h5sclose_f')

end subroutine HDF5_write_pInt2

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt3(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt3: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt3: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt3: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt3: h5sclose_f')

end subroutine HDF5_write_pInt3

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt4(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt4: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt4: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt4: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt4: h5sclose_f')

end subroutine HDF5_write_pInt4

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt5(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt5: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt5: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt5: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt5: h5sclose_f')

end subroutine HDF5_write_pInt5

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt6(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt6: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt6: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt6: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt6: h5sclose_f')

end subroutine HDF5_write_pInt6

!--------------------------------------------------------------------------------------------------
!> @brief subroutine for writing dataset of the type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_pInt7(dataset,loc_id,datasetName)

 implicit none
 integer(pInt),      intent(out), dimension(:,:,:,:,:,:,:) ::    dataset
 integer(HID_T),   intent(in) :: loc_id                                                             !< file or group handle
 character(len=*), intent(in) :: datasetName                                                        !< name of the dataset in the file

 integer(pInt), dimension(:), allocatable :: myShape                                                          !<shape of the dataset
 integer :: hdferr
 integer(HID_T) :: dset_id, space_id

 myShape = shape(dataset)
!--------------------------------------------------------------------------------------------------
! create dataspace
call h5screate_simple_f(size(myShape), int(myShape,HSIZE_T), space_id, hdferr, &
                           int(myShape,HSIZE_T))
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt7: h5screate_simple_f')
!--------------------------------------------------------------------------------------------------
! create dataset
call h5dcreate_f(loc_id, trim(datasetName), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_pInt7: h5dcreate_f')

CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(myShape,HSIZE_T), hdferr)

!--------------------------------------------------------------------------------------------------
!close types, dataspaces
 call h5dclose_f(dset_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt7: h5dclose_f')
call h5sclose_f(space_id, hdferr)
if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_HDF5_write_pInt7: h5sclose_f')

end subroutine HDF5_write_pInt7


end module HDF5_Utilities
