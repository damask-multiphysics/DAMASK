!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!--------------------------------------------------------------------------------------------------
module HDF5_utilities
  use IO
  use HDF5
#ifdef PETSC
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
  use MPI, only: MPI_INFO_NULL_F90 => MPI_INFO_NULL
#endif
#endif

  use prec
  use parallelization
  use misc

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

!--------------------------------------------------------------------------------------------------
!> @brief Read integer or real data of defined shape from file.
!> @details for parallel IO, all dimension except for the last need to match
!--------------------------------------------------------------------------------------------------
  interface HDF5_read
    module procedure HDF5_read_real1
    module procedure HDF5_read_real2
    module procedure HDF5_read_real3
    module procedure HDF5_read_real4
    module procedure HDF5_read_real5
    module procedure HDF5_read_real6
    module procedure HDF5_read_real7

    module procedure HDF5_read_int1
    module procedure HDF5_read_int2
    module procedure HDF5_read_int3
    module procedure HDF5_read_int4
    module procedure HDF5_read_int5
    module procedure HDF5_read_int6
    module procedure HDF5_read_int7
  end interface HDF5_read

!--------------------------------------------------------------------------------------------------
!> @brief Write integer or real data of defined shape to file.
!> @details for parallel IO, all dimension except for the last need to match
!--------------------------------------------------------------------------------------------------
  interface HDF5_write
#if defined(__GFORTRAN__)
    module procedure HDF5_write_real1
    module procedure HDF5_write_real2
    module procedure HDF5_write_real3
    module procedure HDF5_write_real4
    module procedure HDF5_write_real5
    module procedure HDF5_write_real6
    module procedure HDF5_write_real7
    module procedure HDF5_write_int1
    module procedure HDF5_write_int2
    module procedure HDF5_write_int3
    module procedure HDF5_write_int4
    module procedure HDF5_write_int5
    module procedure HDF5_write_int6
    module procedure HDF5_write_int7
#else
    module procedure HDF5_write_real
    module procedure HDF5_write_int
#endif
  end interface HDF5_write

!--------------------------------------------------------------------------------------------------
!> @brief Attach attribute of type char, integer, or real to a file/dataset/group.
!--------------------------------------------------------------------------------------------------
  interface HDF5_addAttribute
    module procedure HDF5_addAttribute_str
    module procedure HDF5_addAttribute_int
    module procedure HDF5_addAttribute_real
    module procedure HDF5_addAttribute_str_array
    module procedure HDF5_addAttribute_int_array
    module procedure HDF5_addAttribute_real_array
  end interface HDF5_addAttribute

#ifdef PETSC
  logical, parameter :: parallel_default = .true.
#else
  logical, parameter :: parallel_default = .false.
#endif
  logical :: compression_possible

  public :: &
    HDF5_utilities_init, &
    HDF5_chkerr, &
    HDF5_read, &
    HDF5_write, &
    HDF5_write_str, &
    HDF5_addAttribute, &
    HDF5_addGroup, &
    HDF5_openGroup, &
    HDF5_closeGroup, &
    HDF5_openFile, &
    HDF5_closeFile, &
    HDF5_objectExists, &
    HDF5_setLink

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize HDF5 libary and perform sanity checks.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_utilities_init()

  integer :: hdferr, HDF5_major, HDF5_minor, HDF5_release, configFlags
  logical :: avail
  integer(SIZE_T) :: typeSize


  print'(/,1x,a)', '<<<+-  HDF5_utilities init  -+>>>'


  call H5Open_f(hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tget_size_f(H5T_NATIVE_INTEGER,typeSize, hdferr)
  call HDF5_chkerr(hdferr)
  if (int(bit_size(0),SIZE_T)/=typeSize*8) &
    error stop 'Default integer size does not match H5T_NATIVE_INTEGER'

  call H5Tget_size_f(H5T_NATIVE_DOUBLE,typeSize, hdferr)
  call HDF5_chkerr(hdferr)
  if (int(storage_size(0.0_pREAL),SIZE_T)/=typeSize*8) &
    error stop 'pREAL does not match H5T_NATIVE_DOUBLE'

  call H5get_libversion_f(HDF5_major,HDF5_minor,HDF5_release,hdferr)
  call HDF5_chkerr(hdferr)
  compression_possible = (HDF5_major == 1 .and. HDF5_minor >= 12)                                   ! https://forum.hdfgroup.org/t/6186

  call H5Zfilter_avail_f(H5Z_FILTER_DEFLATE_F,avail,hdferr)
  call HDF5_chkerr(hdferr)
  compression_possible = compression_possible .and. avail

  if (avail) then
    call H5Zget_filter_info_f(H5Z_FILTER_DEFLATE_F,configFlags,hdferr)
    call HDF5_chkerr(hdferr)
    compression_possible = compression_possible .and. iand(H5Z_FILTER_ENCODE_ENABLED_F,configFlags) > 0
  end if

  call H5Zfilter_avail_f(H5Z_FILTER_SHUFFLE_F,avail,hdferr)
  call HDF5_chkerr(hdferr)
  compression_possible = compression_possible .and. avail

  if (avail) then
    call H5Zget_filter_info_f(H5Z_FILTER_SHUFFLE_F,configFlags,hdferr)
    call HDF5_chkerr(hdferr)
    compression_possible = compression_possible .and. iand(H5Z_FILTER_ENCODE_ENABLED_F,configFlags) > 0
  end if

end subroutine HDF5_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief Check for HDF5 error.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_chkerr(e)

  integer, intent(in) :: e


  if (e<0) error stop 'HDF5 error'

end subroutine HDF5_chkerr


!--------------------------------------------------------------------------------------------------
!> @brief Open and initialize HDF5 file.
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openFile(fileName,mode,parallel)

  character(len=*), intent(in)           :: fileName
  character,        intent(in), optional :: mode
  logical,          intent(in), optional :: parallel

  character      :: m
  integer(HID_T) :: plist_id
  integer :: hdferr
  logical :: exist


  m = misc_optional(mode,'r')

  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
  call HDF5_chkerr(hdferr)

#ifdef PETSC
  if (misc_optional(parallel,.true.)) &
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
    call H5Pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL_F90, hdferr)
#else
    call H5Pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL,     hdferr)
#endif
  call HDF5_chkerr(hdferr)
#endif

  if     (m == 'w') then
    call H5Fcreate_f(fileName,H5F_ACC_TRUNC_F,HDF5_openFile,hdferr,access_prp = plist_id)
    call HDF5_chkerr(hdferr)
  elseif (m == 'a') then
    call H5Fopen_f(fileName,H5F_ACC_RDWR_F,HDF5_openFile,hdferr,access_prp = plist_id)
    call HDF5_chkerr(hdferr)
  elseif (m == 'r') then
    inquire(file=fileName,exist=exist)
    if (.not. exist) call IO_error(100,trim(fileName))
    call H5Fopen_f(fileName,H5F_ACC_RDONLY_F,HDF5_openFile,hdferr,access_prp = plist_id)
    call HDF5_chkerr(hdferr)
  else
    error stop 'unknown access mode'
  end if

  call H5Pclose_f(plist_id, hdferr)
  call HDF5_chkerr(hdferr)

end function HDF5_openFile


!--------------------------------------------------------------------------------------------------
!> @brief Close an open HDF5 file.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeFile(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer     :: hdferr

  call H5Fclose_f(fileHandle,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_closeFile


!--------------------------------------------------------------------------------------------------
!> @brief Add a new group to the file.
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup(fileHandle,groupName)

  integer(HID_T),   intent(in) :: fileHandle
  character(len=*), intent(in) :: groupName

  integer        :: hdferr
  integer(HID_T) :: aplist_id

!-------------------------------------------------------------------------------------------------
! creating a property list for data access properties
  call H5Pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  call HDF5_chkerr(hdferr)

!-------------------------------------------------------------------------------------------------
! setting I/O mode to collective
#ifdef PETSC
  call H5Pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
  call HDF5_chkerr(hdferr)
#endif

!-------------------------------------------------------------------------------------------------
! Create group
  call H5Gcreate_f(fileHandle, trim(groupName), HDF5_addGroup, hdferr, OBJECT_NAMELEN_DEFAULT_F,gapl_id = aplist_id)
  call HDF5_chkerr(hdferr)

  call H5Pclose_f(aplist_id,hdferr)

end function HDF5_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief Open an existing group in a file.
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup(fileHandle,groupName)

  integer(HID_T),   intent(in) :: fileHandle
  character(len=*), intent(in) :: groupName


  integer        :: hdferr
  integer(HID_T) :: aplist_id
  logical        :: is_collective


 !-------------------------------------------------------------------------------------------------
 ! creating a property list for data access properties
  call H5Pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  call HDF5_chkerr(hdferr)

 !-------------------------------------------------------------------------------------------------
 ! setting I/O mode to collective
#ifdef PETSC
  call H5Pget_all_coll_metadata_ops_f(aplist_id, is_collective, hdferr)
  call HDF5_chkerr(hdferr)
#endif

 !-------------------------------------------------------------------------------------------------
 ! opening the group
  call H5Gopen_f(fileHandle, trim(groupName), HDF5_openGroup, hdferr, gapl_id = aplist_id)
  call HDF5_chkerr(hdferr)

  call H5Pclose_f(aplist_id,hdferr)

end function HDF5_openGroup


!--------------------------------------------------------------------------------------------------
!> @brief Close a group.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeGroup(group_id)

  integer(HID_T), intent(in) :: group_id

  integer :: hdferr

  call H5Gclose_f(group_id, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_closeGroup


!--------------------------------------------------------------------------------------------------
!> @brief Check whether a group or a dataset exists.
!--------------------------------------------------------------------------------------------------
logical function HDF5_objectExists(loc_id,path)

  integer(HID_T),   intent(in)            :: loc_id
  character(len=*), intent(in), optional  :: path

  integer :: hdferr
  character(len=:), allocatable :: p


  p = trim(misc_optional(path,'.'))

  call H5Lexists_f(loc_id, p, HDF5_objectExists, hdferr)
  call HDF5_chkerr(hdferr)

  if (HDF5_objectExists) then
    call H5Oexists_by_name_f(loc_id, p, HDF5_objectExists, hdferr)
    call HDF5_chkerr(hdferr)
  end if

end function HDF5_objectExists


!--------------------------------------------------------------------------------------------------
!> @brief Add a string attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_str(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path

  integer(HID_T) :: attr_id, space_id
  logical        :: attrExists
  integer        :: hdferr
  character(len=:), allocatable :: p
  character(len=len_trim(attrValue)+1,kind=C_CHAR), dimension(1), target :: attrValue_
  type(C_PTR), target, dimension(1) :: ptr


  p = trim(misc_optional(path,'.'))

  attrValue_(1) = trim(attrValue)//C_NULL_CHAR
  ptr(1) = c_loc(attrValue_(1))

  call H5Screate_f(H5S_SCALAR_F,space_id,hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_STRING,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_STRING, c_loc(ptr), hdferr)                                          ! ptr instead of c_loc(ptr) works on gfortran, not on ifort
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_int(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel
  integer,          intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path

  integer(HID_T) :: attr_id, space_id
  integer        :: hdferr
  logical        :: attrExists
  character(len=:), allocatable :: p


  p = trim(misc_optional(path,'.'))

  call H5Screate_f(H5S_SCALAR_F,space_id,hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, int([1],HSIZE_T), hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief Add a real attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel
  real(pREAL),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path

  integer(HID_T) :: attr_id, space_id
  integer        :: hdferr
  logical        :: attrExists
  character(len=:), allocatable :: p


  p = trim(misc_optional(path,'.'))

  call H5Screate_f(H5S_SCALAR_F,space_id,hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, int([1],HSIZE_T), hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief Add a string array attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_str_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  character(len=*), intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path

  integer(HID_T)                :: attr_id, space_id
  logical                       :: attrExists
  integer                       :: hdferr,i
  character(len=:), allocatable :: p
  character(len=len(attrValue)+1,kind=C_CHAR), dimension(size(attrValue)), target :: attrValue_
  type(C_PTR), target, dimension(size(attrValue))  :: ptr


  p = trim(misc_optional(path,'.'))

  do i=1,size(attrValue)
    attrValue_(i) = attrValue(i)//C_NULL_CHAR
    ptr(i) = c_loc(attrValue_(i))
  end do

  call H5Screate_simple_f(1,shape(attrValue_,kind=HSIZE_T),space_id,hdferr,shape(attrValue_,kind=HSIZE_T))
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_STRING,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_STRING, c_loc(ptr), hdferr)                                          ! ptr instead of c_loc(ptr) works on gfortran, not on ifort
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_str_array


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer array attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_int_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  integer,          intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path

  integer(HSIZE_T),dimension(1) :: array_size
  integer(HID_T)                :: attr_id, space_id
  integer                       :: hdferr
  logical                       :: attrExists
  character(len=:), allocatable :: p


  p = trim(misc_optional(path,'.'))

  array_size = size(attrValue,kind=HSIZE_T)

  call H5Screate_simple_f(1, array_size, space_id, hdferr, array_size)
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, array_size, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief Add a real array attribute to the path given relative to the location.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  real(pREAL),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path

  integer(HSIZE_T),dimension(1) :: array_size
  integer(HID_T)                :: attr_id, space_id
  integer                       :: hdferr
  logical                       :: attrExists
  character(len=:), allocatable :: p


  p = trim(misc_optional(path,'.'))

  array_size = size(attrValue,kind=HSIZE_T)

  call H5Screate_simple_f(1, array_size, space_id, hdferr, array_size)
  call HDF5_chkerr(hdferr)

  call H5Aexists_by_name_f(loc_id,p,attrLabel,attrExists,hdferr)
  call HDF5_chkerr(hdferr)
  if (attrExists) then
    call H5Adelete_by_name_f(loc_id,p,attrLabel,hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Acreate_by_name_f(loc_id,p,trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, array_size, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Aclose_f(attr_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id,hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief Set link to object in results file.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_setLink(loc_id,target_name,link_name)

  character(len=*), intent(in) :: target_name, link_name
  integer(HID_T),   intent(in) :: loc_id
  integer                      :: hdferr
  logical                      :: linkExists

  call H5Lexists_f(loc_id, link_name,linkExists, hdferr)
  call HDF5_chkerr(hdferr)
  if (linkExists) then
    call H5Ldelete_f(loc_id,link_name, hdferr)
    call HDF5_chkerr(hdferr)
  end if
  call H5Lcreate_soft_f(target_name, loc_id, link_name, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_setLink


!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real1(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:) :: dataset                                            !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id, aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real1

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real2(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:) :: dataset                                          !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real2

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real3(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:,:) :: dataset                                        !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real3

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real4(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:,:,:) :: dataset                                      !< read data
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real4

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real5(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:,:,:,:) :: dataset                                    !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real5

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real6(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:,:,:,:,:) :: dataset                                  !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real6

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real7(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(out), dimension(:,:,:,:,:,:,:) :: dataset                                !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real7


!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int1(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:) :: dataset                                            !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int1

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int2(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:) :: dataset                                          !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int2

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int3(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:,:) :: dataset                                        !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                         myStart,totalShape,loc_id,myShape,datasetName, &
                         misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int3

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer withh 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int4(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:,:,:) :: dataset                                      !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                         myStart,totalShape,loc_id,myShape,datasetName, &
                         misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int4

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int5(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:,:,:,:) :: dataset                                    !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                         myStart,totalShape,loc_id,myShape,datasetName, &
                         misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int5

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int6(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:,:,:,:,:) :: dataset                                  !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int6

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type integer with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int7(dataset,loc_id,datasetName,parallel)

  integer,          intent(out), dimension(:,:,:,:,:,:,:) :: dataset                                !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr


  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_read(dset_id,filespace_id,memspace_id,plist_id,aplist_id, &
                       myStart,totalShape,loc_id,myShape,datasetName, &
                       misc_optional(parallel,parallel_default))

  if (any(totalShape == 0)) return

  call H5Dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  call HDF5_chkerr(hdferr)

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int7

#if defined(__GFORTRAN__)

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real1(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:) :: dataset                                             !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real1

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real2(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:) :: dataset                                           !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real2

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real3(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:,:) :: dataset                                         !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real3

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real4(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:,:,:) :: dataset                                       !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real4


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real5(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:,:,:,:) :: dataset                                     !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real5

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real6(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:,:,:,:,:) :: dataset                                   !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real6

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real7(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(:,:,:,:,:,:,:) :: dataset                                 !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real7

#else

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 1-7 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real(dataset,loc_id,datasetName,parallel)

  real(pREAL),      intent(in), dimension(..) :: dataset                                            !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    select rank(dataset)
      rank (1)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (2)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (3)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (4)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (5)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (6)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank (7)
        call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                        file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    end select
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real
#endif


!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type string (scalar).
!> @details Not collective, must be called by one process at at time.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_str(dataset,loc_id,datasetName)

  character(len=*), intent(in) :: dataset
  integer(HID_T),   intent(in) :: loc_id
  character(len=*), intent(in) :: datasetName                                                      !< name of the dataset in the file

  integer(HID_T)  :: filetype_id, memtype_id, space_id, dataset_id, dcpl
  integer :: hdferr
  character(len=len_trim(dataset,pI64),kind=C_CHAR), target :: dataset_


  dataset_ = trim(dataset)

  call H5Tcopy_f(H5T_C_S1, filetype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tset_size_f(filetype_id, len(dataset_,SIZE_T)+1_SIZE_T, hdferr)                            ! +1 for NULL
  call HDF5_chkerr(hdferr)

  call H5Tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tset_size_f(memtype_id, len(dataset_,SIZE_T), hdferr)
  call HDF5_chkerr(hdferr)

  call H5Pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Pset_chunk_f(dcpl, 1, [1_HSIZE_T], hdferr)
  call HDF5_chkerr(hdferr)
  call H5Pset_Fletcher32_f(dcpl,hdferr)
  call HDF5_chkerr(hdferr)
  if (compression_possible .and. len(dataset,pI64) > 1024_pI64*256_pI64) then
    call H5Pset_shuffle_f(dcpl, hdferr)
    call HDF5_chkerr(hdferr)
    call H5Pset_deflate_f(dcpl, 6, hdferr)
    call HDF5_chkerr(hdferr)
  end if

  call H5Screate_simple_f(1, [1_HSIZE_T], space_id, hdferr)
  call HDF5_chkerr(hdferr)
  CALL H5Dcreate_f(loc_id, datasetName, filetype_id, space_id, dataset_id, hdferr, dcpl)
  call HDF5_chkerr(hdferr)

  call H5Dwrite_f(dataset_id, memtype_id, c_loc(dataset_(1:1)), hdferr)
  call HDF5_chkerr(hdferr)

  call H5Pclose_f(dcpl, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Dclose_f(dataset_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(space_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(memtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(filetype_id, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine HDF5_write_str

#if defined(__GFORTRAN__)

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 1 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int1(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:) :: dataset                                             !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int1

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 2 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int2(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:) :: dataset                                           !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int2

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 3 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int3(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:) :: dataset                                         !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int3

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 4 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int4(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:) :: dataset                                       !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int4

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 5 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int5(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:) :: dataset                                     !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
   call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int5

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 6 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int6(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:,:) :: dataset                                   !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int6

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 7 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int7(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:,:,:) :: dataset                                 !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int7

#else

!--------------------------------------------------------------------------------------------------
!> @brief Write dataset of type integer with 1-7 dimensions.
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(..) :: dataset                                            !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(rank(dataset)) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  call initialize_write(dset_id,filespace_id,memspace_id,plist_id, &
                        myStart,totalShape,loc_id,myShape,datasetName,H5T_NATIVE_INTEGER, &
                        misc_optional(parallel,parallel_default))

  if (product(totalShape) /= 0) then
    select rank(dataset)
      rank(1)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(2)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(3)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(4)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(5)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(6)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
      rank(7)
        call H5Dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                       file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    end select
    call HDF5_chkerr(hdferr)
  end if

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int
#endif

!--------------------------------------------------------------------------------------------------
!> @brief Initialize read handles and determine global shape in case of parallel IO.
!--------------------------------------------------------------------------------------------------
subroutine initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                           myStart, globalShape, &
                           loc_id,localShape,datasetName,parallel)

  integer(HID_T),    intent(in) :: loc_id                                                           !< file or group handle
  character(len=*),  intent(in) :: datasetName                                                      !< name of the dataset in the file
  logical,           intent(in) :: parallel
  integer(HSIZE_T),  intent(in),  dimension(:) :: localShape
  integer(HSIZE_T),  intent(out), dimension(size(localShape)) :: &
    myStart, &
    globalShape                                                                                     !< shape of the dataset (all processes)
  integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id, aplist_id

  integer(MPI_INTEGER_KIND), dimension(worldsize) :: readSize                                       !< contribution of all processes
  integer :: hdferr
  integer(MPI_INTEGER_KIND) :: err_MPI

!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective for MPI)
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
  readSize = 0_MPI_INTEGER_KIND
  readSize(worldrank+1) = int(localShape(ubound(localShape,1)),MPI_INTEGER_KIND)
#ifdef PETSC
  if (parallel) then
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call HDF5_chkerr(hdferr)
    call MPI_Allgather(int(localShape(ubound(localShape,1)),MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                       readSize,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  end if
#endif
  myStart                   = int(0,HSIZE_T)
  myStart(ubound(myStart))  = int(sum(readSize(1:worldrank)),HSIZE_T)
  globalShape = [localShape(1:ubound(localShape,1)-1),int(sum(readSize),HSIZE_T)]

  if (any(globalShape == 0)) then
    call H5Pclose_f(plist_id, hdferr)
    call HDF5_chkerr(hdferr)
    return
  end if

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape)
  call H5Screate_simple_f(size(localShape), localShape, memspace_id, hdferr, localShape)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! creating a property list for IO and set it to collective
  call H5Pcreate_f(H5P_DATASET_ACCESS_F, aplist_id, hdferr)
  call HDF5_chkerr(hdferr)
#ifdef PETSC
  call H5Pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
  call HDF5_chkerr(hdferr)
#endif

!--------------------------------------------------------------------------------------------------
! open the dataset in the file and get the space ID
  call H5Dopen_f(loc_id,datasetName,dset_id,hdferr, dapl_id = aplist_id)
  call HDF5_chkerr(hdferr)
  call H5Dget_space_f(dset_id, filespace_id, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! select a hyperslab (the portion of the current process) in the file
  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, localShape, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine initialize_read


!--------------------------------------------------------------------------------------------------
!> @brief Close read handles.
!--------------------------------------------------------------------------------------------------
subroutine finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

  integer(HID_T), intent(in) :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer :: hdferr

  call H5Pclose_f(plist_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Pclose_f(aplist_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Dclose_f(dset_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(filespace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(memspace_id, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine finalize_read


!--------------------------------------------------------------------------------------------------
!> @brief Initialize write handles and determine global shape in case of parallel IO.
!--------------------------------------------------------------------------------------------------
subroutine initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                            myStart, totalShape, &
                            loc_id,localShape,datasetName,datatype,parallel)

  integer(HID_T),    intent(in) :: loc_id                                                           !< file or group handle
  character(len=*),  intent(in) :: datasetName                                                      !< name of the dataset in the file
  logical,           intent(in) :: parallel
  integer(HID_T),    intent(in) :: datatype
  integer(HSIZE_T),  intent(in),  dimension(:) :: localShape
  integer(HSIZE_T),  intent(out), dimension(size(localShape)) :: &
    myStart, &
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id

  integer(MPI_INTEGER_KIND), dimension(worldsize) :: writeSize                                      !< contribution of all processes
  integer(HID_T) :: dcpl
  integer :: hdferr
  integer(MPI_INTEGER_KIND) :: err_MPI
  integer(HSIZE_T), parameter :: chunkSize = 1024_HSIZE_T**2/8_HSIZE_T


!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective when writing in parallel)
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  call HDF5_chkerr(hdferr)
#ifdef PETSC
  if (parallel) then
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
    call HDF5_chkerr(hdferr)
  end if
#endif

!--------------------------------------------------------------------------------------------------
! determine the global data layout among all processes
  writeSize              = 0_MPI_INTEGER_KIND
  writeSize(worldrank+1) = int(localShape(ubound(localShape,1)),MPI_INTEGER_KIND)
#ifdef PETSC
  if (parallel) then
    call MPI_Allgather(int(localShape(ubound(localShape,1)),MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                       writeSize,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  end if
#endif
  myStart                   = int(0,HSIZE_T)
  myStart(ubound(myStart))  = int(sum(writeSize(1:worldrank)),HSIZE_T)
  totalShape = [localShape(1:ubound(localShape,1)-1),int(sum(writeSize),HSIZE_T)]

!--------------------------------------------------------------------------------------------------
! chunk dataset, enable compression for larger datasets
  call H5Pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
  call HDF5_chkerr(hdferr)

  if (product(totalShape) > 0) then
    call H5Pset_Fletcher32_f(dcpl,hdferr)
    call HDF5_chkerr(hdferr)

    if (product(totalShape) >= chunkSize*2_HSIZE_T) then
      call H5Pset_chunk_f(dcpl, size(totalShape), getChunks(totalShape,chunkSize), hdferr)
      call HDF5_chkerr(hdferr)
      if (compression_possible) then
        call H5Pset_shuffle_f(dcpl, hdferr)
        call HDF5_chkerr(hdferr)
        call H5Pset_deflate_f(dcpl, 6, hdferr)
        call HDF5_chkerr(hdferr)
      end if
    else
      call H5Pset_chunk_f(dcpl, size(totalShape), totalShape, hdferr)
      call HDF5_chkerr(hdferr)
    end if
  end if

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape) and in file (global shape)
  call H5Screate_simple_f(size(localShape), localShape, memspace_id, hdferr, localShape)
  call HDF5_chkerr(hdferr)
  call H5Screate_simple_f(size(totalShape), totalShape, filespace_id, hdferr, totalShape)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! create dataset in the file and select a hyperslab from it (the portion of the current process)
  call H5Dcreate_f(loc_id, trim(datasetName), datatype, filespace_id, dset_id, hdferr, dcpl)
  call HDF5_chkerr(hdferr)
  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, localShape, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Pclose_f(dcpl , hdferr)
  call HDF5_chkerr(hdferr)

  contains
  !------------------------------------------------------------------------------------------------
  !> @brief Determine chunk layout.
  !------------------------------------------------------------------------------------------------
  pure function getChunks(totalShape,chunkSize)

    integer(HSIZE_T), dimension(:), intent(in)    :: totalShape
    integer(HSIZE_T),               intent(in)    :: chunkSize
    integer(HSIZE_T), dimension(size(totalShape)) :: getChunks

    getChunks = [totalShape(1:size(totalShape)-1),&
                 chunkSize/product(totalShape(1:size(totalShape)-1))]

  end function getChunks

end subroutine initialize_write


!--------------------------------------------------------------------------------------------------
!> @brief Closes write handles.
!--------------------------------------------------------------------------------------------------
subroutine finalize_write(plist_id, dset_id, filespace_id, memspace_id)

  integer(HID_T), intent(in) :: dset_id, filespace_id, memspace_id, plist_id
  integer :: hdferr

  call H5Pclose_f(plist_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Dclose_f(dset_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(filespace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(memspace_id, hdferr)
  call HDF5_chkerr(hdferr)

end subroutine finalize_write


end module HDF5_Utilities
