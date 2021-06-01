!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module HDF5_utilities
  use HDF5
#ifdef PETSc
  use PETSC
#endif

  use prec
  use parallelization

  implicit none
  public

!--------------------------------------------------------------------------------------------------
!> @brief reads integer or float data of defined shape from file
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
!> @brief writes integer or real data of defined shape to file
!> @details for parallel IO, all dimension except for the last need to match
!--------------------------------------------------------------------------------------------------
  interface HDF5_write
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
  end interface HDF5_write

!--------------------------------------------------------------------------------------------------
!> @brief attached attributes of type char, integer or real to a file/dataset/group
!--------------------------------------------------------------------------------------------------
  interface HDF5_addAttribute
    module procedure HDF5_addAttribute_str
    module procedure HDF5_addAttribute_int
    module procedure HDF5_addAttribute_real
    module procedure HDF5_addAttribute_int_array
    module procedure HDF5_addAttribute_real_array
  end interface HDF5_addAttribute

#ifdef PETSc
  logical, parameter, private :: parallel_default = .true.
#else
  logical, parameter, private :: parallel_default = .false.
#endif

contains


!--------------------------------------------------------------------------------------------------
!> @brief initialize HDF5 libary and do sanity checks
!--------------------------------------------------------------------------------------------------
subroutine HDF5_utilities_init

  integer         :: hdferr
  integer(SIZE_T) :: typeSize

  print'(/,a)', ' <<<+-  HDF5_Utilities init  -+>>>'

!--------------------------------------------------------------------------------------------------
!initialize HDF5 library and check if integer and float type size match
  call h5open_f(hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tget_size_f(H5T_NATIVE_INTEGER,typeSize, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (int(bit_size(0),SIZE_T)/=typeSize*8) &
    error stop 'Default integer size does not match H5T_NATIVE_INTEGER'

  call h5tget_size_f(H5T_NATIVE_DOUBLE,typeSize, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (int(storage_size(0.0_pReal),SIZE_T)/=typeSize*8) &
    error stop 'pReal does not match H5T_NATIVE_DOUBLE'

end subroutine HDF5_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief open and initializes HDF5 output file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openFile(fileName,mode)

  character(len=*), intent(in)           :: fileName
  character,        intent(in), optional :: mode

  character                              :: m
  integer(HID_T)                         :: plist_id
  integer                 :: hdferr


  if (present(mode)) then
    m = mode
  else
    m = 'r'
  endif

  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

#ifdef PETSc
  call h5pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#endif

  if    (m == 'w') then
    call h5fcreate_f(fileName,H5F_ACC_TRUNC_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  elseif(m == 'a') then
    call h5fopen_f(fileName,H5F_ACC_RDWR_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  elseif(m == 'r') then
    call h5fopen_f(fileName,H5F_ACC_RDONLY_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  else
    error stop 'unknown access mode'
  endif

  call h5pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end function HDF5_openFile


!--------------------------------------------------------------------------------------------------
!> @brief close the opened HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeFile(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer     :: hdferr

  call h5fclose_f(fileHandle,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_closeFile


!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the fileHandle
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup(fileHandle,groupName)

  integer(HID_T),   intent(in) :: fileHandle
  character(len=*), intent(in) :: groupName

  integer        :: hdferr
  integer(HID_T) :: aplist_id

!-------------------------------------------------------------------------------------------------
! creating a property list for data access properties
  call h5pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!-------------------------------------------------------------------------------------------------
! setting I/O mode to collective
#ifdef PETSc
  call h5pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#endif

!-------------------------------------------------------------------------------------------------
! Create group
  call h5gcreate_f(fileHandle, trim(groupName), HDF5_addGroup, hdferr, OBJECT_NAMELEN_DEFAULT_F,gapl_id = aplist_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5pclose_f(aplist_id,hdferr)

end function HDF5_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief open an existing group of a file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup(fileHandle,groupName)

  integer(HID_T),   intent(in) :: fileHandle
  character(len=*), intent(in) :: groupName


  integer        :: hdferr
  integer(HID_T) :: aplist_id
  logical        :: is_collective


 !-------------------------------------------------------------------------------------------------
 ! creating a property list for data access properties
  call h5pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

 !-------------------------------------------------------------------------------------------------
 ! setting I/O mode to collective
#ifdef PETSc
  call h5pget_all_coll_metadata_ops_f(aplist_id, is_collective, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#endif

 !-------------------------------------------------------------------------------------------------
 ! opening the group
  call h5gopen_f(fileHandle, trim(groupName), HDF5_openGroup, hdferr, gapl_id = aplist_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5pclose_f(aplist_id,hdferr)

end function HDF5_openGroup


!--------------------------------------------------------------------------------------------------
!> @brief close a group
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeGroup(group_id)

  integer(HID_T), intent(in) :: group_id

  integer :: hdferr

  call h5gclose_f(group_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_closeGroup


!--------------------------------------------------------------------------------------------------
!> @brief check whether a group or a dataset exists
!--------------------------------------------------------------------------------------------------
logical function HDF5_objectExists(loc_id,path)

  integer(HID_T),   intent(in)            :: loc_id
  character(len=*), intent(in), optional  :: path

  integer :: hdferr
  character(len=:), allocatable :: p

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  call h5lexists_f(loc_id, p, HDF5_objectExists, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  if(HDF5_objectExists) then
    call h5oexists_by_name_f(loc_id, p, HDF5_objectExists, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

end function HDF5_objectExists


!--------------------------------------------------------------------------------------------------
!> @brief adds a string attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_str(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path

  integer(HID_T) :: attr_id, space_id, type_id
  logical        :: attrExists
  integer        :: hdferr
  character(len=:), allocatable :: p
  character(len=:,kind=C_CHAR), allocatable,target :: attrValue_
  type(c_ptr), target, dimension(1) :: ptr

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  attrValue_ = trim(attrValue)//C_NULL_CHAR
  ptr(1) = c_loc(attrValue_)

  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tcopy_f(H5T_STRING, type_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (attrExists) then
  call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),type_id,space_id,attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5awrite_f(attr_id, type_id, c_loc(ptr(1)), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aclose_f(attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(type_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
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

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, int([1],HSIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aclose_f(attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel
  real(pReal),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path

  integer(HID_T) :: attr_id, space_id
  integer        :: hdferr
  logical        :: attrExists
  character(len=:), allocatable :: p

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, int([1],HSIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aclose_f(attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
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

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  array_size = size(attrValue,kind=HSIZE_T)

  call h5screate_simple_f(1, array_size, space_id, hdferr, array_size)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, array_size, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aclose_f(attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief adds a real attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  real(pReal),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path

  integer(HSIZE_T),dimension(1) :: array_size
  integer(HID_T)                :: attr_id, space_id
  integer                       :: hdferr
  logical                       :: attrExists
  character(len=:), allocatable :: p

  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif

  array_size = size(attrValue,kind=HSIZE_T)

  call h5screate_simple_f(1, array_size, space_id, hdferr, array_size)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, array_size, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5aclose_f(attr_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(space_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief set link to object in results file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_setLink(loc_id,target_name,link_name)

  character(len=*), intent(in) :: target_name, link_name
  integer(HID_T),   intent(in) :: loc_id
  integer                      :: hdferr
  logical                      :: linkExists

  call h5lexists_f(loc_id, link_name,linkExists, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if (linkExists) then
    call h5ldelete_f(loc_id,link_name, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
  call h5lcreate_soft_f(target_name, loc_id, link_name, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine HDF5_setLink


!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real1(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:) :: dataset                                            !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &                                            ! ToDo: Fortran 2018 size(shape(A)) = rank(A)
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real1

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real2(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:) :: dataset                                          !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real2

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real3(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:,:) :: dataset                                        !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real3

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real4(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:,:,:) :: dataset                                      !< read data
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real4

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real5(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:,:,:,:) :: dataset                                    !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real5

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real6(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:,:,:,:,:) :: dataset                                  !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real6

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type real with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real7(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(out), dimension(:,:,:,:,:,:,:) :: dataset                                !< data read from file
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes

  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
   call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                        myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

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
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer :: hdferr

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

!---------------------------------------------------------------------------------------------------
! initialize HDF5 data structures
  if (present(parallel)) then
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel)
  else
    call initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                         myStart, totalShape, loc_id,myShape,datasetName,parallel_default)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if(hdferr < 0) error stop 'HDF5 error'

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int7


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real1(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:) :: dataset                                             !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real1

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real2(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:) :: dataset                                           !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real2

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real3(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:,:) :: dataset                                         !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real3

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real4(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:,:,:) :: dataset                                       !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real4


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real5(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:,:,:,:) :: dataset                                     !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real5

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real6(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:,:,:,:,:) :: dataset                                   !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real6

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type real with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real7(dataset,loc_id,datasetName,parallel)

  real(pReal),      intent(in), dimension(:,:,:,:,:,:,:) :: dataset                                 !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real7


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int1(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:) :: dataset                                             !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int1

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int2(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:) :: dataset                                           !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int2

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int3(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:) :: dataset                                         !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int3

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int4(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:) :: dataset                                       !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int4

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int5(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:) :: dataset                                     !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
   if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int5

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int6(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:,:) :: dataset                                   !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int6

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type integer with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int7(dataset,loc_id,datasetName,parallel)

  integer,          intent(in), dimension(:,:,:,:,:,:,:) :: dataset                                 !< data written to file
  integer(HID_T),   intent(in)  :: loc_id                                                           !< file or group handle
  character(len=*), intent(in)  :: datasetName                                                      !< name of the dataset in the file
  logical, intent(in), optional :: parallel                                                         !< dataset is distributed over multiple processes


  integer :: hdferr
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
  myShape = int(shape(dataset),HSIZE_T)
  if (any(myShape(1:size(myShape)-1) == 0)) return                                                  !< empty dataset (last dimension can be empty)

  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel)
  else
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,parallel_default)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if(hdferr < 0) error stop 'HDF5 error'
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int7


!--------------------------------------------------------------------------------------------------
!> @brief initialize HDF5 handles, determines global shape and start for parallel read
!--------------------------------------------------------------------------------------------------
subroutine initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                           myStart, globalShape, &
                           loc_id,localShape,datasetName,parallel)

  integer(HID_T),    intent(in) :: loc_id                                                           !< file or group handle
  character(len=*),  intent(in) :: datasetName                                                      !< name of the dataset in the file
  logical,           intent(in) :: parallel
  integer(HSIZE_T),  intent(in),   dimension(:) :: &
    localShape
  integer(HSIZE_T),  intent(out), dimension(size(localShape,1)):: &
    myStart, &
    globalShape                                                                                    !< shape of the dataset (all processes)
  integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id, aplist_id

  integer, dimension(worldsize) :: &
    readSize                                                                                        !< contribution of all processes
  integer :: ierr
  integer :: hdferr

!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective for MPI)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
  readSize = 0
  readSize(worldrank+1) = int(localShape(ubound(localShape,1)))
#ifdef PETSc
  if (parallel) then
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
    call MPI_allreduce(MPI_IN_PLACE,readSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)       ! get total output size over each process
    if (ierr /= 0) error stop 'MPI error'
  endif
#endif
  myStart                   = int(0,HSIZE_T)
  myStart(ubound(myStart))  = int(sum(readSize(1:worldrank)),HSIZE_T)
  globalShape = [localShape(1:ubound(localShape,1)-1),int(sum(readSize),HSIZE_T)]

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape)
  call h5screate_simple_f(size(localShape), localShape, memspace_id, hdferr, localShape)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! creating a property list for IO and set it to collective
  call h5pcreate_f(H5P_DATASET_ACCESS_F, aplist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#ifdef PETSc
  call h5pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#endif

!--------------------------------------------------------------------------------------------------
! open the dataset in the file and get the space ID
  call h5dopen_f(loc_id,datasetName,dset_id,hdferr, dapl_id = aplist_id)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dget_space_f(dset_id, filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! select a hyperslab (the portion of the current process) in the file
  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, localShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine initialize_read


!--------------------------------------------------------------------------------------------------
!> @brief closes HDF5 handles
!--------------------------------------------------------------------------------------------------
subroutine finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

  integer(HID_T), intent(in) :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
  integer :: hdferr

  call h5pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5pclose_f(aplist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine finalize_read


!--------------------------------------------------------------------------------------------------
!> @brief initialize HDF5 handles, determines global shape and start for parallel write
!--------------------------------------------------------------------------------------------------
subroutine initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                            myStart, totalShape, &
                            loc_id,myShape,datasetName,datatype,parallel)

  integer(HID_T),    intent(in) :: loc_id                                                           !< file or group handle
  character(len=*),  intent(in) :: datasetName                                                      !< name of the dataset in the file
  logical,           intent(in) :: parallel
  integer(HID_T),    intent(in) :: datatype
  integer(HSIZE_T),  intent(in),   dimension(:) :: &
    myShape
  integer(HSIZE_T),  intent(out), dimension(size(myShape,1)):: &
    myStart, &
    totalShape                                                                                      !< shape of the dataset (all processes)
  integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id

  integer,          dimension(worldsize)      :: writeSize                                          !< contribution of all processes
  integer(HID_T) ::  dcpl
  integer :: ierr, hdferr, HDF5_major, HDF5_minor, HDF5_release
  integer(HSIZE_T), parameter :: chunkSize = 1024_HSIZE_T**2/8_HSIZE_T

!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective when reading in parallel)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
#ifdef PETSc
  if (parallel) then
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
    if(hdferr < 0) error stop 'HDF5 error'
  endif
#endif

!--------------------------------------------------------------------------------------------------
! determine the global data layout among all processes
  writeSize              = 0
  writeSize(worldrank+1) = int(myShape(ubound(myShape,1)))
#ifdef PETSc
  if (parallel) then
    call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)      ! get total output size over each process
    if (ierr /= 0) error stop 'MPI error'
  endif
#endif
  myStart                   = int(0,HSIZE_T)
  myStart(ubound(myStart))  = int(sum(writeSize(1:worldrank)),HSIZE_T)
  totalShape = [myShape(1:ubound(myShape,1)-1),int(sum(writeSize),HSIZE_T)]

!--------------------------------------------------------------------------------------------------
! compress (and chunk) larger datasets
  call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  if(product(totalShape) >= chunkSize*2_HSIZE_T) then
    call H5get_libversion_f(HDF5_major,HDF5_minor,HDF5_release,hdferr)
    if (hdferr < 0) error stop 'HDF5 error'
    if (HDF5_major == 1 .and. HDF5_minor >= 12) then                                                ! https://forum.hdfgroup.org/t/6186
      call h5pset_chunk_f(dcpl, size(totalShape), getChunks(totalShape,chunkSize), hdferr)
      if (hdferr < 0) error stop 'HDF5 error'
      call h5pset_shuffle_f(dcpl, hdferr)
      if (hdferr < 0) error stop 'HDF5 error'
      call h5pset_deflate_f(dcpl, 6, hdferr)
      if (hdferr < 0) error stop 'HDF5 error'
      call h5pset_Fletcher32_f(dcpl,hdferr)
      if (hdferr < 0) error stop 'HDF5 error'
    endif
  endif

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape) and in file (global shape)
  call h5screate_simple_f(size(myShape), myShape, memspace_id, hdferr, myShape)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5screate_simple_f(size(totalShape), totalShape, filespace_id, hdferr, totalShape)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create dataset in the file and select a hyperslab from it (the portion of the current process)
  call h5dcreate_f(loc_id, trim(datasetName), datatype, filespace_id, dset_id, hdferr, dcpl)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, myShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5pclose_f(dcpl , hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  contains
  !------------------------------------------------------------------------------------------------
  !> @brief determine chunk layout
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
!> @brief closes HDF5 handles
!--------------------------------------------------------------------------------------------------
subroutine finalize_write(plist_id, dset_id, filespace_id, memspace_id)

  integer(HID_T), intent(in) :: dset_id, filespace_id, memspace_id, plist_id
  integer :: hdferr

  call h5pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

end subroutine finalize_write


end module HDF5_Utilities
