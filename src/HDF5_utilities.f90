!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module HDF5_utilities
  use prec
  use IO
  use HDF5
  use rotations
  use numerics
#ifdef PETSc
  use PETSC
#endif

 implicit none
 public

!--------------------------------------------------------------------------------------------------
!> @brief reads pInt or pReal data of defined shape from file                                       ! ToDo: order of arguments wrong
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
!> @brief writes pInt or pReal data of defined shape to file                                        ! ToDo: order of arguments wrong
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
    
    module procedure HDF5_write_rotation
 
  end interface HDF5_write
 
!--------------------------------------------------------------------------------------------------
!> @brief attached attributes of type char,pInt or pReal to a file/dataset/group
!--------------------------------------------------------------------------------------------------
  interface HDF5_addAttribute
    module procedure HDF5_addAttribute_str
    module procedure HDF5_addAttribute_int
    module procedure HDF5_addAttribute_real
    module procedure HDF5_addAttribute_int_array
    module procedure HDF5_addAttribute_real_array
  end interface HDF5_addAttribute
 
 
!--------------------------------------------------------------------------------------------------
  public :: &
    HDF5_utilities_init, &
    HDF5_openFile, &
    HDF5_closeFile, &
    HDF5_addAttribute, &
    HDF5_closeGroup ,&
    HDF5_openGroup, &
    HDF5_addGroup, &
    HDF5_read, &
    HDF5_write, &
    HDF5_setLink, &
    HDF5_objectExists
contains


!--------------------------------------------------------------------------------------------------
!> @brief open libary and do sanity checks
!--------------------------------------------------------------------------------------------------
subroutine HDF5_utilities_init

  integer :: hdferr
  integer(SIZE_T)        :: typeSize

  write(6,'(/,a)') ' <<<+-  HDF5_Utilities init  -+>>>'

!--------------------------------------------------------------------------------------------------
!initialize HDF5 library and check if integer and float type size match
  call h5open_f(hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg='HDF5_Utilities_init: h5open_f')

  call h5tget_size_f(H5T_NATIVE_INTEGER,typeSize, hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg='HDF5_Utilities_init: h5tget_size_f (int)')
  if (int(bit_size(0),SIZE_T)/=typeSize*8) &
    call IO_error(0_pInt,ext_msg='Default integer size does not match H5T_NATIVE_INTEGER')

  call h5tget_size_f(H5T_NATIVE_DOUBLE,typeSize, hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg='HDF5_Utilities_init: h5tget_size_f (double)')
  if (int(storage_size(0.0_pReal),SIZE_T)/=typeSize*8) &
    call IO_error(0,ext_msg='pReal does not match H5T_NATIVE_DOUBLE')

end subroutine HDF5_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief open and initializes HDF5 output file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openFile(fileName,mode,parallel) ! ToDo: simply "open" is enough

  character(len=*), intent(in)           :: fileName
  character,        intent(in), optional :: mode
  logical,          intent(in), optional :: parallel
 
  character                              :: m
  integer(HID_T)                         :: plist_id
  integer                 :: hdferr
 
  if (present(mode)) then
    m = mode
  else
    m = 'r'
  endif
 
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5pcreate_f')

#ifdef PETSc
  if (present(parallel)) then; if (parallel) then
    call h5pset_fapl_mpio_f(plist_id, PETSC_COMM_WORLD, MPI_INFO_NULL, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5pset_fapl_mpio_f')
  endif; endif
#endif

  if    (m == 'w') then
    call h5fcreate_f(fileName,H5F_ACC_TRUNC_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fcreate_f (w)')
  elseif(m == 'a') then
    call h5fopen_f(fileName,H5F_ACC_RDWR_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (a)')
  elseif(m == 'r') then
    call h5fopen_f(fileName,H5F_ACC_RDONLY_F,HDF5_openFile,hdferr,access_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f (r)')
  else
    call IO_error(1_pInt,ext_msg='HDF5_openFile: h5fopen_f unknown access mode: '//trim(m))
  endif
 
  call h5pclose_f(plist_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_openFile: h5pclose_f')

end function HDF5_openFile


!--------------------------------------------------------------------------------------------------
!> @brief close the opened HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeFile(fileHandle)

  integer(HID_T), intent(in) :: fileHandle
 
  integer     :: hdferr
 
  call h5fclose_f(fileHandle,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_closeFile: h5fclose_f')

end subroutine HDF5_closeFile


!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the fileHandle
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_addGroup(fileHandle,groupName)

  integer(HID_T), intent(in)   :: fileHandle
  character(len=*), intent(in) :: groupName
 
  integer       :: hdferr
  integer(HID_T)               :: aplist_id

!-------------------------------------------------------------------------------------------------
! creating a property list for data access properties
  call h5pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup: h5pcreate_f ('//trim(groupName)//')')

!-------------------------------------------------------------------------------------------------
! setting I/O mode to collective
#ifdef PETSc
  call h5pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup: h5pset_all_coll_metadata_ops_f ('//trim(groupName)//')')
#endif
 
!-------------------------------------------------------------------------------------------------
! Create group
  call h5gcreate_f(fileHandle, trim(groupName), HDF5_addGroup, hdferr, OBJECT_NAMELEN_DEFAULT_F,gapl_id = aplist_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_addGroup: h5gcreate_f ('//trim(groupName)//')')

  call h5pclose_f(aplist_id,hdferr)

end function HDF5_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief open an existing group of a file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function HDF5_openGroup(fileHandle,groupName)

  integer(HID_T),   intent(in) :: fileHandle
  character(len=*), intent(in) :: groupName
 
 
  integer       :: hdferr
  integer(HID_T)   :: aplist_id
  logical          :: is_collective
 
 
 !-------------------------------------------------------------------------------------------------
 ! creating a property list for data access properties
  call h5pcreate_f(H5P_GROUP_ACCESS_F, aplist_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup: h5pcreate_f ('//trim(groupName)//')')

 !-------------------------------------------------------------------------------------------------
 ! setting I/O mode to collective
#ifdef PETSc
  call h5pget_all_coll_metadata_ops_f(aplist_id, is_collective, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup: h5pset_all_coll_metadata_ops_f ('//trim(groupName)//')')
#endif
 
 !-------------------------------------------------------------------------------------------------
 ! opening the group
  call h5gopen_f(fileHandle, trim(groupName), HDF5_openGroup, hdferr, gapl_id = aplist_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_openGroup: h5gopen_f ('//trim(groupName)//')')

  call h5pclose_f(aplist_id,hdferr)

end function HDF5_openGroup


!--------------------------------------------------------------------------------------------------
!> @brief close a group
!--------------------------------------------------------------------------------------------------
subroutine HDF5_closeGroup(group_id)

  integer(HID_T), intent(in) :: group_id
  integer     :: hdferr
 
  call h5gclose_f(group_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_closeGroup: h5gclose_f (el is ID)', el = int(group_id,pInt))

end subroutine HDF5_closeGroup


!--------------------------------------------------------------------------------------------------
!> @brief check whether a group or a dataset exists
!--------------------------------------------------------------------------------------------------
logical function HDF5_objectExists(loc_id,path)

  integer(HID_T),   intent(in)            :: loc_id
  character(len=*), intent(in), optional  :: path
  
  integer            :: hdferr
  character(len=256) :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
 
  call h5lexists_f(loc_id, p, HDF5_objectExists, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_objectExists: h5oexists_by_name_f')
  
  if(HDF5_objectExists) then
    call h5oexists_by_name_f(loc_id, p, HDF5_objectExists, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_objectExists: h5oexists_by_name_f')
  endif

end function HDF5_objectExists


!--------------------------------------------------------------------------------------------------
!> @brief adds a string attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_str(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path
  
  integer            :: hdferr
  integer(HID_T)     :: attr_id, space_id, type_id
  logical            :: attrExists
  character(len=256) :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
 
  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5screate_f')
  call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5tcopy_f')
  call h5tset_size_f(type_id, int(len(trim(attrValue)),HSIZE_T), hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5tset_size_f')
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5aexists_by_name_f')
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5adelete_by_name_f')
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),type_id,space_id,attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5acreate_f')
  call h5awrite_f(attr_id, type_id, trim(attrValue), int([1],HSIZE_T), hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5awrite_f')
  call h5aclose_f(attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5aclose_f')
  call h5tclose_f(type_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5tclose_f')
  call h5sclose_f(space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_str: h5sclose_f')

end subroutine HDF5_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_int(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)            :: loc_id
  character(len=*), intent(in)            :: attrLabel
  integer(pInt),    intent(in)            :: attrValue
  character(len=*), intent(in), optional  :: path
  
  integer            :: hdferr
  integer(HID_T)     :: attr_id, space_id
  logical            :: attrExists
  character(len=256) :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
 
  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5screate_f')
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5aexists_by_name_f')
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5adelete_by_name_f')
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5acreate_f')
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, int([1],HSIZE_T), hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5awrite_f')
  call h5aclose_f(attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5tclose_f')
  call h5sclose_f(space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int: h5sclose_f')

end subroutine HDF5_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)           :: loc_id
  character(len=*), intent(in)           :: attrLabel
  real(pReal),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path
  
  integer            :: hdferr
  integer(HID_T)     :: attr_id, space_id
  logical            :: attrExists
  character(len=256) :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
 
  call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5screate_f')
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5aexists_by_name_f')
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5adelete_by_name_f')
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5acreate_f')
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, int([1],HSIZE_T), hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5awrite_f')
  call h5aclose_f(attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5tclose_f')
  call h5sclose_f(space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_real: h5sclose_f')

end subroutine HDF5_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief adds a integer attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_int_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  integer(pInt),    intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path
  
  integer                       :: hdferr
  integer(HID_T)                :: attr_id, space_id
  integer(HSIZE_T),dimension(1) :: array_size
  logical                       :: attrExists
  character(len=256)            :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
  
  array_size = size(attrValue,kind=HSIZE_T)
 
  call h5screate_simple_f(1, array_size, space_id, hdferr, array_size)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5screate_f')
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5aexists_by_name_f')
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5adelete_by_name_f')
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5acreate_f')
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attrValue, array_size, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5awrite_f')
  call h5aclose_f(attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5tclose_f')
  call h5sclose_f(space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5sclose_f')

end subroutine HDF5_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief adds a real attribute to the path given relative to the location
!--------------------------------------------------------------------------------------------------
subroutine HDF5_addAttribute_real_array(loc_id,attrLabel,attrValue,path)

  integer(HID_T),   intent(in)               :: loc_id
  character(len=*), intent(in)               :: attrLabel
  real(pReal),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path
  
  integer                       :: hdferr
  integer(HID_T)                :: attr_id, space_id
  integer(HSIZE_T),dimension(1) :: array_size
  logical                       :: attrExists
  character(len=256)            :: p
  
  if (present(path)) then
    p = trim(path)
  else
    p = '.'
  endif
  
  array_size = size(attrValue,kind=HSIZE_T)
 
  call h5screate_simple_f(1, array_size, space_id, hdferr, array_size)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5screate_f')
  call h5aexists_by_name_f(loc_id,trim(p),attrLabel,attrExists,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5aexists_by_name_f')
  if (attrExists) then
    call h5adelete_by_name_f(loc_id, trim(p), attrLabel, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5adelete_by_name_f')
  endif
  call h5acreate_by_name_f(loc_id,trim(p),trim(attrLabel),H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5acreate_f')
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attrValue, array_size, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5awrite_f')
  call h5aclose_f(attr_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5tclose_f')
  call h5sclose_f(space_id,hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_addAttribute_int_array: h5sclose_f')

end subroutine HDF5_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief set link to object in results file
!--------------------------------------------------------------------------------------------------
subroutine HDF5_setLink(loc_id,target_name,link_name)

  character(len=*), intent(in) :: target_name, link_name
  integer(HID_T),  intent(in) :: loc_id
  integer       :: hdferr
  logical                      :: linkExists
 
  call h5lexists_f(loc_id, link_name,linkExists, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5lexists_soft_f ('//trim(link_name)//')')
  if (linkExists) then
    call h5ldelete_f(loc_id,link_name, hdferr)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5ldelete_soft_f ('//trim(link_name)//')')
  endif
  call h5lcreate_soft_f(target_name, loc_id, link_name, hdferr)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg = 'HDF5_setLink: h5lcreate_soft_f ('//trim(target_name)//' '//trim(link_name)//')')

end subroutine HDF5_setLink


!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real1(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real1: h5dread_f')
 
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)
   
end subroutine HDF5_read_real1

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real2(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real2: h5dread_f')
 
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real2

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real3(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
  
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real3: h5dread_f')

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real3

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real4(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel

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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real4: h5dread_f')
 
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real4

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real5(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real5: h5dread_f')
 
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real5

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real6(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real6: h5dread_f')

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real6

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_real7(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel

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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_real7: h5dread_f')
 
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_real7


!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int1(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int1: h5dread_f')

  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)
  
end subroutine HDF5_read_int1

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int2(loc_id,dataset,datasetName,parallel)

  integer,           intent(inout), dimension(:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel

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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int2: h5dread_f')
  
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int2

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int3(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel

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
                        myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int3: h5dread_f')
   
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int3

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt withh 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int4(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int4: h5dread_f')
   
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int4

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int5(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int5: h5dread_f')
   
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int5

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int6(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel

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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int6: h5dread_f')
   
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int6

!--------------------------------------------------------------------------------------------------
!> @brief read dataset of type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_read_int7(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:,:,:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel
 
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
                         myStart, totalShape, loc_id,myShape,datasetName,.false.)
  endif
 
  call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dataset,totalShape, hdferr,&
                 file_space_id = filespace_id, xfer_prp = plist_id, mem_space_id = memspace_id)
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_read_int7: h5dread_f')
   
  call finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

end subroutine HDF5_read_int7


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real1(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape,loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real1: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real1

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real2(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real2: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real2

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real3(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real3: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real3

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real4(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif
 
  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real4: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real4


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real5(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif
 
  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real5: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real5

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real6(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real6: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real6

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pReal with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_real7(loc_id,dataset,datasetName,parallel)

  real(pReal),      intent(inout), dimension(:,:,:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_DOUBLE,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_real7: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_real7


!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 1 dimension
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int1(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int1: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int1

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 2 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int2(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int2: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int2

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 3 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int3(loc_id,dataset,datasetName,parallel)

  integer,           intent(inout), dimension(:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int3: h5dwrite_f')
  endif
  
  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int3

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 4 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int4(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int4: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int4

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 5 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int5(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int5: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int5

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 6 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int6(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int6: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int6

!--------------------------------------------------------------------------------------------------
!> @brief write dataset of type pInt with 7 dimensions
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_int7(loc_id,dataset,datasetName,parallel)

  integer,          intent(inout), dimension(:,:,:,:,:,:,:) ::    dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical, intent(in), optional :: parallel


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
                          myStart, totalShape, loc_id,myShape,datasetName,H5T_NATIVE_INTEGER,.false.)
  endif

  if (product(totalShape) /= 0) then
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,dataset,int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_int7: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_int7


!--------------------------------------------------------------------------------------------------
!> @brief writes a scalar orientation dataset
! ToDo: It might be possible to write the dataset as a whole
! ToDo: We could optionally write out other representations (axis angle, euler, ...)
!--------------------------------------------------------------------------------------------------
subroutine HDF5_write_rotation(loc_id,dataset,datasetName,parallel)

  type(rotation),   intent(in), dimension(:) :: dataset
  integer(HID_T),   intent(in) :: loc_id                                                            !< file or group handle
  character(len=*), intent(in) :: datasetName                                                       !< name of the dataset in the file
  logical,          intent(in), optional :: parallel

  integer :: hdferr
  real(pReal),      dimension(4,size(dataset)) :: dataset_asArray
  integer(HID_T)   :: dset_id, filespace_id, memspace_id, plist_id,dtype_id,w_id,x_id,y_id,z_id
  integer(HSIZE_T), dimension(size(shape(dataset))) :: &
    myStart, &
    myShape, &                                                                                      !< shape of the dataset (this process)
    totalShape                                                                                      !< shape of the dataset (all processes)
   integer(SIZE_T) :: type_size_real
   integer :: i

   do i = 1, size(dataset)
     dataset_asArray(1:4,i) = dataset(i)%asQuaternion()
   enddo

!---------------------------------------------------------------------------------------------------
! determine shape of dataset
   myShape = int(shape(dataset),HSIZE_T)
  
!---------------------------------------------------------------------------------------------------
! compound type: name of each quaternion component
  call h5tget_size_f(H5T_NATIVE_DOUBLE, type_size_real, hdferr)
  
  call h5tcreate_f(H5T_COMPOUND_F, type_size_real*4_SIZE_T, dtype_id, hdferr)
  call h5tinsert_f(dtype_id, "w",  type_size_real*0_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
  call h5tinsert_f(dtype_id, "x",  type_size_real*1_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
  call h5tinsert_f(dtype_id, "y",  type_size_real*2_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
  call h5tinsert_f(dtype_id, "z",  type_size_real*3_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
  
  if (present(parallel)) then
    call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,dtype_id,parallel)
  else
   call initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                          myStart, totalShape, loc_id,myShape,datasetName,dtype_id,.false.)
  endif

  call h5pset_preserve_f(plist_id, .TRUE., hdferr)
 
  if (product(totalShape) /= 0) then
    call h5tcreate_f(H5T_COMPOUND_F, type_size_real, x_id, hdferr)
    call h5tinsert_f(x_id, "x", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
    call h5tcreate_f(H5T_COMPOUND_F, type_size_real, w_id, hdferr)
    call h5tinsert_f(w_id, "w", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
    call h5tcreate_f(H5T_COMPOUND_F, type_size_real, y_id, hdferr)
    call h5tinsert_f(y_id, "y", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
    call h5tcreate_f(H5T_COMPOUND_F, type_size_real, z_id, hdferr)
    call h5tinsert_f(z_id, "z", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdferr)
    
    call h5dwrite_f(dset_id, w_id,dataset_asArray(1,:),int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call h5dwrite_f(dset_id, x_id,dataset_asArray(2,:),int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call h5dwrite_f(dset_id, y_id,dataset_asArray(3,:),int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    call h5dwrite_f(dset_id, z_id,dataset_asArray(4,:),int(totalShape,HSIZE_T), hdferr,&
                   file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    if (hdferr < 0) call IO_error(1_pInt,ext_msg='HDF5_write_rotation: h5dwrite_f')
  endif

  call finalize_write(plist_id, dset_id, filespace_id, memspace_id)

end subroutine HDF5_write_rotation


!--------------------------------------------------------------------------------------------------
!> @brief initialize HDF5 handles, determines global shape and start for parallel read
!--------------------------------------------------------------------------------------------------
subroutine initialize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id, &
                           myStart, globalShape, &
                           loc_id,localShape,datasetName,parallel)

 integer(HID_T),    intent(in) :: loc_id                                                             !< file or group handle
 character(len=*),  intent(in) :: datasetName                                                        !< name of the dataset in the file
 logical,           intent(in) :: parallel
 integer(HSIZE_T),  intent(in),   dimension(:) :: &
   localShape   
 integer(HSIZE_T),  intent(out), dimension(size(localShape,1)):: &
   myStart, &
   globalShape                                                                                      !< shape of the dataset (all processes)
 integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
    
 integer(pInt), dimension(worldsize) :: &
   readSize                                                                                         !< contribution of all processes
 integer :: ierr
 integer :: hdferr
 
!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective for MPI)
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5pcreate_f')

!--------------------------------------------------------------------------------------------------
 readSize = 0_pInt
 readSize(worldrank+1) = int(localShape(ubound(localShape,1)),pInt)
#ifdef PETSc
 if (parallel) then
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5pset_dxpl_mpio_f')
   call MPI_allreduce(MPI_IN_PLACE,readSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)       ! get total output size over each process
   if (ierr /= 0) call IO_error(894_pInt,ext_msg='initialize_read: MPI_allreduce')
 endif
#endif
 myStart                   = int(0,HSIZE_T)
 myStart(ubound(myStart))  = int(sum(readSize(1:worldrank)),HSIZE_T)
 globalShape = [localShape(1:ubound(localShape,1)-1),int(sum(readSize),HSIZE_T)]

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape)
 call h5screate_simple_f(size(localShape), localShape, memspace_id, hdferr, localShape)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5screate_simple_f/memspace_id')

!--------------------------------------------------------------------------------------------------
! creating a property list for IO and set it to collective
 call h5pcreate_f(H5P_DATASET_ACCESS_F, aplist_id, hdferr) 
  if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5pcreate_f')
#ifdef PETSc
 call h5pset_all_coll_metadata_ops_f(aplist_id, .true., hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5pset_all_coll_metadata_ops_f')
#endif

!--------------------------------------------------------------------------------------------------
! open the dataset in the file and get the space ID
 call h5dopen_f(loc_id,datasetName,dset_id,hdferr, dapl_id = aplist_id)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5dopen_f')
 call h5dget_space_f(dset_id, filespace_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5dget_space_f')

!--------------------------------------------------------------------------------------------------
! select a hyperslab (the portion of the current process) in the file
 call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, localShape, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_read: h5sselect_hyperslab_f')

end subroutine initialize_read


!--------------------------------------------------------------------------------------------------
!> @brief closes HDF5 handles
!--------------------------------------------------------------------------------------------------
subroutine finalize_read(dset_id, filespace_id, memspace_id, plist_id, aplist_id)

 integer(HID_T), intent(in)   :: dset_id, filespace_id, memspace_id, plist_id, aplist_id
 integer :: hdferr

 call h5pclose_f(plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_read: plist_id')
 call h5pclose_f(aplist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_read: aplist_id')
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_read: h5dclose_f')
 call h5sclose_f(filespace_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_read: h5sclose_f/filespace_id')
 call h5sclose_f(memspace_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_read: h5sclose_f/memspace_id')

end subroutine finalize_read


!--------------------------------------------------------------------------------------------------
!> @brief initialize HDF5 handles, determines global shape and start for parallel write
!--------------------------------------------------------------------------------------------------
subroutine initialize_write(dset_id, filespace_id, memspace_id, plist_id, &
                           myStart, totalShape, &
                           loc_id,myShape,datasetName,datatype,parallel)

 integer(HID_T),    intent(in) :: loc_id                                                             !< file or group handle
 character(len=*),  intent(in) :: datasetName                                                        !< name of the dataset in the file
 logical,           intent(in) :: parallel
 integer(HID_T),    intent(in) :: datatype
 integer(HSIZE_T),  intent(in),   dimension(:) :: &
   myShape   
 integer(HSIZE_T),  intent(out), dimension(size(myShape,1)):: &
   myStart, &
   totalShape                                                                                      !< shape of the dataset (all processes)
 integer(HID_T),    intent(out) :: dset_id, filespace_id, memspace_id, plist_id
    
 integer, dimension(worldsize) :: &
   writeSize                                                                                        !< contribution of all processes
 integer :: ierr
 integer :: hdferr

!-------------------------------------------------------------------------------------------------
! creating a property list for transfer properties (is collective when reading in parallel)
 call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5pcreate_f')
#ifdef PETSc
if (parallel) then
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
   if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5pset_dxpl_mpio_f')
 endif
#endif

!--------------------------------------------------------------------------------------------------
! determine the global data layout among all processes
 writeSize              = 0_pInt
 writeSize(worldrank+1) = int(myShape(ubound(myShape,1)),pInt)
#ifdef PETSc
if (parallel) then
   call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)       ! get total output size over each process
   if (ierr /= 0) call IO_error(894_pInt,ext_msg='initialize_write: MPI_allreduce')
 endif
#endif
 myStart                   = int(0,HSIZE_T)
 myStart(ubound(myStart))  = int(sum(writeSize(1:worldrank)),HSIZE_T)
 totalShape = [myShape(1:ubound(myShape,1)-1),int(sum(writeSize),HSIZE_T)]

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape) and in file (global shape)
 call h5screate_simple_f(size(myShape), myShape, memspace_id, hdferr, myShape)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5dopen_f')
 call h5screate_simple_f(size(totalShape), totalShape, filespace_id, hdferr, totalShape)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5dget_space_f')

!--------------------------------------------------------------------------------------------------
! create dataset in the file and select a hyperslab from it (the portion of the current process)
 call h5dcreate_f(loc_id, trim(datasetName), datatype, filespace_id, dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5dcreate_f')
 call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myStart, myShape, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='initialize_write: h5sselect_hyperslab_f')


end subroutine initialize_write


!--------------------------------------------------------------------------------------------------
!> @brief closes HDF5 handles
!--------------------------------------------------------------------------------------------------
subroutine finalize_write(plist_id, dset_id, filespace_id, memspace_id)

 integer(HID_T), intent(in) :: dset_id, filespace_id, memspace_id, plist_id
 integer :: hdferr

 call h5pclose_f(plist_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_write: plist_id')
 call h5dclose_f(dset_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_write: h5dclose_f')
 call h5sclose_f(filespace_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_write: h5sclose_f/filespace_id')
 call h5sclose_f(memspace_id, hdferr)
 if (hdferr < 0) call IO_error(1_pInt,ext_msg='finalize_write: h5sclose_f/memspace_id')

end subroutine finalize_write

end module HDF5_Utilities
