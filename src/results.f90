!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module results
  use DAMASK_interface
  use rotations
  use numerics
  use HDF5_utilities
#ifdef PETSc
  use PETSC
#endif

  implicit none
  private
  
#if defined(PETSc) || defined(DAMASK_HDF5)
  integer(HID_T) :: resultsFile

  interface results_writeDataset
  
    module procedure results_writeTensorDataset_real
    module procedure results_writeVectorDataset_real
    module procedure results_writeScalarDataset_real
    
    module procedure results_writeTensorDataset_int
    module procedure results_writeVectorDataset_int
    
    module procedure results_writeScalarDataset_rotation
    
  end interface results_writeDataset
  
  interface results_addAttribute
  
    module procedure results_addAttribute_real
    module procedure results_addAttribute_int
    module procedure results_addAttribute_str
    
    module procedure results_addAttribute_int_array
    module procedure results_addAttribute_real_array
    
  end interface results_addAttribute

  public :: &
    results_init, &
    results_openJobFile, &
    results_closeJobFile, &
    results_addIncrement, &
    results_addGroup, &
    results_openGroup, &
    results_writeDataset, &
    results_setLink, &
    results_addAttribute, &
    results_removeLink, &
    results_mapping_constituent, &
    results_mapping_materialpoint
contains

subroutine results_init

  character(len=pStringLen) :: commandLine

  write(6,'(/,a)') ' <<<+-  results init  -+>>>'

  write(6,'(/,a)') ' Diehl et al., Integrating Materials and Manufacturing Innovation 6(1):83–91, 2017'
  write(6,'(a)')   ' https://doi.org/10.1007/s40192-017-0084-5'

  resultsFile = HDF5_openFile(trim(getSolverJobName())//'.hdf5','w',.true.)
  call HDF5_addAttribute(resultsFile,'DADF5-version',0.2)
  call HDF5_addAttribute(resultsFile,'DADF5-major',0)
  call HDF5_addAttribute(resultsFile,'DADF5-minor',2)
  call HDF5_addAttribute(resultsFile,'DAMASK',DAMASKVERSION)
  call get_command(commandLine)
  call HDF5_addAttribute(resultsFile,'call',trim(commandLine))
  call HDF5_closeGroup(results_addGroup('mapping'))
  call HDF5_closeGroup(results_addGroup('mapping/cellResults'))
  call HDF5_closeFile(resultsFile)

end subroutine results_init


!--------------------------------------------------------------------------------------------------
!> @brief opens the results file to append data
!--------------------------------------------------------------------------------------------------
subroutine results_openJobFile

  resultsFile = HDF5_openFile(trim(getSolverJobName())//'.hdf5','a',.true.)
 
end subroutine results_openJobFile


!--------------------------------------------------------------------------------------------------
!> @brief closes the results file
!--------------------------------------------------------------------------------------------------
subroutine results_closeJobFile

  call HDF5_closeFile(resultsFile)

end subroutine results_closeJobFile


!--------------------------------------------------------------------------------------------------
!> @brief creates the group of increment and adds time as attribute to the file
!--------------------------------------------------------------------------------------------------
subroutine results_addIncrement(inc,time)
 
  integer,       intent(in) :: inc
  real(pReal),   intent(in) :: time
  character(len=pStringLen) :: incChar

  write(incChar,'(i5.5)') inc                                                                       ! allow up to 99999 increments
  call HDF5_closeGroup(results_addGroup(trim('inc'//trim(adjustl(incChar)))))
  call results_setLink(trim('inc'//trim(adjustl(incChar))),'current')
  call HDF5_addAttribute(resultsFile,'time/s',time,trim('inc'//trim(adjustl(incChar))))
  
  call HDF5_closeGroup(results_addGroup('current/constituent'))
  call HDF5_closeGroup(results_addGroup('current/materialpoint'))

end subroutine results_addIncrement

!--------------------------------------------------------------------------------------------------
!> @brief open a group from the results file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function results_openGroup(groupName)

  character(len=*), intent(in) :: groupName
  
  results_openGroup = HDF5_openGroup(resultsFile,groupName)

end function results_openGroup


!--------------------------------------------------------------------------------------------------
!> @brief adds a new group to the results file
!--------------------------------------------------------------------------------------------------
integer(HID_T) function results_addGroup(groupName)

  character(len=*), intent(in) :: groupName
  
  results_addGroup = HDF5_addGroup(resultsFile,groupName)

end function results_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief set link to object in results file
!--------------------------------------------------------------------------------------------------
subroutine results_setLink(path,link)

  character(len=*), intent(in) :: path, link

  call HDF5_setLink(resultsFile,path,link)

end subroutine results_setLink


!--------------------------------------------------------------------------------------------------
!> @brief adds a string attribute to an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_str(attrLabel,attrValue,path)

  character(len=*), intent(in) :: attrLabel, attrValue, path

  call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)

end subroutine results_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief adds an integer attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int(attrLabel,attrValue,path)

  character(len=*), intent(in) :: attrLabel, path
  integer,          intent(in) :: attrValue

  call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)

end subroutine results_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief adds a real attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real(attrLabel,attrValue,path)

  character(len=*), intent(in) :: attrLabel, path
  real(pReal),      intent(in) :: attrValue

  call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)

end subroutine results_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief adds an integer array attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel, path
  integer,          intent(in), dimension(:) :: attrValue

  call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)

end subroutine results_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief adds a real array attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel, path
  real(pReal),      intent(in), dimension(:) :: attrValue

  call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)

end subroutine results_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief remove link to an object
!--------------------------------------------------------------------------------------------------
subroutine results_removeLink(link)

  character(len=*), intent(in) :: link
  integer                      :: hdferr

  call h5ldelete_f(resultsFile,link, hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg = 'results_removeLink: h5ldelete_soft_f ('//trim(link)//')')

end subroutine results_removeLink


!--------------------------------------------------------------------------------------------------
!> @brief stores a scalar dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeScalarDataset_real(group,dataset,label,description,SIunit)

  character(len=*), intent(in)                  :: label,group,description
  character(len=*), intent(in),    optional     :: SIunit
  real(pReal),      intent(inout), dimension(:) :: dataset
  
  integer(HID_T) :: groupHandle
 
  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset,label,.true.)
#else
  call HDF5_write(groupHandle,dataset,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(SIunit)) &
    call HDF5_addAttribute(groupHandle,'Unit',SIunit,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeScalarDataset_real

!--------------------------------------------------------------------------------------------------
!> @brief stores a vector dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeVectorDataset_real(group,dataset,label,description,SIunit)

  character(len=*), intent(in)                    :: label,group,description
  character(len=*), intent(in),    optional       :: SIunit
  real(pReal),      intent(inout), dimension(:,:) :: dataset
  
  integer(HID_T) :: groupHandle
 
  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset,label,.true.)
#else
  call HDF5_write(groupHandle,dataset,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(SIunit)) &
    call HDF5_addAttribute(groupHandle,'Unit',SIunit,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeVectorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief stores a tensor dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeTensorDataset_real(group,dataset,label,description,SIunit)

  character(len=*), intent(in)                   :: label,group,description
  character(len=*), intent(in), optional         :: SIunit
  real(pReal),      intent(in), dimension(:,:,:) :: dataset
  
  integer :: i 
  integer(HID_T) :: groupHandle
  real(pReal), dimension(:,:,:), allocatable :: dataset_transposed


  allocate(dataset_transposed,mold=dataset)
  do i=1,size(dataset,3)
    dataset_transposed(1:3,1:3,i) = transpose(dataset(1:3,1:3,i))
  enddo

  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset_transposed,label,.true.)
#else
  call HDF5_write(groupHandle,dataset_transposed,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(SIunit)) &
    call HDF5_addAttribute(groupHandle,'Unit',SIunit,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeTensorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief stores a vector dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeVectorDataset_int(group,dataset,label,description,SIunit)

  character(len=*), intent(in)                :: label,group,description
  character(len=*), intent(in), optional      :: SIunit
  integer,      intent(inout), dimension(:,:) :: dataset
  
  integer(HID_T) :: groupHandle
 
  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset,label,.true.)
#else
  call HDF5_write(groupHandle,dataset,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(SIunit)) &
    call HDF5_addAttribute(groupHandle,'Unit',SIunit,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeVectorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief stores a tensor dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeTensorDataset_int(group,dataset,label,description,SIunit)

  character(len=*), intent(in)                  :: label,group,description
  character(len=*), intent(in), optional        :: SIunit
  integer,      intent(inout), dimension(:,:,:) :: dataset
  
  integer(HID_T) :: groupHandle
 
  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset,label,.true.)
#else
  call HDF5_write(groupHandle,dataset,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(SIunit)) &
    call HDF5_addAttribute(groupHandle,'Unit',SIunit,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeTensorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief stores a scalar dataset in a group
!--------------------------------------------------------------------------------------------------
subroutine results_writeScalarDataset_rotation(group,dataset,label,description,lattice_structure)

  character(len=*), intent(in)                  :: label,group,description
  character(len=*), intent(in), optional        :: lattice_structure
  type(rotation),   intent(inout), dimension(:) :: dataset
  
  integer(HID_T) :: groupHandle
 
  groupHandle = results_openGroup(group)
  
#ifdef PETSc
  call HDF5_write(groupHandle,dataset,label,.true.)
#else
  call HDF5_write(groupHandle,dataset,label,.false.)
#endif
  
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Description',description,label)
  if (HDF5_objectExists(groupHandle,label) .and. present(lattice_structure)) &
    call HDF5_addAttribute(groupHandle,'Lattice',lattice_structure,label)
  if (HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(groupHandle,'Creator','DAMASK '//DAMASKVERSION,label)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeScalarDataset_rotation


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine results_mapping_constituent(phaseAt,memberAt,label)
    
  integer,           dimension(:,:),   intent(in)  :: phaseAt                                       !< phase section at (constituent,element)
  integer,           dimension(:,:,:), intent(in)  :: memberAt                                      !< phase member at (constituent,IP,element)
  character(len=64), dimension(:),     intent(in)  :: label                                         !< label of each phase section
  
  integer, dimension(size(memberAt,1),size(memberAt,2),size(memberAt,3)) :: &
    phaseAt_perIP, &
    memberAt_total
  integer, dimension(size(label),0:worldsize-1) :: memberOffset                                     !< offset in member counting per process
  integer, dimension(0:worldsize-1)             :: writeSize                                        !< amount of data written per process
  integer(HSIZE_T), dimension(2) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)
  
  integer(HID_T) :: &
    loc_id, &                                                                                       !< identifier of group in file
    dtype_id, &                                                                                     !< identifier of compound data type
    name_id, &                                                                                      !< identifier of name (string) in compound data type
    position_id, &                                                                                  !< identifier of position/index (integer) in compound data type
    dset_id, &
    memspace_id, &
    filespace_id, &
    plist_id, &
    dt_id

  
  integer(SIZE_T) :: type_size_string, type_size_int
  integer         :: ierr, i
  
!---------------------------------------------------------------------------------------------------
! compound type: name of phase section + position/index within results array
  call h5tcopy_f(H5T_NATIVE_CHARACTER, dt_id, ierr)
  call h5tset_size_f(dt_id, int(len(label(1)),SIZE_T), ierr)
  call h5tget_size_f(dt_id, type_size_string, ierr)
  
  call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_int, ierr)
  
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, ierr)
  call h5tinsert_f(dtype_id, "Name", 0_SIZE_T, dt_id,ierr)
  call h5tinsert_f(dtype_id, "Position", type_size_string, H5T_NATIVE_INTEGER, ierr)
  
!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string, name_id, ierr)
  call h5tinsert_f(name_id, "Name", 0_SIZE_T, dt_id, ierr)
  
  call h5tcreate_f(H5T_COMPOUND_F, type_size_int, position_id, ierr)
  call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_NATIVE_INTEGER, ierr)
  
  call h5tclose_f(dt_id, ierr)

!--------------------------------------------------------------------------------------------------
! prepare MPI communication (transparent for non-MPI runs)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
  memberOffset = 0
  do i=1, size(label)
    memberOffset(i,worldrank) = count(phaseAt == i)*size(memberAt,2)                                ! number of points/instance of this process
  enddo
  writeSize = 0
  writeSize(worldrank) = size(memberAt(1,:,:))                                                      ! total number of points by this process

!--------------------------------------------------------------------------------------------------
! MPI settings and communication
#ifdef PETSc
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5pset_dxpl_mpio_f')
  
  call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)        ! get output at each process
  if (ierr /= 0) call IO_error(894,ext_msg='results_mapping_constituent: MPI_allreduce/writeSize')
  
  call MPI_allreduce(MPI_IN_PLACE,memberOffset,size(memberOffset),MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)! get offset at each process
  if (ierr /= 0) call IO_error(894,ext_msg='results_mapping_constituent: MPI_allreduce/memberOffset')
#endif

  myShape    = int([size(phaseAt,1),writeSize(worldrank)],  HSIZE_T)
  myOffset   = int([0,sum(writeSize(0:worldrank-1))],       HSIZE_T)
  totalShape = int([size(phaseAt,1),sum(writeSize)],        HSIZE_T)
  
!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call h5screate_simple_f(2,myShape,memspace_id,ierr,myShape)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5screate_simple_f/memspace_id')
  
  call h5screate_simple_f(2,totalShape,filespace_id,ierr,totalShape)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5screate_simple_f/filespace_id')
  
  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5sselect_hyperslab_f')

!---------------------------------------------------------------------------------------------------
! expand phaseAt to consider IPs (is not stored per IP)
  do i = 1, size(phaseAt_perIP,2)
    phaseAt_perIP(:,i,:) = phaseAt
  enddo
  
!---------------------------------------------------------------------------------------------------
! renumber member from my process to all processes
  do i = 1, size(label)
    where(phaseAt_perIP == i) memberAt_total = memberAt + sum(memberOffset(i,0:worldrank-1)) -1     ! convert to 0-based
  enddo

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call h5pset_preserve_f(plist_id, .TRUE., ierr)
  
  loc_id = results_openGroup('/mapping/cellResults')
  call h5dcreate_f(loc_id, 'constituent', dtype_id, filespace_id, dset_id, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5dcreate_f')
  
  call h5dwrite_f(dset_id, name_id, reshape(label(pack(phaseAt_perIP,.true.)),myShape), &
                  myShape, ierr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5dwrite_f/name_id')
  call h5dwrite_f(dset_id, position_id, reshape(pack(memberAt_total,.true.),myShape), &
                  myShape, ierr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_constituent: h5dwrite_f/position_id')

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call h5pclose_f(plist_id, ierr)
  call h5sclose_f(filespace_id, ierr)
  call h5sclose_f(memspace_id, ierr)
  call h5dclose_f(dset_id, ierr)
  call h5tclose_f(dtype_id, ierr)
  call h5tclose_f(name_id, ierr)
  call h5tclose_f(position_id, ierr)

end subroutine results_mapping_constituent


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine results_mapping_materialpoint(homogenizationAt,memberAt,label)
    
  integer,           dimension(:),   intent(in)  :: homogenizationAt                                !< homogenization section at (element)
  integer,           dimension(:,:), intent(in)  :: memberAt                                        !< homogenization member at (IP,element)
  character(len=64), dimension(:),   intent(in)  :: label                                           !< label of each homogenization section
  
  integer, dimension(size(memberAt,1),size(memberAt,2)) :: &
    homogenizationAt_perIP, &
    memberAt_total
  integer, dimension(size(label),0:worldsize-1) :: memberOffset                                     !< offset in member counting per process
  integer, dimension(0:worldsize-1)             :: writeSize                                        !< amount of data written per process
  integer(HSIZE_T), dimension(1) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)
  
  integer(HID_T) :: &
    loc_id, &                                                                                       !< identifier of group in file
    dtype_id, &                                                                                     !< identifier of compound data type
    name_id, &                                                                                      !< identifier of name (string) in compound data type
    position_id, &                                                                                  !< identifier of position/index (integer) in compound data type
    dset_id, &
    memspace_id, &
    filespace_id, &
    plist_id, &
    dt_id

  
  integer(SIZE_T) :: type_size_string, type_size_int
  integer         :: ierr, i
  
!---------------------------------------------------------------------------------------------------
! compound type: name of phase section + position/index within results array
  call h5tcopy_f(H5T_NATIVE_CHARACTER, dt_id, ierr)
  call h5tset_size_f(dt_id, int(len(label(1)),SIZE_T), ierr)
  call h5tget_size_f(dt_id, type_size_string, ierr)
  
  call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_int, ierr)
  
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, ierr)
  call h5tinsert_f(dtype_id, "Name", 0_SIZE_T, dt_id,ierr)
  call h5tinsert_f(dtype_id, "Position", type_size_string, H5T_NATIVE_INTEGER, ierr)
  
!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string, name_id, ierr)
  call h5tinsert_f(name_id, "Name", 0_SIZE_T, dt_id, ierr)
  
  call h5tcreate_f(H5T_COMPOUND_F, type_size_int, position_id, ierr)
  call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_NATIVE_INTEGER, ierr)
  
  call h5tclose_f(dt_id, ierr)

!--------------------------------------------------------------------------------------------------
! prepare MPI communication (transparent for non-MPI runs)
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
  memberOffset = 0
  do i=1, size(label)
    memberOffset(i,worldrank) = count(homogenizationAt == i)*size(memberAt,1)                       ! number of points/instance of this process
  enddo
  writeSize = 0
  writeSize(worldrank) = size(memberAt)                                                             ! total number of points by this process

!--------------------------------------------------------------------------------------------------
! MPI settings and communication
#ifdef PETSc
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5pset_dxpl_mpio_f')
  
  call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)        ! get output at each process
  if (ierr /= 0) call IO_error(894,ext_msg='results_mapping_materialpoint: MPI_allreduce/writeSize')
  
  call MPI_allreduce(MPI_IN_PLACE,memberOffset,size(memberOffset),MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)! get offset at each process
  if (ierr /= 0) call IO_error(894,ext_msg='results_mapping_materialpoint: MPI_allreduce/memberOffset')
#endif

  myShape    = int([writeSize(worldrank)],          HSIZE_T)
  myOffset   = int([sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([sum(writeSize)],                HSIZE_T)
  
!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call h5screate_simple_f(1,myShape,memspace_id,ierr,myShape)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5screate_simple_f/memspace_id')
  
  call h5screate_simple_f(1,totalShape,filespace_id,ierr,totalShape)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5screate_simple_f/filespace_id')
  
  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5sselect_hyperslab_f')

!---------------------------------------------------------------------------------------------------
! expand phaseAt to consider IPs (is not stored per IP)
  do i = 1, size(homogenizationAt_perIP,1)
    homogenizationAt_perIP(i,:) = homogenizationAt
  enddo
  
!---------------------------------------------------------------------------------------------------
! renumber member from my process to all processes
  do i = 1, size(label)
    where(homogenizationAt_perIP == i) memberAt_total = memberAt + sum(memberOffset(i,0:worldrank-1)) - 1  ! convert to 0-based
  enddo

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call h5pset_preserve_f(plist_id, .TRUE., ierr)
  
  loc_id = results_openGroup('/mapping/cellResults')
  call h5dcreate_f(loc_id, 'materialpoint', dtype_id, filespace_id, dset_id, ierr)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5dcreate_f')
  
  call h5dwrite_f(dset_id, name_id, reshape(label(pack(homogenizationAt_perIP,.true.)),myShape), &
                  myShape, ierr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5dwrite_f/name_id')
  call h5dwrite_f(dset_id, position_id, reshape(pack(memberAt_total,.true.),myShape), &
                  myShape, ierr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if (ierr < 0) call IO_error(1,ext_msg='results_mapping_materialpoint: h5dwrite_f/position_id')

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call h5pclose_f(plist_id, ierr)
  call h5sclose_f(filespace_id, ierr)
  call h5sclose_f(memspace_id, ierr)
  call h5dclose_f(dset_id, ierr)
  call h5tclose_f(dtype_id, ierr)
  call h5tclose_f(name_id, ierr)
  call h5tclose_f(position_id, ierr)

end subroutine results_mapping_materialpoint


!!--------------------------------------------------------------------------------------------------
!!> @brief adds the backward mapping from spatial position and constituent ID to results
!!--------------------------------------------------------------------------------------------------
!subroutine HDF5_backwardMappingPhase(material_phase,phasememberat,phase_name,dataspace_size,mpiOffset,mpiOffset_phase)
! use hdf5

! integer(pInt),    intent(in), dimension(:,:,:) :: material_phase, phasememberat
! character(len=*), intent(in), dimension(:)     :: phase_name
! integer(pInt),    intent(in), dimension(:)     :: dataspace_size, mpiOffset_phase
! integer(pInt),    intent(in)                   :: mpiOffset

! integer(pInt)   :: hdferr, NmatPoints, Nconstituents, i, j
! integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, position_id, plist_id, memspace
! integer(SIZE_T) :: type_size

! integer(pInt), dimension(:,:), allocatable :: arr

! integer(HSIZE_T),  dimension(1) :: counter
! integer(HSSIZE_T), dimension(1) :: fileOffset

! character(len=64) :: phaseID

! Nconstituents = size(phasememberat,1)
! NmatPoints = count(material_phase /=0)/Nconstituents

! allocate(arr(2,NmatPoints*Nconstituents))

! do i=1, NmatPoints
!   do j=Nconstituents-1, 0, -1
!     arr(1,Nconstituents*i-j) = i-1
!   enddo
! enddo
! arr(2,:) = pack(material_phase,material_phase/=0)

! do i=1, size(phase_name)
!   write(phaseID, '(i0)') i
!   mapping_ID = results_openGroup('/current/constitutive/'//trim(phaseID)//'_'//phase_name(i))
!   NmatPoints = count(material_phase == i)

!!--------------------------------------------------------------------------------------------------
!  ! create dataspace
!   call h5screate_simple_f(1, int([dataspace_size(i)],HSIZE_T), space_id, hdferr, &
!                              int([dataspace_size(i)],HSIZE_T))
!   if (hdferr < 0) call IO_error(1,ext_msg='HDF5_writeBackwardMapping')

!!--------------------------------------------------------------------------------------------------
!  ! compound type
!   call h5tget_size_f(H5T_STD_I32LE, type_size, hdferr)
!   call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='HDF5_writeBackwardMapping: h5tcreate_f dtype_id')

!   call h5tinsert_f(dtype_id, "Position",          0_SIZE_T, H5T_STD_I32LE, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5tinsert_f 0')

!!--------------------------------------------------------------------------------------------------
!  ! create Dataset
!   call h5dcreate_f(mapping_id, 'mapGeometry', dtype_id, space_id, dset_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase')

!!--------------------------------------------------------------------------------------------------
!  ! Create memory types (one compound datatype for each member)
!   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5tcreate_f position_id')
!   call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5tinsert_f position_id')

!!--------------------------------------------------------------------------------------------------
!  ! Define and select hyperslabs
!   counter = NmatPoints                       ! how big i am
!   fileOffset = mpiOffset_phase(i)            ! where i start to write my data

!   call h5screate_simple_f(1, counter, memspace, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5screate_simple_f')
!   call h5dget_space_f(dset_id, space_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5dget_space_f')
!   call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5sselect_hyperslab_f')

!!--------------------------------------------------------------------------------------------------
! ! Create property list for collective dataset write
!#ifdef PETSc
!   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5pcreate_f')
!   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5pset_dxpl_mpio_f')
!#endif

!!--------------------------------------------------------------------------------------------------
!  ! write data by fields in the datatype. Fields order is not important.
!   call h5dwrite_f(dset_id, position_id, pack(arr(1,:),arr(2,:)==i)+mpiOffset, int([dataspace_size(i)],HSIZE_T),&
!                              hdferr, file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5dwrite_f instance_id')

!!--------------------------------------------------------------------------------------------------
!  !close types, dataspaces
!   call h5tclose_f(dtype_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5tclose_f dtype_id')
!   call h5tclose_f(position_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5tclose_f position_id')
!   call h5dclose_f(dset_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5dclose_f')
!   call h5sclose_f(space_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5sclose_f space_id')
!   call h5sclose_f(memspace, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5sclose_f memspace')
!   call h5pclose_f(plist_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingPhase: h5pclose_f')
!   call HDF5_closeGroup(mapping_ID)

! enddo

!end subroutine HDF5_backwardMappingPhase


!!--------------------------------------------------------------------------------------------------
!!> @brief adds the backward mapping from spatial position and constituent ID to results
!!--------------------------------------------------------------------------------------------------
!subroutine HDF5_backwardMappingHomog(material_homog,homogmemberat,homogenization_name,dataspace_size,mpiOffset,mpiOffset_homog)
! use hdf5

! integer(pInt),    intent(in), dimension(:,:) :: material_homog, homogmemberat
! character(len=*), intent(in), dimension(:)   :: homogenization_name
! integer(pInt),    intent(in), dimension(:)   :: dataspace_size, mpiOffset_homog
! integer(pInt),    intent(in)                 :: mpiOffset

! integer(pInt)   :: hdferr, NmatPoints, i
! integer(HID_T)  :: mapping_id, dtype_id, dset_id, space_id, position_id, plist_id, memspace
! integer(SIZE_T) :: type_size

! integer(pInt),     dimension(:,:), allocatable :: arr

! integer(HSIZE_T),  dimension(1) :: counter
! integer(HSSIZE_T), dimension(1) :: fileOffset

! character(len=64) :: homogID

! NmatPoints = count(material_homog /=0)
! allocate(arr(2,NmatPoints))

! arr(1,:) = (/(i, i=0,NmatPoints-1)/)
! arr(2,:) = pack(material_homog,material_homog/=0)

! do i=1, size(homogenization_name)
!   write(homogID, '(i0)') i
!   mapping_ID = results_openGroup('/current/homogenization/'//trim(homogID)//'_'//homogenization_name(i))

!!--------------------------------------------------------------------------------------------------
!  ! create dataspace
!   call h5screate_simple_f(1, int([dataspace_size(i)],HSIZE_T), space_id, hdferr, &
!                              int([dataspace_size(i)],HSIZE_T))
!   if (hdferr < 0) call IO_error(1,ext_msg='HDF5_writeBackwardMapping')

!!--------------------------------------------------------------------------------------------------
!  ! compound type
!   call h5tget_size_f(H5T_STD_I32LE, type_size, hdferr)
!   call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='HDF5_writeBackwardMapping: h5tcreate_f dtype_id')

!   call h5tinsert_f(dtype_id, "Position",          0_SIZE_T, H5T_STD_I32LE, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5tinsert_f 0')

!!--------------------------------------------------------------------------------------------------
!  ! create Dataset
!   call h5dcreate_f(mapping_id, 'mapGeometry', dtype_id, space_id, dset_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog')

!!--------------------------------------------------------------------------------------------------
!  ! Create memory types (one compound datatype for each member)
!   call h5tcreate_f(H5T_COMPOUND_F, int(pInt,SIZE_T), position_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5tcreate_f position_id')
!   call h5tinsert_f(position_id, "Position", 0_SIZE_T, H5T_STD_I32LE, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5tinsert_f position_id')

!!--------------------------------------------------------------------------------------------------
!  ! Define and select hyperslabs
!   counter = NmatPoints                           ! how big i am
!   fileOffset = mpiOffset_homog(i)                ! where i start to write my data

!   call h5screate_simple_f(1, counter, memspace, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5screate_simple_f')
!   call h5dget_space_f(dset_id, space_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5dget_space_f')
!   call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, fileOffset, counter, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5sselect_hyperslab_f')

!!--------------------------------------------------------------------------------------------------
! ! Create property list for collective dataset write
!#ifdef PETSc
!   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5pcreate_f')
!   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5pset_dxpl_mpio_f')
!#endif

!!--------------------------------------------------------------------------------------------------
!  ! write data by fields in the datatype. Fields order is not important.
!   call h5dwrite_f(dset_id, position_id, pack(arr(1,:),arr(2,:)==i)+mpiOffset,int([dataspace_size(i)],HSIZE_T),&
!                   hdferr, file_space_id = space_id, mem_space_id = memspace, xfer_prp = plist_id)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5dwrite_f instance_id')

!!--------------------------------------------------------------------------------------------------
!  !close types, dataspaces
!   call h5tclose_f(dtype_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5tclose_f dtype_id')
!   call h5tclose_f(position_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5tclose_f position_id')
!   call h5dclose_f(dset_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5dclose_f')
!   call h5sclose_f(space_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5sclose_f space_id')
!   call h5sclose_f(memspace, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5sclose_f memspace')
!   call h5pclose_f(plist_id, hdferr)
!   if (hdferr < 0) call IO_error(1,ext_msg='IO_backwardMappingHomog: h5pclose_f')
!   call HDF5_closeGroup(mapping_ID)

! enddo

!end subroutine HDF5_backwardMappingHomog


!!--------------------------------------------------------------------------------------------------
!!> @brief adds the unique cell to node mapping
!!--------------------------------------------------------------------------------------------------
!subroutine HDF5_mappingCells(mapping)
! use hdf5

! integer(pInt), intent(in), dimension(:) :: mapping

! integer        :: hdferr, Nnodes
! integer(HID_T) :: mapping_id, dset_id, space_id

! Nnodes=size(mapping)
! mapping_ID = results_openGroup("mapping")

!!--------------------------------------------------------------------------------------------------
!! create dataspace
! call h5screate_simple_f(1, int([Nnodes],HSIZE_T), space_id, hdferr, &
!                            int([Nnodes],HSIZE_T))
! if (hdferr < 0) call IO_error(1,ext_msg='IO_mappingCells: h5screate_simple_f')

!!--------------------------------------------------------------------------------------------------
!! create Dataset
! call h5dcreate_f(mapping_id, "Cell",H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
! if (hdferr < 0) call IO_error(1,ext_msg='IO_mappingCells')

!!--------------------------------------------------------------------------------------------------
!! write data by fields in the datatype. Fields order is not important.
! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, mapping, int([Nnodes],HSIZE_T), hdferr)
! if (hdferr < 0) call IO_error(1,ext_msg='IO_mappingCells: h5dwrite_f instance_id')

!!--------------------------------------------------------------------------------------------------
!!close types, dataspaces
! call h5dclose_f(dset_id, hdferr)
! if (hdferr < 0) call IO_error(1,ext_msg='IO_mappingConstitutive: h5dclose_f')
! call h5sclose_f(space_id, hdferr)
! if (hdferr < 0) call IO_error(1,ext_msg='IO_mappingConstitutive: h5sclose_f')
! call HDF5_closeGroup(mapping_ID)

!end subroutine HDF5_mappingCells

#endif
end module results
