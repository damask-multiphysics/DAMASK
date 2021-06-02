!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module results
  use DAMASK_interface
  use parallelization
  use IO
  use HDF5_utilities
#ifdef PETSc
  use PETSC
#endif

  implicit none
  private

  integer(HID_T) :: resultsFile

  interface results_writeDataset
    module procedure results_writeTensorDataset_real
    module procedure results_writeVectorDataset_real
    module procedure results_writeScalarDataset_real

    module procedure results_writeTensorDataset_int
    module procedure results_writeVectorDataset_int
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
    results_finalizeIncrement, &
    results_addGroup, &
    results_openGroup, &
    results_closeGroup, &
    results_writeDataset, &
    results_setLink, &
    results_addAttribute, &
    results_removeLink, &
    results_mapping_phase, &
    results_mapping_homogenization
contains

subroutine results_init(restart)

  logical, intent(in) :: restart

  character(len=pPathLen) :: commandLine

  print'(/,a)', ' <<<+-  results init  -+>>>'; flush(IO_STDOUT)

  print*, 'M. Diehl et al., Integrating Materials and Manufacturing Innovation 6(1):83–91, 2017'
  print*, 'https://doi.org/10.1007/s40192-017-0084-5'//IO_EOL

  if(.not. restart) then
    resultsFile = HDF5_openFile(getSolverJobName()//'.hdf5','w')
    call results_addAttribute('DADF5_version_major',0)
    call results_addAttribute('DADF5_version_minor',13)
    call get_command_argument(0,commandLine)
    call results_addAttribute('creator',trim(commandLine)//' '//DAMASKVERSION)
    call results_addAttribute('created',now())
    call get_command(commandLine)
    call results_addAttribute('call',trim(commandLine))
    call results_closeGroup(results_addGroup('cell_to'))
    call results_addAttribute('description','mappings to place data in space','cell_to')
    call results_closeJobFile
  endif

end subroutine results_init


!--------------------------------------------------------------------------------------------------
!> @brief opens the results file to append data
!--------------------------------------------------------------------------------------------------
subroutine results_openJobFile

  resultsFile = HDF5_openFile(getSolverJobName()//'.hdf5','a')

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


  write(incChar,'(i10)') inc
  call results_closeGroup(results_addGroup(trim('increment_'//trim(adjustl(incChar)))))
  call results_setLink(trim('increment_'//trim(adjustl(incChar))),'current')
  call results_addAttribute('t/s',time,trim('increment_'//trim(adjustl(incChar))))

end subroutine results_addIncrement


!--------------------------------------------------------------------------------------------------
!> @brief finalize increment
!> @details remove soft link
!--------------------------------------------------------------------------------------------------
subroutine results_finalizeIncrement

  call results_removeLink('current')

end subroutine results_finalizeIncrement


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
!> @brief close a group
!--------------------------------------------------------------------------------------------------
subroutine results_closeGroup(group_id)

  integer(HID_T), intent(in) :: group_id


  call HDF5_closeGroup(group_id)

end subroutine results_closeGroup


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

  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  endif

end subroutine results_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief adds an integer attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  integer,          intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  endif

end subroutine results_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief adds a real attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  real(pReal),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  endif

end subroutine results_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief adds an integer array attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  integer,          intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  endif

end subroutine results_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief adds a real array attribute an object in the results file
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  real(pReal),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  endif

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
!> @brief Store real scalar dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine results_writeScalarDataset_real(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                  :: label,group,description
  character(len=*), intent(in),    optional     :: SIunit
  real(pReal),      intent(in),    dimension(:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeScalarDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store real vector dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine results_writeVectorDataset_real(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                    :: label,group,description
  character(len=*), intent(in),    optional       :: SIunit
  real(pReal),      intent(in),    dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeVectorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store real tensor dataset with associated metadata.
!> @details Data is transposed to compenstate transposed storage order.
!--------------------------------------------------------------------------------------------------
subroutine results_writeTensorDataset_real(dataset,group,label,description,SIunit,transposed)

  character(len=*), intent(in)                   :: label,group,description
  character(len=*), intent(in), optional         :: SIunit
  logical,          intent(in), optional         :: transposed
  real(pReal),      intent(in), dimension(:,:,:) :: dataset

  integer :: i
  logical :: transposed_
  integer(HID_T) :: groupHandle
  real(pReal), dimension(:,:,:), allocatable :: dataset_transposed


  if(present(transposed)) then
    transposed_ = transposed
  else
    transposed_ = .true.
  endif

  groupHandle = results_openGroup(group)
  if(transposed_) then
    if(size(dataset,1) /= size(dataset,2)) error stop 'transpose non-symmetric tensor'
    allocate(dataset_transposed,mold=dataset)
    do i=1,size(dataset_transposed,3)
      dataset_transposed(:,:,i) = transpose(dataset(:,:,i))
    enddo
    call HDF5_write(dataset_transposed,groupHandle,label)
  else
    call HDF5_write(dataset,groupHandle,label)
  endif
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeTensorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store integer vector dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine results_writeVectorDataset_int(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                 :: label,group,description
  character(len=*), intent(in), optional       :: SIunit
  integer,          intent(in), dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeVectorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief Store integer tensor dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine results_writeTensorDataset_int(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                   :: label,group,description
  character(len=*), intent(in), optional         :: SIunit
  integer,          intent(in), dimension(:,:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)


end subroutine results_writeTensorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine results_mapping_phase(ID,entry,label)

  integer,          dimension(:,:), intent(in) :: ID                                                !< phase ID at (co,ce)
  integer,          dimension(:,:), intent(in) :: entry                                             !< phase entry at (co,ce)
  character(len=*), dimension(:),   intent(in) :: label                                             !< label of each phase section

  integer, dimension(size(entry,1),size(entry,2)) :: &
    entryGlobal
  integer, dimension(size(label),0:worldsize-1) :: entryOffset                                      !< offset in entry counting per process
  integer, dimension(0:worldsize-1)             :: writeSize                                        !< amount of data written per process
  integer(HSIZE_T), dimension(2) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)

  integer(HID_T) :: &
    loc_id, &                                                                                       !< identifier of group in file
    dtype_id, &                                                                                     !< identifier of compound data type
    label_id, &                                                                                     !< identifier of label (string) in compound data type
    entry_id, &                                                                                     !< identifier of entry (integer) in compound data type
    dset_id, &
    memspace_id, &
    filespace_id, &
    plist_id, &
    dt_id

  integer(SIZE_T) :: type_size_string, type_size_int
  integer         :: hdferr, ierr, ce, co


  writeSize = 0
  writeSize(worldrank) = size(entry(1,:))                                                           ! total number of entries of this process

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

#ifndef PETSc
  entryGlobal = entry -1                                                                            ! 0-based
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)        ! get output at each process
  if(ierr /= 0) error stop 'MPI error'

  entryOffset = 0
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryOffset(ID(co,ce),worldrank) = entryOffset(ID(co,ce),worldrank) +1
    enddo
  enddo
  call MPI_allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)! get offset at each process
  if(ierr /= 0) error stop 'MPI error'
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryGlobal(co,ce) = entry(co,ce) -1 + entryOffset(ID(co,ce),worldrank)
    enddo
  enddo
#endif

  myShape = int([size(ID,1),writeSize(worldrank)], HSIZE_T)
  myOffset = int([0,sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([size(ID,1),sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call h5tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tget_size_f(dt_id, type_size_string, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_int, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(dtype_id, 'entry', type_size_string, H5T_NATIVE_INTEGER, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string, label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(entry_id, 'entry', 0_SIZE_T, H5T_NATIVE_INTEGER, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tclose_f(dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call h5screate_simple_f(2,myShape,memspace_id,hdferr,myShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5screate_simple_f(2,totalShape,filespace_id,hdferr,totalShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call h5pset_preserve_f(plist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  loc_id = results_openGroup('/cell_to')
  call h5dcreate_f(loc_id, 'phase', dtype_id, filespace_id, dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call h5pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(entry_id, hdferr)

  call executionStamp('cell_to/phase','cell ID and constituent ID to phase results')

end subroutine results_mapping_phase


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine results_mapping_homogenization(ID,entry,label)

  integer,          dimension(:), intent(in) :: ID                                                  !< homogenization ID at (ce)
  integer,          dimension(:), intent(in) :: entry                                               !< homogenization entry at (ce)
  character(len=*), dimension(:), intent(in) :: label                                               !< label of each homogenization section

  integer, dimension(size(entry,1)) :: &
    entryGlobal
  integer, dimension(size(label),0:worldsize-1) :: entryOffset                                      !< offset in entry counting per process
  integer, dimension(0:worldsize-1)             :: writeSize                                        !< amount of data written per process
  integer(HSIZE_T), dimension(1) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)

  integer(HID_T) :: &
    loc_id, &                                                                                       !< identifier of group in file
    dtype_id, &                                                                                     !< identifier of compound data type
    label_id, &                                                                                     !< identifier of label (string) in compound data type
    entry_id, &                                                                                     !< identifier of entry (integer) in compound data type
    dset_id, &
    memspace_id, &
    filespace_id, &
    plist_id, &
    dt_id

  integer(SIZE_T) :: type_size_string, type_size_int
  integer         :: hdferr, ierr, ce


  writeSize = 0
  writeSize(worldrank) = size(entry)                                                                ! total number of entries of this process

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

#ifndef PETSc
  entryGlobal = entry -1                                                                            ! 0-based
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call MPI_allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)        ! get output at each process
  if(ierr /= 0) error stop 'MPI error'

  entryOffset = 0
  do ce = 1, size(ID,1)
    entryOffset(ID(ce),worldrank) = entryOffset(ID(ce),worldrank) +1
  enddo
  call MPI_allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr) ! get offset at each process
  if(ierr /= 0) error stop 'MPI error'
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do ce = 1, size(ID,1)
    entryGlobal(ce) = entry(ce) -1 + entryOffset(ID(ce),worldrank)
  enddo
#endif

  myShape = int([writeSize(worldrank)], HSIZE_T)
  myOffset = int([sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call h5tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tget_size_f(dt_id, type_size_string, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tget_size_f(H5T_NATIVE_INTEGER, type_size_int, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(dtype_id, 'entry', type_size_string, H5T_NATIVE_INTEGER, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call h5tcreate_f(H5T_COMPOUND_F, type_size_string, label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tinsert_f(entry_id, 'entry', 0_SIZE_T, H5T_NATIVE_INTEGER, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5tclose_f(dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call h5screate_simple_f(1,myShape,memspace_id,hdferr,myShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5screate_simple_f(1,totalShape,filespace_id,hdferr,totalShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call h5pset_preserve_f(plist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  loc_id = results_openGroup('/cell_to')
  call h5dcreate_f(loc_id, 'homogenization', dtype_id, filespace_id, dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call h5dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call h5pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call h5tclose_f(entry_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call executionStamp('cell_to/homogenization','cell ID to homogenization results')

end subroutine results_mapping_homogenization


!--------------------------------------------------------------------------------------------------
!> @brief Add default information to a dataset.
!--------------------------------------------------------------------------------------------------
subroutine executionStamp(path,description,SIunit)

  character(len=*), intent(in)           :: path,description
  character(len=*), intent(in), optional :: SIunit


  if (HDF5_objectExists(resultsFile,path)) &
    call HDF5_addAttribute(resultsFile,'creator','DAMASK '//DAMASKVERSION,path)
  if (HDF5_objectExists(resultsFile,path)) &
    call HDF5_addAttribute(resultsFile,'created',now(),path)
  if (HDF5_objectExists(resultsFile,path)) &
    call HDF5_addAttribute(resultsFile,'description',description,path)
  if (HDF5_objectExists(resultsFile,path) .and. present(SIunit)) &
    call HDF5_addAttribute(resultsFile,'unit',SIunit,path)

end subroutine executionStamp


!--------------------------------------------------------------------------------------------------
!> @brief Return current date and time (including time zone information).
!--------------------------------------------------------------------------------------------------
character(len=24) function now()

  character(len=5)      :: zone
  integer, dimension(8) :: values


  call date_and_time(values=values,zone=zone)
  write(now,'(i4.4,5(a,i2.2),a)') &
    values(1),'-',values(2),'-',values(3),' ',values(5),':',values(6),':',values(7),zone

end function now


end module results
