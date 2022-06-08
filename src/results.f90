!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module results
  use prec
  use parallelization
  use IO
  use HDF5_utilities
  use HDF5
#ifdef PETSC
  use CLI
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif
#else
  use DAMASK_interface
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
    module procedure results_addAttribute_str
    module procedure results_addAttribute_int
    module procedure results_addAttribute_real

    module procedure results_addAttribute_str_array
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
    results_writeDataset_str, &
    results_setLink, &
    results_addAttribute, &
    results_removeLink, &
    results_mapping_phase, &
    results_mapping_homogenization
contains

subroutine results_init(restart)

  logical, intent(in) :: restart

  character(len=pPathLen) :: commandLine
  integer :: hdferr
  character(len=:), allocatable :: date


  print'(/,1x,a)', '<<<+-  results init  -+>>>'; flush(IO_STDOUT)

  print'(/,1x,a)', 'M. Diehl et al., Integrating Materials and Manufacturing Innovation 6(1):83–91, 2017'
  print'(  1x,a)', 'https://doi.org/10.1007/s40192-017-0084-5'

  if (.not. restart) then
    resultsFile = HDF5_openFile(getSolverJobName()//'.hdf5','w')
    call results_addAttribute('DADF5_version_major',0)
    call results_addAttribute('DADF5_version_minor',14)
    call get_command_argument(0,commandLine)
    call results_addAttribute('creator',trim(commandLine)//' '//DAMASKVERSION)
    call results_addAttribute('created',now())
    call get_command(commandLine)
    call results_addAttribute('call',trim(commandLine))
    call results_closeGroup(results_addGroup('cell_to'))
    call results_addAttribute('description','mappings to place data in space','cell_to')
    call results_closeGroup(results_addGroup('setup'))
    call results_addAttribute('description','input data used to run the simulation','setup')
  else
    date = now()
    call results_openJobFile
    call get_command(commandLine)
    call results_addAttribute('call (restart at '//date//')',trim(commandLine))
    call H5Gmove_f(resultsFile,'setup','tmp',hdferr)
    call results_addAttribute('description','input data used to run the simulation up to restart at '//date,'tmp')
    call results_closeGroup(results_addGroup('setup'))
    call results_addAttribute('description','input data used to run the simulation','setup')
    call H5Gmove_f(resultsFile,'tmp','setup/previous',hdferr)
  end if

  call results_closeJobFile

end subroutine results_init


!--------------------------------------------------------------------------------------------------
!> @brief opens the results file to append data
!--------------------------------------------------------------------------------------------------
subroutine results_openJobFile(parallel)

  logical, intent(in), optional :: parallel


  resultsFile = HDF5_openFile(getSolverJobName()//'.hdf5','a',parallel)

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
!> @brief Add a string attribute to an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_str(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer attribute an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  integer,          intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief Add a real attribute an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  real(pReal),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief Add a string array attribute an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_str_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  character(len=*), intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_str_array


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer array attribute an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_int_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  integer,          intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief Add a real array attribute an object in the results file.
!--------------------------------------------------------------------------------------------------
subroutine results_addAttribute_real_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  real(pReal),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  if (present(path)) then
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue, path)
  else
    call HDF5_addAttribute(resultsFile,attrLabel, attrValue)
  end if

end subroutine results_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief remove link to an object
!--------------------------------------------------------------------------------------------------
subroutine results_removeLink(link)

  character(len=*), intent(in) :: link
  integer                      :: hdferr


  call H5Ldelete_f(resultsFile,link, hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg = 'results_removeLink: H5Ldelete_soft_f ('//trim(link)//')')

end subroutine results_removeLink


!--------------------------------------------------------------------------------------------------
!> @brief Store string dataset.
!> @details Not collective, must be called by one process at at time.
!--------------------------------------------------------------------------------------------------
subroutine results_writeDataset_str(dataset,group,label,description)

  character(len=*), intent(in) :: label,group,description,dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write_str(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeDataset_str

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
subroutine results_writeVectorDataset_real(dataset,group,label,description,SIunit,systems)

  character(len=*), intent(in)                    :: label,group,description
  character(len=*), intent(in),    optional       :: SIunit
  character(len=*), intent(in),    dimension(:), optional :: systems
  real(pReal),      intent(in),    dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  if (present(systems) .and. HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(resultsFile,'systems',systems,group//'/'//label)
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
  end if

  groupHandle = results_openGroup(group)
  if(transposed_) then
    if(size(dataset,1) /= size(dataset,2)) error stop 'transpose non-symmetric tensor'
    allocate(dataset_transposed,mold=dataset)
    do i=1,size(dataset_transposed,3)
      dataset_transposed(:,:,i) = transpose(dataset(:,:,i))
    end do
    call HDF5_write(dataset_transposed,groupHandle,label)
  else
    call HDF5_write(dataset,groupHandle,label)
  end if
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine results_writeTensorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store integer vector dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine results_writeVectorDataset_int(dataset,group,label,description,SIunit,systems)

  character(len=*), intent(in)                 :: label,group,description
  character(len=*), intent(in), optional       :: SIunit
  character(len=*), intent(in),    dimension(:), optional :: systems
  integer,          intent(in), dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = results_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  if (present(systems) .and. HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(resultsFile,'systems',systems,group//'/'//label)
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

  integer(pI64), dimension(size(entry,1),size(entry,2)) :: &
    entryGlobal
  integer(pI64), dimension(size(label),0:worldsize-1) :: entryOffset                                !< offset in entry counting per process
  integer, dimension(0:worldsize-1) :: writeSize                                                    !< amount of data written per process
  integer(HSIZE_T), dimension(2) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)

  integer(HID_T) :: &
    pI64_t, &                                                                                       !< HDF5 type for pI64 (8 bit integer)
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
  integer         :: hdferr, ce, co
  integer(MPI_INTEGER_KIND) :: err_MPI


  writeSize = 0
  writeSize(worldrank) = size(entry(1,:))                                                           ! total number of entries of this process

  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

#ifndef PETSC
  entryGlobal = int(entry -1,pI64)                                                                  ! 0-based
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call MPI_Allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err_MPI)   ! get output at each process
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  entryOffset = 0_pI64
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryOffset(ID(co,ce),worldrank) = entryOffset(ID(co,ce),worldrank) +1_pI64
    end do
  end do
  call MPI_Allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,err_MPI)! get offset at each process
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryGlobal(co,ce) = int(entry(co,ce),pI64) -1_pI64 + entryOffset(ID(co,ce),worldrank)
    end do
  end do
#endif

  myShape = int([size(ID,1),writeSize(worldrank)], HSIZE_T)
  myOffset = int([0,sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([size(ID,1),sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call H5Tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tget_size_f(dt_id, type_size_string, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  pI64_t = h5kind_to_type(kind(entryGlobal),H5_INTEGER_KIND)
  call H5Tget_size_f(pI64_t, type_size_int, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(dtype_id, 'entry', type_size_string, pI64_t, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call H5Tcreate_f(H5T_COMPOUND_F, type_size_string, label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(entry_id, 'entry', 0_SIZE_T, pI64_t, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tclose_f(dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call H5Screate_simple_f(2,myShape,memspace_id,hdferr,myShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Screate_simple_f(2,totalShape,filespace_id,hdferr,totalShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call H5Pset_preserve_f(plist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  loc_id = results_openGroup('/cell_to')
  call H5Dcreate_f(loc_id, 'phase', dtype_id, filespace_id, dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call H5Pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(entry_id, hdferr)

  call executionStamp('cell_to/phase','cell ID and constituent ID to phase results')

end subroutine results_mapping_phase


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine results_mapping_homogenization(ID,entry,label)

  integer,          dimension(:), intent(in) :: ID                                                  !< homogenization ID at (ce)
  integer,          dimension(:), intent(in) :: entry                                               !< homogenization entry at (ce)
  character(len=*), dimension(:), intent(in) :: label                                               !< label of each homogenization section

  integer(pI64), dimension(size(entry,1)) :: &
    entryGlobal
  integer(pI64), dimension(size(label),0:worldsize-1) :: entryOffset                                !< offset in entry counting per process
  integer, dimension(0:worldsize-1) :: writeSize                                                    !< amount of data written per process
  integer(HSIZE_T), dimension(1) :: &
    myShape, &                                                                                      !< shape of the dataset (this process)
    myOffset, &
    totalShape                                                                                      !< shape of the dataset (all processes)

  integer(HID_T) :: &
    pI64_t, &                                                                                       !< HDF5 type for pI64 (8 bit integer)
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
  integer         :: hdferr, ce
  integer(MPI_INTEGER_KIND) :: err_MPI


  writeSize = 0
  writeSize(worldrank) = size(entry)                                                                ! total number of entries of this process

  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

#ifndef PETSC
  entryGlobal = int(entry -1,pI64)                                                                  ! 0-based
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call MPI_Allreduce(MPI_IN_PLACE,writeSize,worldsize,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err_MPI)   ! get output at each process
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  entryOffset = 0_pI64
  do ce = 1, size(ID,1)
    entryOffset(ID(ce),worldrank) = entryOffset(ID(ce),worldrank) +1_pI64
  end do
  call MPI_Allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,err_MPI)! get offset at each process
  if(err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do ce = 1, size(ID,1)
    entryGlobal(ce) = int(entry(ce),pI64) -1_pI64 + entryOffset(ID(ce),worldrank)
  end do
#endif

  myShape = int([writeSize(worldrank)], HSIZE_T)
  myOffset = int([sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call H5Tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tget_size_f(dt_id, type_size_string, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  pI64_t = h5kind_to_type(kind(entryGlobal),H5_INTEGER_KIND)
  call H5Tget_size_f(pI64_t, type_size_int, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_string + type_size_int, dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(dtype_id, 'entry', type_size_string, pI64_t, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call H5Tcreate_f(H5T_COMPOUND_F, type_size_string, label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tinsert_f(entry_id, 'entry', 0_SIZE_T, pI64_t, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Tclose_f(dt_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call H5Screate_simple_f(1,myShape,memspace_id,hdferr,myShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Screate_simple_f(1,totalShape,filespace_id,hdferr,totalShape)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call H5Pset_preserve_f(plist_id, .true., hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  loc_id = results_openGroup('/cell_to')
  call H5Dcreate_f(loc_id, 'homogenization', dtype_id, filespace_id, dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'

  call H5Dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  if(hdferr < 0) error stop 'HDF5 error'

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call H5Pclose_f(plist_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Sclose_f(filespace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Sclose_f(memspace_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Dclose_f(dset_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(dtype_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(label_id, hdferr)
  if(hdferr < 0) error stop 'HDF5 error'
  call H5Tclose_f(entry_id, hdferr)
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
