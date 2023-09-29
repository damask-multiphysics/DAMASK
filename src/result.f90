!--------------------------------------------------------------------------------------------------
!> @author Vitesh Shah, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi-Chin Yang, Max-Planck-Institut für Eisenforschung GmbH
!> @author Jennifer Nastola, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
module result
  use prec
  use misc
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

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

  integer(HID_T) :: resultFile

  interface result_writeDataset
    module procedure result_writeTensorDataset_real
    module procedure result_writeVectorDataset_real
    module procedure result_writeScalarDataset_real

    module procedure result_writeTensorDataset_int
    module procedure result_writeVectorDataset_int
  end interface result_writeDataset

  interface result_addAttribute
    module procedure result_addAttribute_str
    module procedure result_addAttribute_int
    module procedure result_addAttribute_real

    module procedure result_addAttribute_str_array
    module procedure result_addAttribute_int_array
    module procedure result_addAttribute_real_array
  end interface result_addAttribute

  public :: &
    result_init, &
    result_openJobFile, &
    result_closeJobFile, &
    result_addIncrement, &
    result_finalizeIncrement, &
    result_addGroup, &
    result_openGroup, &
    result_closeGroup, &
    result_writeDataset, &
    result_writeDataset_str, &
    result_setLink, &
    result_addSetupFile, &
    result_addAttribute, &
    result_removeLink, &
    result_mapping_phase, &
    result_mapping_homogenization
contains

subroutine result_init(restart)

  logical, intent(in) :: restart

  character(len=pPathLen) :: commandLine
  integer :: hdferr
  character(len=:), allocatable :: date


  print'(/,1x,a)', '<<<+-  result init  -+>>>'; flush(IO_STDOUT)

  print'(/,1x,a)', 'M. Diehl et al., Integrating Materials and Manufacturing Innovation 6(1):83–91, 2017'
  print'(  1x,a)', 'https://doi.org/10.1007/s40192-017-0084-5'

  if (.not. restart) then
    resultFile = HDF5_openFile(getSolverJobName()//'.hdf5','w')
    call result_addAttribute('DADF5_version_major',1)
    call result_addAttribute('DADF5_version_minor',0)
    call get_command_argument(0,commandLine)
    call result_addAttribute('creator',trim(commandLine)//' '//DAMASKVERSION)
    call result_addAttribute('created',now())
    call get_command(commandLine)
    call result_addAttribute('call',trim(commandLine))
    call result_closeGroup(result_addGroup('cell_to'))
    call result_addAttribute('description','mappings to place data in space','cell_to')
    call result_closeGroup(result_addGroup('setup'))
    call result_addAttribute('description','input data used to run the simulation','setup')
  else
    date = now()
    call result_openJobFile()
    call get_command(commandLine)
    call result_addAttribute('call (restart at '//date//')',trim(commandLine))
    call H5Gmove_f(resultFile,'setup','tmp',hdferr)
    call result_addAttribute('description','input data used to run the simulation up to restart at '//date,'tmp')
    call result_closeGroup(result_addGroup('setup'))
    call result_addAttribute('description','input data used to run the simulation','setup')
    call H5Gmove_f(resultFile,'tmp','setup/previous',hdferr)
  end if

  call result_closeJobFile()

end subroutine result_init


!--------------------------------------------------------------------------------------------------
!> @brief opens the result file to append data
!--------------------------------------------------------------------------------------------------
subroutine result_openJobFile(parallel)

  logical, intent(in), optional :: parallel


  resultFile = HDF5_openFile(getSolverJobName()//'.hdf5','a',parallel)

end subroutine result_openJobFile


!--------------------------------------------------------------------------------------------------
!> @brief closes the result file
!--------------------------------------------------------------------------------------------------
subroutine result_closeJobFile

  call HDF5_closeFile(resultFile)

end subroutine result_closeJobFile


!--------------------------------------------------------------------------------------------------
!> @brief creates the group of increment and adds time as attribute to the file
!--------------------------------------------------------------------------------------------------
subroutine result_addIncrement(inc,time)

  integer,       intent(in) :: inc
  real(pREAL),   intent(in) :: time

  character(len=pSTRLEN) :: incChar


  write(incChar,'(i10)') inc
  call result_closeGroup(result_addGroup(trim('increment_'//trim(adjustl(incChar)))))
  call result_setLink(trim('increment_'//trim(adjustl(incChar))),'current')
  call result_addAttribute('t/s',time,trim('increment_'//trim(adjustl(incChar))))

end subroutine result_addIncrement


!--------------------------------------------------------------------------------------------------
!> @brief finalize increment
!> @details remove soft link
!--------------------------------------------------------------------------------------------------
subroutine result_finalizeIncrement

  call result_removeLink('current')

end subroutine result_finalizeIncrement


!--------------------------------------------------------------------------------------------------
!> @brief Open a group from the result file.
!--------------------------------------------------------------------------------------------------
integer(HID_T) function result_openGroup(groupName)

  character(len=*), intent(in) :: groupName


  result_openGroup = HDF5_openGroup(resultFile,groupName)

end function result_openGroup


!--------------------------------------------------------------------------------------------------
!> @brief Add a new group to the result file.
!--------------------------------------------------------------------------------------------------
integer(HID_T) function result_addGroup(groupName)

  character(len=*), intent(in) :: groupName


  result_addGroup = HDF5_addGroup(resultFile,groupName)

end function result_addGroup


!--------------------------------------------------------------------------------------------------
!> @brief Close a group.
!--------------------------------------------------------------------------------------------------
subroutine result_closeGroup(group_id)

  integer(HID_T), intent(in) :: group_id


  call HDF5_closeGroup(group_id)

end subroutine result_closeGroup


!--------------------------------------------------------------------------------------------------
!> @brief Set link to object in result file.
!--------------------------------------------------------------------------------------------------
subroutine result_setLink(path,link)

  character(len=*), intent(in) :: path, link


  call HDF5_setLink(resultFile,path,link)

end subroutine result_setLink


!--------------------------------------------------------------------------------------------------
!> @brief Add file to setup folder and ensure unique name.
!--------------------------------------------------------------------------------------------------
subroutine result_addSetupFile(content,fname,description)

  character(len=*), intent(in) :: content, fname, description

  integer(HID_T) :: groupHandle
  character(len=:), allocatable :: name,suffix
  integer :: i

  groupHandle = result_openGroup('setup')
  name = fname(scan(fname,'/',.true.)+1:)
  suffix = ''
  i = 0

  do while (HDF5_objectExists(groupHandle,name//suffix))
      i = i+1
      suffix = '.'//IO_intAsStr(i)
  end do
  call result_writeDataset_str(content,'setup',name//suffix,description)
  call result_closeGroup(groupHandle)

end subroutine result_addSetupFile


!--------------------------------------------------------------------------------------------------
!> @brief Add a string attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_str(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel, attrValue
  character(len=*), intent(in), optional :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_str


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_int(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  integer,          intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_int


!--------------------------------------------------------------------------------------------------
!> @brief Add a real (float) attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_real(attrLabel,attrValue,path)

  character(len=*), intent(in)           :: attrLabel
  real(pREAL),      intent(in)           :: attrValue
  character(len=*), intent(in), optional :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_real


!--------------------------------------------------------------------------------------------------
!> @brief Add a string array attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_str_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  character(len=*), intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_str_array


!--------------------------------------------------------------------------------------------------
!> @brief Add an integer array attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_int_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  integer,          intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_int_array


!--------------------------------------------------------------------------------------------------
!> @brief Add a real (float) array attribute to an object in the result file.
!--------------------------------------------------------------------------------------------------
subroutine result_addAttribute_real_array(attrLabel,attrValue,path)

  character(len=*), intent(in)               :: attrLabel
  real(pREAL),      intent(in), dimension(:) :: attrValue
  character(len=*), intent(in), optional     :: path


  call HDF5_addAttribute(resultFile,attrLabel, attrValue, path)

end subroutine result_addAttribute_real_array


!--------------------------------------------------------------------------------------------------
!> @brief remove link to an object
!--------------------------------------------------------------------------------------------------
subroutine result_removeLink(link)

  character(len=*), intent(in) :: link
  integer                      :: hdferr


  call H5Ldelete_f(resultFile,link, hdferr)
  if (hdferr < 0) call IO_error(1,ext_msg = 'result_removeLink: H5Ldelete_soft_f ('//trim(link)//')')

end subroutine result_removeLink


!--------------------------------------------------------------------------------------------------
!> @brief Store string dataset.
!> @details Not collective, must be called by one process at at time.
!--------------------------------------------------------------------------------------------------
subroutine result_writeDataset_str(dataset,group,label,description)

  character(len=*), intent(in) :: label,group,description,dataset

  integer(HID_T) :: groupHandle


  groupHandle = result_openGroup(group)
  call HDF5_write_str(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description)
  call HDF5_closeGroup(groupHandle)

end subroutine result_writeDataset_str

!--------------------------------------------------------------------------------------------------
!> @brief Store real scalar dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine result_writeScalarDataset_real(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                  :: label,group,description
  character(len=*), intent(in),    optional     :: SIunit
  real(pREAL),      intent(in),    dimension(:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = result_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)

end subroutine result_writeScalarDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store real vector dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine result_writeVectorDataset_real(dataset,group,label,description,SIunit,systems)

  character(len=*), intent(in)                    :: label,group,description
  character(len=*), intent(in),    optional       :: SIunit
  character(len=*), intent(in),    dimension(:), optional :: systems
  real(pREAL),      intent(in),    dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = result_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  if (present(systems) .and. HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(resultFile,'systems',systems,group//'/'//label)
  call HDF5_closeGroup(groupHandle)

end subroutine result_writeVectorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store real tensor dataset with associated metadata.
!> @details Data is transposed to compenstate transposed storage order.
!--------------------------------------------------------------------------------------------------
subroutine result_writeTensorDataset_real(dataset,group,label,description,SIunit,transposed)

  character(len=*), intent(in)                   :: label,group,description
  character(len=*), intent(in), optional         :: SIunit
  logical,          intent(in), optional         :: transposed
  real(pREAL),      intent(in), dimension(:,:,:) :: dataset

  integer :: i
  integer(HID_T) :: groupHandle
  real(pREAL), dimension(:,:,:), allocatable :: dataset_transposed


  groupHandle = result_openGroup(group)
  if (misc_optional(transposed,.true.)) then
    if (size(dataset,1) /= size(dataset,2)) error stop 'transpose non-symmetric tensor'
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

end subroutine result_writeTensorDataset_real


!--------------------------------------------------------------------------------------------------
!> @brief Store integer vector dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine result_writeVectorDataset_int(dataset,group,label,description,SIunit,systems)

  character(len=*), intent(in)                 :: label,group,description
  character(len=*), intent(in), optional       :: SIunit
  character(len=*), intent(in),    dimension(:), optional :: systems
  integer,          intent(in), dimension(:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = result_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  if (present(systems) .and. HDF5_objectExists(groupHandle,label)) &
    call HDF5_addAttribute(resultFile,'systems',systems,group//'/'//label)
  call HDF5_closeGroup(groupHandle)

end subroutine result_writeVectorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief Store integer tensor dataset with associated metadata.
!--------------------------------------------------------------------------------------------------
subroutine result_writeTensorDataset_int(dataset,group,label,description,SIunit)

  character(len=*), intent(in)                   :: label,group,description
  character(len=*), intent(in), optional         :: SIunit
  integer,          intent(in), dimension(:,:,:) :: dataset

  integer(HID_T) :: groupHandle


  groupHandle = result_openGroup(group)
  call HDF5_write(dataset,groupHandle,label)
  call executionStamp(group//'/'//label,description,SIunit)
  call HDF5_closeGroup(groupHandle)


end subroutine result_writeTensorDataset_int


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine result_mapping_phase(ID,entry,label)

  integer,          dimension(:,:), intent(in) :: ID                                                !< phase ID at (co,ce)
  integer,          dimension(:,:), intent(in) :: entry                                             !< phase entry at (co,ce)
  character(len=*), dimension(:),   intent(in) :: label                                             !< label of each phase section

  integer(pI64), dimension(size(entry,1),size(entry,2)) :: entryGlobal
  integer(pI64), dimension(size(label),0:worldsize-1) :: entryOffset                                !< offset in entry counting per process
  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: writeSize                                  !< amount of data written per process
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

  integer(SIZE_T) :: type_size_str, type_size_int
  integer         :: hdferr, ce, co
  integer(MPI_INTEGER_KIND) :: err_MPI


  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  call HDF5_chkerr(hdferr)

#ifndef PETSC
  entryGlobal = int(entry-1,pI64)                                                                   ! 0-based
  writeSize(0) = size(entry,dim=2,kind=MPI_INTEGER_KIND)                                                         ! total number of entries of this process
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  call HDF5_chkerr(hdferr)
  call MPI_Allgather(size(entry,dim=2,kind=MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     writeSize,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  entryOffset = 0_pI64
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryOffset(ID(co,ce),worldrank) = entryOffset(ID(co,ce),worldrank) +1_pI64
    end do
  end do
  call MPI_Allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,err_MPI)! get offset at each process
  call parallelization_chkerr(err_MPI)
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do co = 1, size(ID,1)
    do ce = 1, size(ID,2)
      entryGlobal(co,ce) = int(entry(co,ce),pI64) -1_pI64 + entryOffset(ID(co,ce),worldrank)
    end do
  end do
#endif

  myShape = int([size(ID,1,MPI_INTEGER_KIND),writeSize(worldrank)], HSIZE_T)
  myOffset = int([0_MPI_INTEGER_KIND,sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([size(ID,1,MPI_INTEGER_KIND),sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call H5Tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tget_size_f(dt_id, type_size_str, hdferr)
  call HDF5_chkerr(hdferr)

  pI64_t = h5kind_to_type(kind(entryGlobal),H5_INTEGER_KIND)
  call H5Tget_size_f(pI64_t, type_size_int, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_str + type_size_int, dtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(dtype_id, 'entry', type_size_str, pI64_t, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call H5Tcreate_f(H5T_COMPOUND_F, type_size_str, label_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(entry_id, 'entry', 0_SIZE_T, pI64_t, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tclose_f(dt_id, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call H5Screate_simple_f(2,myShape,memspace_id,hdferr,myShape)
  call HDF5_chkerr(hdferr)

  call H5Screate_simple_f(2,totalShape,filespace_id,hdferr,totalShape)
  call HDF5_chkerr(hdferr)

  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call H5Pset_preserve_f(plist_id, .true., hdferr)
  call HDF5_chkerr(hdferr)

  loc_id = result_openGroup('/cell_to')
  call H5Dcreate_f(loc_id, 'phase', dtype_id, filespace_id, dset_id, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  call HDF5_chkerr(hdferr)
  call H5Dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call H5Pclose_f(plist_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(filespace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(memspace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Dclose_f(dset_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(dtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(label_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(entry_id, hdferr)

  call executionStamp('cell_to/phase','cell ID and constituent ID to phase results')

end subroutine result_mapping_phase


!--------------------------------------------------------------------------------------------------
!> @brief adds the unique mapping from spatial position and constituent ID to results
!--------------------------------------------------------------------------------------------------
subroutine result_mapping_homogenization(ID,entry,label)

  integer,          dimension(:), intent(in) :: ID                                                  !< homogenization ID at (ce)
  integer,          dimension(:), intent(in) :: entry                                               !< homogenization entry at (ce)
  character(len=*), dimension(:), intent(in) :: label                                               !< label of each homogenization section

  integer(pI64), dimension(size(entry,1)) :: entryGlobal
  integer(pI64), dimension(size(label),0:worldsize-1) :: entryOffset                                !< offset in entry counting per process
  integer(MPI_INTEGER_KIND), dimension(0:worldsize-1) :: writeSize                                  !< amount of data written per process
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

  integer(SIZE_T) :: type_size_str, type_size_int
  integer         :: hdferr, ce
  integer(MPI_INTEGER_KIND) :: err_MPI


  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  call HDF5_chkerr(hdferr)

#ifndef PETSC
  entryGlobal = int(entry-1,pI64)
  writeSize(0) = size(entry,kind=MPI_INTEGER_KIND)                                                         ! total number of entries of this process                                                            ! 0-based
#else
!--------------------------------------------------------------------------------------------------
! MPI settings and communication
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  call HDF5_chkerr(hdferr)
  call MPI_Allgather(size(entry,kind=MPI_INTEGER_KIND),1_MPI_INTEGER_KIND,MPI_INTEGER,&
                     writeSize,1_MPI_INTEGER_KIND,MPI_INTEGER,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  entryOffset = 0_pI64
  do ce = 1, size(ID)
    entryOffset(ID(ce),worldrank) = entryOffset(ID(ce),worldrank) +1_pI64
  end do
  call MPI_Allreduce(MPI_IN_PLACE,entryOffset,size(entryOffset),MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,err_MPI)! get offset at each process
  call parallelization_chkerr(err_MPI)
  entryOffset(:,worldrank) = sum(entryOffset(:,0:worldrank-1),2)
  do ce = 1, size(ID)
    entryGlobal(ce) = int(entry(ce),pI64) -1_pI64 + entryOffset(ID(ce),worldrank)
  end do
#endif

  myShape = int([writeSize(worldrank)], HSIZE_T)
  myOffset = int([sum(writeSize(0:worldrank-1))], HSIZE_T)
  totalShape = int([sum(writeSize)], HSIZE_T)

!---------------------------------------------------------------------------------------------------
! compound type: label(ID) + entry
  call H5Tcopy_f(H5T_NATIVE_CHARACTER, dt_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tset_size_f(dt_id, int(len(label(1)),SIZE_T), hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tget_size_f(dt_id, type_size_str, hdferr)
  call HDF5_chkerr(hdferr)

  pI64_t = h5kind_to_type(kind(entryGlobal),H5_INTEGER_KIND)
  call H5Tget_size_f(pI64_t, type_size_int, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_str + type_size_int, dtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(dtype_id, 'label', 0_SIZE_T, dt_id,hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(dtype_id, 'entry', type_size_str, pI64_t, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! create memory types for each component of the compound type
  call H5Tcreate_f(H5T_COMPOUND_F, type_size_str, label_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(label_id, 'label', 0_SIZE_T, dt_id, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tcreate_f(H5T_COMPOUND_F, type_size_int, entry_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tinsert_f(entry_id, 'entry', 0_SIZE_T, pI64_t, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Tclose_f(dt_id, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! create dataspace in memory (local shape = hyperslab) and in file (global shape)
  call H5Screate_simple_f(1,myShape,memspace_id,hdferr,myShape)
  call HDF5_chkerr(hdferr)

  call H5Screate_simple_f(1,totalShape,filespace_id,hdferr,totalShape)
  call HDF5_chkerr(hdferr)

  call H5Sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, myOffset, myShape, hdferr)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! write the components of the compound type individually
  call H5Pset_preserve_f(plist_id, .true., hdferr)
  call HDF5_chkerr(hdferr)

  loc_id = result_openGroup('/cell_to')
  call H5Dcreate_f(loc_id, 'homogenization', dtype_id, filespace_id, dset_id, hdferr)
  call HDF5_chkerr(hdferr)

  call H5Dwrite_f(dset_id, label_id, reshape(label(pack(ID,.true.)),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  call HDF5_chkerr(hdferr)
  call H5Dwrite_f(dset_id, entry_id, reshape(pack(entryGlobal,.true.),myShape), &
                  myShape, hdferr, file_space_id = filespace_id, mem_space_id = memspace_id, xfer_prp = plist_id)
  call HDF5_chkerr(hdferr)

!--------------------------------------------------------------------------------------------------
! close all
  call HDF5_closeGroup(loc_id)
  call H5Pclose_f(plist_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(filespace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Sclose_f(memspace_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Dclose_f(dset_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(dtype_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(label_id, hdferr)
  call HDF5_chkerr(hdferr)
  call H5Tclose_f(entry_id, hdferr)
  call HDF5_chkerr(hdferr)

  call executionStamp('cell_to/homogenization','cell ID to homogenization results')

end subroutine result_mapping_homogenization


!--------------------------------------------------------------------------------------------------
!> @brief Add default information to a dataset.
!--------------------------------------------------------------------------------------------------
subroutine executionStamp(path,description,SIunit)

  character(len=*), intent(in)           :: path,description
  character(len=*), intent(in), optional :: SIunit


  if (HDF5_objectExists(resultFile,path)) &
    call HDF5_addAttribute(resultFile,'creator','DAMASK '//DAMASKVERSION,path)
  if (HDF5_objectExists(resultFile,path)) &
    call HDF5_addAttribute(resultFile,'created',now(),path)
  if (HDF5_objectExists(resultFile,path)) &
    call HDF5_addAttribute(resultFile,'description',description,path)
  if (HDF5_objectExists(resultFile,path) .and. present(SIunit)) &
    call HDF5_addAttribute(resultFile,'unit',SIunit,path)

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


end module result
