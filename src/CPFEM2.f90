!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module CPFEM2
  use prec
  use numerics
  use debug
  use config
  use FEsolving
  use math
  use rotations
  use material
  use lattice
  use IO
  use DAMASK_interface
  use results
  use discretization
  use HDF5
  use HDF5_utilities
  use homogenization
  use constitutive
  use crystallite
#ifdef FEM
  use FEM_Zoo
  use mesh
#else
  use mesh_grid
#endif

  implicit none
  private

  public :: &
    CPFEM_forward, &
    CPFEM_initAll, &
    CPFEM_results, &
    CPFEM_restartWrite

contains


!--------------------------------------------------------------------------------------------------
!> @brief call all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll

  call DAMASK_interface_init                                                                        ! Spectral and FEM interface to commandline
  call prec_init
  call IO_init
#ifdef FEM
  call FEM_Zoo_init
#endif
  call numerics_init
  call debug_init
  call config_init
  call math_init
  call rotations_init
  call lattice_init
  call HDF5_utilities_init
  call results_init
  call mesh_init
  call material_init
  call constitutive_init
  call crystallite_init
  call homogenization_init
  call materialpoint_postResults
  call CPFEM_init

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init
 
  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'; flush(6)

  if (interface_restartInc > 0) then
    write(6,'(/,a,i0,a)') ' reading restart information of increment ', interface_restartInc, ' from file'

    write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
    fileHandle = HDF5_openFile(fileName)
    
    call HDF5_read(fileHandle,crystallite_F0, 'F')
    call HDF5_read(fileHandle,crystallite_Fp0,'Fp')
    call HDF5_read(fileHandle,crystallite_Fi0,'Fi')
    call HDF5_read(fileHandle,crystallite_Lp0,'Lp')
    call HDF5_read(fileHandle,crystallite_Li0,'Li')
    call HDF5_read(fileHandle,crystallite_S0, 'S')
    
    groupHandle = HDF5_openGroup(fileHandle,'constituent')
    do i = 1,size(phase_plasticity)
      write(datasetName,'(i0,a)') i,'_omega_plastic'
      call HDF5_read(groupHandle,plasticState(i)%state0,datasetName)
    enddo
    call HDF5_closeGroup(groupHandle)

    groupHandle = HDF5_openGroup(fileHandle,'materialpoint')
    do i = 1, material_Nhomogenization
      write(datasetName,'(i0,a)') i,'_omega_homogenization'
      call HDF5_read(groupHandle,homogState(i)%state0,datasetName)
    enddo
    call HDF5_closeGroup(groupHandle)

    call HDF5_closeFile(fileHandle)
  endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_forward

  integer :: i, j

  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) &
    write(6,'(a)') '<< CPFEM >> aging states'

  crystallite_F0  = crystallite_partionedF
  crystallite_Fp0 = crystallite_Fp
  crystallite_Lp0 = crystallite_Lp
  crystallite_Fi0 = crystallite_Fi
  crystallite_Li0 = crystallite_Li
  crystallite_S0  = crystallite_S

  do i = 1, size(plasticState)
    plasticState(i)%state0 = plasticState(i)%state
  enddo 
  do i = 1, size(sourceState)
    do j = 1,phase_Nsources(i)
      sourceState(i)%p(j)%state0 = sourceState(i)%p(j)%state
  enddo; enddo
  do i = 1, material_Nhomogenization
    homogState  (i)%state0 = homogState  (i)%state
    thermalState(i)%state0 = thermalState(i)%state
    damageState (i)%state0 = damageState (i)%state
  enddo

end subroutine CPFEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write current  restart information (Field and constitutive data) to file.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_restartWrite

  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  write(6,'(a)') ' writing field and constitutive data required for restart to file';flush(6)
    
  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName,'a')
    
  call HDF5_write(fileHandle,crystallite_partionedF,'F')
  call HDF5_write(fileHandle,crystallite_Fp,        'Fp')
  call HDF5_write(fileHandle,crystallite_Fi,        'Fi')
  call HDF5_write(fileHandle,crystallite_Lp,        'Lp')
  call HDF5_write(fileHandle,crystallite_Li,        'Li')
  call HDF5_write(fileHandle,crystallite_S,         'S')
    
  groupHandle = HDF5_addGroup(fileHandle,'constituent')
  do i = 1,size(phase_plasticity)
    write(datasetName,'(i0,a)') i,'_omega_plastic'
    call HDF5_write(groupHandle,plasticState(i)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_addGroup(fileHandle,'materialpoint')
  do i = 1, material_Nhomogenization
    write(datasetName,'(i0,a)') i,'_omega_homogenization'
    call HDF5_write(groupHandle,homogState(i)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)
    
  call HDF5_closeFile(fileHandle)

end subroutine CPFEM_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)
 
  integer,     intent(in) :: inc
  real(pReal), intent(in) :: time
 
  call results_openJobFile
  call results_addIncrement(inc,time)
  call constitutive_results
  call crystallite_results
  call homogenization_results
  call discretization_results
  call results_removeLink('current') ! ToDo: put this into closeJobFile?
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM2
