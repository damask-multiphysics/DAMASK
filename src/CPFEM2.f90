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
  use HDF5
  use DAMASK_interface
  use results
  use discretization
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
    CPFEM_age, &
    CPFEM_initAll, &
    CPFEM_results, &
    CPFEM_restartWrite

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
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
  call mesh_init
  call lattice_init
  call HDF5_utilities_init
  call results_init
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
 
  integer :: ph,homog
  character(len=1024) :: rankStr, PlasticItem, HomogItem
  integer(HID_T) :: fileHandle, groupPlasticID, groupHomogID

  write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
  flush(6)

  ! *** restore the last converged values of each essential variable
  if (interface_restartInc > 0) then
    if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0) then
      write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from hdf5 file'
      flush(6)
    endif
    
    write(rankStr,'(a1,i0)')'_',worldrank

    fileHandle = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5')
   
    call HDF5_read(fileHandle,crystallite_F0, 'convergedF')
    call HDF5_read(fileHandle,crystallite_Fp0,'convergedFp')
    call HDF5_read(fileHandle,crystallite_Fi0,'convergedFi')
    call HDF5_read(fileHandle,crystallite_Lp0,'convergedLp')
    call HDF5_read(fileHandle,crystallite_Li0,'convergedLi')
    call HDF5_read(fileHandle,crystallite_S0, 'convergedS')
    
    groupPlasticID = HDF5_openGroup(fileHandle,'PlasticPhases')
    do ph = 1,size(phase_plasticity)
      write(PlasticItem,*) ph,'_'
      call HDF5_read(groupPlasticID,plasticState(ph)%state0,trim(PlasticItem)//'convergedStateConst')
    enddo
    call HDF5_closeGroup(groupPlasticID)
    
    groupHomogID = HDF5_openGroup(fileHandle,'HomogStates')
    do homog = 1, material_Nhomogenization
      write(HomogItem,*) homog,'_'
      call HDF5_read(groupHomogID,homogState(homog)%state0, trim(HomogItem)//'convergedStateHomog')
    enddo
    call HDF5_closeGroup(groupHomogID)

    call HDF5_closeFile(fileHandle)
  endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_age

  integer :: i, homog, mySource

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
    do mySource = 1,phase_Nsources(i)
      sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state
  enddo; enddo
  do homog = 1, material_Nhomogenization
    homogState       (homog)%state0 =  homogState       (homog)%state
    thermalState     (homog)%state0 =  thermalState     (homog)%state
    damageState      (homog)%state0 =  damageState      (homog)%state
  enddo

end subroutine CPFEM_age


!--------------------------------------------------------------------------------------------------
!> @brief Store DAMASK restart data.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_restartWrite

  integer           :: ph, homog
  character(len=32) :: rankStr, PlasticItem, HomogItem
  integer(HID_T)    :: fileHandle, groupPlastic, groupHomog

  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) &
    write(6,'(a)') '<< CPFEM >> writing restart variables of last converged step to hdf5 file'
    
  write(rankStr,'(a1,i0)')'_',worldrank
  fileHandle = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5','a')
    
  call HDF5_write(fileHandle,crystallite_F0,  'convergedF')
  call HDF5_write(fileHandle,crystallite_Fp0, 'convergedFp')
  call HDF5_write(fileHandle,crystallite_Fi0, 'convergedFi')
  call HDF5_write(fileHandle,crystallite_Lp0, 'convergedLp')
  call HDF5_write(fileHandle,crystallite_Li0, 'convergedLi')
  call HDF5_write(fileHandle,crystallite_S0,  'convergedS')
    
  groupPlastic = HDF5_addGroup(fileHandle,'PlasticPhases')
  do ph = 1,size(phase_plasticity)
    write(PlasticItem,*) ph,'_'
    call HDF5_write(groupPlastic,plasticState(ph)%state0,trim(PlasticItem)//'convergedStateConst')
  enddo
  call HDF5_closeGroup(groupPlastic)

  groupHomog = HDF5_addGroup(fileHandle,'HomogStates')
  do homog = 1, material_Nhomogenization
    write(HomogItem,*) homog,'_'
    call HDF5_write(groupHomog,homogState(homog)%state0,trim(HomogItem)//'convergedStateHomog')
  enddo
  call HDF5_closeGroup(groupHomog)
    
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
  call results_removeLink('current') ! ToDo: put this into closeJobFile
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM2
