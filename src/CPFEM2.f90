!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief needs a good name and description
!--------------------------------------------------------------------------------------------------
module CPFEM2

 implicit none
 private

 public :: &
   CPFEM_age, &
   CPFEM_initAll, &
   CPFEM_results
contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll()
 use prec, only: &
   pInt, &
   prec_init
 use numerics, only: &
   numerics_init
 use debug, only: &
   debug_init
 use config, only: &
   config_init
 use FEsolving, only: &
   FE_init
 use math, only: &
   math_init
 use mesh, only: &
   mesh_init
 use material, only: &
   material_init
 use HDF5_utilities, only: &
   HDF5_utilities_init
 use results, only: &
   results_init
 use lattice, only: &
   lattice_init
 use constitutive, only: &
   constitutive_init
 use crystallite, only: &
   crystallite_init
 use homogenization, only: &
   homogenization_init, &
   materialpoint_postResults
 use IO, only: &
   IO_init
 use DAMASK_interface
#ifdef FEM
 use FEM_Zoo, only: &
   FEM_Zoo_init
#endif

 implicit none

 call DAMASK_interface_init                                                                         ! Spectral and FEM interface to commandline
 call prec_init
 call IO_init
#ifdef FEM
 call FEM_Zoo_init
#endif
 call numerics_init
 call debug_init
 call config_init
 call math_init
 call FE_init
 call mesh_init
 call lattice_init
 call material_init
 call HDF5_utilities_init
 call results_init
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
 use prec, only: &
   pInt, pReal
 use IO, only: &
   IO_error
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive
 use FEsolving, only: &
   restartRead
 use material, only: &
   material_phase, &
   homogState, &
   phase_plasticity, &
   plasticState
 use config, only: &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Lp0, &
   crystallite_Fi0, &
   crystallite_Li0, &
   crystallite_S0
 use hdf5
 use HDF5_utilities, only: &
   HDF5_openFile, &
   HDF5_closeFile, &
   HDF5_openGroup, &
   HDF5_closeGroup, &
   HDF5_read
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer(pInt) :: ph,homog
 character(len=1024) :: rankStr, PlasticItem, HomogItem
 integer(HID_T) :: fileHandle, groupPlasticID, groupHomogID

 write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
 flush(6)

 ! *** restore the last converged values of each essential variable from the binary file
 if (restartRead) then
   if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from hdf5 file'
     flush(6)
   endif
   
   write(rankStr,'(a1,i0)')'_',worldrank

   fileHandle = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5')
  
   call HDF5_read(fileHandle,material_phase,      'recordedPhase')
   call HDF5_read(fileHandle,crystallite_F0,      'convergedF')
   call HDF5_read(fileHandle,crystallite_Fp0,     'convergedFp')
   call HDF5_read(fileHandle,crystallite_Fi0,     'convergedFi')
   call HDF5_read(fileHandle,crystallite_Lp0,     'convergedLp')
   call HDF5_read(fileHandle,crystallite_Li0,     'convergedLi')
   call HDF5_read(fileHandle,crystallite_S0,      'convergedS')
   
   groupPlasticID = HDF5_openGroup(fileHandle,'PlasticPhases')
   do ph = 1_pInt,size(phase_plasticity)
    write(PlasticItem,*) ph,'_'
    call HDF5_read(groupPlasticID,plasticState(ph)%state0,trim(PlasticItem)//'convergedStateConst')
   enddo
   call HDF5_closeGroup(groupPlasticID)
   
   groupHomogID = HDF5_openGroup(fileHandle,'HomogStates')
   do homog = 1_pInt, material_Nhomogenization
    write(HomogItem,*) homog,'_'
    call HDF5_read(groupHomogID,homogState(homog)%state0, trim(HomogItem)//'convergedStateHomog')
   enddo
   call HDF5_closeGroup(groupHomogID)

   call HDF5_closeFile(fileHandle)
   
   restartRead = .false.
 endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief forwards data after successful increment
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_age()
 use prec, only: &
   pReal, &
   pInt
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective
 use FEsolving, only: &
   restartWrite
 use material, only: &
   plasticState, &
   sourceState, &
   homogState, &
   thermalState, &
   damageState, &
   material_phase, &
   phase_plasticity, &
   phase_Nsources
 use config, only: &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_partionedF,&
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Fi0, &
   crystallite_Fi, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_Li0, &
   crystallite_Li, &
   crystallite_S0, &
   crystallite_S
 use HDF5_utilities, only: &
   HDF5_openFile, &
   HDF5_closeFile, &
   HDF5_addGroup, &
   HDF5_closeGroup, &
   HDF5_write
 use hdf5
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer(pInt)    ::  i, ph, homog, mySource
 character(len=32) :: rankStr, PlasticItem, HomogItem
 integer(HID_T) :: fileHandle, groupPlastic, groupHomog

 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
   write(6,'(a)') '<< CPFEM >> aging states'

 crystallite_F0  = crystallite_partionedF
 crystallite_Fp0 = crystallite_Fp
 crystallite_Lp0 = crystallite_Lp
 crystallite_Fi0 = crystallite_Fi
 crystallite_Li0 = crystallite_Li
 crystallite_S0 = crystallite_S
 
 forall (i = 1:size(plasticState)) plasticState(i)%state0 = plasticState(i)%state            ! copy state in this lengthy way because: A component cannot be an array if the encompassing structure is an array
 
 do i = 1, size(sourceState)
   do mySource = 1,phase_Nsources(i)
     sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state                    ! copy state in this lengthy way because: A component cannot be an array if the encompassing structure is an array
 enddo; enddo
 
 do homog = 1_pInt, material_Nhomogenization
   homogState       (homog)%state0 =  homogState       (homog)%state
   thermalState     (homog)%state0 =  thermalState     (homog)%state
   damageState      (homog)%state0 =  damageState      (homog)%state
 enddo

 if (restartWrite) then
   if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
     write(6,'(a)') '<< CPFEM >> writing restart variables of last converged step to hdf5 file'
   
   write(rankStr,'(a1,i0)')'_',worldrank
   fileHandle = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5','w')
   
   call HDF5_write(fileHandle,material_phase,      'recordedPhase')
   call HDF5_write(fileHandle,crystallite_F0,      'convergedF')
   call HDF5_write(fileHandle,crystallite_Fp0,     'convergedFp')
   call HDF5_write(fileHandle,crystallite_Fi0,     'convergedFi')
   call HDF5_write(fileHandle,crystallite_Lp0,     'convergedLp')
   call HDF5_write(fileHandle,crystallite_Li0,     'convergedLi')
   call HDF5_write(fileHandle,crystallite_S0,      'convergedS')
   
   groupPlastic = HDF5_addGroup(fileHandle,'PlasticPhases')
   do ph = 1_pInt,size(phase_plasticity)
     write(PlasticItem,*) ph,'_'
     call HDF5_write(groupPlastic,plasticState(ph)%state0,trim(PlasticItem)//'convergedStateConst')
   enddo
   call HDF5_closeGroup(groupPlastic)

   groupHomog = HDF5_addGroup(fileHandle,'HomogStates')
   do homog = 1_pInt, material_Nhomogenization
     write(HomogItem,*) homog,'_'
     call HDF5_write(groupHomog,homogState(homog)%state0,trim(HomogItem)//'convergedStateHomog')
   enddo
   call HDF5_closeGroup(groupHomog)
   
   call HDF5_closeFile(fileHandle)
   restartWrite = .false.
 endif

 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
   write(6,'(a)') '<< CPFEM >> done aging states'

end subroutine CPFEM_age

!--------------------------------------------------------------------------------------------------
!> @brief triggers writing of the results
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)
 use prec, only: &
   pInt
 use results
 use HDF5_utilities
 use constitutive, only: &
   constitutive_results

 implicit none
 integer(pInt), intent(in) :: inc
 real(pReal),   intent(in) :: time

 call results_openJobFile
 call results_addIncrement(inc,time)
 call constitutive_results()
 call results_removeLink('current') ! ToDo: put this into closeJobFile
 call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM2
