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
   CPFEM_initAll

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll(el,ip)
 use prec, only: &
   pInt
 use prec, only: &
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
 integer(pInt), intent(in) ::                        el, &                                          !< FE el number
                                                     ip                                             !< FE integration point number

 call DAMASK_interface_init                                                                    ! Spectral and FEM interface to commandline
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
 call mesh_init(ip, el)                                                                        ! pass on coordinates to alter calcMode of first ip
 call lattice_init
 call material_init
 call HDF5_utilities_init
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
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pInt, pReal, pLongInt
 use IO, only: &
   IO_read_realFile,&
   IO_read_intFile, &
   IO_timeStamp, &
   IO_error
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive
 use FEsolving, only: &
   restartRead, &
   modelName
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
   crystallite_dPdF0, &
   crystallite_Tstar0_v
 use hdf5
 use HDF5_utilities, only: &
   HDF5_openFile, &
   HDF5_openGroup2, &
   HDF5_read
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none
 integer(pInt) :: k,l,m,ph,homog
 character(len=1024) :: rankStr
 integer(HID_T) :: fileReadID, groupPlasticID, groupHomogID
 integer        :: hdferr

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
   flush(6)
 endif mainProcess

 ! *** restore the last converged values of each essential variable from the binary file
 if (restartRead) then
   if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from hdf5 file'
     flush(6)
   endif
   
   write(rankStr,'(a1,i0)')'_',worldrank

   fileReadID = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5')
  
   call HDF5_read(material_phase,      fileReadID,'recordedPhase')
   write(6,*) material_phase
   call HDF5_read(crystallite_F0,      fileReadID,'convergedF')
   write(6,*) crystallite_F0
   call HDF5_read(crystallite_Fp0,     fileReadID,'convergedFp')
   call HDF5_read(crystallite_Fi0,     fileReadID,'convergedFi')
   call HDF5_read(crystallite_Lp0,     fileReadID,'convergedLp')
   call HDF5_read(crystallite_Li0,     fileReadID,'convergedLi')
   call HDF5_read(crystallite_dPdF0,   fileReadID,'convergeddPdF')
   call HDF5_read(crystallite_Tstar0_v,fileReadID,'convergedTstar')
   
   groupPlasticID = HDF5_openGroup2(fileReadID,'PlasticPhases')
   do ph = 1_pInt,size(phase_plasticity)
    call HDF5_read(plasticState(ph)%state0,groupPlasticID,'convergedStateConst')
    write(6,*) plasticState(ph)%state0
   enddo
   
   groupHomogID = HDF5_openGroup2(fileReadID,'material_Nhomogenization')
   do homog = 1_pInt, material_Nhomogenization
    call HDF5_read(homogState(homog)%state0, groupHomogID,'convergedStateHomog')
   enddo
   
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
   vacancyfluxState, &
   hydrogenfluxState, &
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
   crystallite_dPdF0, &
   crystallite_dPdF, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v
 use IO, only: &
   IO_write_jobRealFile, &
   IO_warning
 use HDF5_utilities, only: &
   HDF5_openFile, &
   HDF5_closeFile, &
   HDF5_closeGroup, &
   HDF5_addGroup2, &
   HDF5_write
 use hdf5
 use DAMASK_interface, only: &
   getSolverJobName

 implicit none

 integer(pInt) ::                                    i, k, l, m, ph, homog, mySource
 character(len=32) :: rankStr, groupItem
 integer(HID_T) :: fileHandle, groupPlastic, groupHomog
 integer        :: hdferr
 integer(HSIZE_T)  :: hdfsize
 

if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
  write(6,'(a)') '<< CPFEM >> aging states'

  crystallite_F0  = crystallite_partionedF                                                    ! crystallite deformation (_subF is perturbed...)
  crystallite_Fp0 = crystallite_Fp                                                            ! crystallite plastic deformation
  crystallite_Lp0 = crystallite_Lp                                                            ! crystallite plastic velocity
  crystallite_Fi0 = crystallite_Fi                                                            ! crystallite intermediate deformation
  crystallite_Li0 = crystallite_Li                                                            ! crystallite intermediate velocity
  crystallite_dPdF0 = crystallite_dPdF                                                        ! crystallite stiffness
  crystallite_Tstar0_v = crystallite_Tstar_v                                                  ! crystallite 2nd Piola Kirchhoff stress
  
  forall (i = 1:size(plasticState)) plasticState(i)%state0 = plasticState(i)%state            ! copy state in this lengthy way because: A component cannot be an array if the encompassing structure is an array
  
  do i = 1, size(sourceState)
    do mySource = 1,phase_Nsources(i)
      sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state                    ! copy state in this lengthy way because: A component cannot be an array if the encompassing structure is an array
  enddo; enddo
  
  do homog = 1_pInt, material_Nhomogenization
    homogState       (homog)%state0 =  homogState       (homog)%state
    thermalState     (homog)%state0 =  thermalState     (homog)%state
    damageState      (homog)%state0 =  damageState      (homog)%state
    vacancyfluxState (homog)%state0 =  vacancyfluxState (homog)%state
    hydrogenfluxState(homog)%state0 =  hydrogenfluxState(homog)%state
  enddo

if (restartWrite) then
  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
    write(6,'(a)') '<< CPFEM >> writing restart variables of last converged step to hdf5 file'
  write(rankStr,'(a1,i0)')'_',worldrank

  fileHandle = HDF5_openFile(trim(getSolverJobName())//trim(rankStr)//'.hdf5','w')
  
  call HDF5_write(material_phase,      fileHandle,'recordedPhase')
  call HDF5_write(crystallite_F0,      fileHandle,'convergedF')
  call HDF5_write(crystallite_Fp0,     fileHandle,'convergedFp')
  call HDF5_write(crystallite_Fi0,     fileHandle,'convergedFi')
  call HDF5_write(crystallite_Lp0,     fileHandle,'convergedLp')
  call HDF5_write(crystallite_Li0,     fileHandle,'convergedLi')
  call HDF5_write(crystallite_dPdF0,   fileHandle,'convergeddPdF')
  call HDF5_write(crystallite_Tstar0_v,fileHandle,'convergedTstar')
  
  groupPlastic = HDF5_addGroup2(fileHandle,'PlasticPhases')
  do ph = 1_pInt,size(phase_plasticity)
    !write(groupItem,'(a1,i0)') ph
    !write(6,*) groupItem
    call HDF5_write(plasticState(ph)%state0,groupPlastic,'convergedStateConst')
  enddo
  call HDF5_closeGroup(groupPlastic)

  groupHomog = HDF5_addGroup2(fileHandle,'material_Nhomogenization')
  do homog = 1_pInt, material_Nhomogenization
    call HDF5_write(homogState(homog)%state0,groupHomog,'convergedStateHomog')
  enddo
  call HDF5_closeGroup(groupHomog)
  
  call HDF5_closeFile(fileHandle)
  restartWrite = .false.
endif

if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
  write(6,'(a)') '<< CPFEM >> done aging states'

end subroutine CPFEM_age

end module CPFEM2
