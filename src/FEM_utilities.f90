!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_utilities
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdmda.h>
 use prec, only: pReal, pInt

use PETScdmda
use PETScis

 implicit none
 private
#include <petsc/finclude/petsc.h>
!--------------------------------------------------------------------------------------------------
! 
 logical,       public             :: cutBack = .false.                                             !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
 integer(pInt), public, parameter  :: maxFields = 6_pInt
 integer(pInt), public             :: nActiveFields = 0_pInt
 
!--------------------------------------------------------------------------------------------------
! grid related information information
 real(pReal),   public                :: wgt                                                        !< weighting factor 1/Nelems
 real(pReal),   public                :: wgtDof                                                     !< weighting factor 1/Nelems
 real(pReal),   public                :: C_volAvg(3,3,3,3)
 
!--------------------------------------------------------------------------------------------------
! output data
 PetscViewer,                    public :: resUnit
 Vec,                            public :: coordinatesVec
 Vec,               allocatable, public :: homogenizationResultsVec(:), &
                                           crystalliteResultsVec(:,:), &
                                           phaseResultsVec(:,:) 

!--------------------------------------------------------------------------------------------------
! field labels information
 character(len=*),                         parameter,            public :: &
   FIELD_MECH_label     = 'mechanical', &
   FIELD_THERMAL_label  = 'thermal', &
   FIELD_DAMAGE_label   = 'damage', &
   FIELD_SOLUTE_label   = 'solute', &
   FIELD_MGTWIN_label   = 'mgtwin'
 
 enum, bind(c)
   enumerator :: FIELD_UNDEFINED_ID, &
                 FIELD_MECH_ID, &
                 FIELD_THERMAL_ID, &
                 FIELD_DAMAGE_ID, &
                 FIELD_SOLUTE_ID, &
                 FIELD_MGTWIN_ID
 end enum
 enum, bind(c)
   enumerator :: COMPONENT_UNDEFINED_ID, &
                 COMPONENT_MECH_X_ID, &
                 COMPONENT_MECH_Y_ID, &
                 COMPONENT_MECH_Z_ID, &
                 COMPONENT_THERMAL_T_ID, &
                 COMPONENT_DAMAGE_PHI_ID, &
                 COMPONENT_SOLUTE_CV_ID, &
                 COMPONENT_SOLUTE_CVPOT_ID, &
                 COMPONENT_SOLUTE_CH_ID, &
                 COMPONENT_SOLUTE_CHPOT_ID, &
                 COMPONENT_SOLUTE_CVaH_ID, &
                 COMPONENT_SOLUTE_CVaHPOT_ID, &
                 COMPONENT_MGTWIN_PHI_ID
 end enum
 
!--------------------------------------------------------------------------------------------------
! variables controlling debugging
 logical, private :: &
   debugGeneral, &                                                                                  !< general debugging of FEM solver
   debugRotation, &                                                                                 !< also printing out results in lab frame
   debugPETSc                                                                                       !< use some in debug defined options for more verbose PETSc solution

!--------------------------------------------------------------------------------------------------
! derived types
 type, public :: tSolutionState                                                                     !< return type of solution from FEM solver variants
   logical       :: converged         = .true.   
   logical       :: stagConverged     = .true.   
   logical       :: regrid            = .false.   
   integer(pInt) :: iterationsNeeded  = 0_pInt
 end type tSolutionState

 type, public :: tComponentBC
   integer(kind(COMPONENT_UNDEFINED_ID))          :: ID
   real(pReal),                       allocatable :: Value(:)
   logical,                           allocatable :: Mask(:)   
 end type tComponentBC

 type, public :: tFieldBC
   integer(kind(FIELD_UNDEFINED_ID))              :: ID
   integer(pInt)                                  :: nComponents = 0_pInt
   type(tComponentBC),                allocatable :: componentBC(:)
 end type tFieldBC

 type, public :: tLoadCase
   real(pReal)                       :: time                   = 0.0_pReal                               !< length of increment
   integer(pInt)                     :: incs                   = 0_pInt, &                               !< number of increments
                                        outputfrequency        = 1_pInt, &                               !< frequency of result writes
                                        restartfrequency       = 0_pInt, &                               !< frequency of restart writes
                                        logscale               = 0_pInt                                  !< linear/logarithmic time inc flag
   logical                           :: followFormerTrajectory = .true.                                  !< follow trajectory of former loadcase
   integer(pInt),        allocatable :: faceID(:)
   type(tFieldBC),       allocatable :: fieldBC(:)
 end type tLoadCase

 type, public :: tFEMInterpolation
   integer(pInt)                                :: n                             
   real(pReal),   dimension(:,:)  , allocatable :: shapeFunc, shapeDerivReal, geomShapeDerivIso
   real(pReal),   dimension(:,:,:), allocatable :: shapeDerivIso
 end type tFEMInterpolation
 
 type, public :: tQuadrature
   integer(pInt)                            :: n
   real(pReal), dimension(:)  , allocatable :: Weights
   real(pReal), dimension(:,:), allocatable :: Points
 end type tQuadrature   
 
 public :: &
   utilities_init, &
   utilities_constitutiveResponse, &
   utilities_indexBoundaryDofs, &
   utilities_projectBCValues, &
   utilities_indexActiveSet, &
   utilities_destroy, &
   FIELD_MECH_ID, &
   FIELD_THERMAL_ID, &
   FIELD_DAMAGE_ID, &
   FIELD_SOLUTE_ID, &
   FIELD_MGTWIN_ID, &
   COMPONENT_MECH_X_ID, &
   COMPONENT_MECH_Y_ID, &
   COMPONENT_MECH_Z_ID, &
   COMPONENT_THERMAL_T_ID, &
   COMPONENT_DAMAGE_PHI_ID, &
   COMPONENT_SOLUTE_CV_ID, &
   COMPONENT_SOLUTE_CVPOT_ID, &
   COMPONENT_SOLUTE_CH_ID, &
   COMPONENT_SOLUTE_CHPOT_ID, &
   COMPONENT_SOLUTE_CVaH_ID, &
   COMPONENT_SOLUTE_CVaHPOT_ID, &
   COMPONENT_MGTWIN_PHI_ID

 external :: &
   MPI_Allreduce, &
   PetscOptionsInsertString, &
   PetscObjectSetName, &
   DMPlexGetHeightStratum, &
   DMGetLabelIdIS, &
   DMPlexGetChart, &
   DMPlexLabelComplete, &
   PetscViewerHDF5Open, &
   PetscViewerHDF5PushGroup, &
   PetscViewerHDF5PopGroup, &
   PetscViewerDestroy

contains 

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags
!--------------------------------------------------------------------------------------------------
subroutine utilities_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use DAMASK_interface, only: &
  getSolverJobName
 use IO, only: &
   IO_error, &
   IO_warning, &
   IO_timeStamp, &
   IO_open_file
 use numerics, only: &                        
   integrationOrder, &
   worldsize, & 
   worldrank, &
   petsc_defaultOptions, &
   petsc_options
 use debug, only: &
   debug_level, &
   debug_SPECTRAL, &
   debug_LEVELBASIC, &
   debug_SPECTRALPETSC, &
   debug_SPECTRALROTATION
 use debug, only: &
   PETSCDEBUG
 use math                                                                                           ! must use the whole module for use of FFTW
 use mesh, only: &
   mesh_NcpElemsGlobal, &
   mesh_maxNips, &
   geomMesh, &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   material_homog, &
   material_phase, &
   microstructure_crystallite

 implicit none

 character(len=1024)                :: petsc_optionsPhysics, grainStr
 integer(pInt)                      :: dimPlex
 integer(pInt)                      :: headerID = 205_pInt
 PetscInt,    dimension(:), pointer :: points
 PetscInt,    allocatable           :: nEntities(:), nOutputCells(:), nOutputNodes(:), mappingCells(:)
 PetscInt                           :: cellStart, cellEnd, cell, ip, dim, ctr, qPt
 PetscInt                           :: homog, cryst, grain, phase 
 PetscInt,              allocatable :: connectivity(:,:)
 Vec                                :: connectivityVec
 PetscScalar, dimension(:), pointer :: results
 PetscErrorCode                     :: ierr

 if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  DAMASK_FEM_utilities init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif
 
!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_SPECTRAL),debug_LEVELBASIC)         /= 0
 debugRotation   = iand(debug_level(debug_SPECTRAL),debug_SPECTRALROTATION)   /= 0
 debugPETSc      = iand(debug_level(debug_SPECTRAL),debug_SPECTRALPETSC)      /= 0
 if(debugPETSc) write(6,'(3(/,a),/)') &
                ' Initializing PETSc with debug options: ', &
                trim(PETScDebug), &
                ' add more using the PETSc_Options keyword in numerics.config '
 flush(6)
 call PetscOptionsClear(PETSC_NULL_OPTIONS,ierr)
 CHKERRQ(ierr)
 if(debugPETSc) call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(PETSCDEBUG),ierr)
 CHKERRQ(ierr)
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_defaultOptions),ierr)
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_options),ierr)
 CHKERRQ(ierr)
 !write(petsc_optionsPhysics,'(a,i0)') '-mechFE_petscspace_order '   , structOrder
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(petsc_optionsPhysics),ierr)
 CHKERRQ(ierr)
 
 wgt = 1.0/real(mesh_maxNips*mesh_NcpElemsGlobal,pReal)

 call PetscViewerHDF5Open(PETSC_COMM_WORLD, trim(getSolverJobName())//'.h5', &
                          FILE_MODE_WRITE, resUnit, ierr); CHKERRQ(ierr)
 call PetscViewerHDF5PushGroup(resUnit, '/', ierr); CHKERRQ(ierr)  
 call DMGetDimension(geomMesh,dimPlex,ierr); CHKERRQ(ierr) 
 allocate(nEntities(dimPlex+1), source=0)
 allocate(nOutputNodes(worldsize), source = 0)
 allocate(nOutputCells(worldsize), source = 0)
 do dim = 0, dimPlex
   call DMGetStratumSize(geomMesh,'depth',dim,nEntities(dim+1),ierr)
   CHKERRQ(ierr)
 enddo
 select case (integrationOrder)
   case(1_pInt)
     nOutputNodes(worldrank+1) = nEntities(1)
   case(2_pInt)  
     nOutputNodes(worldrank+1) = sum(nEntities)
   case default
     nOutputNodes(worldrank+1) = mesh_maxNips*nEntities(dimPlex+1)
 end select    
 nOutputCells(worldrank+1) = count(material_homog > 0_pInt)
 call MPI_Allreduce(MPI_IN_PLACE,nOutputNodes,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,nOutputCells,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
 if (worldrank == 0_pInt) then
   open(unit=headerID, file=trim(getSolverJobName())//'.header', &
        form='FORMATTED', status='REPLACE')
   write(headerID, '(a,i0)') 'dimension : ',   dimPlex                                             
   write(headerID, '(a,i0)') 'number of nodes : ', sum(nOutputNodes)
   write(headerID, '(a,i0)') 'number of cells : ', sum(nOutputCells)
 endif 

 allocate(connectivity(2**dimPlex,nOutputCells(worldrank+1)))
 call DMPlexGetHeightStratum(geomMesh,0,cellStart,cellEnd,ierr)
 CHKERRQ(ierr)
 ctr = 0
 select case (integrationOrder)
   case(1_pInt)
     do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
       call DMPlexGetTransitiveClosure(geomMesh,cell,PETSC_TRUE,points,ierr)
       CHKERRQ(ierr)
       if (dimPlex == 2) then
         connectivity(:,ctr+1) =  [points( 9), points(11), points(13), points(13)]  - nEntities(dimPlex+1)
         ctr = ctr + 1
       else
         connectivity(:,ctr+1) =  [points(23), points(25), points(27), points(27), &
                                   points(29), points(29), points(29), points(29)] - nEntities(dimPlex+1)
         ctr = ctr + 1
       endif 
     enddo 
   
   case(2_pInt)  
     do cell = cellStart, cellEnd-1                                                                     !< loop over all elements 
       call DMPlexGetTransitiveClosure(geomMesh,cell,PETSC_TRUE,points,ierr)
       CHKERRQ(ierr)
       if (dimPlex == 2) then
         connectivity(:,ctr+1) = [points(9 ), points(3), points(1), points(7)]
         connectivity(:,ctr+2) = [points(11), points(5), points(1), points(3)]
         connectivity(:,ctr+3) = [points(13), points(7), points(1), points(5)]
         ctr = ctr + 3
       else
         connectivity(:,ctr+1) = [points(23), points(11), points(3), points(15), points(17), points(5), points(1), points(7)]
         connectivity(:,ctr+2) = [points(25), points(13), points(3), points(11), points(19), points(9), points(1), points(5)]
         connectivity(:,ctr+3) = [points(27), points(15), points(3), points(13), points(21), points(7), points(1), points(9)]
         connectivity(:,ctr+4) = [points(29), points(17), points(7), points(21), points(19), points(5), points(1), points(9)]
         ctr = ctr + 4_pInt
       endif  
     enddo 
   
   case default
     do cell = cellStart, cellEnd-1; do ip = 0, mesh_maxNips-1
       connectivity(:,ctr+1) = cell*mesh_maxNips + ip
       ctr = ctr + 1
     enddo; enddo        
   
 end select 
 connectivity = connectivity + sum(nOutputNodes(1:worldrank))

 call VecCreateMPI(PETSC_COMM_WORLD,dimPlex*nOutputNodes(worldrank+1),dimPlex*sum(nOutputNodes), &
                   coordinatesVec,ierr);CHKERRQ(ierr)
 call PetscObjectSetName(coordinatesVec, 'NodalCoordinates',ierr)
 call VecSetFromOptions(coordinatesVec, ierr); CHKERRQ(ierr)
 
 !allocate(mappingCells(worldsize), source = 0)
 !do homog = 1, material_Nhomogenization
 !  mappingCells = 0_pInt; mappingCells(worldrank+1) = homogOutput(homog)%sizeIpCells
 !  call MPI_Allreduce(MPI_IN_PLACE,mappingCells,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
 !  call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1),sum(mappingCells), &
 !                    homogenizationResultsVec(homog),ierr);CHKERRQ(ierr)
 !  if (sum(mappingCells) > 0) then
 !    call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1)*2**dimPlex,sum(mappingCells)*2**dimPlex, &
 !                      connectivityVec,ierr);CHKERRQ(ierr)
 !    call PetscObjectSetName(connectivityVec,'mapping_'//trim(homogenization_name(homog)),ierr)
 !    CHKERRQ(ierr)
 !    call VecGetArrayF90(connectivityVec,results,ierr); CHKERRQ(ierr) 
 !    results = 0.0_pReal; ctr = 1_pInt
 !    do cell = cellStart, cellEnd-1; do qPt = 1, mesh_maxNips
 !      if (material_homog(qPt,cell+1) == homog) then
 !        results(ctr:ctr+2**dimPlex-1) = real(reshape(connectivity(1:2**dimPlex,mesh_maxNips*cell+qPt), &
 !                                                     shape=[2**dimPlex]))
 !        ctr = ctr + 2**dimPlex
 !      endif  
 !    enddo; enddo
 !    call VecRestoreArrayF90(connectivityVec, results, ierr); CHKERRQ(ierr)
 !    call VecAssemblyBegin(connectivityVec, ierr); CHKERRQ(ierr)
 !    call VecAssemblyEnd  (connectivityVec, ierr); CHKERRQ(ierr)
 !    call VecView(connectivityVec, resUnit, ierr); CHKERRQ(ierr)
 !    call VecDestroy(connectivityVec, ierr); CHKERRQ(ierr)
 !  endif
 !enddo  
 !do cryst = 1, material_Ncrystallite; do grain = 1, homogenization_maxNgrains
 !  mappingCells              = 0_pInt
 !  mappingCells(worldrank+1) = crystalliteOutput(cryst,grain)%sizeIpCells
 !  call MPI_Allreduce(MPI_IN_PLACE,mappingCells,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
 !  call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1),sum(mappingCells), &
 !                    crystalliteResultsVec(cryst,grain),ierr);CHKERRQ(ierr)
 !  if (sum(mappingCells) > 0) then
 !    call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1)*2**dimPlex,sum(mappingCells)*2**dimPlex, &
 !                      connectivityVec,ierr);CHKERRQ(ierr)
 !    write(grainStr,'(a,i0)') 'Grain',grain
 !    call PetscObjectSetName(connectivityVec,'mapping_'// &
 !                                            trim(crystallite_name(cryst))//'_'// &
 !                                            trim(grainStr),ierr)
 !    CHKERRQ(ierr)
 !    call VecGetArrayF90(connectivityVec, results, ierr); CHKERRQ(ierr) 
 !    results = 0.0_pReal; ctr = 1_pInt
 !    do cell = cellStart, cellEnd-1; do qPt = 1, mesh_maxNips
 !      if (homogenization_Ngrains    (mesh_element(3,cell+1)) >= grain .and. &
 !          microstructure_crystallite(mesh_element(4,cell+1)) == cryst) then
 !        results(ctr:ctr+2**dimPlex-1) = real(reshape(connectivity(1:2**dimPlex,mesh_maxNips*cell+qPt), &
 !                                                     shape=[2**dimPlex]))
 !        ctr = ctr + 2**dimPlex
 !      endif  
 !    enddo; enddo
 !    call VecRestoreArrayF90(connectivityVec, results, ierr); CHKERRQ(ierr)
 !    call VecAssemblyBegin(connectivityVec, ierr); CHKERRQ(ierr)
 !    call VecAssemblyEnd  (connectivityVec, ierr); CHKERRQ(ierr)
 !    call VecView(connectivityVec, resUnit, ierr); CHKERRQ(ierr)
 !    call VecDestroy(connectivityVec, ierr); CHKERRQ(ierr)
 !  endif
 !enddo; enddo
 !do phase = 1, material_Nphase; do grain = 1, homogenization_maxNgrains
 !  mappingCells              = 0_pInt
 !  mappingCells(worldrank+1) = phaseOutput(phase,grain)%sizeIpCells
 !  call MPI_Allreduce(MPI_IN_PLACE,mappingCells,worldsize,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
 !  call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1),sum(mappingCells), &
 !                    phaseResultsVec(phase,grain),ierr);CHKERRQ(ierr)
 !  if (sum(mappingCells) > 0) then
 !    call VecCreateMPI(PETSC_COMM_WORLD,mappingCells(worldrank+1)*2**dimPlex,sum(mappingCells)*2**dimPlex, &
 !                      connectivityVec,ierr);CHKERRQ(ierr)
 !    write(grainStr,'(a,i0)') 'Grain',grain
 !    call PetscObjectSetName(connectivityVec,&
 !                            'mapping_'//trim(phase_name(phase))//'_'// &
 !                            trim(grainStr),ierr)
 !    CHKERRQ(ierr)
 !    call VecGetArrayF90(connectivityVec, results, ierr)
 !    CHKERRQ(ierr) 
 !    results = 0.0_pReal; ctr = 1_pInt
 !    do cell = cellStart, cellEnd-1; do qPt = 1, mesh_maxNips
 !      if (material_phase(grain,qPt,cell+1) == phase) then
 !        results(ctr:ctr+2**dimPlex-1) = real(reshape(connectivity(1:2**dimPlex,mesh_maxNips*cell+qPt), &
 !                                                     shape=[2**dimPlex]))
 !        ctr = ctr + 2**dimPlex
 !      endif  
 !    enddo; enddo
 !    call VecRestoreArrayF90(connectivityVec, results, ierr)
 !    CHKERRQ(ierr)
 !    call VecAssemblyBegin(connectivityVec, ierr);CHKERRQ(ierr)
 !    call VecAssemblyEnd  (connectivityVec, ierr);CHKERRQ(ierr)
 !    call VecView(connectivityVec, resUnit, ierr);CHKERRQ(ierr)
 !    call VecDestroy(connectivityVec, ierr); CHKERRQ(ierr)
 !  endif
 !enddo; enddo  
 !if (worldrank == 0_pInt) then
 !  do homog = 1, material_Nhomogenization
 !    call VecGetSize(homogenizationResultsVec(homog),mappingCells(1),ierr)
 !    CHKERRQ(ierr)
 !    if (mappingCells(1) > 0) &
 !      write(headerID, '(a,i0)') 'number of homog_'// &
 !                                trim(homogenization_name(homog))//'_'// &
 !                                'cells : ', mappingCells(1)
 !  enddo
 !  do cryst = 1, material_Ncrystallite; do grain = 1, homogenization_maxNgrains
 !    call VecGetSize(crystalliteResultsVec(cryst,grain),mappingCells(1),ierr)
 !    CHKERRQ(ierr)
 !    write(grainStr,'(a,i0)') 'Grain',grain
 !    if (mappingCells(1) > 0) &
 !      write(headerID, '(a,i0)') 'number of cryst_'// &
 !                                trim(crystallite_name(cryst))//'_'// &
 !                                trim(grainStr)//'_'// &
 !                                'cells : ', mappingCells(1)
 !  enddo; enddo
 !  do phase = 1, material_Nphase; do grain = 1, homogenization_maxNgrains
 !    call VecGetSize(phaseResultsVec(phase,grain),mappingCells(1),ierr)
 !    CHKERRQ(ierr)
 !    write(grainStr,'(a,i0)') 'Grain',grain
 !    if (mappingCells(1) > 0) &
 !      write(headerID, '(a,i0)') 'number of phase_'// &
 !                                trim(phase_name(phase))//'_'//trim(grainStr)//'_'// &
 !                                'cells : ', mappingCells(1)
 !  enddo; enddo  
 !  close(headerID) 
 !endif 

end subroutine utilities_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(timeinc,P_av,forwardData)
 use debug, only: &
   debug_reset, &
   debug_info
 use numerics, only: &
   worldrank
 use math, only: &
   math_transpose33, &
   math_rotate_forward33, &
   math_det33
 use FEsolving, only: &
   restartWrite
 use homogenization, only: &
   materialpoint_F0, &
   materialpoint_F, &
   materialpoint_P, &
   materialpoint_dPdF, &
   materialpoint_stressAndItsTangent
 use mesh, only: &
   mesh_NcpElems  
 
 implicit none
 real(pReal), intent(in)                                         :: timeinc                         !< loading time
 logical,     intent(in)                                         :: forwardData                     !< age results
 
 real(pReal),intent(out), dimension(3,3)                         :: P_av                            !< average PK stress
 
 logical :: &
   age

 integer(pInt) :: &
   j
 real(pReal)   :: defgradDetMin, defgradDetMax, defgradDet
 PetscErrorCode :: ierr

 if (worldrank == 0) &
   write(6,'(/,a)') ' ... evaluating constitutive response ......................................'

 age = .False.
 if (forwardData) then                                                                              ! aging results
   age = .True.
 endif
 if (cutBack) then                                                                                  ! restore saved variables
   age = .False.
 endif
 call debug_reset()

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
 if(debugGeneral) then
   defgradDetMax = -huge(1.0_pReal)
   defgradDetMin = +huge(1.0_pReal)
   do j = 1_pInt, mesh_NcpElems
     defgradDet = math_det33(materialpoint_F(1:3,1:3,1,j))
     defgradDetMax = max(defgradDetMax,defgradDet)
     defgradDetMin = min(defgradDetMin,defgradDet) 
   end do
   write(6,'(a,1x,es11.4)') ' max determinant of deformation =', defgradDetMax
   write(6,'(a,1x,es11.4)') ' min determinant of deformation =', defgradDetMin
   flush(6)
 endif
  
 call materialpoint_stressAndItsTangent(.true.,timeinc)                                             ! calculate P field

 call debug_info()
 
 restartWrite = .false.                                                                             ! reset restartWrite status
 cutBack = .false.                                                                                  ! reset cutBack status
 
 P_av = sum(sum(materialpoint_P,dim=4),dim=3) * wgt                                                    ! average of P 
 C_volAvg = sum(sum(materialpoint_dPdF,dim=6),dim=5) * wgt
 call MPI_Allreduce(MPI_IN_PLACE,P_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,C_volAvg,81,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD, ierr)

end subroutine utilities_constitutiveResponse

!--------------------------------------------------------------------------------------------------
!> @brief Create index sets of boundary dofs (in local and global numbering)
!--------------------------------------------------------------------------------------------------
subroutine utilities_indexBoundaryDofs(dm_local,nFaceSets,numFields,local2global,section,localIS,globalIS)

 implicit none
 
 DM                                    :: dm_local
 ISLocalToGlobalMapping                :: local2global
 PetscSection                          :: section
 PetscInt                              :: nFaceSets, numFields, nDof
 IS,    dimension(nFaceSets,numFields) :: localIS, globalIS
 PetscInt                              :: field, faceSet, point, dof, offset
 PetscInt                              :: localSize, storageSize, ISSize
 PetscInt, dimension(:)  , allocatable :: localIndices
 IS                                    :: faceSetIS, BC_IS, dummyIS
 PetscInt, dimension(:)  ,     pointer :: pFaceSets, pBCvertex, pBCvertexlc
 DMLabel                               :: BCLabel
 PetscErrorCode                        :: ierr

 call DMGetLabel(dm_local,'Face Sets',BCLabel,ierr); CHKERRQ(ierr)
 call DMPlexLabelComplete(dm_local,BCLabel,ierr); CHKERRQ(ierr)
 call PetscSectionGetStorageSize(section,storageSize,ierr); CHKERRQ(ierr)
 call DMGetLabelIdIS(dm_local,'Face Sets',faceSetIS,ierr); CHKERRQ(ierr)
 call ISGetIndicesF90(faceSetIS,pFaceSets,ierr); CHKERRQ(ierr)
 allocate(localIndices (storageSize))
 do faceSet = 1, nFaceSets
   call DMGetStratumSize(dm_local,'Face Sets',pFaceSets(faceSet),ISSize,ierr)
   CHKERRQ(ierr)
   call DMGetStratumIS(dm_local,'Face Sets',pFaceSets(faceSet),BC_IS,ierr)
   CHKERRQ(ierr)
   if (ISSize > 0) call ISGetIndicesF90(BC_IS,pBCvertex,ierr)
   do field = 1, numFields
     localSize  = 0
     do point = 1, ISSize
       call PetscSectionGetFieldDof(section,pBCvertex(point),field-1,nDof,ierr)
       CHKERRQ(ierr)
       call PetscSectionGetFieldOffset(section,pBCvertex(point),field-1,offset,ierr)
       CHKERRQ(ierr)
       do dof = 1, nDof
         localSize = localSize + 1
         localIndices(localSize) = offset + dof - 1
       enddo  
     enddo
     call ISCreateGeneral(PETSC_COMM_SELF,localSize,localIndices,PETSC_COPY_VALUES, &
                          localIS(faceSet,field),ierr)
     CHKERRQ(ierr)
     call ISLocalToGlobalMappingApplyIS(local2global,localIS(faceSet,field), &
                                        globalIS(faceSet,field),ierr)
     CHKERRQ(ierr)
   enddo
   if (ISSize > 0) call ISRestoreIndicesF90(BC_IS,pBCvertex,ierr)
   call ISDestroy(BC_IS,ierr); CHKERRQ(ierr)
 enddo  
 call ISRestoreIndicesF90(faceSetIS,pFaceSets,ierr); CHKERRQ(ierr)
 call ISDestroy(faceSetIS,ierr); CHKERRQ(ierr)
 
 do faceSet = 1, nFaceSets; do field = 1, numFields
   call ISGetSize(globalIS(faceSet,field),ISSize,ierr); CHKERRQ(ierr)
   if (ISSize > 0) then
     call ISGetIndicesF90(localIS(faceSet,field),pBCvertexlc,ierr); CHKERRQ(ierr)
     call ISGetIndicesF90(globalIS(faceSet,field),pBCvertex,ierr); CHKERRQ(ierr)
   endif  
   localSize = 0
   do point = 1, ISSize
     if (pBCvertex(point) >= 0) then
       localSize = localSize + 1
       localIndices(localSize) = pBCvertexlc(point)
     endif  
   enddo
   if (ISSize > 0) then
     call ISRestoreIndicesF90(localIS(faceSet,field),pBCvertexlc,ierr); CHKERRQ(ierr)
     call ISRestoreIndicesF90(globalIS(faceSet,field),pBCvertex,ierr); CHKERRQ(ierr)
   endif  
   call ISDestroy(globalIS(faceSet,field),ierr); CHKERRQ(ierr) 
   call ISCreateGeneral(PETSC_COMM_SELF,localSize,localIndices,PETSC_COPY_VALUES, &
                        globalIS(faceSet,field),ierr) 
   CHKERRQ(ierr)
   if (ISSize > 0) then
     call ISDuplicate(localIS(faceSet,field),dummyIS,ierr); CHKERRQ(ierr)
     call ISDestroy(localIS(faceSet,field),ierr); CHKERRQ(ierr)
     call ISDifference(dummyIS,globalIS(faceSet,field),localIS(faceSet,field),ierr)
     CHKERRQ(ierr)
     call ISDestroy(dummyIS,ierr); CHKERRQ(ierr)
   endif
 enddo; enddo                       
 deallocate(localIndices)
 
end subroutine utilities_indexBoundaryDofs

!--------------------------------------------------------------------------------------------------
!> @brief Project BC values to local vector
!--------------------------------------------------------------------------------------------------
subroutine utilities_projectBCValues(localVec,section,field,comp,bcPointsIS,BCValue,BCDotValue,timeinc)

 implicit none
 
 Vec                  :: localVec
 PetscInt             :: field, comp, nBcPoints, point, dof, numDof, numComp, offset
 PetscSection         :: section
 IS                   :: bcPointsIS
 PetscInt,    pointer :: bcPoints(:)
 PetscScalar, pointer :: localArray(:)
 PetscScalar          :: BCValue,BCDotValue,timeinc
 PetscErrorCode       :: ierr

 call PetscSectionGetFieldComponents(section,field,numComp,ierr); CHKERRQ(ierr)
 call ISGetSize(bcPointsIS,nBcPoints,ierr); CHKERRQ(ierr)
 if (nBcPoints > 0) call ISGetIndicesF90(bcPointsIS,bcPoints,ierr)
 call VecGetArrayF90(localVec,localArray,ierr); CHKERRQ(ierr)
 do point = 1, nBcPoints
   call PetscSectionGetFieldDof(section,bcPoints(point),field,numDof,ierr)
   CHKERRQ(ierr)
   call PetscSectionGetFieldOffset(section,bcPoints(point),field,offset,ierr)
   CHKERRQ(ierr)
   do dof = offset+comp+1, offset+numDof, numComp
     localArray(dof) = localArray(dof) + BCValue + BCDotValue*timeinc
   enddo
 enddo    
 call VecRestoreArrayF90(localVec,localArray,ierr); CHKERRQ(ierr)
 call VecAssemblyBegin(localVec, ierr); CHKERRQ(ierr)
 call VecAssemblyEnd  (localVec, ierr); CHKERRQ(ierr)
 if (nBcPoints > 0) call ISRestoreIndicesF90(bcPointsIS,bcPoints,ierr)
 
end subroutine utilities_projectBCValues

!--------------------------------------------------------------------------------------------------
!> @brief Create index sets of boundary dofs (in local and global numbering)
!--------------------------------------------------------------------------------------------------
subroutine utilities_indexActiveSet(field,section,x_local,f_local,localIS,globalIS)
 use mesh, only: &
   geomMesh

 implicit none
 
 ISLocalToGlobalMapping                   :: local2global
 PetscSection                             :: section
 Vec                                      :: x_local, f_local
 PetscInt                                 :: field
 IS                                       :: localIS, globalIS, dummyIS
 PetscScalar, dimension(:)  ,     pointer :: x_scal, f_scal
 PetscInt                                 :: ISSize
 PetscInt                                 :: chart, chartStart, chartEnd, nDof, dof, offset
 PetscInt                                 :: localSize
 PetscInt,    dimension(:)  , allocatable :: localIndices
 PetscInt,    dimension(:)  ,     pointer :: pBCvertex, pBCvertexlc
 PetscErrorCode                           :: ierr

 call DMGetLocalToGlobalMapping(geomMesh,local2global,ierr)
 CHKERRQ(ierr)
 call DMPlexGetChart(geomMesh,chartStart,chartEnd,ierr)
 CHKERRQ(ierr)
 call VecGetArrayF90(x_local,x_scal,ierr); CHKERRQ(ierr)
 call VecGetArrayF90(f_local,f_scal,ierr); CHKERRQ(ierr)
 localSize = 0
 do chart = chartStart, chartEnd-1
   call PetscSectionGetFieldDof(section,chart,field-1,nDof,ierr); CHKERRQ(ierr)
   call PetscSectionGetFieldOffset(section,chart,field-1,offset,ierr); CHKERRQ(ierr)
   do dof = offset+1, offset+nDof
     if (((x_scal(dof) <       1.0e-8) .and. (f_scal(dof) > 0.0)) .or. &
         ((x_scal(dof) > 1.0 - 1.0e-8) .and. (f_scal(dof) < 0.0))) localSize = localSize + 1
   enddo
 enddo
 allocate(localIndices(localSize))        
 localSize = 0
 do chart = chartStart, chartEnd-1
   call PetscSectionGetFieldDof(section,chart,field-1,nDof,ierr); CHKERRQ(ierr)
   call PetscSectionGetFieldOffset(section,chart,field-1,offset,ierr); CHKERRQ(ierr)
   do dof = offset+1, offset+nDof
     if (((x_scal(dof) <       1.0e-8) .and. (f_scal(dof) > 0.0)) .or. &
         ((x_scal(dof) > 1.0 - 1.0e-8) .and. (f_scal(dof) < 0.0))) then
       localSize = localSize + 1
       localIndices(localSize) = dof-1
     endif  
   enddo
 enddo
 call VecRestoreArrayF90(x_local,x_scal,ierr); CHKERRQ(ierr)
 call VecRestoreArrayF90(f_local,f_scal,ierr); CHKERRQ(ierr)
 call ISCreateGeneral(PETSC_COMM_SELF,localSize,localIndices,PETSC_COPY_VALUES,localIS,ierr)
 CHKERRQ(ierr)
 call ISLocalToGlobalMappingApplyIS(local2global,localIS,globalIS,ierr)
 CHKERRQ(ierr)
 call ISGetSize(globalIS,ISSize,ierr); CHKERRQ(ierr)
 if (ISSize > 0) then
   call ISGetIndicesF90(localIS,pBCvertexlc,ierr); CHKERRQ(ierr)
   call ISGetIndicesF90(globalIS,pBCvertex,ierr); CHKERRQ(ierr)
 endif  
 localSize = 0
 do chart = 1, ISSize
   if (pBCvertex(chart) >= 0) then
     localSize = localSize + 1
     localIndices(localSize) = pBCvertexlc(chart)
   endif  
 enddo
 if (ISSize > 0) then
   call ISRestoreIndicesF90(localIS,pBCvertexlc,ierr); CHKERRQ(ierr)
   call ISRestoreIndicesF90(globalIS,pBCvertex,ierr); CHKERRQ(ierr)
 endif  
 call ISDestroy(globalIS,ierr); CHKERRQ(ierr) 
 call ISCreateGeneral(PETSC_COMM_SELF,localSize,localIndices,PETSC_COPY_VALUES,globalIS,ierr) 
 CHKERRQ(ierr)
 if (ISSize > 0) then
   call ISDuplicate(localIS,dummyIS,ierr); CHKERRQ(ierr)
   call ISDestroy(localIS,ierr); CHKERRQ(ierr)
   call ISDifference(dummyIS,globalIS,localIS,ierr)
   CHKERRQ(ierr)
   call ISDestroy(dummyIS,ierr); CHKERRQ(ierr)
 endif                      
 deallocate(localIndices)
 
end subroutine utilities_indexActiveSet

!--------------------------------------------------------------------------------------------------
!> @brief cleans up
!--------------------------------------------------------------------------------------------------
subroutine utilities_destroy()
 use material, only: &
   homogenization_Ngrains

 !implicit none
 !PetscInt       :: homog, cryst, grain, phase 
 !PetscErrorCode :: ierr

 !call PetscViewerHDF5PopGroup(resUnit, ierr); CHKERRQ(ierr)
 !call VecDestroy(coordinatesVec,ierr); CHKERRQ(ierr)
 !do homog = 1, material_Nhomogenization
 !  call VecDestroy(homogenizationResultsVec(homog),ierr);CHKERRQ(ierr)
 !  do cryst = 1, material_Ncrystallite; do grain = 1, homogenization_Ngrains(homog)
 !    call VecDestroy(crystalliteResultsVec(cryst,grain),ierr);CHKERRQ(ierr)
 !  enddo; enddo
 !  do phase = 1, material_Nphase; do grain = 1, homogenization_Ngrains(homog)
 !    call VecDestroy(phaseResultsVec(phase,grain),ierr);CHKERRQ(ierr)
 !  enddo; enddo  
 !enddo      
 !call PetscViewerDestroy(resUnit, ierr); CHKERRQ(ierr)

end subroutine utilities_destroy


end module FEM_utilities
