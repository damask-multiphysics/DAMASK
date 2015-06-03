!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the different spectral solver variants
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_utilities
 use, intrinsic :: iso_c_binding
 use prec, only: &
   pReal, &
   pInt
 use math, only: &
  math_I3 

 implicit none
 private
#ifdef PETSc
#include <petsc-finclude/petscsys.h>
#endif
 include 'fftw3-mpi.f03'
 
 logical,       public             :: cutBack =.false.                                              !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
 integer(pInt), public, parameter  :: maxPhaseFields = 2_pInt
 integer(pInt), public             :: nActiveFields = 0_pInt
 
!--------------------------------------------------------------------------------------------------
! field labels information
 enum, bind(c)
   enumerator :: FIELD_UNDEFINED_ID, &
                 FIELD_MECH_ID, &
                 FIELD_THERMAL_ID, &
                 FIELD_DAMAGE_ID, &
                 FIELD_VACANCYDIFFUSION_ID
 end enum
 
!--------------------------------------------------------------------------------------------------
! grid related information information
 real(pReal),   public                :: wgt                                                        !< weighting factor 1/Nelems
 
!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 integer(pInt), public                                                    :: grid1Red               !< grid(1)/2
 real   (C_DOUBLE),        public,  dimension(:,:,:,:,:),     pointer     :: tensorField_realMPI    !< real representation (some stress or deformation) of field_fourier
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:,:),     pointer     :: tensorField_fourierMPI !< field on which the Fourier transform operates
 real(C_DOUBLE),           public,  dimension(:,:,:,:),       pointer     :: vectorField_realMPI    !< vector field real representation for fftw
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:),       pointer     :: vectorField_fourierMPI !< vector field fourier representation for fftw
 real(C_DOUBLE),           public,  dimension(:,:,:),         pointer     :: scalarField_realMPI    !< scalar field real representation for fftw
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:),         pointer     :: scalarField_fourierMPI !< scalar field fourier representation for fftw
 real(pReal),              private, dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat              !< gamma operator (field) for spectral method
 real(pReal),              private, dimension(:,:,:,:),       allocatable :: xi                     !< wave vector field for divergence and for gamma operator
 real(pReal),              private, dimension(3,3,3,3)                    :: C_ref                  !< mechanic reference stiffness
 real(pReal), protected,   public,  dimension(3)                          :: scaledGeomSize         !< scaled geometry size for calculation of divergence (Basic, Basic PETSc)
  
!--------------------------------------------------------------------------------------------------
! plans for FFTW
 type(C_PTR),   private :: &
   planTensorForthMPI, &                                                                            !< FFTW MPI plan P(x) to P(k)
   planTensorBackMPI, &                                                                             !< FFTW MPI plan F(k) to F(x)
   planVectorForthMPI, &                                                                            !< FFTW MPI plan P(x) to P(k)
   planVectorBackMPI, &                                                                             !< FFTW MPI plan F(k) to F(x)
   planScalarForthMPI, &                                                                            !< FFTW MPI plan P(x) to P(k)
   planScalarBackMPI, &                                                                             !< FFTW MPI plan F(k) to F(x)
   planDebugForthMPI, &                                                                             !< FFTW MPI plan for scalar field
   planDebugBackMPI, &                                                                              !< FFTW MPI plan for scalar field inverse
   planDivMPI                                                                                       !< FFTW MPI plan for FFTW in case of debugging divergence calculation

!--------------------------------------------------------------------------------------------------
! variables controlling debugging
 logical, private :: &
   debugGeneral, &                                                                                  !< general debugging of spectral solver
   debugDivergence, &                                                                               !< debugging of divergence calculation (comparison to function used for post processing)
   debugFFTW, &                                                                                     !< doing additional FFT on scalar field and compare to results of strided 3D FFT
   debugRotation, &                                                                                 !< also printing out results in lab frame
   debugPETSc                                                                                       !< use some in debug defined options for more verbose PETSc solution

!--------------------------------------------------------------------------------------------------
! derived types
 type, public :: tSolutionState                                                                     !< return type of solution from spectral solver variants
   logical       :: converged         = .true.   
   logical       :: regrid            = .false.   
   logical       :: stagConverged     = .true.   
   logical       :: termIll           = .false.   
   integer(pInt) :: iterationsNeeded  = 0_pInt
 end type tSolutionState

 type, public :: tBoundaryCondition                                                                 !< set of parameters defining a boundary condition
   real(pReal), dimension(3,3) :: values      = 0.0_pReal
   real(pReal), dimension(3,3) :: maskFloat   = 0.0_pReal
   logical,     dimension(3,3) :: maskLogical = .false.
   character(len=64)           :: myType      = 'None'
 end type tBoundaryCondition

 type, public :: tLoadCase
   real(pReal), dimension (3,3) :: rotation               = math_I3                                 !< rotation of BC
   type(tBoundaryCondition) ::     P, &                                                             !< stress BC
                                   deformation                                                      !< deformation BC (Fdot or L)
   real(pReal) ::                  time                   = 0.0_pReal, &                            !< length of increment
                                   temperature            = 300.0_pReal                             !< isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               !< number of increments
                                   outputfrequency        = 1_pInt, &                               !< frequency of result writes
                                   restartfrequency       = 0_pInt, &                               !< frequency of restart writes
                                   logscale               = 0_pInt                                  !< linear/logarithmic time inc flag
   logical ::                      followFormerTrajectory = .true.                                  !< follow trajectory of former loadcase 
   integer(kind(FIELD_UNDEFINED_ID)), allocatable :: ID(:)
 end type tLoadCase
 
 type, public :: tSolutionParams                                                                    !< @todo use here the type definition for a full loadcase including mask
   real(pReal), dimension(3,3) :: P_BC, rotation_BC
   real(pReal) :: timeinc
   real(pReal) :: timeincOld
   real(pReal) :: temperature
   real(pReal) :: density
 end type tSolutionParams
 
 type(tSolutionParams),   private :: params
 
 type, public :: phaseFieldDataBin                                                                  !< set of parameters defining a phase field
   real(pReal)       :: diffusion      = 0.0_pReal, &                                               !< thermal conductivity
                        mobility       = 0.0_pReal, &                                               !< thermal mobility
                        phaseField0    = 0.0_pReal                                                  !< homogeneous damage field starting condition 
   logical           :: active         = .false.
   character(len=64) :: label      = ''
 end type phaseFieldDataBin
 
 public :: &
   utilities_init, &
   utilities_updateGamma, &
   utilities_FFTtensorForward, &
   utilities_FFTtensorBackward, &
   utilities_FFTvectorForward, &
   utilities_FFTvectorBackward, &
   utilities_FFTscalarForward, &
   utilities_FFTscalarBackward, &
   utilities_fourierGammaConvolution, &
   utilities_fourierGreenConvolution, &
   utilities_inverseLaplace, &
   utilities_divergenceRMS, &
   utilities_curlRMS, &
   utilities_fourierScalarGradient, &
   utilities_fourierVectorDivergence, &
   utilities_maskedCompliance, &
   utilities_constitutiveResponse, &
   utilities_calculateRate, &
   utilities_forwardField, &
   utilities_destroy, &
   utilities_updateIPcoords, &
   FIELD_UNDEFINED_ID, &
   FIELD_MECH_ID, &
   FIELD_THERMAL_ID, &
   FIELD_DAMAGE_ID
 private :: &
   utilities_getFilter

contains 

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags, create plans for FFTW
!> @details Sets the debug levels for general, divergence, restart and FFTW from the biwise coding 
!> provided by the debug module to logicals.
!> Allocates all fields used by FFTW and create the corresponding plans depending on the debug
!> level chosen.
!> Initializes FFTW.
!--------------------------------------------------------------------------------------------------
subroutine utilities_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use IO, only: &
   IO_error, &
   IO_warning, &
   IO_timeStamp, &
   IO_open_file
 use numerics, only: &
   fftw_planner_flag, &
   fftw_timelimit, &
   memory_efficient, &
   petsc_options, &
   divergence_correction, &
   worldrank
 use debug, only: &
   debug_level, &
   debug_SPECTRAL, &
   debug_LEVELBASIC, &
   debug_SPECTRALDIVERGENCE, &
   debug_SPECTRALFFTW, &
   debug_SPECTRALPETSC, &
   debug_SPECTRALROTATION
#ifdef PETSc
 use debug, only: &
   PETSCDEBUG
#endif
 use math                                                                                    
 use mesh, only: &
   gridGlobal, &
   gridLocal, &
   gridOffset, &
   geomSizeGlobal

 implicit none
#ifdef PETSc
 external :: &
   PETScOptionsClear, &
   PETScOptionsInsertString, &
   MPI_Abort
 PetscErrorCode :: ierr
#endif  
 integer(pInt)               :: i, j, k
 integer(pInt), dimension(3) :: k_s
 type(C_PTR) :: &
   tensorFieldMPI, &                                                                                !< field cotaining data for FFTW in real and fourier space (in place)
   vectorFieldMPI, &                                                                                !< field cotaining data for FFTW in real space when debugging FFTW (no in place)
   scalarFieldMPI                                                                                   !< field cotaining data for FFTW in real space when debugging FFTW (no in place)
 integer(C_INTPTR_T) :: gridFFTW(3), alloc_local, local_K, local_K_offset
 integer(C_INTPTR_T), parameter :: &
   scalarSize = 1_C_INTPTR_T, &
   vecSize = 3_C_INTPTR_T, &
   tensorSize = 9_C_INTPTR_T
 
 mainProcess: if (worldrank == 0) then
   write(6,'(/,a)')   ' <<<+-  DAMASK_spectral_utilities init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_SPECTRAL),debug_LEVELBASIC)         /= 0
 debugDivergence = iand(debug_level(debug_SPECTRAL),debug_SPECTRALDIVERGENCE) /= 0
 debugFFTW       = iand(debug_level(debug_SPECTRAL),debug_SPECTRALFFTW)       /= 0
 debugRotation   = iand(debug_level(debug_SPECTRAL),debug_SPECTRALROTATION)   /= 0
 debugPETSc      = iand(debug_level(debug_SPECTRAL),debug_SPECTRALPETSC)      /= 0

 if(debugPETSc .and. worldrank == 0_pInt) write(6,'(3(/,a),/)') &
                ' Initializing PETSc with debug options: ', &
                trim(PETScDebug), &
                ' add more using the PETSc_Options keyword in numerics.config '
 flush(6)
 call PetscOptionsClear(ierr); CHKERRQ(ierr)
 if(debugPETSc) call PetscOptionsInsertString(trim(PETSCDEBUG),ierr); CHKERRQ(ierr)
 call PetscOptionsInsertString(trim(petsc_options),ierr); CHKERRQ(ierr)

 grid1Red = gridLocal(1)/2_pInt + 1_pInt
 wgt = 1.0/real(product(gridGlobal),pReal)

 if (worldrank == 0) then
   write(6,'(a,3(i12  ))')  ' grid     a b c: ', gridGlobal
   write(6,'(a,3(es12.5))') ' size     x y z: ', geomSizeGlobal
 endif
 
!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and reso-
! lution-independent divergence
 if (divergence_correction == 1_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomSizeGlobal,1) .and. j /= maxloc(geomSizeGlobal,1)) &
      scaledGeomSize = geomSizeGlobal/geomSizeGlobal(j)
   enddo
 elseif (divergence_correction == 2_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomSizeGlobal/gridGlobal,1) .and. j /= maxloc(geomSizeGlobal/gridGlobal,1)) &
      scaledGeomSize = geomSizeGlobal/geomSizeGlobal(j)*gridGlobal(j)
   enddo
 else
   scaledGeomSize = geomSizeGlobal
 endif


!--------------------------------------------------------------------------------------------------
! MPI allocation
 gridFFTW = int(gridGlobal,C_INTPTR_T)
 alloc_local = fftw_mpi_local_size_3d(gridFFTW(3), gridFFTW(2), gridFFTW(1)/2 +1, &
                                      MPI_COMM_WORLD, local_K, local_K_offset)
 
 tensorFieldMPI = fftw_alloc_complex(tensorSize*alloc_local)
 call c_f_pointer(tensorFieldMPI, tensorField_realMPI,    [3_C_INTPTR_T,3_C_INTPTR_T, &
                  2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])      ! place a pointer for a real tensor representation
 call c_f_pointer(tensorFieldMPI, tensorField_fourierMPI, [3_C_INTPTR_T,3_C_INTPTR_T, &
                  gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T ,              gridFFTW(2),local_K])      ! place a pointer for a fourier tensor representation

 vectorFieldMPI = fftw_alloc_complex(vecSize*alloc_local)
 call c_f_pointer(vectorFieldMPI, vectorField_realMPI,   [3_C_INTPTR_T,&
                  2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])      ! place a pointer for a real vector representation
 call c_f_pointer(vectorFieldMPI, vectorField_fourierMPI,[3_C_INTPTR_T,&
                  gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T,               gridFFTW(2),local_K])      ! place a pointer for a fourier vector representation

 scalarFieldMPI = fftw_alloc_complex(scalarSize*alloc_local)                                        ! allocate data for real representation (no in place transform)
 call c_f_pointer(scalarFieldMPI,    scalarField_realMPI, &
                  [2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1),gridFFTW(2),local_K])                ! place a pointer for a real scalar representation
 call c_f_pointer(scalarFieldMPI, scalarField_fourierMPI, &
                   [             gridFFTW(1)/2_C_INTPTR_T + 1 ,gridFFTW(2),local_K])                ! place a pointer for a fourier scarlar representation
 allocate (xi(3,grid1Red,gridLocal(2),gridLocal(3)),source = 0.0_pReal)                             ! frequencies, only half the size for first dimension
 
!--------------------------------------------------------------------------------------------------
! tensor MPI fftw plans
 planTensorForthMPI = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order
                                       tensorSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                      tensorField_realMPI, tensorField_fourierMPI, &! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! use all processors, planer precision
 planTensorBackMPI  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order
                                      tensorSize,  FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                       tensorField_fourierMPI,tensorField_realMPI, &! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! all processors, planer precision
   
!--------------------------------------------------------------------------------------------------
! vector MPI fftw plans
 planVectorForthMPI = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order
                                          vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                      vectorField_realMPI, vectorField_fourierMPI, &! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! use all processors, planer precision
 planVectorBackMPI  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)],  &       ! dimension, logical length in each dimension in reversed order
                                        vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &  ! no. of transforms, default iblock and oblock
                                                     vectorField_fourierMPI,vectorField_realMPI, &  ! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! all processors, planer precision
   
!--------------------------------------------------------------------------------------------------
! scalar MPI fftw plans                            
 planScalarForthMPI = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order 
                                      scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, & ! no. of transforms, default iblock and oblock
                                                      scalarField_realMPI, scalarField_fourierMPI, & ! input data, output data
                                                               MPI_COMM_WORLD, fftw_planner_flag)   ! use all processors, planer precision
 planScalarBackMPI  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order, no. of transforms
                                      scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, & ! no. of transforms, default iblock and oblock
                                                      scalarField_fourierMPI,scalarField_realMPI, & ! input data, output data
                                                               MPI_COMM_WORLD, fftw_planner_flag)   ! use all processors, planer precision
!--------------------------------------------------------------------------------------------------
! depending on debug options, allocate more memory and create additional plans 
 if (debugDivergence) then
   planDivMPI  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2) ,gridFFTW(1)],vecSize, &
                                            FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                               vectorField_fourierMPI, vectorField_realMPI, &
                                                         MPI_COMM_WORLD, fftw_planner_flag)
 endif

 if (debugFFTW) then
   planDebugForthMPI = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &
                                scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                               scalarField_realMPI, scalarField_fourierMPI, &
                                                         MPI_COMM_WORLD, fftw_planner_flag)
   planDebugBackMPI  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &
                                scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                                scalarField_fourierMPI,scalarField_realMPI, &
                                                         MPI_COMM_WORLD, fftw_planner_flag)
 endif 
!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(0_pInt,ext_msg='Fortran to C')             ! check for correct precision in C
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

 if (debugGeneral .and. worldrank == 0_pInt) write(6,'(/,a)') ' FFTW initialized'
   flush(6)
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 do k = gridOffset+1_pInt, gridOffset+gridLocal(3)
   k_s(3) = k - 1_pInt
   if(k > gridGlobal(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - gridGlobal(3)                            ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
     do j = 1_pInt, gridLocal(2)
       k_s(2) = j - 1_pInt
       if(j > gridGlobal(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - gridGlobal(2)                        ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
         do i = 1_pInt, grid1Red
           k_s(1) = i - 1_pInt                                                                      ! symmetry, junst running from 0,1,...,N/2,N/2+1
           xi(1:3,i,j,k-gridOffset) = real(k_s, pReal)/scaledGeomSize                               ! if divergence_correction is set, frequencies are calculated on unit length
 enddo; enddo; enddo
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(3,3,3,3,1,1,1), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(3,3,3,3,grid1Red,gridLocal(2),gridLocal(3)), source = 0.0_pReal)
 endif

end subroutine utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief updates references stiffness and potentially precalculated gamma operator
!> @details Sets the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of a on-the-fly calculation, only the reference stiffness is updated.
!> The gamma operator is filtered depening on the filter selected in numerics.
!> Also writes out the current reference stiffness for restart.
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateGamma(C,saveReference)
 use IO, only: &
  IO_write_jobRealFile
 use numerics, only: &
   memory_efficient, &
   worldrank
 use mesh, only: &
   gridOffset, &
   gridLocal
 use math, only: &
   math_inv33

 implicit none
 real(pReal), intent(in), dimension(3,3,3,3) :: C                                                   !< input stiffness to store as reference stiffness
 logical    , intent(in)                     :: saveReference                                       !< save reference stiffness to file for restart
 real(pReal),                 dimension(3,3) :: temp33_Real, xiDyad
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o
 
 C_ref = C
 if (saveReference) then
   if (worldrank == 0_pInt) then
     write(6,'(/,a)') ' writing reference stiffness to file'
     flush(6)
     call IO_write_jobRealFile(777,'C_ref',size(C_ref))
     write (777,rec=1) C_ref
     close(777)
   endif
 endif
 
 if(.not. memory_efficient) then                                                   
   do k = gridOffset+1_pInt, gridOffset+gridLocal(3); do j = 1_pInt, gridLocal(2); do i = 1_pInt, grid1Red
     if (any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k-gridOffset)*xi(m, i,j,k-gridOffset)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o,i,j,k-gridOffset) =  temp33_Real(l,n)*xiDyad(m,o)
     endif  
   enddo; enddo; enddo
 endif  

end subroutine utilities_updateGamma

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorForward()
 use math
 use numerics, only: &
   worldrank
 use mesh, only: &
   gridLocal

 implicit none
 integer(pInt)  :: row, column                                                                      ! if debug FFTW, compare 3D array field of row and column
 real(pReal), dimension(2) :: myRand, maxScalarField                                                ! random numbers
 integer(pInt) ::  i, j, k
 PetscErrorCode :: ierr

!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
 if (debugFFTW) then
   if (worldrank == 0_pInt) then
     call random_number(myRand)                                                                     ! two numbers: 0 <= x < 1
     row    = nint(myRand(1)*2_pReal + 1_pReal,pInt)
     column = nint(myRand(2)*2_pReal + 1_pReal,pInt)
   endif
   call MPI_Bcast(row   ,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
   call MPI_Bcast(column,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
   scalarField_realMPI = 0.0_pReal
   scalarField_realMPI(1:gridLocal(1),1:gridLocal(2),1:gridLocal(3)) = &
     tensorField_realMPI(row,column,1:gridLocal(1),1:gridLocal(2),1:gridLocal(3))                   ! store the selected component 
 endif

!--------------------------------------------------------------------------------------------------
! doing the FFT
 call fftw_mpi_execute_dft_r2c(planTensorForthMPI,tensorField_realMPI,tensorField_fourierMPI)
  
!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
 if (debugFFTW) then
   call fftw_mpi_execute_dft_r2c(planDebugForthMPI,scalarField_realMPI,scalarField_fourierMPI)
   where(abs(scalarField_fourierMPI(1:grid1Red,1:gridLocal(2),1:gridLocal(3))) > tiny(1.0_pReal))   ! avoid division by zero
     scalarField_fourierMPI(1:grid1Red,1:gridLocal(2),1:gridLocal(3)) = &  
                   (scalarField_fourierMPI(1:grid1Red,1:gridLocal(2),1:gridLocal(3))-&
                    tensorField_fourierMPI(row,column,1:grid1Red,1:gridLocal(2),1:gridLocal(3)))/&
                    scalarField_fourierMPI(1:grid1Red,1:gridLocal(2),1:gridLocal(3))
   else where
     scalarField_realMPI = cmplx(0.0,0.0,pReal)
   end where
   maxScalarField(1) = maxval(real (scalarField_fourierMPI(1:grid1Red,1:gridLocal(2), &
                                                                     1:gridLocal(3))))
   maxScalarField(2) = maxval(aimag(scalarField_fourierMPI(1:grid1Red,1:gridLocal(2), &
                                                                     1:gridLocal(3))))
   call MPI_reduce(MPI_IN_PLACE,maxScalarField,2,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD,ierr)
   if (worldrank == 0_pInt) then
     write(6,'(/,a,i1,1x,i1,a)') ' .. checking FT results of compontent ', row, column, ' ..'
     write(6,'(/,a,2(es11.4,1x))')  ' max FT relative error = ',&                                   ! print real and imaginary part seperately
       maxScalarField(1),maxScalarField(2)
     flush(6)
   endif
 endif

!--------------------------------------------------------------------------------------------------
! applying filter
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   tensorField_fourierMPI(1:3,1:3,i,j,k) = utilities_getFilter(xi(1:3,i,j,k))* &
                                           tensorField_fourierMPI(1:3,1:3,i,j,k)
 enddo; enddo; enddo
 
end subroutine utilities_FFTtensorForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorBackward()
 use math
 use numerics, only: &
   worldrank          
 use mesh, only: &
   gridLocal

 implicit none
 integer(pInt) :: row, column                                                                       !< if debug FFTW, compare 3D array field of row and column
 real(pReal), dimension(2) :: myRand 
 real(pReal) :: maxScalarField
 PetscErrorCode :: ierr

!--------------------------------------------------------------------------------------------------
! unpack FFT data for conj complex symmetric part. This data is not transformed when using c2r
 if (debugFFTW) then
  if (worldrank == 0_pInt) then
     call random_number(myRand)                                                                     ! two numbers: 0 <= x < 1
     row    = nint(myRand(1)*2_pReal + 1_pReal,pInt)
     column = nint(myRand(2)*2_pReal + 1_pReal,pInt)
   endif
   call MPI_Bcast(row   ,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
   call MPI_Bcast(column,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr)
   scalarField_fourierMPI(1:grid1Red,1:gridLocal(2),1:gridLocal(3)) = &
     tensorField_fourierMPI(row,column,1:grid1Red,1:gridLocal(2),1:gridLocal(3))
 endif 
 
!--------------------------------------------------------------------------------------------------
! doing the iFFT
 call fftw_mpi_execute_dft_c2r(planTensorBackMPI,tensorField_fourierMPI,tensorField_realMPI)        ! back transform of fluct deformation gradient

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
 if (debugFFTW) then
   call fftw_mpi_execute_dft_c2r(planDebugBackMPI,scalarField_fourierMPI,scalarField_realMPI)
   where(abs(real(scalarField_realMPI,pReal)) > tiny(1.0_pReal))                                    ! avoid division by zero
     scalarField_realMPI(1:gridLocal(1),1:gridLocal(2),1:gridLocal(3)) = &
   (scalarField_realMPI(1:gridLocal(1),1:gridLocal(2),1:gridLocal(3)) &
     -  tensorField_realMPI      (row,column,1:gridLocal(1),1:gridLocal(2),1:gridLocal(3)))/ &
       scalarField_realMPI(1:gridLocal(1),1:gridLocal(2),1:gridLocal(3))
   else where
     scalarField_realMPI = cmplx(0.0,0.0,pReal)
   end where
   maxScalarField = maxval(real (scalarField_realMPI(1:gridLocal(1),1:gridLocal(2),1:gridLocal(3))))
   call MPI_reduce(MPI_IN_PLACE,maxScalarField,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD,ierr)
   if (worldrank == 0_pInt) then
     write(6,'(/,a,i1,1x,i1,a)') ' ... checking iFT results of compontent ', row, column, ' ..'
     write(6,'(/,a,es11.4)')     ' max iFT relative error = ', maxScalarField
     flush(6)
   endif
 endif
 
 tensorField_realMPI = tensorField_realMPI * wgt                                                    ! normalize the result by number of elements

end subroutine utilities_FFTtensorBackward

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the scalar
!> is independetly transformed complex to complex and compared to the whole scalar transform
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarForward()
 use math
 use mesh, only: &
   gridLocal
   
 integer(pInt) ::  i, j, k
   
! doing the scalar FFT
 call fftw_mpi_execute_dft_r2c(planScalarForthMPI,scalarField_realMPI,scalarField_fourierMPI)

! applying filter
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   scalarField_fourierMPI(i,j,k) = utilities_getFilter(xi(1:3,i,j,k))* &
                                   scalarField_fourierMPI(i,j,k)
 enddo; enddo; enddo
 
end subroutine utilities_FFTscalarForward

!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the scalar
!> is independetly transformed complex to complex and compared to the whole scalar transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarBackward()
 use math
 
! doing the scalar iFFT
 call fftw_mpi_execute_dft_c2r(planScalarBackMPI,scalarField_fourierMPI,scalarField_realMPI)

 scalarField_realMPI = scalarField_realMPI * wgt                                                    ! normalize the result by number of elements

end subroutine utilities_FFTscalarBackward

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the vector
!> is independetly transformed complex to complex and compared to the whole vector transform
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorForward()
 use math
 use mesh, only: &
   gridLocal
   
 integer(pInt) ::  i, j, k
   
! doing the vecotr FFT
 call fftw_mpi_execute_dft_r2c(planVectorForthMPI,vectorField_realMPI,vectorField_fourierMPI)

! applying filter
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   vectorField_fourierMPI(1:3,i,j,k) = utilities_getFilter(xi(1:3,i,j,k))* &
                                       vectorField_fourierMPI(1:3,i,j,k)
 enddo; enddo; enddo
 
end subroutine utilities_FFTvectorForward

!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the vector
!> is independetly transformed complex to complex and compared to the whole vector transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorBackward()
 use math
 
! doing the vector iFFT
 call fftw_mpi_execute_dft_c2r(planVectorBackMPI,vectorField_fourierMPI,vectorField_realMPI)

 vectorField_realMPI = vectorField_realMPI * wgt                                                    ! normalize the result by number of elements

end subroutine utilities_FFTvectorBackward

!--------------------------------------------------------------------------------------------------
!> @brief doing convolution with inverse laplace kernel
!--------------------------------------------------------------------------------------------------
subroutine utilities_inverseLaplace()
 use math, only: &
   math_inv33, &
   PI
 use numerics, only: &
   worldrank 
 use mesh, only: &
   gridLocal, &
   gridOffset, &
   geomSizeGlobal

 implicit none  
 integer(pInt) :: i, j, k
 real(pReal), dimension(3) :: k_s

 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... doing inverse laplace .................................................'
   flush(6)
 endif
 
 do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2); do i = 1_pInt, grid1Red
   k_s = xi(1:3,i,j,k)*scaledGeomSize
   if (any(k_s /= 0_pInt)) tensorField_fourierMPI(1:3,1:3,i,j,k-gridOffset) =  &
                             tensorField_fourierMPI(1:3,1:3,i,j,k-gridOffset)/ &
                           cmplx(-sum((2.0_pReal*PI*k_s/geomSizeGlobal)* &
                    (2.0_pReal*PI*k_s/geomSizeGlobal)),0.0_pReal,pReal)
 enddo; enddo; enddo
 
 if (gridOffset == 0_pInt) &
  tensorField_fourierMPI(1:3,1:3,1,1,1) = cmplx(0.0_pReal,0.0_pReal,pReal)

end subroutine utilities_inverseLaplace
 

!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real, ensuring that average value = fieldAim
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierGammaConvolution(fieldAim)
 use numerics, only: &
   memory_efficient
 use math, only: &
   math_inv33
 use numerics, only: &
   worldrank 
 use mesh, only: &
   gridLocal, &
   gridOffset

 implicit none  
 real(pReal), intent(in), dimension(3,3) :: fieldAim                                                !< desired average value of the field after convolution
 real(pReal),             dimension(3,3) :: xiDyad, temp33_Real
 complex(pReal),          dimension(3,3) :: temp33_complex
 
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o

 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... doing gamma convolution ................................................'
   flush(6)
 endif
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation (mechanical equilibrium)
 if(memory_efficient) then                                                                          ! memory saving version, on-the-fly calculation of gamma_hat
   do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2) ;do i = 1_pInt, grid1Red
     if(any([i,j,k+gridOffset] /= 1_pInt)) then                                                     ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o, 1,1,1) =  temp33_Real(l,n)*xiDyad(m,o)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, 1,1,1) * &
                               tensorField_fourierMPI(1:3,1:3,i,j,k))
       tensorField_fourierMPI(1:3,1:3,i,j,k) = temp33_Complex 
     endif             
   enddo; enddo; enddo
 else                                                                                               ! use precalculated gamma-operator
   do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
     forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
       temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3,i,j,k) * &
                             tensorField_fourierMPI(1:3,1:3,i,j,k))
     tensorField_fourierMPI(1:3,1:3,i,j,k) = temp33_Complex
   enddo; enddo; enddo
 endif
 
 if (gridOffset == 0_pInt) &
   tensorField_fourierMPI(1:3,1:3,1,1,1) = cmplx(fieldAim/wgt,0.0_pReal,pReal)                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1  

end subroutine utilities_fourierGammaConvolution
 
!--------------------------------------------------------------------------------------------------
!> @brief doing convolution DamageGreenOp_hat * field_real
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierGreenConvolution(D_ref, mobility_ref, deltaT)

 use numerics, only: &
   memory_efficient
 use math, only: &
   math_mul33x3, &
   PI
 use numerics, only: &
   worldrank 
 use mesh, only: &
   gridLocal, &
   geomSizeGlobal

 implicit none  
 real(pReal), dimension(3,3), intent(in) :: D_ref                                                   !< desired average value of the field after convolution
 real(pReal),                 intent(in) :: mobility_ref, deltaT                                    !< desired average value of the field after convolution
 real(pReal), dimension(3)               :: k_s
 real(pReal)                             :: GreenOp_hat
 integer(pInt)                           :: i, j, k
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation
 do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2) ;do i = 1_pInt, grid1Red
   k_s = xi(1:3,i,j,k)*scaledGeomSize
   GreenOp_hat =  1.0_pReal/ &
                  (mobility_ref + deltaT*sum((2.0_pReal*PI*k_s/geomSizeGlobal)* &
                                             math_mul33x3(D_ref,(2.0_pReal*PI*k_s/geomSizeGlobal))))!< GreenOp_hat = iK^{T} * D_ref * iK, K is frequency
   scalarField_fourierMPI(i,j,k) = scalarField_fourierMPI(i,j,k)*GreenOp_hat                 
 enddo; enddo; enddo

end subroutine utilities_fourierGreenConvolution

!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS()
 use math    
 use numerics, only: &
   worldrank                                   
 use mesh, only: &
   gridLocal, &
   gridGlobal

 implicit none
 integer(pInt) :: i, j, k 
 real(pReal) :: &
   err_real_div_RMS, &                                                                              !< RMS of divergence in real space
   err_div_max, &                                                                                   !< maximum value of divergence in Fourier space
   err_real_div_max                                                                                 !< maximum value of divergence in real space
 complex(pReal), dimension(3) ::  temp3_complex
 PetscErrorCode :: ierr

 
 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... calculating divergence ................................................'
   flush(6)
 endif

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
 utilities_divergenceRMS = 0.0_pReal
 do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2)
   do i = 2_pInt, grid1Red -1_pInt                                                                  ! Has somewhere a conj. complex counterpart. Therefore count it twice.
     utilities_divergenceRMS = utilities_divergenceRMS &
           + 2.0_pReal*(sum (real(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,i,j,k),&      ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                       xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&        ! --> sum squared L_2 norm of vector 
                       +sum(aimag(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,i,j,k),& 
                                                       xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
   enddo
   utilities_divergenceRMS = utilities_divergenceRMS &                                              ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart (if grid(1) /= 1)
              + sum( real(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,1       ,j,k), &
                                               xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum(aimag(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,1       ,j,k), &
                                               xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum( real(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,grid1Red,j,k), &
                                               xi(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum(aimag(math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,grid1Red,j,k), &
                                               xi(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal)
 enddo; enddo
 if(gridGlobal(1) == 1_pInt) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pReal          ! counted twice in case of grid(1) == 1
 call MPI_Allreduce(MPI_IN_PLACE,utilities_divergenceRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                      ! RMS in real space calculated with Parsevals theorem from Fourier space

!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
 if (debugDivergence) then                                                                          ! calculate divergence again
   err_div_max = 0.0_pReal
   do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2); do i = 1_pInt, grid1Red
     temp3_Complex = math_mul33x3_complex(tensorField_fourierMPI(1:3,1:3,i,j,k)*wgt,&               ! weighting P_fourier
                                             xi(1:3,i,j,k))*TWOPIIMG
     err_div_max = max(err_div_max,sum(abs(temp3_Complex)**2.0_pReal))
     vectorField_fourierMPI(1:3,i,j,k) = temp3_Complex                                              ! need divergence NOT squared
   enddo; enddo; enddo
   
   call fftw_mpi_execute_dft_c2r(planDivMPI,vectorField_fourierMPI,vectorField_realMPI)             ! already weighted

   err_real_div_RMS = sum(vectorField_realMPI**2.0_pReal)
   call MPI_reduce(MPI_IN_PLACE,err_real_div_RMS,1,MPI_DOUBLE,MPI_SUM,0,PETSC_COMM_WORLD,ierr)
   err_real_div_RMS = sqrt(wgt*err_real_div_RMS)                                                    ! RMS in real space
   
   err_real_div_max = maxval(sum(vectorField_realMPI**2.0_pReal,dim=4))                             ! max in real space
   call MPI_reduce(MPI_IN_PLACE,err_real_div_max,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD,ierr)
   err_real_div_max = sqrt(err_real_div_max)

   call MPI_reduce(MPI_IN_PLACE,err_div_max,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD,ierr)
   err_div_max      = sqrt(    err_div_max)                                                         ! max in Fourier space
   
   if (worldrank == 0_pInt) then 
     write(6,'(/,1x,a,es11.4)')      'error divergence  FT  RMS = ',utilities_divergenceRMS
     write(6,'(1x,a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
     write(6,'(1x,a,es11.4)')        'error divergence  FT  max = ',err_div_max
     write(6,'(1x,a,es11.4)')        'error divergence Real max = ',err_real_div_max
     flush(6)
   endif
 endif

end function utilities_divergenceRMS
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate max of curl of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_curlRMS()
 use math       
 use numerics, only: &
   worldrank                                            
 use mesh, only: &
   gridLocal, &
   gridGlobal

 implicit none
 integer(pInt)  ::  i, j, k, l 
 complex(pReal), dimension(3,3) ::  curl_fourier
 PetscErrorCode :: ierr

 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... calculating curl ......................................................'
   flush(6)
 endif

 !--------------------------------------------------------------------------------------------------
! calculating max curl criterion in Fourier space
 utilities_curlRMS = 0.0_pReal
 
 do k = 1_pInt, gridLocal(3); do j = 1_pInt, gridLocal(2); 
 do i = 2_pInt, grid1Red - 1_pInt
   do l = 1_pInt, 3_pInt
     curl_fourier(l,1) = (+tensorField_fourierMPI(l,3,i,j,k)*xi(2,i,j,k)&
                          -tensorField_fourierMPI(l,2,i,j,k)*xi(3,i,j,k))*TWOPIIMG
     curl_fourier(l,2) = (+tensorField_fourierMPI(l,1,i,j,k)*xi(3,i,j,k)&
                          -tensorField_fourierMPI(l,3,i,j,k)*xi(1,i,j,k))*TWOPIIMG
     curl_fourier(l,3) = (+tensorField_fourierMPI(l,2,i,j,k)*xi(1,i,j,k)&
                          -tensorField_fourierMPI(l,1,i,j,k)*xi(2,i,j,k))*TWOPIIMG
   enddo
   utilities_curlRMS = utilities_curlRMS + &
                      2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)
 enddo 
 do l = 1_pInt, 3_pInt
   curl_fourier = (+tensorField_fourierMPI(l,3,1,j,k)*xi(2,1,j,k)&
                   -tensorField_fourierMPI(l,2,1,j,k)*xi(3,1,j,k))*TWOPIIMG
   curl_fourier = (+tensorField_fourierMPI(l,1,1,j,k)*xi(3,1,j,k)&
                   -tensorField_fourierMPI(l,3,1,j,k)*xi(1,1,j,k))*TWOPIIMG
   curl_fourier = (+tensorField_fourierMPI(l,2,1,j,k)*xi(1,1,j,k)&
                   -tensorField_fourierMPI(l,1,1,j,k)*xi(2,1,j,k))*TWOPIIMG
 enddo
 utilities_curlRMS = utilities_curlRMS + &
                     2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)   
 do l = 1_pInt, 3_pInt  
   curl_fourier = (+tensorField_fourierMPI(l,3,grid1Red,j,k)*xi(2,grid1Red,j,k)&
                   -tensorField_fourierMPI(l,2,grid1Red,j,k)*xi(3,grid1Red,j,k))*TWOPIIMG
   curl_fourier = (+tensorField_fourierMPI(l,1,grid1Red,j,k)*xi(3,grid1Red,j,k)&
                   -tensorField_fourierMPI(l,3,grid1Red,j,k)*xi(1,grid1Red,j,k))*TWOPIIMG
   curl_fourier = (+tensorField_fourierMPI(l,2,grid1Red,j,k)*xi(1,grid1Red,j,k)&
                   -tensorField_fourierMPI(l,1,grid1Red,j,k)*xi(2,grid1Red,j,k))*TWOPIIMG
 enddo
 utilities_curlRMS = utilities_curlRMS + &
                     2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal) 
 enddo; enddo
 call MPI_Allreduce(MPI_IN_PLACE,utilities_curlRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 utilities_curlRMS = sqrt(utilities_curlRMS) * wgt
 if(gridGlobal(1) == 1_pInt) utilities_curlRMS = utilities_curlRMS * 0.5_pReal                      ! counted twice in case of grid(1) == 1

end function utilities_curlRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculates mask compliance tensor used to adjust F to fullfill stress BC
!--------------------------------------------------------------------------------------------------
function utilities_maskedCompliance(rot_BC,mask_stress,C)
 use IO, only: &
   IO_error
 use numerics, only: &
   worldrank
 use math, only: &
   math_Plain3333to99, &
   math_plain99to3333, &
   math_rotate_forward3333, &
   math_rotate_forward33, &
   math_invert

 implicit none
 real(pReal),              dimension(3,3,3,3) :: utilities_maskedCompliance                         !< masked compliance
 real(pReal), intent(in) , dimension(3,3,3,3) :: C                                                  !< current average stiffness
 real(pReal), intent(in) , dimension(3,3)     :: rot_BC                                             !< rotation of load frame
 logical,     intent(in),  dimension(3,3)     :: mask_stress                                        !< mask of stress BC
 integer(pInt) :: j, k, m, n 
 logical, dimension(9) :: mask_stressVector
 real(pReal), dimension(9,9) :: temp99_Real   
 integer(pInt) :: size_reduced = 0_pInt 
 real(pReal),              dimension(:,:), allocatable ::  &
   s_reduced, &                                                                                     !< reduced compliance matrix (depending on number of stress BC)
   c_reduced, &                                                                                     !< reduced stiffness (depending on number of stress BC) 
   sTimesC                                                                                          !< temp variable to check inversion
 logical :: errmatinv
 character(len=1024):: formatString
 
 mask_stressVector = reshape(transpose(mask_stress), [9])
 size_reduced = int(count(mask_stressVector), pInt)
 if(size_reduced > 0_pInt )then
   allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (sTimesC(size_reduced,size_reduced),   source =0.0_pReal)

   temp99_Real = math_Plain3333to99(math_rotate_forward3333(C,rot_BC))

   if(debugGeneral .and. worldrank == 0_pInt) then 
     write(6,'(/,a)') ' ... updating masked compliance ............................................'
     write(6,'(/,a,/,9(9(2x,f12.7,1x)/))',advance='no') ' Stiffness C (load) / GPa =',&
                                                  transpose(temp99_Real)/1.e9_pReal
     flush(6)
   endif
   k = 0_pInt                                                                                       ! calculate reduced stiffness
   do n = 1_pInt,9_pInt
     if(mask_stressVector(n)) then
       k = k + 1_pInt
       j = 0_pInt
       do m = 1_pInt,9_pInt
         if(mask_stressVector(m)) then
           j = j + 1_pInt
           c_reduced(k,j) = temp99_Real(n,m)
   endif; enddo; endif; enddo
   call math_invert(size_reduced, c_reduced, s_reduced, errmatinv)                                  ! invert reduced stiffness
   if(errmatinv) call IO_error(error_ID=400_pInt,ext_msg='utilities_maskedCompliance')
   temp99_Real = 0.0_pReal                                                                          ! fill up compliance with zeros
    k = 0_pInt
    do n = 1_pInt,9_pInt
      if(mask_stressVector(n)) then
        k = k + 1_pInt
        j = 0_pInt
        do m = 1_pInt,9_pInt
          if(mask_stressVector(m)) then
            j = j + 1_pInt
            temp99_Real(n,m) = s_reduced(k,j)
   endif; enddo; endif; enddo
   
!--------------------------------------------------------------------------------------------------
! check if inversion was successful
   sTimesC = matmul(c_reduced,s_reduced)
   do m=1_pInt, size_reduced
     do n=1_pInt, size_reduced
       if(m==n .and. abs(sTimesC(m,n)) > (1.0_pReal + 10.0e-12_pReal)) errmatinv = .true.           ! diagonal elements of S*C should be 1
       if(m/=n .and. abs(sTimesC(m,n)) > (0.0_pReal + 10.0e-12_pReal)) errmatinv = .true.           ! off diagonal elements of S*C should be 0
     enddo
   enddo
   if((debugGeneral .or. errmatinv) .and. (worldrank == 0_pInt)) then                               ! report
     write(formatString, '(I16.16)') size_reduced
     formatString = '(/,a,/,'//trim(formatString)//'('//trim(formatString)//'(2x,es9.2,1x)/))'
     write(6,trim(formatString),advance='no') ' C * S (load) ', &
                                                            transpose(matmul(c_reduced,s_reduced))
     write(6,trim(formatString),advance='no') ' S (load) ', transpose(s_reduced)
   endif
   if(errmatinv) call IO_error(error_ID=400_pInt,ext_msg='utilities_maskedCompliance')
   deallocate(c_reduced)
   deallocate(s_reduced)
   deallocate(sTimesC)
 else
   temp99_real = 0.0_pReal
 endif
 if(debugGeneral .and. worldrank == 0_pInt) &                                                       ! report
   write(6,'(/,a,/,9(9(2x,f12.7,1x)/),/)',advance='no') ' Masked Compliance (load) * GPa =', &
                                                    transpose(temp99_Real*1.e9_pReal)
 flush(6)
 utilities_maskedCompliance = math_Plain99to3333(temp99_Real)

end function utilities_maskedCompliance 


!--------------------------------------------------------------------------------------------------
!> @brief calculate scalar gradient in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierScalarGradient()
 use math, only: &
   PI
 use mesh, only: &
   gridLocal, &
   geomSizeGlobal
 integer(pInt)                :: i, j, k

 vectorField_fourierMPI = cmplx(0.0_pReal,0.0_pReal,pReal)
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   vectorField_fourierMPI(1:3,i,j,k) = scalarField_fourierMPI(i,j,k)* &
                                       cmplx(0.0_pReal,2.0_pReal*PI*xi(1:3,i,j,k)* &
                                       scaledGeomSize/geomSizeGlobal,pReal)
 enddo; enddo; enddo
end subroutine utilities_fourierScalarGradient


!--------------------------------------------------------------------------------------------------
!> @brief calculate vector divergence in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierVectorDivergence()
 use math, only: &
   PI
 use mesh, only: &
   gridLocal, &
   geomSizeGlobal
 integer(pInt)                :: i, j, k, m

 scalarField_fourierMPI = cmplx(0.0_pReal,0.0_pReal,pReal)
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   do m = 1_pInt, 3_pInt
     scalarField_fourierMPI(i,j,k) = &
       scalarField_fourierMPI(i,j,k) + &
       vectorField_fourierMPI(m,i,j,k)* &
       cmplx(0.0_pReal,2.0_pReal*PI*xi(m,i,j,k)*scaledGeomSize(m)/geomSizeGlobal(m),pReal)
   enddo
 enddo; enddo; enddo
end subroutine utilities_fourierVectorDivergence

!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(F_lastInc,F,temperature,timeinc,&
                                          P,C_volAvg,C_minmaxAvg,P_av,forwardData,rotation_BC)
 use debug, only: &
   debug_reset, &
   debug_info
 use numerics, only: &
   worldrank
 use math, only: &
   math_transpose33, &
   math_rotate_forward33, &
   math_det33
 use mesh, only: &
   gridLocal
 use FEsolving, only: &
   restartWrite
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_COLLECT, &
   CPFEM_CALCRESULTS, &
   CPFEM_AGERESULTS
 use homogenization, only: &
   materialpoint_F0, &
   materialpoint_F, &
   materialpoint_P, &
   materialpoint_dPdF
 
 implicit none
 real(pReal), intent(in)                                         :: temperature                     !< temperature (no field)
 real(pReal), intent(in), dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: &
   F_lastInc, &                                                                                     !< target deformation gradient
   F                                                                                                !< previous deformation gradient
 real(pReal), intent(in)                                         :: timeinc                         !< loading time
 logical,     intent(in)                                         :: forwardData                     !< age results
 real(pReal), intent(in), dimension(3,3)                         :: rotation_BC                     !< rotation of load frame
 
 real(pReal),intent(out), dimension(3,3,3,3)                     :: C_volAvg, C_minmaxAvg           !< average stiffness
 real(pReal),intent(out), dimension(3,3)                         :: P_av                            !< average PK stress
 real(pReal),intent(out), dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: P                !< PK stress
 
 integer(pInt) :: &
   calcMode, &                                                                                      !< CPFEM mode for calculation
   j,k
 real(pReal), dimension(3,3,3,3) :: max_dPdF, min_dPdF
 real(pReal)   :: max_dPdF_norm, min_dPdF_norm, defgradDetMin, defgradDetMax, defgradDet
 PetscErrorCode :: ierr

  if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... evaluating constitutive response ......................................'
   flush(6)
 endif
 calcMode    = CPFEM_CALCRESULTS

 if (forwardData) then                                                                              ! aging results
   calcMode    = ior(calcMode,    CPFEM_AGERESULTS)             
   materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(gridLocal)])
 endif
 if (cutBack) then                                                                                  ! restore saved variables
  calcMode    = iand(calcMode,    not(CPFEM_AGERESULTS)) 
 endif

 call CPFEM_general(CPFEM_COLLECT,F_lastInc(1:3,1:3,1,1,1),F(1:3,1:3,1,1,1), &
                   temperature,timeinc,1_pInt,1_pInt)

 materialpoint_F  = reshape(F,[3,3,1,product(gridLocal)])
 call debug_reset()

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
 if(debugGeneral) then
   defgradDetMax = -huge(1.0_pReal)
   defgradDetMin = +huge(1.0_pReal)
   do j = 1_pInt, product(gridLocal)
     defgradDet = math_det33(materialpoint_F(1:3,1:3,1,j))
     defgradDetMax = max(defgradDetMax,defgradDet)
     defgradDetMin = min(defgradDetMin,defgradDet) 
   end do
   call MPI_reduce(MPI_IN_PLACE,defgradDetMax,1,MPI_DOUBLE,MPI_MAX,0,PETSC_COMM_WORLD,ierr)
   call MPI_reduce(MPI_IN_PLACE,defgradDetMin,1,MPI_DOUBLE,MPI_MIN,0,PETSC_COMM_WORLD,ierr)
   if (worldrank == 0_pInt) then 
     write(6,'(a,1x,es11.4)') ' max determinant of deformation =', defgradDetMax
     write(6,'(a,1x,es11.4)') ' min determinant of deformation =', defgradDetMin
     flush(6)
   endif
 endif

 call CPFEM_general(calcMode,F_lastInc(1:3,1:3,1,1,1), F(1:3,1:3,1,1,1), &                          ! first call calculates everything
                    temperature,timeinc,1_pInt,1_pInt)

 max_dPdF = 0.0_pReal
 max_dPdF_norm = 0.0_pReal
 min_dPdF = huge(1.0_pReal)
 min_dPdF_norm = huge(1.0_pReal)
 do k = 1_pInt, product(gridLocal)
   if (max_dPdF_norm < sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)) then
     max_dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)
     max_dPdF_norm = sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)
   endif  
   if (min_dPdF_norm > sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)) then
     min_dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)
     min_dPdF_norm = sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)
   endif  
 end do
 call MPI_Allreduce(MPI_IN_PLACE,max_dPdF,81,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD,ierr)
 call MPI_Allreduce(MPI_IN_PLACE,min_dPdF,81,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD,ierr)
 C_minmaxAvg = 0.5_pReal*(max_dPdF + min_dPdF)

 C_volAvg = sum(sum(materialpoint_dPdF,dim=6),dim=5) * wgt
 call MPI_Allreduce(MPI_IN_PLACE,C_volAvg,81,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 
 call debug_info()
 
 restartWrite = .false.                                                                             ! reset restartWrite status
 cutBack = .false.                                                                                  ! reset cutBack status
 
 P = reshape(materialpoint_P, [3,3,gridLocal(1),gridLocal(2),gridLocal(3)])
 P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt                                                    ! average of P 
 call MPI_Allreduce(MPI_IN_PLACE,P_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 if (debugRotation .and. worldrank == 0_pInt) &
 write(6,'(/,a,/,3(3(2x,f12.4,1x)/))',advance='no') ' Piola--Kirchhoff stress (lab) / MPa =',&
                                                     math_transpose33(P_av)*1.e-6_pReal
 P_av = math_rotate_forward33(P_av,rotation_BC)
 if (worldrank == 0_pInt) then
   write(6,'(/,a,/,3(3(2x,f12.4,1x)/))',advance='no') ' Piola--Kirchhoff stress / MPa =',&
                                                     math_transpose33(P_av)*1.e-6_pReal
   flush(6)
 endif

end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief calculates forward rate, either guessing or just add delta/timeinc
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(avRate,timeinc_old,guess,field_lastInc,field)
 use mesh, only: &
   gridLocal
   
 implicit none
 real(pReal), intent(in), dimension(3,3)                      :: avRate                             !< homogeneous addon
 real(pReal), intent(in) :: &
   timeinc_old                                                                                      !< timeinc of last step
 logical, intent(in) :: &
   guess                                                                                            !< guess along former trajectory
 real(pReal), intent(in), dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: &
   field_lastInc, &                                                                                 !< data of previous step
   field                                                                                            !< data of current step
 real(pReal),             dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: &
   utilities_calculateRate
 
 if (guess) then
   utilities_calculateRate = (field-field_lastInc) / timeinc_old
 else
   utilities_calculateRate = spread(spread(spread(avRate,3,gridLocal(1)),4,gridLocal(2)), &
                                                                           5,gridLocal(3))
 endif

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given, 
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
function utilities_forwardField(timeinc,field_lastInc,rate,aim)
 use mesh, only: &
   gridLocal

 implicit none
 real(pReal), intent(in) :: & 
   timeinc                                                                                          !< timeinc of current step
 real(pReal), intent(in),           dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: &
   field_lastInc, &                                                                                 !< initial field
   rate                                                                                             !< rate by which to forward
 real(pReal), intent(in), optional, dimension(3,3) :: &
   aim                                                                                              !< average field value aim
 real(pReal),                       dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)) :: &
   utilities_forwardField
 real(pReal),                       dimension(3,3)                      :: fieldDiff                !< <a + adot*t> - aim
 PetscErrorCode :: ierr
 
 utilities_forwardField = field_lastInc + rate*timeinc
 if (present(aim)) then                                                                             !< correct to match average
   fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt
   call MPI_Allreduce(MPI_IN_PLACE,fieldDiff,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
   fieldDiff = fieldDiff - aim
   utilities_forwardField = utilities_forwardField - &
                    spread(spread(spread(fieldDiff,3,gridLocal(1)),4,gridLocal(2)),5,gridLocal(3))
 endif

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief calculates filter for fourier convolution depending on type given in numerics.config
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_getFilter(k)
 use IO, only: &
   IO_error
 use numerics, only: &                        
   spectral_filter
 use math, only: &
   PI
 use mesh, only: &
   gridGlobal
  
 implicit none
 real(pReal),intent(in), dimension(3) :: k                                                          !< indices of frequency
  
 utilities_getFilter = 1.0_pReal

 select case (spectral_filter)
   case ('none')                                                                                    ! default, no weighting
   case ('cosine')                                                                                  ! cosine curve with 1 for avg and zero for highest freq
     utilities_getFilter = product(1.0_pReal + cos(PI*k*scaledGeomSize/gridGlobal))/8.0_pReal
   case ('gradient')                                                                                ! gradient, might need grid scaling as for cosine filter
     utilities_getFilter = 1.0_pReal/(1.0_pReal + &
                                      (k(1)*k(1) + k(2)*k(2) + k(3)*k(3)))
   case default
     call IO_error(error_ID = 892_pInt, ext_msg = trim(spectral_filter))
 end select 

 if (gridGlobal(1) /= 1_pInt .and. k(1) == real(grid1Red - 1_pInt,    pReal)/scaledGeomSize(1)) &
   utilities_getFilter = 0.0_pReal
 if (gridGlobal(2) /= 1_pInt .and. k(2) == real(gridGlobal(2)/2_pInt, pReal)/scaledGeomSize(2)) &
   utilities_getFilter = 0.0_pReal                                                                  ! do not delete the whole slice in case of 2D calculation
 if (gridGlobal(2) /= 1_pInt .and. &
       k(2) == real(gridGlobal(2)/2_pInt + mod(gridGlobal(2),2_pInt), pReal)/scaledGeomSize(2)) &
   utilities_getFilter = 0.0_pReal                                                                  ! do not delete the whole slice in case of 2D calculation
 if (gridGlobal(3) /= 1_pInt .and. k(3) == real(gridGlobal(3)/2_pInt, pReal)/scaledGeomSize(3)) &
   utilities_getFilter = 0.0_pReal                                                                  ! do not delete the whole slice in case of 2D calculation
 if (gridGlobal(3) /= 1_pInt .and. &
       k(3) == real(gridGlobal(3)/2_pInt + mod(gridGlobal(3),2_pInt), pReal)/scaledGeomSize(3)) &
   utilities_getFilter = 0.0_pReal                                                                  ! do not delete the whole slice in case of 2D calculation

end function utilities_getFilter


!--------------------------------------------------------------------------------------------------
!> @brief cleans up
!--------------------------------------------------------------------------------------------------
subroutine utilities_destroy()
 use math

 implicit none

 if (debugDivergence) call fftw_destroy_plan(planDivMPI)
 if (debugFFTW) call fftw_destroy_plan(planDebugForthMPI)
 if (debugFFTW) call fftw_destroy_plan(planDebugBackMPI)
 call fftw_destroy_plan(planTensorForthMPI)
 call fftw_destroy_plan(planTensorBackMPI)
 call fftw_destroy_plan(planVectorForthMPI)
 call fftw_destroy_plan(planVectorBackMPI)
 call fftw_destroy_plan(planScalarForthMPI)
 call fftw_destroy_plan(planScalarBackMPI)

end subroutine utilities_destroy


!--------------------------------------------------------------------------------------------------
!> @brief calculate coordinates in current configuration for given defgrad field
! using integration in Fourier space. Similar as in mesh.f90, but using data already defined for
! convolution
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateIPcoords(F)
 use math
 use mesh, only: &
   gridGlobal, &
   gridLocal, &
   gridOffset, &
   geomSizeGlobal, &
   mesh_ipCoordinates 
 implicit none

 real(pReal),   dimension(3,3,gridLocal(1),gridLocal(2),gridLocal(3)), intent(in) :: F
 integer(pInt) :: i, j, k, m
 real(pReal),   dimension(3) :: step, offset_coords, integrator
 real(pReal),   dimension(3,3) :: Favg
 PetscErrorCode :: ierr

 tensorField_realMPI = 0.0_pReal
 tensorField_realMPI(1:3,1:3,1:gridLocal(1),1:gridLocal(2),1:gridLocal(3)) = F
 call utilities_FFTtensorForward()

 integrator = geomSizeGlobal * 0.5_pReal / PI
 step = geomSizeGlobal/real(gridGlobal, pReal)
 
!--------------------------------------------------------------------------------------------------
! average F
 if (gridOffset == 0_pInt) Favg = real(tensorField_fourierMPI(1:3,1:3,1,1,1),pReal)*wgt
 call MPI_Bcast(Favg,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! integration in Fourier space
 vectorField_fourierMPI = cmplx(0.0_pReal, 0.0_pReal, pReal)
 do k = 1_pInt, gridLocal(3);  do j = 1_pInt, gridLocal(2);  do i = 1_pInt,grid1Red
   do m = 1_pInt,3_pInt
     vectorField_fourierMPI(m,i,j,k) = sum(tensorField_fourierMPI(m,1:3,i,j,k)*&
                                           cmplx(0.0_pReal,xi(1:3,i,j,k)*scaledGeomSize*integrator,pReal))
   enddo
   if (any(xi(1:3,i,j,k) /= 0.0_pReal)) &
     vectorField_fourierMPI(1:3,i,j,k) = &
     vectorField_fourierMPI(1:3,i,j,k)/cmplx(-sum(xi(1:3,i,j,k)*scaledGeomSize*xi(1:3,i,j,k)* &
                                                     scaledGeomSize),0.0_pReal,pReal)
 enddo; enddo; enddo
 call fftw_mpi_execute_dft_c2r(planVectorBackMPI,vectorField_fourierMPI,vectorField_realMPI)

!--------------------------------------------------------------------------------------------------
! add average to fluctuation and put (0,0,0) on (0,0,0)
 if (gridOffset == 0_pInt) offset_coords = vectorField_realMPI(1:3,1,1,1)
 call MPI_Bcast(offset_coords,3,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
 offset_coords = math_mul33x3(Favg,step/2.0_pReal) - offset_coords
 m = 1_pInt
 do k = 1_pInt,gridLocal(3); do j = 1_pInt,gridLocal(2); do i = 1_pInt,gridLocal(1)
   mesh_ipCoordinates(1:3,1,m) = vectorField_realMPI(1:3,i,j,k) &
                               + offset_coords &
                               + math_mul33x3(Favg,step*real([i,j,k+gridOffset]-1_pInt,pReal))
   m = m+1_pInt
 enddo; enddo; enddo

end subroutine utilities_updateIPcoords


end module DAMASK_spectral_utilities
