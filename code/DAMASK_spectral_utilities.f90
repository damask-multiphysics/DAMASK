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
 use numerics, only: &                        
   spectral_filter

 implicit none
 private
#ifdef PETSc
#include <petsc/finclude/petscsys.h>
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
 real   (C_DOUBLE),        public,  dimension(:,:,:,:,:),     pointer     :: tensorField_real       !< real representation (some stress or deformation) of field_fourier
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:,:),     pointer     :: tensorField_fourier    !< field on which the Fourier transform operates
 real(C_DOUBLE),           public,  dimension(:,:,:,:),       pointer     :: vectorField_real       !< vector field real representation for fftw
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:),       pointer     :: vectorField_fourier    !< vector field fourier representation for fftw
 real(C_DOUBLE),           public,  dimension(:,:,:),         pointer     :: scalarField_real       !< scalar field real representation for fftw
 complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:),         pointer     :: scalarField_fourier    !< scalar field fourier representation for fftw
 real(pReal),              private, dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat              !< gamma operator (field) for spectral method
 real(pReal),              private, dimension(:,:,:,:),       allocatable :: xi1st                  !< wave vector field for first derivatives
 real(pReal),              private, dimension(:,:,:,:),       allocatable :: xi2nd                  !< wave vector field for second derivatives
 real(pReal),              private, dimension(3,3,3,3)                    :: C_ref                  !< mechanic reference stiffness
 real(pReal), protected,   public,  dimension(3)                          :: scaledGeomSize         !< scaled geometry size for calculation of divergence (Basic, Basic PETSc)
  
!--------------------------------------------------------------------------------------------------
! plans for FFTW
 type(C_PTR),   private :: &
   planTensorForth, &                                                                               !< FFTW MPI plan P(x) to P(k)
   planTensorBack, &                                                                                !< FFTW MPI plan F(k) to F(x)
   planVectorForth, &                                                                               !< FFTW MPI plan v(x) to v(k)
   planVectorBack, &                                                                                !< FFTW MPI plan v(k) to v(x)
   planScalarForth, &                                                                               !< FFTW MPI plan s(x) to s(k)
   planScalarBack                                                                                   !< FFTW MPI plan s(k) to s(x)

!--------------------------------------------------------------------------------------------------
! variables controlling debugging
 logical, private :: &
   debugGeneral, &                                                                                  !< general debugging of spectral solver
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
   real(pReal) ::                  time                   = 0.0_pReal                               !< length of increment
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

 enum, bind(c)                                                                                      
   enumerator :: FILTER_NONE_ID, &                                                            
                 FILTER_GRADIENT_ID, &                                                                  
                 FILTER_COSINE_ID
 end enum
 integer(kind(FILTER_NONE_ID)) :: &
   spectral_filter_ID
 
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
   petsc_defaultOptions, &
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
   grid, &
   grid3, &
   grid3Offset, &
   geomSize

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
   tensorField, &                                                                                   !< field containing data for FFTW in real and fourier space (in place)
   vectorField, &                                                                                   !< field containing data for FFTW in real space when debugging FFTW (no in place)
   scalarField                                                                                      !< field containing data for FFTW in real space when debugging FFTW (no in place)
 integer(C_INTPTR_T), dimension(3) :: gridFFTW
 integer(C_INTPTR_T) :: alloc_local, local_K, local_K_offset
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
 debugRotation   = iand(debug_level(debug_SPECTRAL),debug_SPECTRALROTATION)   /= 0
 debugPETSc      = iand(debug_level(debug_SPECTRAL),debug_SPECTRALPETSC)      /= 0

 if(debugPETSc .and. worldrank == 0_pInt) write(6,'(3(/,a),/)') &
                ' Initializing PETSc with debug options: ', &
                trim(PETScDebug), &
                ' add more using the PETSc_Options keyword in numerics.config '
 flush(6)
 call PetscOptionsClear(ierr); CHKERRQ(ierr)
 if(debugPETSc) call PetscOptionsInsertString(trim(PETSCDEBUG),ierr); CHKERRQ(ierr)
 call PetscOptionsInsertString(trim(petsc_defaultOptions),ierr); CHKERRQ(ierr)
 call PetscOptionsInsertString(trim(petsc_options),ierr); CHKERRQ(ierr)

 grid1Red = grid(1)/2_pInt + 1_pInt
 wgt = 1.0/real(product(grid),pReal)

 if (worldrank == 0) then
   write(6,'(a,3(i12  ))')  ' grid     a b c: ', grid
   write(6,'(a,3(es12.5))') ' size     x y z: ', geomSize
 endif
 
!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and
! resolution-independent divergence
 if (divergence_correction == 1_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomSize,1) .and. j /= maxloc(geomSize,1)) &
      scaledGeomSize = geomSize/geomSize(j)
   enddo
 elseif (divergence_correction == 2_pInt) then
   do j = 1_pInt, 3_pInt
    if (j /= minloc(geomSize/grid,1) .and. j /= maxloc(geomSize/grid,1)) &
      scaledGeomSize = geomSize/geomSize(j)*grid(j)
   enddo
 else
   scaledGeomSize = geomSize
 endif


!--------------------------------------------------------------------------------------------------
! MPI allocation
 gridFFTW = int(grid,C_INTPTR_T)
 alloc_local = fftw_mpi_local_size_3d(gridFFTW(3), gridFFTW(2), gridFFTW(1)/2 +1, &
                                      MPI_COMM_WORLD, local_K, local_K_offset)
 allocate (xi1st(3,grid1Red,grid(2),grid3),source = 0.0_pReal)                                      ! frequencies, only half the size for first dimension
 allocate (xi2nd(3,grid1Red,grid(2),grid3),source = 0.0_pReal)                                      ! frequencies, only half the size for first dimension
 
 tensorField = fftw_alloc_complex(tensorSize*alloc_local)
 call c_f_pointer(tensorField, tensorField_real,    [3_C_INTPTR_T,3_C_INTPTR_T, &
                  2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])      ! place a pointer for a real tensor representation
 call c_f_pointer(tensorField, tensorField_fourier, [3_C_INTPTR_T,3_C_INTPTR_T, &
                  gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T ,              gridFFTW(2),local_K])      ! place a pointer for a fourier tensor representation

 vectorField = fftw_alloc_complex(vecSize*alloc_local)
 call c_f_pointer(vectorField, vectorField_real,   [3_C_INTPTR_T,&
                  2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])      ! place a pointer for a real vector representation
 call c_f_pointer(vectorField, vectorField_fourier,[3_C_INTPTR_T,&
                  gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T,               gridFFTW(2),local_K])      ! place a pointer for a fourier vector representation

 scalarField = fftw_alloc_complex(scalarSize*alloc_local)                                           ! allocate data for real representation (no in place transform)
 call c_f_pointer(scalarField,    scalarField_real, &
                  [2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1),gridFFTW(2),local_K])                ! place a pointer for a real scalar representation
 call c_f_pointer(scalarField, scalarField_fourier, &
                   [             gridFFTW(1)/2_C_INTPTR_T + 1 ,gridFFTW(2),local_K])                ! place a pointer for a fourier scarlar representation
 
!--------------------------------------------------------------------------------------------------
! tensor MPI fftw plans
 planTensorForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &           ! dimension, logical length in each dimension in reversed order
                                       tensorSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                      tensorField_real, tensorField_fourier, &      ! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! use all processors, planer precision
 if (.not. C_ASSOCIATED(planTensorForth)) call IO_error(810, ext_msg='planTensorForth')
 planTensorBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &           ! dimension, logical length in each dimension in reversed order
                                      tensorSize,  FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                       tensorField_fourier,tensorField_real, &      ! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! all processors, planer precision
 if (.not. C_ASSOCIATED(planTensorBack)) call IO_error(810, ext_msg='planTensorBack')
   
!--------------------------------------------------------------------------------------------------
! vector MPI fftw plans
 planVectorForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &           ! dimension, logical length in each dimension in reversed order
                                          vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                      vectorField_real, vectorField_fourier, &      ! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! use all processors, planer precision
 if (.not. C_ASSOCIATED(planVectorForth)) call IO_error(810, ext_msg='planVectorForth')
 planVectorBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)],  &          ! dimension, logical length in each dimension in reversed order
                                        vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &  ! no. of transforms, default iblock and oblock
                                                     vectorField_fourier,vectorField_real, &        ! input data, output data
                                                                MPI_COMM_WORLD, fftw_planner_flag)  ! all processors, planer precision
 if (.not. C_ASSOCIATED(planVectorBack)) call IO_error(810, ext_msg='planVectorBack')
   
!--------------------------------------------------------------------------------------------------
! scalar MPI fftw plans                            
 planScalarForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &           ! dimension, logical length in each dimension in reversed order 
                                      scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, & ! no. of transforms, default iblock and oblock
                                                      scalarField_real, scalarField_fourier, &      ! input data, output data
                                                               MPI_COMM_WORLD, fftw_planner_flag)   ! use all processors, planer precision
 if (.not. C_ASSOCIATED(planScalarForth)) call IO_error(810, ext_msg='planScalarForth')
 planScalarBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &        ! dimension, logical length in each dimension in reversed order, no. of transforms
                                      scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, & ! no. of transforms, default iblock and oblock
                                                      scalarField_fourier,scalarField_real, &        ! input data, output data
                                                               MPI_COMM_WORLD, fftw_planner_flag)   ! use all processors, planer precision
 if (.not. C_ASSOCIATED(planScalarBack)) call IO_error(810, ext_msg='planScalarBack')

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(0_pInt,ext_msg='Fortran to C')             ! check for correct precision in C
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

 if (debugGeneral .and. worldrank == 0_pInt) write(6,'(/,a)') ' FFTW initialized'
   flush(6)
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 do k = grid3Offset+1_pInt, grid3Offset+grid3
   k_s(3) = k - 1_pInt
   if(k > grid(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - grid(3)                                        ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
     do j = 1_pInt, grid(2)
       k_s(2) = j - 1_pInt
       if(j > grid(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - grid(2)                                    ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
         do i = 1_pInt, grid1Red
           k_s(1) = i - 1_pInt                                                                      ! symmetry, junst running from 0,1,...,N/2,N/2+1
           xi2nd(1:3,i,j,k-grid3Offset) = real(k_s, pReal)/scaledGeomSize                           ! if divergence_correction is set, frequencies are calculated on unit length
           where(mod(grid,2)==0 .and. [i,j,k] == grid/2+1)                                          ! for even grids, set the Nyquist Freq component to 0.0
             xi1st(1:3,i,j,k-grid3Offset) = 0.0_pReal
           elsewhere
             xi1st(1:3,i,j,k-grid3Offset) = xi2nd(1:3,i,j,k-grid3Offset)
           endwhere
 enddo; enddo; enddo
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(3,3,3,3,1,1,1), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(3,3,3,3,grid1Red,grid(2),grid3), source = 0.0_pReal)
 endif

 select case (spectral_filter)
   case ('none')                                                                                     ! default, no weighting
     spectral_filter_ID = FILTER_NONE_ID
   case ('cosine')                                                                                   ! cosine curve with 1 for avg and zero for highest freq
     spectral_filter_ID = FILTER_COSINE_ID
   case ('gradient')                                                                                 ! gradient, might need grid scaling as for cosine filter
     spectral_filter_ID = FILTER_GRADIENT_ID
   case default
     call IO_error(892_pInt,ext_msg=trim(spectral_filter))
 end select

end subroutine utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief updates references stiffness and potentially precalculated gamma operator
!> @details Sets the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of a on-the-fly calculation, only the reference stiffness is updated.
!> Also writes out the current reference stiffness for restart.
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateGamma(C,saveReference)
 use IO, only: &
  IO_write_jobRealFile
 use numerics, only: &
   memory_efficient, &
   worldrank
 use mesh, only: &
   grid3Offset, &
   grid3,&
   grid
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
   do k = grid3Offset+1_pInt, grid3Offset+grid3; do j = 1_pInt, grid(2); do i = 1_pInt, grid1Red
     if (any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi1st(l, i,j,k-grid3Offset)*xi1st(m, i,j,k-grid3Offset)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o,i,j,k-grid3Offset) =  temp33_Real(l,n)*xiDyad(m,o)
     endif  
   enddo; enddo; enddo
 endif  

end subroutine utilities_updateGamma

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier
!> @details Does an unweighted filtered FFT transform from real to complex
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorForward()
 use mesh, only: &
   grid3, &
   grid

 implicit none
 integer(pInt) ::  i, j, k

!--------------------------------------------------------------------------------------------------
! doing the tensor FFT
 call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)
  
!--------------------------------------------------------------------------------------------------
! applying filter
 forall(k = 1_pInt:grid3, j = 1_pInt:grid(2), i = 1_pInt:grid1Red) &
   tensorField_fourier(1:3,1:3,i,j,k) = utilities_getFilter(xi2nd(1:3,i,j,k))* &
                                        tensorField_fourier(1:3,1:3,i,j,k)
 
end subroutine utilities_FFTtensorForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorBackward()
 implicit none

 call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
 tensorField_real = tensorField_real * wgt                                                          ! normalize the result by number of elements

end subroutine utilities_FFTtensorBackward

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in scalarField_real to scalarField_fourier
!> @details Does an unweighted filtered FFT transform from real to complex
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarForward()
 use mesh, only: &
   grid3, &
   grid
 
 implicit none
 integer(pInt) ::  i, j, k
   
!--------------------------------------------------------------------------------------------------
! doing the scalar FFT
 call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)

!--------------------------------------------------------------------------------------------------
! applying filter
 forall(k = 1_pInt:grid3, j = 1_pInt:grid(2), i = 1_pInt:grid1Red) &
   scalarField_fourier(i,j,k) = utilities_getFilter(xi2nd(1:3,i,j,k))* &
                                   scalarField_fourier(i,j,k)
 
end subroutine utilities_FFTscalarForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in scalarField_fourier to scalarField_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarBackward()
 implicit none

 call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
 scalarField_real = scalarField_real * wgt                                                    ! normalize the result by number of elements

end subroutine utilities_FFTscalarBackward


!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted filtered FFT transform from real to complex.
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorForward()
 use mesh, only: &
   grid3, &
   grid
 
 implicit none
 integer(pInt) ::  i, j, k

!--------------------------------------------------------------------------------------------------
! doing the vector FFT
 call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)

!--------------------------------------------------------------------------------------------------
! applying filter
 forall(k = 1_pInt:grid3, j = 1_pInt:grid(2), i = 1_pInt:grid1Red) &
   vectorField_fourier(1:3,i,j,k) = utilities_getFilter(xi2nd(1:3,i,j,k))* &
                                       vectorField_fourier(1:3,i,j,k)
 
end subroutine utilities_FFTvectorForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorBackward()
 implicit none

 call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
 vectorField_real = vectorField_real * wgt                                                    ! normalize the result by number of elements

end subroutine utilities_FFTvectorBackward


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
   grid3, &
   grid, &
   grid3Offset

 implicit none  
 real(pReal), intent(in), dimension(3,3) :: fieldAim                                                !< desired average value of the field after convolution
 real(pReal),             dimension(3,3) :: xiDyad, temp33_Real
 complex(pReal),          dimension(3,3) :: temp33_complex
 
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o

 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... doing gamma convolution ...............................................'
   flush(6)
 endif
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation (mechanical equilibrium)
 memoryEfficient: if(memory_efficient) then
   do k = 1_pInt, grid3; do j = 1_pInt, grid(2) ;do i = 1_pInt, grid1Red
     if(any([i,j,k+grid3Offset] /= 1_pInt)) then                                                     ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi1st(l, i,j,k)*xi1st(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o, 1,1,1) =  temp33_Real(l,n)*xiDyad(m,o)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, 1,1,1) * &
                               tensorField_fourier(1:3,1:3,i,j,k))
       tensorField_fourier(1:3,1:3,i,j,k) = temp33_Complex 
     endif             
   enddo; enddo; enddo
 else memoryEfficient
   do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid1Red
     forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
       temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3,i,j,k) * &
                             tensorField_fourier(1:3,1:3,i,j,k))
     tensorField_fourier(1:3,1:3,i,j,k) = temp33_Complex
   enddo; enddo; enddo
 endif memoryEfficient
 
 if (grid3Offset == 0_pInt) &
   tensorField_fourier(1:3,1:3,1,1,1) = cmplx(fieldAim/wgt,0.0_pReal,pReal)                        ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1  

end subroutine utilities_fourierGammaConvolution
 

!--------------------------------------------------------------------------------------------------
!> @brief doing convolution DamageGreenOp_hat * field_real
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierGreenConvolution(D_ref, mobility_ref, deltaT)

 use math, only: &
   math_mul33x3, &
   PI
 use mesh, only: &
   grid, &
   grid3, &
   geomSize

 implicit none  
 real(pReal), dimension(3,3), intent(in) :: D_ref                                                   !< desired average value of the field after convolution
 real(pReal),                 intent(in) :: mobility_ref, deltaT                                    !< desired average value of the field after convolution
 real(pReal), dimension(3)               :: k_s
 real(pReal)                             :: GreenOp_hat
 integer(pInt)                           :: i, j, k
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2) ;do i = 1_pInt, grid1Red
   k_s = xi2nd(1:3,i,j,k)*scaledGeomSize
   GreenOp_hat =  1.0_pReal/ &
                  (mobility_ref + deltaT*sum((2.0_pReal*PI*k_s/geomSize)* &
                                             math_mul33x3(D_ref,(2.0_pReal*PI*k_s/geomSize))))      !< GreenOp_hat = iK^{T} * D_ref * iK, K is frequency
   scalarField_fourier(i,j,k) = scalarField_fourier(i,j,k)*GreenOp_hat                 
 enddo; enddo; enddo

end subroutine utilities_fourierGreenConvolution


!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS()   
 use math, only: &
   TWOPIIMG, &
   math_mul33x3_complex
 use numerics, only: &
   worldrank                                   
 use mesh, only: &
   grid, &
   grid3

 implicit none
 integer(pInt) :: i, j, k 
 PetscErrorCode :: ierr

 
 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... calculating divergence ................................................'
   flush(6)
 endif

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
 utilities_divergenceRMS = 0.0_pReal
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2)
   do i = 2_pInt, grid1Red -1_pInt                                                                  ! Has somewhere a conj. complex counterpart. Therefore count it twice.
     utilities_divergenceRMS = utilities_divergenceRMS &
           + 2.0_pReal*(sum (real(math_mul33x3_complex(tensorField_fourier(1:3,1:3,i,j,k),&         ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                       xi1st(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&     ! --> sum squared L_2 norm of vector 
                       +sum(aimag(math_mul33x3_complex(tensorField_fourier(1:3,1:3,i,j,k),& 
                                                       xi1st(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
   enddo
   utilities_divergenceRMS = utilities_divergenceRMS &                                              ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart (if grid(1) /= 1)
              + sum( real(math_mul33x3_complex(tensorField_fourier(1:3,1:3,1       ,j,k), &
                                               xi1st(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum(aimag(math_mul33x3_complex(tensorField_fourier(1:3,1:3,1       ,j,k), &
                                               xi1st(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum( real(math_mul33x3_complex(tensorField_fourier(1:3,1:3,grid1Red,j,k), &
                                               xi1st(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal) &
              + sum(aimag(math_mul33x3_complex(tensorField_fourier(1:3,1:3,grid1Red,j,k), &
                                               xi1st(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal)
 enddo; enddo
 if(grid(1) == 1_pInt) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pReal                ! counted twice in case of grid(1) == 1
 call MPI_Allreduce(MPI_IN_PLACE,utilities_divergenceRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                      ! RMS in real space calculated with Parsevals theorem from Fourier space


end function utilities_divergenceRMS
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate max of curl of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_curlRMS()
 use math       
 use numerics, only: &
   worldrank                                            
 use mesh, only: &
   grid, &
   grid3

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
 
 do k = 1_pInt, grid3; do j = 1_pInt, grid(2); 
   do i = 2_pInt, grid1Red - 1_pInt                                                                  
     do l = 1_pInt, 3_pInt
       curl_fourier(l,1) = (+tensorField_fourier(l,3,i,j,k)*xi1st(2,i,j,k)&
                            -tensorField_fourier(l,2,i,j,k)*xi1st(3,i,j,k))*TWOPIIMG
       curl_fourier(l,2) = (+tensorField_fourier(l,1,i,j,k)*xi1st(3,i,j,k)&
                            -tensorField_fourier(l,3,i,j,k)*xi1st(1,i,j,k))*TWOPIIMG
       curl_fourier(l,3) = (+tensorField_fourier(l,2,i,j,k)*xi1st(1,i,j,k)&
                            -tensorField_fourier(l,1,i,j,k)*xi1st(2,i,j,k))*TWOPIIMG
     enddo
     utilities_curlRMS = utilities_curlRMS + &
                        2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)! Has somewhere a conj. complex counterpart. Therefore count it twice.
   enddo 
   do l = 1_pInt, 3_pInt                                            
      curl_fourier = (+tensorField_fourier(l,3,1,j,k)*xi1st(2,1,j,k)&
                      -tensorField_fourier(l,2,1,j,k)*xi1st(3,1,j,k))*TWOPIIMG
      curl_fourier = (+tensorField_fourier(l,1,1,j,k)*xi1st(3,1,j,k)&
                      -tensorField_fourier(l,3,1,j,k)*xi1st(1,1,j,k))*TWOPIIMG
      curl_fourier = (+tensorField_fourier(l,2,1,j,k)*xi1st(1,1,j,k)&
                      -tensorField_fourier(l,1,1,j,k)*xi1st(2,1,j,k))*TWOPIIMG
   enddo
   utilities_curlRMS = utilities_curlRMS + &
                                  sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)! this layer (DC) does not have a conjugate complex counterpart (if grid(1) /= 1)
   do l = 1_pInt, 3_pInt  
     curl_fourier = (+tensorField_fourier(l,3,grid1Red,j,k)*xi1st(2,grid1Red,j,k)&
                     -tensorField_fourier(l,2,grid1Red,j,k)*xi1st(3,grid1Red,j,k))*TWOPIIMG
     curl_fourier = (+tensorField_fourier(l,1,grid1Red,j,k)*xi1st(3,grid1Red,j,k)&
                     -tensorField_fourier(l,3,grid1Red,j,k)*xi1st(1,grid1Red,j,k))*TWOPIIMG
     curl_fourier = (+tensorField_fourier(l,2,grid1Red,j,k)*xi1st(1,grid1Red,j,k)&
                     -tensorField_fourier(l,1,grid1Red,j,k)*xi1st(2,grid1Red,j,k))*TWOPIIMG
   enddo
   utilities_curlRMS = utilities_curlRMS + &
                                  sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)! this layer (Nyquist) does not have a conjugate complex counterpart (if grid(1) /= 1) 
 enddo; enddo

 call MPI_Allreduce(MPI_IN_PLACE,utilities_curlRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
 utilities_curlRMS = sqrt(utilities_curlRMS) * wgt
 if(grid(1) == 1_pInt) utilities_curlRMS = utilities_curlRMS * 0.5_pReal                             ! counted twice in case of grid(1) == 1

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
   grid3, &
   grid, &
   geomSize

 implicit none
 integer(pInt)                :: i, j, k

 vectorField_fourier = cmplx(0.0_pReal,0.0_pReal,pReal)
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid1Red
   vectorField_fourier(1:3,i,j,k) = scalarField_fourier(i,j,k)* &
                                       cmplx(0.0_pReal,2.0_pReal*PI*xi1st(1:3,i,j,k)* &
                                       scaledGeomSize/geomSize,pReal)
 enddo; enddo; enddo
end subroutine utilities_fourierScalarGradient


!--------------------------------------------------------------------------------------------------
!> @brief calculate vector divergence in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierVectorDivergence()
 use math, only: &
   PI
 use mesh, only: &
   grid3, &
   grid, &
   geomSize

 implicit none
 integer(pInt)                :: i, j, k, m

 scalarField_fourier = cmplx(0.0_pReal,0.0_pReal,pReal)
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid1Red
   do m = 1_pInt, 3_pInt
     scalarField_fourier(i,j,k) = &
       scalarField_fourier(i,j,k) + &
       vectorField_fourier(m,i,j,k)* &
       cmplx(0.0_pReal,2.0_pReal*PI*xi1st(m,i,j,k)*scaledGeomSize(m)/geomSize(m),pReal)
   enddo
 enddo; enddo; enddo
end subroutine utilities_fourierVectorDivergence


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(F_lastInc,F,timeinc,&
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
   grid,&
   grid3
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
 real(pReal), intent(in), dimension(3,3,grid(1),grid(2),grid3) :: &
   F_lastInc, &                                                                                     !< target deformation gradient
   F                                                                                                !< previous deformation gradient
 real(pReal), intent(in)                                         :: timeinc                         !< loading time
 logical,     intent(in)                                         :: forwardData                     !< age results
 real(pReal), intent(in), dimension(3,3)                         :: rotation_BC                     !< rotation of load frame
 
 real(pReal),intent(out), dimension(3,3,3,3)                     :: C_volAvg, C_minmaxAvg           !< average stiffness
 real(pReal),intent(out), dimension(3,3)                         :: P_av                            !< average PK stress
 real(pReal),intent(out), dimension(3,3,grid(1),grid(2),grid3) :: P                !< PK stress
 
 integer(pInt) :: &
   calcMode, &                                                                                      !< CPFEM mode for calculation
   j,k
 real(pReal), dimension(3,3,3,3) :: max_dPdF, min_dPdF
 real(pReal)   :: max_dPdF_norm, min_dPdF_norm, defgradDetMin, defgradDetMax, defgradDet
 PetscErrorCode :: ierr
 
 external :: &
   MPI_Allreduce

 if (worldrank == 0_pInt) then 
   write(6,'(/,a)') ' ... evaluating constitutive response ......................................'
   flush(6)
 endif
 calcMode    = CPFEM_CALCRESULTS

 if (forwardData) then                                                                              ! aging results
   calcMode    = ior(calcMode,    CPFEM_AGERESULTS)             
   materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid(1:2))*grid3])
 endif
 if (cutBack) then                                                                                  ! restore saved variables
  calcMode    = iand(calcMode,    not(CPFEM_AGERESULTS)) 
 endif

 call CPFEM_general(CPFEM_COLLECT,F_lastInc(1:3,1:3,1,1,1),F(1:3,1:3,1,1,1), &
                   timeinc,1_pInt,1_pInt)

 materialpoint_F  = reshape(F,[3,3,1,product(grid(1:2))*grid3])
 call debug_reset()

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
 if(debugGeneral) then
   defgradDetMax = -huge(1.0_pReal)
   defgradDetMin = +huge(1.0_pReal)
   do j = 1_pInt, product(grid(1:2))*grid3
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
                    timeinc,1_pInt,1_pInt)

 max_dPdF = 0.0_pReal
 max_dPdF_norm = 0.0_pReal
 min_dPdF = huge(1.0_pReal)
 min_dPdF_norm = huge(1.0_pReal)
 do k = 1_pInt, product(grid(1:2))*grid3
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
 
 P = reshape(materialpoint_P, [3,3,grid(1),grid(2),grid3])
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
   grid3, &
   grid
   
 implicit none
 real(pReal), intent(in), dimension(3,3)                      :: avRate                             !< homogeneous addon
 real(pReal), intent(in) :: &
   timeinc_old                                                                                      !< timeinc of last step
 logical, intent(in) :: &
   guess                                                                                            !< guess along former trajectory
 real(pReal), intent(in), dimension(3,3,grid(1),grid(2),grid3) :: &
   field_lastInc, &                                                                                 !< data of previous step
   field                                                                                            !< data of current step
 real(pReal),             dimension(3,3,grid(1),grid(2),grid3) :: &
   utilities_calculateRate
 
 if (guess) then
   utilities_calculateRate = (field-field_lastInc) / timeinc_old
 else
   utilities_calculateRate = spread(spread(spread(avRate,3,grid(1)),4,grid(2)),5,grid3)
 endif

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given, 
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
function utilities_forwardField(timeinc,field_lastInc,rate,aim)
 use mesh, only: &
   grid3, &
   grid
  
 implicit none
 real(pReal), intent(in) :: & 
   timeinc                                                                                          !< timeinc of current step
 real(pReal), intent(in),           dimension(3,3,grid(1),grid(2),grid3) :: &
   field_lastInc, &                                                                                 !< initial field
   rate                                                                                             !< rate by which to forward
 real(pReal), intent(in), optional, dimension(3,3) :: &
   aim                                                                                              !< average field value aim
 real(pReal),                       dimension(3,3,grid(1),grid(2),grid3) :: &
   utilities_forwardField
 real(pReal),                       dimension(3,3)                       :: fieldDiff               !< <a + adot*t> - aim
 PetscErrorCode :: ierr
 
 external :: &
  MPI_Allreduce
 
 utilities_forwardField = field_lastInc + rate*timeinc
 if (present(aim)) then                                                                             !< correct to match average
   fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt
   call MPI_Allreduce(MPI_IN_PLACE,fieldDiff,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
   fieldDiff = fieldDiff - aim
   utilities_forwardField = utilities_forwardField - &
                    spread(spread(spread(fieldDiff,3,grid(1)),4,grid(2)),5,grid3)
 endif

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief calculates filter for fourier convolution depending on type given in numerics.config
!--------------------------------------------------------------------------------------------------
pure function utilities_getFilter(k)
 use math, only: &
   PI
 use mesh, only: &
   grid

 implicit none
 real(pReal),   intent(in), dimension(3) :: k                                                          !< indices of frequency
 complex(pReal)                          :: utilities_getFilter

 select case (spectral_filter_ID)
   case (FILTER_NONE_ID)                                                                            ! default, no weighting
     utilities_getFilter = (1.0_pReal,0.0_pReal)
   case (FILTER_COSINE_ID)                                                                          ! cosine curve with 1 for avg and zero for highest freq
     utilities_getFilter = cmplx(product(1.0_pReal + cos(PI*k*scaledGeomSize/grid))/8.0_pReal,&
                                                                                         0.0_pReal)
   case (FILTER_GRADIENT_ID)                                                                        ! gradient, might need grid scaling as for cosine filter
     utilities_getFilter = cmplx(1.0_pReal/(1.0_pReal + sum(k**2)),0.0_pReal)
   case default
     utilities_getFilter = (0.0_pReal,0.0_pReal)
 end select

end function


!--------------------------------------------------------------------------------------------------
!> @brief calculate coordinates in current configuration for given defgrad field
! using integration in Fourier space. Similar as in mesh.f90, but using data already defined for
! convolution
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateIPcoords(F)
 use math, only: &
   PI, &
   math_mul33x3
 use mesh, only: &
   grid, &
   grid3, &
   grid3Offset, &
   geomSize, &
   mesh_ipCoordinates 
 implicit none

 real(pReal),   dimension(3,3,grid(1),grid(2),grid3), intent(in) :: F
 integer(pInt) :: i, j, k, m
 real(pReal),   dimension(3) :: step, offset_coords, integrator
 real(pReal),   dimension(3,3) :: Favg
 PetscErrorCode :: ierr
 external &
   MPI_Bcast

 tensorField_real = 0.0_pReal
 tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = F
 call utilities_FFTtensorForward()

 integrator = geomSize * 0.5_pReal / PI
 step = geomSize/real(grid, pReal)
 
!--------------------------------------------------------------------------------------------------
! average F
 if (grid3Offset == 0_pInt) Favg = real(tensorField_fourier(1:3,1:3,1,1,1),pReal)*wgt
 call MPI_Bcast(Favg,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)

!--------------------------------------------------------------------------------------------------
! integration in Fourier space
 vectorField_fourier = cmplx(0.0_pReal, 0.0_pReal, pReal)
 do k = 1_pInt, grid3;  do j = 1_pInt, grid(2);  do i = 1_pInt,grid1Red
   do m = 1_pInt,3_pInt
     vectorField_fourier(m,i,j,k) = sum(tensorField_fourier(m,1:3,i,j,k)*&
                                           cmplx(0.0_pReal,xi2nd(1:3,i,j,k)*scaledGeomSize*integrator,pReal))
   enddo
   if (any(abs(xi2nd(1:3,i,j,k)) > tiny(0.0_pReal))) &
     vectorField_fourier(1:3,i,j,k) = &
     vectorField_fourier(1:3,i,j,k)/cmplx(-sum(xi2nd(1:3,i,j,k)*scaledGeomSize*xi2nd(1:3,i,j,k)* &
                                                     scaledGeomSize),0.0_pReal,pReal)
 enddo; enddo; enddo
 call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)

!--------------------------------------------------------------------------------------------------
! add average to fluctuation and put (0,0,0) on (0,0,0)
 if (grid3Offset == 0_pInt) offset_coords = vectorField_real(1:3,1,1,1)
 call MPI_Bcast(offset_coords,3,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
 offset_coords = math_mul33x3(Favg,step/2.0_pReal) - offset_coords
 m = 1_pInt
 do k = 1_pInt,grid3; do j = 1_pInt,grid(2); do i = 1_pInt,grid(1)
   mesh_ipCoordinates(1:3,1,m) = vectorField_real(1:3,i,j,k) &
                               + offset_coords &
                               + math_mul33x3(Favg,step*real([i,j,k+grid3Offset]-1_pInt,pReal))
   m = m+1_pInt
 enddo; enddo; enddo

end subroutine utilities_updateIPcoords


!--------------------------------------------------------------------------------------------------
!> @brief cleans up
!--------------------------------------------------------------------------------------------------
subroutine utilities_destroy()
 implicit none

 call fftw_destroy_plan(planTensorForth)
 call fftw_destroy_plan(planTensorBack)
 call fftw_destroy_plan(planVectorForth)
 call fftw_destroy_plan(planVectorBack)
 call fftw_destroy_plan(planScalarForth)
 call fftw_destroy_plan(planScalarBack)

end subroutine utilities_destroy


end module DAMASK_spectral_utilities
