!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the different spectral solver variants
!--------------------------------------------------------------------------------------------------
module spectral_utilities
  use, intrinsic :: iso_c_binding

#include <petsc/finclude/petscsys.h>
  use PETScSys

  use prec
  use DAMASK_interface
  use parallelization
  use math
  use rotations
  use IO
  use config
  use discretization_grid
  use discretization
  use homogenization

  implicit none
  private

  include 'fftw3-mpi.f03'

!--------------------------------------------------------------------------------------------------
! field labels information
  enum, bind(c); enumerator :: &
    FIELD_UNDEFINED_ID, &
    FIELD_MECH_ID, &
    FIELD_THERMAL_ID, &
    FIELD_DAMAGE_ID
  end enum

!--------------------------------------------------------------------------------------------------
! grid related information information
  real(pReal), protected,  public                :: wgt                                             !< weighting factor 1/Nelems
  integer,     protected,  public                :: grid1Red                                        !< grid(1)/2
  real(pReal), protected,  public,  dimension(3) :: scaledGeomSize                                  !< scaled geometry size for calculation of divergence

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW

  real   (C_DOUBLE),        public,  dimension(:,:,:,:,:),     pointer     :: tensorField_real      !< real representation (some stress or deformation) of field_fourier
  complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:,:),     pointer     :: tensorField_fourier   !< field on which the Fourier transform operates
  real(C_DOUBLE),           public,  dimension(:,:,:,:),       pointer     :: vectorField_real      !< vector field real representation for fftw
  complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:,:),       pointer     :: vectorField_fourier   !< vector field fourier representation for fftw
  real(C_DOUBLE),           public,  dimension(:,:,:),         pointer     :: scalarField_real      !< scalar field real representation for fftw
  complex(C_DOUBLE_COMPLEX),public,  dimension(:,:,:),         pointer     :: scalarField_fourier   !< scalar field fourier representation for fftw
  complex(pReal),                    dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat             !< gamma operator (field) for spectral method
  complex(pReal),                    dimension(:,:,:,:),       allocatable :: xi1st                 !< wave vector field for first derivatives
  complex(pReal),                    dimension(:,:,:,:),       allocatable :: xi2nd                 !< wave vector field for second derivatives
  real(pReal),                       dimension(3,3,3,3)                    :: C_ref                 !< mechanic reference stiffness


!--------------------------------------------------------------------------------------------------
! plans for FFTW
  type(C_PTR) :: &
    planTensorForth, &                                                                              !< FFTW MPI plan P(x) to P(k)
    planTensorBack, &                                                                               !< FFTW MPI plan F(k) to F(x)
    planVectorForth, &                                                                              !< FFTW MPI plan v(x) to v(k)
    planVectorBack, &                                                                               !< FFTW MPI plan v(k) to v(x)
    planScalarForth, &                                                                              !< FFTW MPI plan s(x) to s(k)
    planScalarBack                                                                                  !< FFTW MPI plan s(k) to s(x)

!--------------------------------------------------------------------------------------------------
! variables controlling debugging
  logical :: &
    debugGeneral, &                                                                                 !< general debugging of spectral solver
    debugRotation, &                                                                                !< also printing out results in lab frame
    debugPETSc                                                                                      !< use some in debug defined options for more verbose PETSc solution

!--------------------------------------------------------------------------------------------------
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from spectral solver variants
    integer :: &
       iterationsNeeded  = 0
    logical :: &
       converged         = .true., &
       stagConverged     = .true., &
       termIll           = .false.
  end type tSolutionState

  type, public :: tBoundaryCondition                                                                !< set of parameters defining a boundary condition
    real(pReal), dimension(3,3)   :: values = 0.0_pReal
    logical,     dimension(3,3)   :: mask   = .false.
    character(len=:), allocatable :: myType
  end type tBoundaryCondition

  type, public :: tSolutionParams
    real(pReal), dimension(3,3) :: stress_BC
    logical, dimension(3,3)     :: stress_mask
    type(rotation)              :: rotation_BC
    real(pReal) :: timeinc
  end type tSolutionParams

  type :: tNumerics
    integer :: &
      divergence_correction                                                                         !< scale divergence/curl calculation: [0: no correction, 1: size scaled to 1, 2: size scaled to Npoints]
    logical :: &
      memory_efficient                                                                              !< calculate gamma operator on the fly
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  enum, bind(c); enumerator :: &
    DERIVATIVE_CONTINUOUS_ID, &
    DERIVATIVE_CENTRAL_DIFF_ID, &
    DERIVATIVE_FWBW_DIFF_ID
  end enum

  integer(kind(DERIVATIVE_CONTINUOUS_ID)) :: &
    spectral_derivative_ID

  public :: &
    spectral_utilities_init, &
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
    utilities_fourierVectorGradient, &
    utilities_fourierTensorDivergence, &
    utilities_maskedCompliance, &
    utilities_constitutiveResponse, &
    utilities_calculateRate, &
    utilities_forwardField, &
    utilities_updateCoords, &
    utilities_saveReferenceStiffness, &
    FIELD_UNDEFINED_ID, &
    FIELD_MECH_ID, &
    FIELD_THERMAL_ID, &
    FIELD_DAMAGE_ID

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags, create plans for FFTW
!> @details Sets the debug levels for general, divergence, restart, and FFTW from the bitwise coding
!> provided by the debug module to logicals.
!> Allocate all fields used by FFTW and create the corresponding plans depending on the debug
!> level chosen.
!> Initializes FFTW.
!--------------------------------------------------------------------------------------------------
subroutine spectral_utilities_init

  PetscErrorCode :: ierr
  integer        :: i, j, k, &
    FFTW_planner_flag
  integer, dimension(3) :: k_s
  type(C_PTR) :: &
    tensorField, &                                                                                  !< field containing data for FFTW in real and fourier space (in place)
    vectorField, &                                                                                  !< field containing data for FFTW in real space when debugging FFTW (no in place)
    scalarField                                                                                     !< field containing data for FFTW in real space when debugging FFTW (no in place)
  integer(C_INTPTR_T), dimension(3) :: gridFFTW
  integer(C_INTPTR_T) :: alloc_local, local_K, local_K_offset
  integer(C_INTPTR_T), parameter :: &
    scalarSize = 1_C_INTPTR_T, &
    vecSize    = 3_C_INTPTR_T, &
    tensorSize = 9_C_INTPTR_T
  character(len=*), parameter :: &
    PETSCDEBUG = ' -snes_view -snes_monitor '
  class(tNode) , pointer :: &
    num_grid, &
    debug_grid                                                                                      ! pointer to grid  debug options

  print'(/,a)', ' <<<+-  spectral_utilities init  -+>>>'

  print*, 'M. Diehl, Diploma Thesis TU München, 2010'
  print*, 'https://doi.org/10.13140/2.1.3234.3840'//IO_EOL

  print*, 'P. Eisenlohr et al., International Journal of Plasticity 46:37–53, 2013'
  print*, 'https://doi.org/10.1016/j.ijplas.2012.09.012'//IO_EOL

  print*, 'P. Shanthraj et al., International Journal of Plasticity 66:31–45, 2015'
  print*, 'https://doi.org/10.1016/j.ijplas.2014.02.006'//IO_EOL

  print*, 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print*, 'https://doi.org/10.1007/978-981-10-6855-3_80'

!--------------------------------------------------------------------------------------------------
! set debugging parameters
  num_grid        => config_numerics%get('grid',defaultVal=emptyDict)

  debug_grid      => config_debug%get('grid',defaultVal=emptyList)
  debugGeneral    =  debug_grid%contains('basic')
  debugRotation   =  debug_grid%contains('rotation')
  debugPETSc      =  debug_grid%contains('PETSc')

  if(debugPETSc) print'(3(/,a),/)', &
                 ' Initializing PETSc with debug options: ', &
                 trim(PETScDebug), &
                 ' add more using the "PETSc_options" keyword in numerics.yaml'
  flush(IO_STDOUT)

  call PetscOptionsClear(PETSC_NULL_OPTIONS,ierr)
  CHKERRQ(ierr)
  if(debugPETSc) call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(PETSCDEBUG),ierr)
  CHKERRQ(ierr)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                                num_grid%get_asString('PETSc_options',defaultVal=''),ierr)
  CHKERRQ(ierr)

  grid1Red = grid(1)/2 + 1
  wgt = 1.0/real(product(grid),pReal)

  num%memory_efficient      = num_grid%get_asInt('memory_efficient',      defaultVal=1) > 0         ! ToDo: should be logical in YAML file
  num%divergence_correction = num_grid%get_asInt('divergence_correction', defaultVal=2)

  if (num%divergence_correction < 0 .or. num%divergence_correction > 2) &
    call IO_error(301,ext_msg='divergence_correction')

  select case (num_grid%get_asString('derivative',defaultVal='continuous'))
    case ('continuous')
      spectral_derivative_ID = DERIVATIVE_CONTINUOUS_ID
    case ('central_difference')
      spectral_derivative_ID = DERIVATIVE_CENTRAL_DIFF_ID
    case ('FWBW_difference')
      spectral_derivative_ID = DERIVATIVE_FWBW_DIFF_ID
    case default
      call IO_error(892,ext_msg=trim(num_grid%get_asString('derivative')))
  end select

!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and
! resolution-independent divergence
  if (num%divergence_correction == 1) then
    do j = 1, 3
     if (j /= minloc(geomSize,1) .and. j /= maxloc(geomSize,1)) &
       scaledGeomSize = geomSize/geomSize(j)
    enddo
  elseif (num%divergence_correction == 2) then
    do j = 1, 3
     if (      j /= int(minloc(geomSize/real(grid,pReal),1)) &
         .and. j /= int(maxloc(geomSize/real(grid,pReal),1))) &
       scaledGeomSize = geomSize/geomSize(j)*real(grid(j),pReal)
    enddo
  else
    scaledGeomSize = geomSize
  endif

  select case(IO_lc(num_grid%get_asString('fftw_plan_mode',defaultVal='FFTW_MEASURE')))
    case('fftw_estimate')                                                                           ! ordered from slow execution (but fast plan creation) to fast execution
      FFTW_planner_flag = FFTW_ESTIMATE
    case('fftw_measure')
      FFTW_planner_flag = FFTW_MEASURE
    case('fftw_patient')
      FFTW_planner_flag = FFTW_PATIENT
    case('fftw_exhaustive')
      FFTW_planner_flag = FFTW_EXHAUSTIVE
    case default
      call IO_warning(warning_ID=47,ext_msg=trim(IO_lc(num_grid%get_asString('fftw_plan_mode'))))
      FFTW_planner_flag = FFTW_MEASURE
  end select

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
  if (pReal /= C_DOUBLE .or. kind(1) /= C_INT) error stop 'C and Fortran datatypes do not match'
  call fftw_set_timelimit(num_grid%get_asFloat('fftw_timelimit',defaultVal=-1.0_pReal))

  print*, 'FFTW initialized'; flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! MPI allocation
  gridFFTW = int(grid,C_INTPTR_T)
  alloc_local = fftw_mpi_local_size_3d(gridFFTW(3), gridFFTW(2), gridFFTW(1)/2 +1, &
                                       PETSC_COMM_WORLD, local_K, local_K_offset)
  allocate (xi1st (3,grid1Red,grid(2),grid3),source = cmplx(0.0_pReal,0.0_pReal,pReal))             ! frequencies for first derivatives, only half the size for first dimension
  allocate (xi2nd (3,grid1Red,grid(2),grid3),source = cmplx(0.0_pReal,0.0_pReal,pReal))             ! frequencies for second derivatives, only half the size for first dimension

  tensorField = fftw_alloc_complex(tensorSize*alloc_local)
  call c_f_pointer(tensorField, tensorField_real,    [3_C_INTPTR_T,3_C_INTPTR_T, &
                   2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])     ! place a pointer for a real tensor representation
  call c_f_pointer(tensorField, tensorField_fourier, [3_C_INTPTR_T,3_C_INTPTR_T, &
                   gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T ,              gridFFTW(2),local_K])     ! place a pointer for a fourier tensor representation

  vectorField = fftw_alloc_complex(vecSize*alloc_local)
  call c_f_pointer(vectorField, vectorField_real,   [3_C_INTPTR_T,&
                   2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T),gridFFTW(2),local_K])     ! place a pointer for a real vector representation
  call c_f_pointer(vectorField, vectorField_fourier,[3_C_INTPTR_T,&
                   gridFFTW(1)/2_C_INTPTR_T + 1_C_INTPTR_T,               gridFFTW(2),local_K])     ! place a pointer for a fourier vector representation

  scalarField = fftw_alloc_complex(scalarSize*alloc_local)                                          ! allocate data for real representation (no in place transform)
  call c_f_pointer(scalarField,    scalarField_real, &
                   [2_C_INTPTR_T*(gridFFTW(1)/2_C_INTPTR_T + 1),gridFFTW(2),local_K])               ! place a pointer for a real scalar representation
  call c_f_pointer(scalarField, scalarField_fourier, &
                    [             gridFFTW(1)/2_C_INTPTR_T + 1 ,gridFFTW(2),local_K])               ! place a pointer for a fourier scarlar representation

!--------------------------------------------------------------------------------------------------
! tensor MPI fftw plans
  planTensorForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &          ! dimension, logical length in each dimension in reversed order
                                        tensorSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                       tensorField_real, tensorField_fourier, &     ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! use all processors, planer precision
  if (.not. C_ASSOCIATED(planTensorForth)) error stop 'FFTW error'
  planTensorBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &          ! dimension, logical length in each dimension in reversed order
                                       tensorSize,  FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                        tensorField_fourier,tensorField_real, &     ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! all processors, planer precision
  if (.not. C_ASSOCIATED(planTensorBack)) error stop 'FFTW error'

!--------------------------------------------------------------------------------------------------
! vector MPI fftw plans
  planVectorForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &          ! dimension, logical length in each dimension in reversed order
                                           vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK,&! no. of transforms, default iblock and oblock
                                                       vectorField_real, vectorField_fourier, &     ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! use all processors, planer precision
  if (.not. C_ASSOCIATED(planVectorForth)) error stop 'FFTW error'
  planVectorBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)],  &         ! dimension, logical length in each dimension in reversed order
                                         vecSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, & ! no. of transforms, default iblock and oblock
                                                      vectorField_fourier,vectorField_real, &       ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! all processors, planer precision
  if (.not. C_ASSOCIATED(planVectorBack)) error stop 'FFTW error'

!--------------------------------------------------------------------------------------------------
! scalar MPI fftw plans
  planScalarForth = fftw_mpi_plan_many_dft_r2c(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &          ! dimension, logical length in each dimension in reversed order
                                       scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                       scalarField_real, scalarField_fourier, &     ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! use all processors, planer precision
  if (.not. C_ASSOCIATED(planScalarForth)) error stop 'FFTW error'
  planScalarBack  = fftw_mpi_plan_many_dft_c2r(3, [gridFFTW(3),gridFFTW(2),gridFFTW(1)], &          ! dimension, logical length in each dimension in reversed order, no. of transforms
                                       scalarSize, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &! no. of transforms, default iblock and oblock
                                                       scalarField_fourier,scalarField_real, &      ! input data, output data
                                                               PETSC_COMM_WORLD, FFTW_planner_flag) ! use all processors, planer precision
  if (.not. C_ASSOCIATED(planScalarBack)) error stop 'FFTW error'

!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
  do k = grid3Offset+1, grid3Offset+grid3
    k_s(3) = k - 1
    if(k > grid(3)/2 + 1) k_s(3) = k_s(3) - grid(3)                                                 ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
      do j = 1, grid(2)
        k_s(2) = j - 1
        if(j > grid(2)/2 + 1) k_s(2) = k_s(2) - grid(2)                                             ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
          do i = 1, grid1Red
            k_s(1) = i - 1                                                                          ! symmetry, junst running from 0,1,...,N/2,N/2+1
            xi2nd(1:3,i,j,k-grid3Offset) = utilities_getFreqDerivative(k_s)
            where(mod(grid,2)==0 .and. [i,j,k] == grid/2+1 .and. &
                  spectral_derivative_ID == DERIVATIVE_CONTINUOUS_ID)                               ! for even grids, set the Nyquist Freq component to 0.0
              xi1st(1:3,i,j,k-grid3Offset) = cmplx(0.0_pReal,0.0_pReal,pReal)
            elsewhere
              xi1st(1:3,i,j,k-grid3Offset) = xi2nd(1:3,i,j,k-grid3Offset)
            endwhere
  enddo; enddo; enddo

  if(num%memory_efficient) then                                                                     ! allocate just single fourth order tensor
    allocate (gamma_hat(3,3,3,3,1,1,1), source = cmplx(0.0_pReal,0.0_pReal,pReal))
  else                                                                                              ! precalculation of gamma_hat field
    allocate (gamma_hat(3,3,3,3,grid1Red,grid(2),grid3), source = cmplx(0.0_pReal,0.0_pReal,pReal))
  endif

end subroutine spectral_utilities_init


!---------------------------------------------------------------------------------------------------
!> @brief updates reference stiffness and potentially precalculated gamma operator
!> @details Sets the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of an on-the-fly calculation, only the reference stiffness is updated.
!---------------------------------------------------------------------------------------------------
subroutine utilities_updateGamma(C)

  real(pReal), intent(in), dimension(3,3,3,3) :: C                                                  !< input stiffness to store as reference stiffness
  complex(pReal),              dimension(3,3) :: temp33_complex, xiDyad_cmplx
  real(pReal),                 dimension(6,6) :: A, A_inv
  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err

  C_ref = C

  if(.not. num%memory_efficient) then
    gamma_hat =  cmplx(0.0_pReal,0.0_pReal,pReal)                                                   ! for the singular point and any non invertible A
    do k = grid3Offset+1, grid3Offset+grid3; do j = 1, grid(2); do i = 1, grid1Red
      if (any([i,j,k] /= 1)) then                                                                   ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,j,k-grid3Offset))*xi1st(m,i,j,k-grid3Offset)
        forall(l = 1:3, m = 1:3) &
          temp33_complex(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal)*xiDyad_cmplx)
        A(1:3,1:3) = real(temp33_complex);  A(4:6,4:6) =   real(temp33_complex)
        A(1:3,4:6) = aimag(temp33_complex); A(4:6,1:3) = -aimag(temp33_complex)
        if (abs(math_det33(A(1:3,1:3))) > 1e-16) then
          call math_invert(A_inv, err, A)
          temp33_complex = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pReal)
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            gamma_hat(l,m,n,o,i,j,k-grid3Offset) = temp33_complex(l,n)* &
                                                   conjg(-xi1st(o,i,j,k-grid3Offset))*xi1st(m,i,j,k-grid3Offset)
        endif
      endif
    enddo; enddo; enddo
  endif

end subroutine utilities_updateGamma


!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier
!> @details Does an unweighted FFT transform from real to complex. Extra padding entries are set
! to 0.0
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorForward

  tensorField_real(1:3,1:3,grid(1)+1:grid1Red*2,:,:) = 0.0_pReal
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

end subroutine utilities_FFTtensorForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTtensorBackward

  call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
  tensorField_real = tensorField_real * wgt                                                         ! normalize the result by number of elements

end subroutine utilities_FFTtensorBackward

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in scalarField_real to scalarField_fourier
!> @details Does an unweighted FFT transform from real to complex. Extra padding entries are set
! to 0.0
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarForward

  scalarField_real(grid(1)+1:grid1Red*2,:,:) = 0.0_pReal
  call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)

end subroutine utilities_FFTscalarForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in scalarField_fourier to scalarField_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTscalarBackward

  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  scalarField_real = scalarField_real * wgt                                                         ! normalize the result by number of elements

end subroutine utilities_FFTscalarBackward


!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted FFT transform from real to complex. Extra padding entries are set
! to 0.0
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorForward

  vectorField_real(1:3,grid(1)+1:grid1Red*2,:,:) = 0.0_pReal
  call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)

end subroutine utilities_FFTvectorForward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorBackward

  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  vectorField_real = vectorField_real * wgt                                                         ! normalize the result by number of elements

end subroutine utilities_FFTvectorBackward


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real, ensuring that average value = fieldAim
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierGammaConvolution(fieldAim)

  real(pReal), intent(in), dimension(3,3) :: fieldAim                                               !< desired average value of the field after convolution
  complex(pReal),          dimension(3,3) :: temp33_complex, xiDyad_cmplx
  real(pReal),             dimension(6,6) :: A, A_inv

  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err


  print'(/,a)', ' ... doing gamma convolution ...............................................'
  flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation (mechanical equilibrium)
  memoryEfficient: if(num%memory_efficient) then
    do k = 1, grid3; do j = 1, grid(2); do i = 1, grid1Red
      if (any([i,j,k+grid3Offset] /= 1)) then                                                       ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,j,k))*xi1st(m,i,j,k)
        forall(l = 1:3, m = 1:3) &
          temp33_complex(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal)*xiDyad_cmplx)
        A(1:3,1:3) =  real(temp33_complex); A(4:6,4:6) =   real(temp33_complex)
        A(1:3,4:6) = aimag(temp33_complex); A(4:6,1:3) = -aimag(temp33_complex)
        if (abs(math_det33(A(1:3,1:3))) > 1e-16) then
          call math_invert(A_inv, err, A)
          temp33_complex = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pReal)
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            gamma_hat(l,m,n,o,1,1,1) =  temp33_complex(l,n)*conjg(-xi1st(o,i,j,k))*xi1st(m,i,j,k)
        else
          gamma_hat(1:3,1:3,1:3,1:3,1,1,1) = cmplx(0.0_pReal,0.0_pReal,pReal)
        endif
        forall(l = 1:3, m = 1:3) &
          temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3,1,1,1)*tensorField_fourier(1:3,1:3,i,j,k))
        tensorField_fourier(1:3,1:3,i,j,k) = temp33_Complex
      endif
    enddo; enddo; enddo
  else memoryEfficient
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid1Red
      forall(l = 1:3, m = 1:3) &
        temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3,i,j,k) * tensorField_fourier(1:3,1:3,i,j,k))
      tensorField_fourier(1:3,1:3,i,j,k) = temp33_Complex
    enddo; enddo; enddo
  endif memoryEfficient

  if (grid3Offset == 0) tensorField_fourier(1:3,1:3,1,1,1) = cmplx(fieldAim/wgt,0.0_pReal,pReal)

end subroutine utilities_fourierGammaConvolution


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution DamageGreenOp_hat * field_real
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierGreenConvolution(D_ref, mobility_ref, deltaT)

  real(pReal), dimension(3,3), intent(in) :: D_ref
  real(pReal),                 intent(in) :: mobility_ref, deltaT
  complex(pReal)                          :: GreenOp_hat
  integer                                 :: i, j, k

!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation
  do k = 1, grid3; do j = 1, grid(2) ;do i = 1, grid1Red
    GreenOp_hat =  cmplx(1.0_pReal,0.0_pReal,pReal)/ &
                   (cmplx(mobility_ref,0.0_pReal,pReal) + cmplx(deltaT,0.0_pReal)*&
                    sum(conjg(xi1st(1:3,i,j,k))* matmul(cmplx(D_ref,0.0_pReal),xi1st(1:3,i,j,k))))
    scalarField_fourier(i,j,k) = scalarField_fourier(i,j,k)*GreenOp_hat
  enddo; enddo; enddo

end subroutine utilities_fourierGreenConvolution


!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS()

  integer :: i, j, k, ierr
  complex(pReal), dimension(3)   :: rescaledGeom

  print'(/,a)', ' ... calculating divergence ................................................'
  flush(IO_STDOUT)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pReal)

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
  utilities_divergenceRMS = 0.0_pReal
  do k = 1, grid3; do j = 1, grid(2)
    do i = 2, grid1Red -1                                                                           ! Has somewhere a conj. complex counterpart. Therefore count it twice.
      utilities_divergenceRMS = utilities_divergenceRMS &
            + 2.0_pReal*(sum (real(matmul(tensorField_fourier(1:3,1:3,i,j,k),&                      ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                          conjg(-xi1st(1:3,i,j,k))*rescaledGeom))**2.0_pReal)&      ! --> sum squared L_2 norm of vector
                        +sum(aimag(matmul(tensorField_fourier(1:3,1:3,i,j,k),&
                                          conjg(-xi1st(1:3,i,j,k))*rescaledGeom))**2.0_pReal))
    enddo
    utilities_divergenceRMS = utilities_divergenceRMS &                                             ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart (if grid(1) /= 1)
               + sum( real(matmul(tensorField_fourier(1:3,1:3,1       ,j,k), &
                                  conjg(-xi1st(1:3,1,j,k))*rescaledGeom))**2.0_pReal) &
               + sum(aimag(matmul(tensorField_fourier(1:3,1:3,1       ,j,k), &
                                  conjg(-xi1st(1:3,1,j,k))*rescaledGeom))**2.0_pReal) &
               + sum( real(matmul(tensorField_fourier(1:3,1:3,grid1Red,j,k), &
                                  conjg(-xi1st(1:3,grid1Red,j,k))*rescaledGeom))**2.0_pReal) &
               + sum(aimag(matmul(tensorField_fourier(1:3,1:3,grid1Red,j,k), &
                                  conjg(-xi1st(1:3,grid1Red,j,k))*rescaledGeom))**2.0_pReal)
  enddo; enddo
  if(grid(1) == 1) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pReal                    ! counted twice in case of grid(1) == 1
  call MPI_Allreduce(MPI_IN_PLACE,utilities_divergenceRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  if(ierr /=0) error stop 'MPI error'
  utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                     ! RMS in real space calculated with Parsevals theorem from Fourier space

end function utilities_divergenceRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculate max of curl of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_curlRMS()

  integer  ::  i, j, k, l, ierr
  complex(pReal), dimension(3,3) :: curl_fourier
  complex(pReal), dimension(3)   :: rescaledGeom

  print'(/,a)', ' ... calculating curl ......................................................'
  flush(IO_STDOUT)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pReal)

!--------------------------------------------------------------------------------------------------
! calculating max curl criterion in Fourier space
  utilities_curlRMS = 0.0_pReal

  do k = 1, grid3; do j = 1, grid(2);
    do i = 2, grid1Red - 1
      do l = 1, 3
        curl_fourier(l,1) = (+tensorField_fourier(l,3,i,j,k)*xi1st(2,i,j,k)*rescaledGeom(2) &
                             -tensorField_fourier(l,2,i,j,k)*xi1st(3,i,j,k)*rescaledGeom(3))
        curl_fourier(l,2) = (+tensorField_fourier(l,1,i,j,k)*xi1st(3,i,j,k)*rescaledGeom(3) &
                             -tensorField_fourier(l,3,i,j,k)*xi1st(1,i,j,k)*rescaledGeom(1))
        curl_fourier(l,3) = (+tensorField_fourier(l,2,i,j,k)*xi1st(1,i,j,k)*rescaledGeom(1) &
                             -tensorField_fourier(l,1,i,j,k)*xi1st(2,i,j,k)*rescaledGeom(2))
      enddo
      utilities_curlRMS = utilities_curlRMS &
                        +2.0_pReal*sum(real(curl_fourier)**2.0_pReal+aimag(curl_fourier)**2.0_pReal)! Has somewhere a conj. complex counterpart. Therefore count it twice.
    enddo
    do l = 1, 3
       curl_fourier = (+tensorField_fourier(l,3,1,j,k)*xi1st(2,1,j,k)*rescaledGeom(2) &
                       -tensorField_fourier(l,2,1,j,k)*xi1st(3,1,j,k)*rescaledGeom(3))
       curl_fourier = (+tensorField_fourier(l,1,1,j,k)*xi1st(3,1,j,k)*rescaledGeom(3) &
                       -tensorField_fourier(l,3,1,j,k)*xi1st(1,1,j,k)*rescaledGeom(1))
       curl_fourier = (+tensorField_fourier(l,2,1,j,k)*xi1st(1,1,j,k)*rescaledGeom(1) &
                       -tensorField_fourier(l,1,1,j,k)*xi1st(2,1,j,k)*rescaledGeom(2))
    enddo
    utilities_curlRMS = utilities_curlRMS &
                      + sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)         ! this layer (DC) does not have a conjugate complex counterpart (if grid(1) /= 1)
    do l = 1, 3
      curl_fourier = (+tensorField_fourier(l,3,grid1Red,j,k)*xi1st(2,grid1Red,j,k)*rescaledGeom(2) &
                      -tensorField_fourier(l,2,grid1Red,j,k)*xi1st(3,grid1Red,j,k)*rescaledGeom(3))
      curl_fourier = (+tensorField_fourier(l,1,grid1Red,j,k)*xi1st(3,grid1Red,j,k)*rescaledGeom(3) &
                      -tensorField_fourier(l,3,grid1Red,j,k)*xi1st(1,grid1Red,j,k)*rescaledGeom(1))
      curl_fourier = (+tensorField_fourier(l,2,grid1Red,j,k)*xi1st(1,grid1Red,j,k)*rescaledGeom(1) &
                      -tensorField_fourier(l,1,grid1Red,j,k)*xi1st(2,grid1Red,j,k)*rescaledGeom(2))
    enddo
    utilities_curlRMS = utilities_curlRMS &
                      + sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)         ! this layer (Nyquist) does not have a conjugate complex counterpart (if grid(1) /= 1)
  enddo; enddo

  call MPI_Allreduce(MPI_IN_PLACE,utilities_curlRMS,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  if(ierr /=0) error stop 'MPI error'
  utilities_curlRMS = sqrt(utilities_curlRMS) * wgt
  if(grid(1) == 1) utilities_curlRMS = utilities_curlRMS * 0.5_pReal                                ! counted twice in case of grid(1) == 1

end function utilities_curlRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculates mask compliance tensor used to adjust F to fullfill stress BC
!--------------------------------------------------------------------------------------------------
function utilities_maskedCompliance(rot_BC,mask_stress,C)

  real(pReal),                dimension(3,3,3,3) :: utilities_maskedCompliance                      !< masked compliance
  real(pReal),    intent(in), dimension(3,3,3,3) :: C                                               !< current average stiffness
  type(rotation), intent(in)                     :: rot_BC                                          !< rotation of load frame
  logical,        intent(in), dimension(3,3)     :: mask_stress                                     !< mask of stress BC

  integer :: i, j
  logical, dimension(9)   :: mask_stressVector
  logical, dimension(9,9) :: mask
  real(pReal), dimension(9,9) :: temp99_real
  integer :: size_reduced = 0
  real(pReal),              dimension(:,:), allocatable ::  &
    s_reduced, &                                                                                    !< reduced compliance matrix (depending on number of stress BC)
    c_reduced, &                                                                                    !< reduced stiffness (depending on number of stress BC)
    sTimesC                                                                                         !< temp variable to check inversion
  logical :: errmatinv
  character(len=pStringLen):: formatString

  mask_stressVector = reshape(transpose(mask_stress), [9])
  size_reduced = count(mask_stressVector)
  if(size_reduced > 0) then
    temp99_real = math_3333to99(rot_BC%rotate(C))

    if(debugGeneral) then
      print'(/,a)', ' ... updating masked compliance ............................................'
      print'(/,a,/,8(9(2x,f12.7,1x)/),9(2x,f12.7,1x))', &
        ' Stiffness C (load) / GPa =', transpose(temp99_Real)*1.0e-9_pReal
      flush(IO_STDOUT)
    endif

    do i = 1,9; do j = 1,9
      mask(i,j) = mask_stressVector(i) .and. mask_stressVector(j)
    enddo; enddo
    c_reduced = reshape(pack(temp99_Real,mask),[size_reduced,size_reduced])

    allocate(s_reduced,mold = c_reduced)
    call math_invert(s_reduced, errmatinv, c_reduced)                                               ! invert reduced stiffness
    if (any(IEEE_is_NaN(s_reduced))) errmatinv = .true.

!--------------------------------------------------------------------------------------------------
! check if inversion was successful
    sTimesC = matmul(c_reduced,s_reduced)
    errmatinv = errmatinv .or. any(dNeq(sTimesC,math_eye(size_reduced),1.0e-12_pReal))
    if (debugGeneral .or. errmatinv) then
      write(formatString, '(i2)') size_reduced
      formatString = '(/,a,/,'//trim(formatString)//'('//trim(formatString)//'(2x,es9.2,1x)/))'
      print trim(formatString), ' C * S (load) ', transpose(matmul(c_reduced,s_reduced))
      print trim(formatString), ' S (load) ', transpose(s_reduced)
      if(errmatinv) error stop 'matrix inversion error'
    endif
    temp99_real = reshape(unpack(reshape(s_reduced,[size_reduced**2]),reshape(mask,[81]),0.0_pReal),[9,9])
  else
    temp99_real = 0.0_pReal
  endif

  utilities_maskedCompliance = math_99to3333(temp99_Real)

  if(debugGeneral) then
    print'(/,a,/,9(9(2x,f10.5,1x)/),9(2x,f10.5,1x))', &
      ' Masked Compliance (load) * GPa =', transpose(temp99_Real)*1.0e9_pReal
    flush(IO_STDOUT)
  endif

end function utilities_maskedCompliance


!--------------------------------------------------------------------------------------------------
!> @brief calculate scalar gradient in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierScalarGradient()

  integer :: i, j, k

  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid1Red
    vectorField_fourier(1:3,i,j,k) = scalarField_fourier(i,j,k)*xi1st(1:3,i,j,k)                    ! ToDo: no -conjg?
  enddo; enddo; enddo

end subroutine utilities_fourierScalarGradient


!--------------------------------------------------------------------------------------------------
!> @brief calculate vector divergence in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierVectorDivergence()

  integer :: i, j, k

  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid1Red
    scalarField_fourier(i,j,k) = sum(vectorField_fourier(1:3,i,j,k)*conjg(-xi1st(1:3,i,j,k)))
  enddo; enddo; enddo

end subroutine utilities_fourierVectorDivergence


!--------------------------------------------------------------------------------------------------
!> @brief calculate vector gradient in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierVectorGradient()

  integer :: i, j, k, m, n

  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid1Red
    do m = 1, 3; do n = 1, 3
      tensorField_fourier(m,n,i,j,k) = vectorField_fourier(m,i,j,k)*xi1st(n,i,j,k)
    enddo; enddo
  enddo; enddo; enddo

end subroutine utilities_fourierVectorGradient


!--------------------------------------------------------------------------------------------------
!> @brief calculate tensor divergence in fourier field
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierTensorDivergence()

  integer :: i, j, k

  do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid1Red
    vectorField_fourier(:,i,j,k) = matmul(tensorField_fourier(:,:,i,j,k),conjg(-xi1st(:,i,j,k)))
  enddo; enddo; enddo

end subroutine utilities_fourierTensorDivergence


!--------------------------------------------------------------------------------------------------
!> @brief calculate constitutive response from homogenization_F0 to F during timeinc
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(P,P_av,C_volAvg,C_minmaxAvg,&
                                          F,timeinc,rotation_BC)

  real(pReal),    intent(out), dimension(3,3,3,3)                   :: C_volAvg, C_minmaxAvg        !< average stiffness
  real(pReal),    intent(out), dimension(3,3)                       :: P_av                         !< average PK stress
  real(pReal),    intent(out), dimension(3,3,grid(1),grid(2),grid3) :: P                            !< PK stress
  real(pReal),    intent(in),  dimension(3,3,grid(1),grid(2),grid3) :: F                            !< deformation gradient target
  real(pReal),    intent(in)                                        :: timeinc                      !< loading time
  type(rotation), intent(in),  optional                             :: rotation_BC                  !< rotation of load frame


  integer :: &
    i,ierr
  real(pReal), dimension(3,3,3,3) :: dPdF_max,      dPdF_min
  real(pReal)                     :: dPdF_norm_max, dPdF_norm_min
  real(pReal), dimension(2) :: valueAndRank                                                         !< pair of min/max norm of dPdF to synchronize min/max of dPdF

  print'(/,a)', ' ... evaluating constitutive response ......................................'
  flush(IO_STDOUT)

  homogenization_F  = reshape(F,[3,3,product(grid(1:2))*grid3])                                     ! set materialpoint target F to estimated field

  call materialpoint_stressAndItsTangent(timeinc,[1,1],[1,product(grid(1:2))*grid3])                ! calculate P field

  P = reshape(homogenization_P, [3,3,grid(1),grid(2),grid3])
  P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt                                                   ! average of P
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  if (debugRotation) print'(/,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    ' Piola--Kirchhoff stress (lab) / MPa =', transpose(P_av)*1.e-6_pReal
  if(present(rotation_BC)) P_av = rotation_BC%rotate(P_av)
  print'(/,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    ' Piola--Kirchhoff stress       / MPa =', transpose(P_av)*1.e-6_pReal
  flush(IO_STDOUT)

  dPdF_max = 0.0_pReal
  dPdF_norm_max = 0.0_pReal
  dPdF_min = huge(1.0_pReal)
  dPdF_norm_min = huge(1.0_pReal)
  do i = 1, product(grid(1:2))*grid3
    if (dPdF_norm_max < sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2.0_pReal)) then
      dPdF_max = homogenization_dPdF(1:3,1:3,1:3,1:3,i)
      dPdF_norm_max = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2.0_pReal)
    endif
    if (dPdF_norm_min > sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2.0_pReal)) then
      dPdF_min = homogenization_dPdF(1:3,1:3,1:3,1:3,i)
      dPdF_norm_min = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2.0_pReal)
    endif
  end do

  valueAndRank = [dPdF_norm_max,real(worldrank,pReal)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'
  call MPI_Bcast(dPdF_max,81,MPI_DOUBLE,int(valueAndRank(2)),PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'

  valueAndRank = [dPdF_norm_min,real(worldrank,pReal)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'
  call MPI_Bcast(dPdF_min,81,MPI_DOUBLE,int(valueAndRank(2)),PETSC_COMM_WORLD, ierr)
  if (ierr /= 0) error stop 'MPI error'

  C_minmaxAvg = 0.5_pReal*(dPdF_max + dPdF_min)

  C_volAvg = sum(homogenization_dPdF,dim=5)
  call MPI_Allreduce(MPI_IN_PLACE,C_volAvg,81,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
  if (ierr /= 0) error stop 'MPI error'
  C_volAvg = C_volAvg * wgt


end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief calculates forward rate, either guessing or just add delta/timeinc
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(heterogeneous,field0,field,dt,avRate)

  real(pReal), intent(in), dimension(3,3) :: &
    avRate                                                                                          !< homogeneous addon
  real(pReal), intent(in) :: &
    dt                                                                                              !< timeinc between field0 and field
  logical, intent(in) :: &
    heterogeneous                                                                                   !< calculate field of rates
  real(pReal), intent(in), dimension(3,3,grid(1),grid(2),grid3) :: &
    field0, &                                                                                       !< data of previous step
    field                                                                                           !< data of current step
  real(pReal),             dimension(3,3,grid(1),grid(2),grid3) :: &
    utilities_calculateRate

  if (heterogeneous) then
    utilities_calculateRate = (field-field0) / dt
  else
    utilities_calculateRate = spread(spread(spread(avRate,3,grid(1)),4,grid(2)),5,grid3)
  endif

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given,
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
function utilities_forwardField(timeinc,field_lastInc,rate,aim)

  real(pReal), intent(in) :: &
    timeinc                                                                                         !< timeinc of current step
  real(pReal), intent(in),           dimension(3,3,grid(1),grid(2),grid3) :: &
    field_lastInc, &                                                                                !< initial field
    rate                                                                                            !< rate by which to forward
  real(pReal), intent(in), optional, dimension(3,3) :: &
    aim                                                                                             !< average field value aim
  real(pReal),                       dimension(3,3,grid(1),grid(2),grid3) :: &
    utilities_forwardField
  real(pReal),                       dimension(3,3)                       :: fieldDiff              !< <a + adot*t> - aim
  PetscErrorCode :: ierr

  utilities_forwardField = field_lastInc + rate*timeinc
  if (present(aim)) then                                                                            !< correct to match average
    fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt
    call MPI_Allreduce(MPI_IN_PLACE,fieldDiff,9,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD,ierr)
    fieldDiff = fieldDiff - aim
    utilities_forwardField = utilities_forwardField - &
                     spread(spread(spread(fieldDiff,3,grid(1)),4,grid(2)),5,grid3)
  endif

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief calculates filter for fourier convolution depending on type given in numerics.config
!> @details this is the full operator to calculate derivatives, i.e. 2 \pi i k for the
! standard approach
!--------------------------------------------------------------------------------------------------
pure function utilities_getFreqDerivative(k_s)

  integer, intent(in),  dimension(3) :: k_s                                                         !< indices of frequency
  complex(pReal),       dimension(3) :: utilities_getFreqDerivative

  select case (spectral_derivative_ID)
    case (DERIVATIVE_CONTINUOUS_ID)
      utilities_getFreqDerivative = cmplx(0.0_pReal, 2.0_pReal*PI*real(k_s,pReal)/geomSize,pReal)

    case (DERIVATIVE_CENTRAL_DIFF_ID)
      utilities_getFreqDerivative = cmplx(0.0_pReal, sin(2.0_pReal*PI*real(k_s,pReal)/real(grid,pReal)), pReal)/ &
                                    cmplx(2.0_pReal*geomSize/real(grid,pReal), 0.0_pReal, pReal)

    case (DERIVATIVE_FWBW_DIFF_ID)
      utilities_getFreqDerivative(1) = &
                               cmplx(cos(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)) - 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(1)/real(grid(1),pReal), 0.0_pReal, pReal)
      utilities_getFreqDerivative(2) = &
                               cmplx(cos(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)) - 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(2)/real(grid(2),pReal), 0.0_pReal, pReal)
      utilities_getFreqDerivative(3) = &
                               cmplx(cos(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(1),pReal)/real(grid(1),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)) + 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(2),pReal)/real(grid(2),pReal)), pReal)* &
                               cmplx(cos(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)) - 1.0_pReal, &
                                     sin(2.0_pReal*PI*real(k_s(3),pReal)/real(grid(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(3)/real(grid(3),pReal), 0.0_pReal, pReal)
  end select

end function utilities_getFreqDerivative


!--------------------------------------------------------------------------------------------------
!> @brief calculate coordinates in current configuration for given defgrad field
! using integration in Fourier space. Similar as in mesh.f90, but using data already defined for
! convolution
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateCoords(F)

  real(pReal),   dimension(3,3,grid(1),grid(2),grid3), intent(in) :: F
  real(pReal),   dimension(3,  grid(1),grid(2),grid3)             :: IPcoords
  real(pReal),   dimension(3,  grid(1),grid(2),grid3+2)           :: IPfluct_padded                  ! Fluctuations of cell center displacement (padded along z for MPI)
  real(pReal),   dimension(3,  grid(1)+1,grid(2)+1,grid3+1)       :: nodeCoords
  integer :: &
    i,j,k,n, &
    rank_t, rank_b, &
    c, &
    ierr
  integer, dimension(4) :: request
  integer, dimension(MPI_STATUS_SIZE,4) :: status
  real(pReal),   dimension(3)   :: step
  real(pReal),   dimension(3,3) :: Favg
  integer,       dimension(3) :: me
  integer, dimension(3,8) :: &
    neighbor = reshape([ &
                        0, 0, 0, &
                        1, 0, 0, &
                        1, 1, 0, &
                        0, 1, 0, &
                        0, 0, 1, &
                        1, 0, 1, &
                        1, 1, 1, &
                        0, 1, 1  ], [3,8])

  step = geomSize/real(grid, pReal)
 !--------------------------------------------------------------------------------------------------
 ! integration in Fourier space to get fluctuations of cell center discplacements
  tensorField_real(1:3,1:3,1:grid(1),1:grid(2),1:grid3) = F
  call utilities_FFTtensorForward()

  do k = 1, grid3; do j = 1, grid(2); do i = 1, grid1Red
    if(any([i,j,k+grid3Offset] /= 1)) then
      vectorField_fourier(1:3,i,j,k) = matmul(tensorField_fourier(1:3,1:3,i,j,k),xi2nd(1:3,i,j,k)) &
                                     / sum(conjg(-xi2nd(1:3,i,j,k))*xi2nd(1:3,i,j,k)) * cmplx(wgt,0.0,pReal)
    else
      vectorField_fourier(1:3,i,j,k) = cmplx(0.0,0.0,pReal)
    endif
  enddo; enddo; enddo

  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)

 !--------------------------------------------------------------------------------------------------
 ! average F
  if (grid3Offset == 0) Favg = real(tensorField_fourier(1:3,1:3,1,1,1),pReal)*wgt
  call MPI_Bcast(Favg,9,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
  if(ierr /=0) error stop 'MPI error'

 !--------------------------------------------------------------------------------------------------
 ! pad cell center fluctuations along z-direction (needed when running MPI simulation)
  IPfluct_padded(1:3,1:grid(1),1:grid(2),2:grid3+1) = vectorField_real(1:3,1:grid(1),1:grid(2),1:grid3)
  c = product(shape(IPfluct_padded(:,:,:,1)))                                                        !< amount of data to transfer
  rank_t = modulo(worldrank+1,worldsize)
  rank_b = modulo(worldrank-1,worldsize)

  ! send bottom layer to process below
  call MPI_Isend(IPfluct_padded(:,:,:,2),      c,MPI_DOUBLE,rank_b,0,PETSC_COMM_WORLD,request(1),ierr)
  if(ierr /=0) error stop 'MPI error'
  call MPI_Irecv(IPfluct_padded(:,:,:,grid3+2),c,MPI_DOUBLE,rank_t,0,PETSC_COMM_WORLD,request(2),ierr)
  if(ierr /=0) error stop 'MPI error'

  ! send top layer to process above
  call MPI_Isend(IPfluct_padded(:,:,:,grid3+1),c,MPI_DOUBLE,rank_t,1,PETSC_COMM_WORLD,request(3),ierr)
  if(ierr /=0) error stop 'MPI error'
  call MPI_Irecv(IPfluct_padded(:,:,:,1),      c,MPI_DOUBLE,rank_b,1,PETSC_COMM_WORLD,request(4),ierr)
  if(ierr /=0) error stop 'MPI error'

  call MPI_Waitall(4,request,status,ierr)
  if(ierr /=0) error stop 'MPI error'
  if(any(status(MPI_ERROR,:) /= 0)) error stop 'MPI error'

 !--------------------------------------------------------------------------------------------------
 ! calculate nodal displacements
  nodeCoords = 0.0_pReal
  do k = 0,grid3; do j = 0,grid(2); do i = 0,grid(1)
    nodeCoords(1:3,i+1,j+1,k+1) = matmul(Favg,step*(real([i,j,k+grid3Offset],pReal)))
    averageFluct: do n = 1,8
      me = [i+neighbor(1,n),j+neighbor(2,n),k+neighbor(3,n)]
      nodeCoords(1:3,i+1,j+1,k+1) = nodeCoords(1:3,i+1,j+1,k+1) &
                                  + IPfluct_padded(1:3,modulo(me(1)-1,grid(1))+1,modulo(me(2)-1,grid(2))+1,me(3)+1)*0.125_pReal
    enddo averageFluct
  enddo; enddo; enddo

 !--------------------------------------------------------------------------------------------------
 ! calculate cell center displacements
  do k = 1,grid3; do j = 1,grid(2); do i = 1,grid(1)
    IPcoords(1:3,i,j,k) = vectorField_real(1:3,i,j,k) &
                        + matmul(Favg,step*(real([i,j,k+grid3Offset],pReal)-0.5_pReal))
  enddo; enddo; enddo

  call discretization_setNodeCoords(reshape(NodeCoords,[3,(grid(1)+1)*(grid(2)+1)*(grid3+1)]))
  call discretization_setIPcoords  (reshape(IPcoords,  [3,grid(1)*grid(2)*grid3]))

end subroutine utilities_updateCoords


!---------------------------------------------------------------------------------------------------
!> @brief Write out the current reference stiffness for restart.
!---------------------------------------------------------------------------------------------------
subroutine utilities_saveReferenceStiffness

  integer :: &
    fileUnit,ierr

  if (worldrank == 0) then
    print'(a)', ' writing reference stiffness data required for restart to file'; flush(IO_STDOUT)
    open(newunit=fileUnit, file=getSolverJobName()//'.C_ref',&
         status='replace',access='stream',action='write',iostat=ierr)
    if(ierr /=0) call IO_error(100,ext_msg='could not open file '//getSolverJobName()//'.C_ref')
    write(fileUnit) C_ref
    close(fileUnit)
  endif

end subroutine utilities_saveReferenceStiffness

end module spectral_utilities
