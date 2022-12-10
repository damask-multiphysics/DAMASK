!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the different spectral solver variants
!--------------------------------------------------------------------------------------------------
module spectral_utilities
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif
  use FFTW3

  use prec
  use CLI
  use parallelization
  use math
  use rotations
  use IO
  use config
  use discretization_grid
  use discretization
  use homogenization


#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  implicit none(type,external)
#else
  implicit none
#endif
  private

!--------------------------------------------------------------------------------------------------
! grid related information
  real(pReal), protected, public               :: wgt                                               !< weighting factor 1/Nelems
  real(pReal), protected, public, dimension(3) :: scaledGeomSize                                    !< scaled geometry size for calculation of divergence
  integer :: &
    cells1Red, &                                                                                    !< cells(1)/2+1
    cells2, &                                                                                       !< (local) cells in 2nd direction
    cells2Offset                                                                                    !< (local) cells offset in 2nd direction

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW

  real(C_DOUBLE),            dimension(:,:,:,:,:),     pointer     :: tensorField_real              !< tensor field in real space
  real(C_DOUBLE),            dimension(:,:,:,:),       pointer     :: vectorField_real              !< vector field in real space
  real(C_DOUBLE),            dimension(:,:,:),         pointer     :: scalarField_real              !< scalar field in real space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:,:),     pointer     :: tensorField_fourier           !< tensor field in Fourier space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:),       pointer     :: vectorField_fourier           !< vector field in Fourier space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:),         pointer     :: scalarField_fourier           !< scalar field in Fourier space
  complex(pReal),            dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat                     !< gamma operator (field) for spectral method
  complex(pReal),            dimension(:,:,:,:),       allocatable :: xi1st                         !< wave vector field for first derivatives
  complex(pReal),            dimension(:,:,:,:),       allocatable :: xi2nd                         !< wave vector field for second derivatives
  real(pReal),               dimension(3,3,3,3)                    :: C_ref                         !< mechanic reference stiffness


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
    logical,     dimension(3,3)   :: mask   = .true.
    character(len=:), allocatable :: myType
  end type tBoundaryCondition

  type, public :: tSolutionParams
    real(pReal), dimension(3,3) :: stress_BC
    logical, dimension(3,3)     :: stress_mask
    type(tRotation)             :: rotation_BC
    real(pReal) :: Delta_t
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
    utilities_GammaConvolution, &
    utilities_GreenConvolution, &
    utilities_divergenceRMS, &
    utilities_curlRMS, &
    utilities_scalarGradient, &
    utilities_vectorDivergence, &
    utilities_maskedCompliance, &
    utilities_constitutiveResponse, &
    utilities_calculateRate, &
    utilities_forwardField, &
    utilities_updateCoords, &
    utilities_saveReferenceStiffness

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags, create plans for FFTW
!> @details Sets the debug levels for general, divergence, restart, and FFTW from the bitwise coding
!> provided by the debug module to logicals.
!> Allocate all fields used by FFTW and create the corresponding plans depending on the debug
!> level chosen.
!> Initializes FFTW.
!--------------------------------------------------------------------------------------------------
subroutine spectral_utilities_init()

  PetscErrorCode :: err_PETSc
  integer        :: i, j, k, &
    FFTW_planner_flag
  integer, dimension(3) :: k_s
  type(C_PTR) :: &
    tensorField, &                                                                                  !< tensor data for FFTW in real and Fourier space (in-place)
    vectorField, &                                                                                  !< vector data for FFTW in real and Fourier space (in-place)
    scalarField                                                                                     !< scalar data for FFTW in real and Fourier space (in-place)
  integer(C_INTPTR_T), dimension(3) :: cellsFFTW
  integer(C_INTPTR_T) :: N, &
    cells3FFTW, &                                                                                   !< # of cells in 3. dim on current process in real space
    cells3_offset, &                                                                                !< offset for cells in 3. dim on current process in real space
    cells2FFTW, &                                                                                   !< # of cells in 2. dim on current process in Fourier space
    cells2_offset                                                                                   !< offset for cells in 2. dim on curren process in Fourier space
  integer(C_INTPTR_T), parameter :: &
    vectorSize = 3_C_INTPTR_T, &
    tensorSize = 9_C_INTPTR_T
  character(len=*), parameter :: &
    PETSCDEBUG = ' -snes_view -snes_monitor '
  type(tDict) , pointer :: &
    num_grid
  type(tList) , pointer :: &
    debug_grid


  print'(/,1x,a)', '<<<+-  spectral_utilities init  -+>>>'

  print'(/,1x,a)', 'M. Diehl, Diploma Thesis TU München, 2010'
  print'(  1x,a)', 'https://doi.org/10.13140/2.1.3234.3840'//IO_EOL

  print'(  1x,a)', 'P. Eisenlohr et al., International Journal of Plasticity 46:37–53, 2013'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2012.09.012'//IO_EOL

  print'(  1x,a)', 'P. Shanthraj et al., International Journal of Plasticity 66:31–45, 2015'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2014.02.006'//IO_EOL

  print'(  1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'

!--------------------------------------------------------------------------------------------------
! set debugging parameters
  num_grid        => config_numerics%get_dict('grid',defaultVal=emptyDict)

  debug_grid      => config_debug%get_List('grid',defaultVal=emptyList)
  debugGeneral    =  debug_grid%contains('basic')
  debugRotation   =  debug_grid%contains('rotation')
  debugPETSc      =  debug_grid%contains('PETSc')

  if (debugPETSc) print'(3(/,1x,a),/)', &
                 'Initializing PETSc with debug options: ', &
                 trim(PETScDebug), &
                 'add more using the "PETSc_options" keyword in numerics.yaml'
  flush(IO_STDOUT)

  call PetscOptionsClear(PETSC_NULL_OPTIONS,err_PETSc)
  CHKERRQ(err_PETSc)
  if (debugPETSc) call PetscOptionsInsertString(PETSC_NULL_OPTIONS,trim(PETSCDEBUG),err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                                num_grid%get_asString('PETSc_options',defaultVal=''),err_PETSc)
  CHKERRQ(err_PETSc)

  cells1Red = cells(1)/2 + 1
  wgt = real(product(cells),pReal)**(-1)

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
    end do
  elseif (num%divergence_correction == 2) then
    do j = 1, 3
     if (      j /= int(minloc(geomSize/real(cells,pReal),1)) &
         .and. j /= int(maxloc(geomSize/real(cells,pReal),1))) &
       scaledGeomSize = geomSize/geomSize(j)*real(cells(j),pReal)
    end do
  else
    scaledGeomSize = geomSize
  end if

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
      call IO_warning(47,'using default FFTW_MEASURE instead of "'//trim(num_grid%get_asString('fftw_plan_mode'))//'"')
      FFTW_planner_flag = FFTW_MEASURE
  end select

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
  if (pReal /= C_DOUBLE .or. kind(1) /= C_INT) error stop 'C and Fortran datatypes do not match'
  call fftw_set_timelimit(num_grid%get_asFloat('fftw_timelimit',defaultVal=300.0_pReal))

  print'(/,1x,a)', 'FFTW initialized'; flush(IO_STDOUT)

  cellsFFTW = int(cells,C_INTPTR_T)

  N = fftw_mpi_local_size_many_transposed(3,[cellsFFTW(3),cellsFFTW(2),int(cells1Red,C_INTPTR_T)], &
                                          tensorSize,FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,PETSC_COMM_WORLD, &
                                          cells3FFTW,cells3_offset,cells2FFTW,cells2_offset)
  cells2 = int(cells2FFTW)
  cells2Offset = int(cells2_offset)
  if (int(cells3FFTW) /= cells3) error stop 'domain decomposition mismatch (tensor, real space)'
  tensorField = fftw_alloc_complex(N)
  call c_f_pointer(tensorField,tensorField_real, &
                   [3_C_INTPTR_T,3_C_INTPTR_T,int(cells1Red*2,C_INTPTR_T),cellsFFTW(2),cells3FFTW])
  call c_f_pointer(tensorField,tensorField_fourier, &
                   [3_C_INTPTR_T,3_C_INTPTR_T,int(cells1Red,  C_INTPTR_T),cellsFFTW(3),cells2FFTW])


  N = fftw_mpi_local_size_many_transposed(3,[cellsFFTW(3),cellsFFTW(2),int(cells1Red,C_INTPTR_T)], &
                                          vectorSize,FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK,PETSC_COMM_WORLD, &
                                          cells3FFTW,cells3_offset,cells2FFTW,cells2_offset)
  if (int(cells3FFTW) /= cells3) error stop 'domain decomposition mismatch (vector, real space)'
  if (int(cells2FFTW) /= cells2) error stop 'domain decomposition mismatch (vector, Fourier space)'
  vectorField = fftw_alloc_complex(N)
  call c_f_pointer(vectorField,vectorField_real, &
                   [3_C_INTPTR_T,int(cells1Red*2,C_INTPTR_T),cellsFFTW(2),cells3FFTW])
  call c_f_pointer(vectorField,vectorField_fourier, &
                   [3_C_INTPTR_T,int(cells1Red,  C_INTPTR_T),cellsFFTW(3),cells2FFTW])

  N = fftw_mpi_local_size_3d_transposed(cellsFFTW(3),cellsFFTW(2),int(cells1Red,C_INTPTR_T), &
                                        PETSC_COMM_WORLD,cells3FFTW,cells3_offset,cells2FFTW,cells2_offset)
  if (int(cells3FFTW) /= cells3) error stop 'domain decomposition mismatch (scalar, real space)'
  if (int(cells2FFTW) /= cells2) error stop 'domain decomposition mismatch (scalar, Fourier space)'
  scalarField = fftw_alloc_complex(N)
  call c_f_pointer(scalarField,scalarField_real, &
                   [int(cells1Red*2,C_INTPTR_T),cellsFFTW(2),cells3FFTW])
  call c_f_pointer(scalarField,scalarField_fourier, &
                   [int(cells1Red,  C_INTPTR_T),cellsFFTW(3),cells2FFTW])

!--------------------------------------------------------------------------------------------------
! allocation
  allocate (xi1st (3,cells1Red,cells(3),cells2),source = cmplx(0.0_pReal,0.0_pReal,pReal))          ! frequencies for first derivatives, only half the size for first dimension
  allocate (xi2nd (3,cells1Red,cells(3),cells2),source = cmplx(0.0_pReal,0.0_pReal,pReal))          ! frequencies for second derivatives, only half the size for first dimension

!--------------------------------------------------------------------------------------------------
! tensor MPI fftw plans
  planTensorForth = fftw_mpi_plan_many_dft_r2c(3,cellsFFTW(3:1:-1),tensorSize, &
                                               FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                                               tensorField_real,tensorField_fourier, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_OUT)
  if (.not. c_associated(planTensorForth)) error stop 'FFTW error r2c tensor'
  planTensorBack  = fftw_mpi_plan_many_dft_c2r(3,cellsFFTW(3:1:-1),tensorSize, &
                                               FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                                               tensorField_fourier,tensorField_real, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_IN)
  if (.not. c_associated(planTensorBack))  error stop 'FFTW error c2r tensor'

!--------------------------------------------------------------------------------------------------
! vector MPI fftw plans
  planVectorForth = fftw_mpi_plan_many_dft_r2c(3,cellsFFTW(3:1:-1),vectorSize, &
                                               FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                                               vectorField_real,vectorField_fourier, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_OUT)
  if (.not. c_associated(planVectorForth)) error stop 'FFTW error r2c vector'
  planVectorBack  = fftw_mpi_plan_many_dft_c2r(3,cellsFFTW(3:1:-1),vectorSize, &
                                               FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                                               vectorField_fourier,vectorField_real, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_IN)
  if (.not. c_associated(planVectorBack))  error stop 'FFTW error c2r vector'

!--------------------------------------------------------------------------------------------------
! scalar MPI fftw plans
  planScalarForth = fftw_mpi_plan_dft_r2c_3d(cellsFFTW(3),cellsFFTW(2),cellsFFTW(1), &
                                             scalarField_real,scalarField_fourier, &
                                             PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_OUT)
  if (.not. c_associated(planScalarForth)) error stop 'FFTW error r2c scalar'
  planScalarBack  = fftw_mpi_plan_dft_c2r_3d(cellsFFTW(3),cellsFFTW(2),cellsFFTW(1), &
                                             scalarField_fourier,scalarField_real, &
                                             PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_IN)
  if (.not. c_associated(planScalarBack))  error stop 'FFTW error c2r scalar'

!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
  do j = cells2Offset+1, cells2Offset+cells2
    k_s(2) = j - 1
    if (j > cells(2)/2 + 1) k_s(2) = k_s(2) - cells(2)                                              ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
      do k = 1, cells(3)
        k_s(3) = k - 1
        if (k > cells(3)/2 + 1) k_s(3) = k_s(3) - cells(3)                                          ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
          do i = 1, cells1Red
            k_s(1) = i - 1                                                                          ! symmetry, junst running from 0,1,...,N/2,N/2+1
            xi2nd(1:3,i,k,j-cells2Offset) = utilities_getFreqDerivative(k_s)
            where(mod(cells,2)==0 .and. [i,j,k] == cells/2+1 .and. &
                  spectral_derivative_ID == DERIVATIVE_CONTINUOUS_ID)                               ! for even grids, set the Nyquist Freq component to 0.0
              xi1st(1:3,i,k,j-cells2Offset) = cmplx(0.0_pReal,0.0_pReal,pReal)
            elsewhere
              xi1st(1:3,i,k,j-cells2Offset) = xi2nd(1:3,i,k,j-cells2Offset)
            endwhere
  end do; end do; end do

  if (num%memory_efficient) then                                                                    ! allocate just single fourth order tensor
    allocate (gamma_hat(3,3,3,3,1,1,1), source = cmplx(0.0_pReal,0.0_pReal,pReal))
  else                                                                                              ! precalculation of gamma_hat field
    allocate (gamma_hat(3,3,3,3,cells1Red,cells(3),cells2), source = cmplx(0.0_pReal,0.0_pReal,pReal))
  end if

  call selfTest()

end subroutine spectral_utilities_init


!---------------------------------------------------------------------------------------------------
!> @brief updates reference stiffness and potentially precalculated gamma operator
!> @details Sets the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of an on-the-fly calculation, only the reference stiffness is updated.
!---------------------------------------------------------------------------------------------------
subroutine utilities_updateGamma(C)

  real(pReal), intent(in), dimension(3,3,3,3) :: C                                                  !< input stiffness to store as reference stiffness

  complex(pReal),              dimension(3,3) :: temp33_cmplx, xiDyad_cmplx
  real(pReal),                 dimension(6,6) :: A, A_inv
  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err


  C_ref = C/wgt

  if (.not. num%memory_efficient) then
    gamma_hat = cmplx(0.0_pReal,0.0_pReal,pReal)                                                    ! for the singular point and any non invertible A
    !$OMP PARALLEL DO PRIVATE(l,m,n,o,temp33_cmplx,xiDyad_cmplx,A,A_inv,err)
    do j = cells2Offset+1, cells2Offset+cells2; do k = 1, cells(3); do i = 1, cells1Red
      if (any([i,j,k] /= 1)) then                                                                   ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
#ifndef __INTEL_COMPILER
        do concurrent(l = 1:3, m = 1:3)
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j-cells2Offset))*xi1st(m,i,k,j-cells2Offset)
        end do
        do concurrent(l = 1:3, m = 1:3)
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal,pReal)*xiDyad_cmplx)
        end do
#else
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j-cells2Offset))*xi1st(m,i,k,j-cells2Offset)
        forall(l = 1:3, m = 1:3) &
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal,pReal)*xiDyad_cmplx)
#endif
        A(1:3,1:3) = temp33_cmplx%re; A(4:6,4:6) =  temp33_cmplx%re
        A(1:3,4:6) = temp33_cmplx%im; A(4:6,1:3) = -temp33_cmplx%im
        if (abs(math_det33(A(1:3,1:3))) > 1.e-16_pReal) then
          call math_invert(A_inv, err, A)
          temp33_cmplx = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pReal)
#ifndef __INTEL_COMPILER
          do concurrent(l=1:3, m=1:3, n=1:3, o=1:3)
            gamma_hat(l,m,n,o,i,k,j-cells2Offset) = temp33_cmplx(l,n) * xiDyad_cmplx(o,m)
          end do
#else
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            gamma_hat(l,m,n,o,i,k,j-cells2Offset) = temp33_cmplx(l,n) * xiDyad_cmplx(o,m)
#endif
        end if
      end if
    end do; end do; end do
    !$OMP END PARALLEL DO
  end if

end subroutine utilities_updateGamma


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an weighted inverse FFT transform from complex to real
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTvectorBackward()

  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  vectorField_real = vectorField_real * wgt                                                         ! normalize the result by number of elements

end subroutine utilities_FFTvectorBackward


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real, ensuring that average value = fieldAim
!--------------------------------------------------------------------------------------------------
function utilities_GammaConvolution(field, fieldAim) result(gammaField)

  real(pReal), intent(in), dimension(3,3,cells(1),cells(2),cells3) :: field
  real(pReal), intent(in), dimension(3,3) :: fieldAim                                               !< desired average value of the field after convolution
  real(pReal),             dimension(3,3,cells(1),cells(2),cells3) :: gammaField

  complex(pReal), dimension(3,3) :: temp33_cmplx, xiDyad_cmplx
  real(pReal),    dimension(6,6) :: A, A_inv
  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err


  print'(/,1x,a)', '... doing gamma convolution ...............................................'
  flush(IO_STDOUT)

  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  memoryEfficient: if (num%memory_efficient) then
    !$OMP PARALLEL DO PRIVATE(l,m,n,o,temp33_cmplx,xiDyad_cmplx,A,A_inv,err,gamma_hat)
    do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
      if (any([i,j+cells2Offset,k] /= 1)) then                                                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
#ifndef __INTEL_COMPILER
        do concurrent(l = 1:3, m = 1:3)
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j))*xi1st(m,i,k,j)
        end do
        do concurrent(l = 1:3, m = 1:3)
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal,pReal)*xiDyad_cmplx)
        end do
#else
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j))*xi1st(m,i,k,j)
        forall(l = 1:3, m = 1:3) &
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pReal,pReal)*xiDyad_cmplx)
#endif
        A(1:3,1:3) = temp33_cmplx%re; A(4:6,4:6) =  temp33_cmplx%re
        A(1:3,4:6) = temp33_cmplx%im; A(4:6,1:3) = -temp33_cmplx%im
        if (abs(math_det33(A(1:3,1:3))) > 1.e-16_pReal) then
          call math_invert(A_inv, err, A)
          temp33_cmplx = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pReal)
#ifndef __INTEL_COMPILER
          do concurrent(l=1:3, m=1:3, n=1:3, o=1:3)
            gamma_hat(l,m,n,o,1,1,1) = temp33_cmplx(l,n)*xiDyad_cmplx(o,m)
          end do
          do concurrent(l = 1:3, m = 1:3)
            temp33_cmplx(l,m) = sum(gamma_hat(l,m,1:3,1:3,1,1,1)*tensorField_fourier(1:3,1:3,i,k,j))
          end do
#else
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            gamma_hat(l,m,n,o,1,1,1) = temp33_cmplx(l,n)*xiDyad_cmplx(o,m)
          forall(l = 1:3, m = 1:3) &
            temp33_cmplx(l,m) = sum(gamma_hat(l,m,1:3,1:3,1,1,1)*tensorField_fourier(1:3,1:3,i,k,j))
#endif
          tensorField_fourier(1:3,1:3,i,k,j) = temp33_cmplx
        else
          tensorField_fourier(1:3,1:3,i,k,j) = cmplx(0.0_pReal,0.0_pReal,pReal)
        end if
      end if
    end do; end do; end do
    !$OMP END PARALLEL DO
  else memoryEfficient
    !$OMP PARALLEL DO PRIVATE(l,m,temp33_cmplx)
    do j = 1, cells2;  do k = 1, cells(3);  do i = 1,cells1Red
#ifndef __INTEL_COMPILER
      do concurrent(l = 1:3, m = 1:3)
        temp33_cmplx(l,m) = sum(gamma_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
      end do
#else
      forall(l = 1:3, m = 1:3) &
        temp33_cmplx(l,m) = sum(gamma_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
#endif
      tensorField_fourier(1:3,1:3,i,k,j) = temp33_cmplx
    end do; end do; end do
    !$OMP END PARALLEL DO
  end if memoryEfficient

  if (cells3Offset == 0) tensorField_fourier(1:3,1:3,1,1,1) = cmplx(fieldAim,0.0_pReal,pReal)

  call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
  gammaField = tensorField_real(1:3,1:3,1:cells(1),1:cells(2),1:cells3)

end function utilities_GammaConvolution


!--------------------------------------------------------------------------------------------------
!> @brief Convolution of Greens' operator for damage/thermal.
!--------------------------------------------------------------------------------------------------
function utilities_GreenConvolution(field, D_ref, mu_ref, Delta_t) result(greenField)

  real(pReal), intent(in), dimension(cells(1),cells(2),cells3) :: field
  real(pReal), dimension(3,3), intent(in) :: D_ref
  real(pReal),                 intent(in) :: mu_ref, Delta_t
  real(pReal), dimension(cells(1),cells(2),cells3) :: greenField

  complex(pReal)                          :: GreenOp_hat
  integer                                 :: i, j, k


  scalarField_real(cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  scalarField_real(1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)

  !$OMP PARALLEL DO PRIVATE(GreenOp_hat)
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    GreenOp_hat = cmplx(wgt,0.0_pReal,pReal) &
                / (cmplx(mu_ref,0.0_pReal,pReal) + cmplx(Delta_t,0.0_pReal,pReal) &
                   * sum(conjg(xi1st(1:3,i,k,j))* matmul(cmplx(D_ref,0.0_pReal,pReal),xi1st(1:3,i,k,j))))
    scalarField_fourier(i,k,j) = scalarField_fourier(i,k,j)*GreenOp_hat
  end do; end do; end do
  !$OMP END PARALLEL DO

  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  greenField = scalarField_real(1:cells(1),1:cells(2),1:cells3)

end function utilities_GreenConvolution


!--------------------------------------------------------------------------------------------------
!> @brief Calculate root mean square of divergence.
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS(tensorField)

  real(pReal), dimension(3,3,cells(1),cells(2),cells3), intent(in) :: tensorField

  integer :: i, j, k
  integer(MPI_INTEGER_KIND) :: err_MPI
  complex(pReal), dimension(3) :: rescaledGeom


  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = tensorField
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
  utilities_divergenceRMS = 0.0_pReal
  do j = 1, cells2; do k = 1, cells(3)
    do i = 2, cells1Red -1                                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
      utilities_divergenceRMS = utilities_divergenceRMS &
            + 2.0_pReal*(sum (real(matmul(tensorField_fourier(1:3,1:3,i,k,j), &                     ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2, i.e. do not take square root and square again
                                          conjg(-xi1st(1:3,i,k,j))*rescaledGeom))**2) &             ! --> sum squared L_2 norm of vector
                        +sum(aimag(matmul(tensorField_fourier(1:3,1:3,i,k,j),&
                                          conjg(-xi1st(1:3,i,k,j))*rescaledGeom))**2))
    end do
    utilities_divergenceRMS = utilities_divergenceRMS &                                             ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart (if cells(1) /= 1)
               + sum( real(matmul(tensorField_fourier(1:3,1:3,1       ,k,j), &
                                  conjg(-xi1st(1:3,1,k,j))*rescaledGeom))**2) &
               + sum(aimag(matmul(tensorField_fourier(1:3,1:3,1       ,k,j), &
                                  conjg(-xi1st(1:3,1,k,j))*rescaledGeom))**2) &
               + sum( real(matmul(tensorField_fourier(1:3,1:3,cells1Red,k,j), &
                                  conjg(-xi1st(1:3,cells1Red,k,j))*rescaledGeom))**2) &
               + sum(aimag(matmul(tensorField_fourier(1:3,1:3,cells1Red,k,j), &
                                  conjg(-xi1st(1:3,cells1Red,k,j))*rescaledGeom))**2)
  end do; end do
  call MPI_Allreduce(MPI_IN_PLACE,utilities_divergenceRMS,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                     ! RMS in real space calculated with Parsevals theorem from Fourier space
  if (cells(1) == 1) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pReal                  ! counted twice in case of cells(1) == 1

end function utilities_divergenceRMS


!--------------------------------------------------------------------------------------------------
!> @brief Calculate root mean square of curl.
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_curlRMS(tensorField)

  real(pReal), dimension(3,3,cells(1),cells(2),cells3), intent(in) :: tensorField

  integer  ::  i, j, k, l
  integer(MPI_INTEGER_KIND) :: err_MPI
  complex(pReal), dimension(3,3) :: curl_fourier
  complex(pReal), dimension(3)   :: rescaledGeom


  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = tensorField
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! calculating max curl criterion in Fourier space
  utilities_curlRMS = 0.0_pReal

  do j = 1, cells2; do k = 1, cells(3);
    do i = 2, cells1Red - 1
      do l = 1, 3
        curl_fourier(l,1) = (+tensorField_fourier(l,3,i,k,j)*xi1st(2,i,k,j)*rescaledGeom(2) &
                             -tensorField_fourier(l,2,i,k,j)*xi1st(3,i,k,j)*rescaledGeom(3))
        curl_fourier(l,2) = (+tensorField_fourier(l,1,i,k,j)*xi1st(3,i,k,j)*rescaledGeom(3) &
                             -tensorField_fourier(l,3,i,k,j)*xi1st(1,i,k,j)*rescaledGeom(1))
        curl_fourier(l,3) = (+tensorField_fourier(l,2,i,k,j)*xi1st(1,i,k,j)*rescaledGeom(1) &
                             -tensorField_fourier(l,1,i,k,j)*xi1st(2,i,k,j)*rescaledGeom(2))
      end do
      utilities_curlRMS = utilities_curlRMS &
                        +2.0_pReal*sum(curl_fourier%re**2+curl_fourier%im**2)                       ! Has somewhere a conj. complex counterpart. Therefore count it twice.
    end do
    do l = 1, 3
       curl_fourier = (+tensorField_fourier(l,3,1,k,j)*xi1st(2,1,k,j)*rescaledGeom(2) &
                       -tensorField_fourier(l,2,1,k,j)*xi1st(3,1,k,j)*rescaledGeom(3))
       curl_fourier = (+tensorField_fourier(l,1,1,k,j)*xi1st(3,1,k,j)*rescaledGeom(3) &
                       -tensorField_fourier(l,3,1,k,j)*xi1st(1,1,k,j)*rescaledGeom(1))
       curl_fourier = (+tensorField_fourier(l,2,1,k,j)*xi1st(1,1,k,j)*rescaledGeom(1) &
                       -tensorField_fourier(l,1,1,k,j)*xi1st(2,1,k,j)*rescaledGeom(2))
    end do
    utilities_curlRMS = utilities_curlRMS &
                      + sum(curl_fourier%re**2 + curl_fourier%im**2)                                ! this layer (DC) does not have a conjugate complex counterpart (if cells(1) /= 1)
    do l = 1, 3
      curl_fourier = (+tensorField_fourier(l,3,cells1Red,k,j)*xi1st(2,cells1Red,k,j)*rescaledGeom(2) &
                      -tensorField_fourier(l,2,cells1Red,k,j)*xi1st(3,cells1Red,k,j)*rescaledGeom(3))
      curl_fourier = (+tensorField_fourier(l,1,cells1Red,k,j)*xi1st(3,cells1Red,k,j)*rescaledGeom(3) &
                      -tensorField_fourier(l,3,cells1Red,k,j)*xi1st(1,cells1Red,k,j)*rescaledGeom(1))
      curl_fourier = (+tensorField_fourier(l,2,cells1Red,k,j)*xi1st(1,cells1Red,k,j)*rescaledGeom(1) &
                      -tensorField_fourier(l,1,cells1Red,k,j)*xi1st(2,cells1Red,k,j)*rescaledGeom(2))
    end do
    utilities_curlRMS = utilities_curlRMS &
                      + sum(curl_fourier%re**2 + curl_fourier%im**2)                                ! this layer (Nyquist) does not have a conjugate complex counterpart (if cells(1) /= 1)
  end do; end do

  call MPI_Allreduce(MPI_IN_PLACE,utilities_curlRMS,1_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  utilities_curlRMS = sqrt(utilities_curlRMS) * wgt                                                 ! RMS in real space calculated with Parsevals theorem from Fourier space
  if (cells(1) == 1) utilities_curlRMS = utilities_curlRMS * 0.5_pReal                              ! counted twice in case of cells(1) == 1

end function utilities_curlRMS


!--------------------------------------------------------------------------------------------------
!> @brief Calculate masked compliance tensor used to adjust F to fullfill stress BC.
!--------------------------------------------------------------------------------------------------
function utilities_maskedCompliance(rot_BC,mask_stress,C)

  real(pReal),                dimension(3,3,3,3) :: utilities_maskedCompliance                      !< masked compliance
  real(pReal),    intent(in), dimension(3,3,3,3) :: C                                               !< current average stiffness
  type(tRotation), intent(in)                    :: rot_BC                                          !< rotation of load frame
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

  mask_stressVector = .not. reshape(transpose(mask_stress), [9])
  size_reduced = count(mask_stressVector)
  if (size_reduced > 0) then
    temp99_real = math_3333to99(rot_BC%rotate(C))

    if (debugGeneral) then
      print'(/,1x,a)', '... updating masked compliance ............................................'
      print'(/,1x,a,/,8(9(2x,f12.7,1x)/),9(2x,f12.7,1x))', &
        'Stiffness C (load) / GPa =', transpose(temp99_Real)*1.0e-9_pReal
      flush(IO_STDOUT)
    end if

    do i = 1,9; do j = 1,9
      mask(i,j) = mask_stressVector(i) .and. mask_stressVector(j)
    end do; end do
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
      formatString = '(/,1x,a,/,'//trim(formatString)//'('//trim(formatString)//'(2x,es9.2,1x)/))'
      print trim(formatString), 'C * S (load) ', transpose(matmul(c_reduced,s_reduced))
      print trim(formatString), 'S (load) ', transpose(s_reduced)
      if (errmatinv) error stop 'matrix inversion error'
    end if
    temp99_real = reshape(unpack(reshape(s_reduced,[size_reduced**2]),reshape(mask,[81]),0.0_pReal),[9,9])
  else
    temp99_real = 0.0_pReal
  end if

  utilities_maskedCompliance = math_99to3333(temp99_Real)

  if (debugGeneral) then
    print'(/,1x,a,/,9(9(2x,f10.5,1x)/),9(2x,f10.5,1x))', &
      'Masked Compliance (load) * GPa =', transpose(temp99_Real)*1.0e9_pReal
    flush(IO_STDOUT)
  end if

end function utilities_maskedCompliance


!--------------------------------------------------------------------------------------------------
!> @brief Calculate gradient of scalar field.
!--------------------------------------------------------------------------------------------------
function utilities_scalarGradient(field) result(grad)

  real(pReal), intent(in), dimension(  cells(1),cells(2),cells3) :: field
  real(pReal),             dimension(3,cells(1),cells(2),cells3) :: grad

  integer :: i, j, k


  scalarField_real(cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  scalarField_real(1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)
  do j = 1, cells2;  do k = 1, cells(3);  do i = 1,cells1Red
    vectorField_fourier(1:3,i,k,j) = scalarField_fourier(i,k,j)*xi1st(1:3,i,k,j)
  end do; end do; end do
  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  grad = vectorField_real(1:3,1:cells(1),1:cells(2),1:cells3)*wgt

end function utilities_scalarGradient


!--------------------------------------------------------------------------------------------------
!> @brief Calculate divergence of vector field.
!--------------------------------------------------------------------------------------------------
function utilities_vectorDivergence(field) result(div)

  real(pReal), intent(in), dimension(3,cells(1),cells(2),cells3) :: field
  real(pReal),             dimension(  cells(1),cells(2),cells3) :: div


  vectorField_real(1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  vectorField_real(1:3,1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)
  scalarField_fourier(1:cells1Red,1:cells(3),1:cells2) = sum(vectorField_fourier(1:3,1:cells1Red,1:cells(3),1:cells2) &
                                                             *conjg(-xi1st),1)                      ! ToDo: use "xi1st" instead of "conjg(-xi1st)"?
  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  div = scalarField_real(1:cells(1),1:cells(2),1:cells3)*wgt

end function utilities_vectorDivergence


!--------------------------------------------------------------------------------------------------
!> @brief calculate constitutive response from homogenization_F0 to F during Delta_t
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(P,P_av,C_volAvg,C_minmaxAvg,&
                                          F,Delta_t,rotation_BC)

  real(pReal),    intent(out), dimension(3,3,3,3)                   :: C_volAvg, C_minmaxAvg        !< average stiffness
  real(pReal),    intent(out), dimension(3,3)                       :: P_av                         !< average PK stress
  real(pReal),    intent(out), dimension(3,3,cells(1),cells(2),cells3) :: P                         !< PK stress
  real(pReal),    intent(in),  dimension(3,3,cells(1),cells(2),cells3) :: F                         !< deformation gradient target
  real(pReal),    intent(in)                                        :: Delta_t                      !< loading time
  type(tRotation), intent(in),  optional                            :: rotation_BC                  !< rotation of load frame


  integer :: i
  integer(MPI_INTEGER_KIND) :: err_MPI
  real(pReal), dimension(3,3,3,3) :: dPdF_max,      dPdF_min
  real(pReal)                     :: dPdF_norm_max, dPdF_norm_min
  real(pReal), dimension(2) :: valueAndRank                                                         !< pair of min/max norm of dPdF to synchronize min/max of dPdF

  print'(/,1x,a)', '... evaluating constitutive response ......................................'
  flush(IO_STDOUT)

  homogenization_F  = reshape(F,[3,3,product(cells(1:2))*cells3])                                   ! set materialpoint target F to estimated field

  call homogenization_mechanical_response(Delta_t,1,product(cells(1:2))*cells3)                     ! calculate P field
  if (.not. terminallyIll) &
    call homogenization_thermal_response(Delta_t,1,product(cells(1:2))*cells3)
  if (.not. terminallyIll) &
    call homogenization_mechanical_response2(Delta_t,[1,1],[1,product(cells(1:2))*cells3])

  P = reshape(homogenization_P, [3,3,cells(1),cells(2),cells3])
  P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt
  call MPI_Allreduce(MPI_IN_PLACE,P_av,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  if (debugRotation) print'(/,1x,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    'Piola--Kirchhoff stress (lab) / MPa =', transpose(P_av)*1.e-6_pReal
  if (present(rotation_BC)) P_av = rotation_BC%rotate(P_av)
  print'(/,1x,a,/,2(3(2x,f12.4,1x)/),3(2x,f12.4,1x))', &
    'Piola--Kirchhoff stress       / MPa =', transpose(P_av)*1.e-6_pReal
  flush(IO_STDOUT)

  dPdF_max = 0.0_pReal
  dPdF_norm_max = 0.0_pReal
  dPdF_min = huge(1.0_pReal)
  dPdF_norm_min = huge(1.0_pReal)
  do i = 1, product(cells(1:2))*cells3
    if (dPdF_norm_max < sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2)) then
      dPdF_max = homogenization_dPdF(1:3,1:3,1:3,1:3,i)
      dPdF_norm_max = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2)
    end if
    if (dPdF_norm_min > sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2)) then
      dPdF_min = homogenization_dPdF(1:3,1:3,1:3,1:3,i)
      dPdF_norm_min = sum(homogenization_dPdF(1:3,1:3,1:3,1:3,i)**2)
    end if
  end do

  valueAndRank = [dPdF_norm_max,real(worldrank,pReal)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1_MPI_INTEGER_KIND,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Bcast(dPdF_max,81_MPI_INTEGER_KIND,MPI_DOUBLE,int(valueAndRank(2),MPI_INTEGER_KIND),MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  valueAndRank = [dPdF_norm_min,real(worldrank,pReal)]
  call MPI_Allreduce(MPI_IN_PLACE,valueAndRank,1_MPI_INTEGER_KIND,MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Bcast(dPdF_min,81_MPI_INTEGER_KIND,MPI_DOUBLE,int(valueAndRank(2),MPI_INTEGER_KIND),MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  C_minmaxAvg = 0.5_pReal*(dPdF_max + dPdF_min)

  C_volAvg = sum(homogenization_dPdF,dim=5)
  call MPI_Allreduce(MPI_IN_PLACE,C_volAvg,81_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  C_volAvg = C_volAvg * wgt


end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief calculates forward rate, either guessing or just add delta/Delta_t
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(heterogeneous,field0,field,dt,avRate)

  real(pReal), intent(in), dimension(3,3) :: &
    avRate                                                                                          !< homogeneous addon
  real(pReal), intent(in) :: &
    dt                                                                                              !< Delta_t between field0 and field
  logical, intent(in) :: &
    heterogeneous                                                                                   !< calculate field of rates
  real(pReal), intent(in), dimension(3,3,cells(1),cells(2),cells3) :: &
    field0, &                                                                                       !< data of previous step
    field                                                                                           !< data of current step
  real(pReal),             dimension(3,3,cells(1),cells(2),cells3) :: &
    utilities_calculateRate


  utilities_calculateRate = merge((field-field0) / dt, &
                                  spread(spread(spread(avRate,3,cells(1)),4,cells(2)),5,cells3), &
                                  heterogeneous)

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given,
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
function utilities_forwardField(Delta_t,field_lastInc,rate,aim)

  real(pReal), intent(in) :: &
    Delta_t                                                                                         !< Delta_t of current step
  real(pReal), intent(in),           dimension(3,3,cells(1),cells(2),cells3) :: &
    field_lastInc, &                                                                                !< initial field
    rate                                                                                            !< rate by which to forward
  real(pReal), intent(in), optional, dimension(3,3) :: &
    aim                                                                                             !< average field value aim

  real(pReal),                       dimension(3,3,cells(1),cells(2),cells3) :: &
    utilities_forwardField
  real(pReal),                       dimension(3,3) :: fieldDiff                                    !< <a + adot*t> - aim
  integer(MPI_INTEGER_KIND) :: err_MPI


  utilities_forwardField = field_lastInc + rate*Delta_t
  if (present(aim)) then                                                                            !< correct to match average
    fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt
    call MPI_Allreduce(MPI_IN_PLACE,fieldDiff,9_MPI_INTEGER_KIND,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
    if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
    fieldDiff = fieldDiff - aim
    utilities_forwardField = utilities_forwardField &
                           - spread(spread(spread(fieldDiff,3,cells(1)),4,cells(2)),5,cells3)
  end if

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief Calculate Filter for Fourier convolution.
!> @details this is the full operator to calculate derivatives, i.e. 2 \pi i k for the
! standard approach
!--------------------------------------------------------------------------------------------------
pure function utilities_getFreqDerivative(k_s)

  integer, intent(in),  dimension(3) :: k_s                                                         !< indices of frequency

  complex(pReal),       dimension(3) :: utilities_getFreqDerivative


  select case (spectral_derivative_ID)
    case (DERIVATIVE_CONTINUOUS_ID)
      utilities_getFreqDerivative = cmplx(0.0_pReal, TAU*real(k_s,pReal)/geomSize,pReal)

    case (DERIVATIVE_CENTRAL_DIFF_ID)
      utilities_getFreqDerivative = cmplx(0.0_pReal, sin(TAU*real(k_s,pReal)/real(cells,pReal)), pReal)/ &
                                    cmplx(2.0_pReal*geomSize/real(cells,pReal), 0.0_pReal, pReal)

    case (DERIVATIVE_FWBW_DIFF_ID)
      utilities_getFreqDerivative(1) = &
                               cmplx(cos(TAU*real(k_s(1),pReal)/real(cells(1),pReal)) - 1.0_pReal, &
                                     sin(TAU*real(k_s(1),pReal)/real(cells(1),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(2),pReal)/real(cells(2),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(2),pReal)/real(cells(2),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(3),pReal)/real(cells(3),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(3),pReal)/real(cells(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(1)/real(cells(1),pReal), 0.0_pReal, pReal)
      utilities_getFreqDerivative(2) = &
                               cmplx(cos(TAU*real(k_s(1),pReal)/real(cells(1),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(1),pReal)/real(cells(1),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(2),pReal)/real(cells(2),pReal)) - 1.0_pReal, &
                                     sin(TAU*real(k_s(2),pReal)/real(cells(2),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(3),pReal)/real(cells(3),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(3),pReal)/real(cells(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(2)/real(cells(2),pReal), 0.0_pReal, pReal)
      utilities_getFreqDerivative(3) = &
                               cmplx(cos(TAU*real(k_s(1),pReal)/real(cells(1),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(1),pReal)/real(cells(1),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(2),pReal)/real(cells(2),pReal)) + 1.0_pReal, &
                                     sin(TAU*real(k_s(2),pReal)/real(cells(2),pReal)), pReal)* &
                               cmplx(cos(TAU*real(k_s(3),pReal)/real(cells(3),pReal)) - 1.0_pReal, &
                                     sin(TAU*real(k_s(3),pReal)/real(cells(3),pReal)), pReal)/ &
                               cmplx(4.0_pReal*geomSize(3)/real(cells(3),pReal), 0.0_pReal, pReal)
  end select

end function utilities_getFreqDerivative


!--------------------------------------------------------------------------------------------------
!> @brief Calculate coordinates in current configuration for given defgrad field
! using integration in Fourier space.
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateCoords(F)

  real(pReal),   dimension(3,3,cells(1),cells(2),cells3), intent(in) :: F

  real(pReal),   dimension(3,  cells(1),cells(2),cells3)             :: IPcoords
  real(pReal),   dimension(3,  cells(1),cells(2),cells3+2)           :: IPfluct_padded              ! Fluctuations of cell center displacement (padded along z for MPI)
  real(pReal),   dimension(3,  cells(1)+1,cells(2)+1,cells3+1)       :: nodeCoords
  integer :: &
    i,j,k,n, &
    c
  integer(MPI_INTEGER_KIND) :: &
    rank_t, rank_b
  integer(MPI_INTEGER_KIND) :: err_MPI
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  type(MPI_Request), dimension(4) :: request
  type(MPI_Status),  dimension(4) :: status
#else
  integer, dimension(4) :: request
  integer, dimension(MPI_STATUS_SIZE,4) :: status
#endif
  real(pReal),   dimension(3)   :: step
  real(pReal),   dimension(3,3) :: Favg
  integer,       dimension(3)   :: me
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


  step = geomSize/real(cells, pReal)

  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = F
  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pReal
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

 !--------------------------------------------------------------------------------------------------
 ! average F
  if (cells3Offset == 0) Favg = tensorField_fourier(1:3,1:3,1,1,1)%re*wgt
  call MPI_Bcast(Favg,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

 !--------------------------------------------------------------------------------------------------
 ! integration in Fourier space to get fluctuations of cell center discplacements
  !$OMP PARALLEL DO
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    if (any([i,j+cells2Offset,k] /= 1)) then
      vectorField_fourier(1:3,i,k,j) = matmul(tensorField_fourier(1:3,1:3,i,k,j),xi2nd(1:3,i,k,j)) &
                                     / sum(conjg(-xi2nd(1:3,i,k,j))*xi2nd(1:3,i,k,j))
    else
      vectorField_fourier(1:3,i,k,j) = cmplx(0.0,0.0,pReal)
    end if
  end do; end do; end do
  !$OMP END PARALLEL DO

  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  vectorField_real = vectorField_real * wgt                                                         ! normalize the result by number of elements

 !--------------------------------------------------------------------------------------------------
 ! pad cell center fluctuations along z-direction (needed when running MPI simulation)
  IPfluct_padded(1:3,1:cells(1),1:cells(2),2:cells3+1) = vectorField_real(1:3,1:cells(1),1:cells(2),1:cells3)
  c = product(shape(IPfluct_padded(:,:,:,1)))                                                       !< amount of data to transfer
  rank_t = modulo(worldrank+1_MPI_INTEGER_KIND,worldsize)
  rank_b = modulo(worldrank-1_MPI_INTEGER_KIND,worldsize)

  ! send bottom layer to process below
  call MPI_Isend(IPfluct_padded(:,:,:,2),       c,MPI_DOUBLE,rank_b,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(1),err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Irecv(IPfluct_padded(:,:,:,cells3+2),c,MPI_DOUBLE,rank_t,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(2),err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  ! send top layer to process above
  call MPI_Isend(IPfluct_padded(:,:,:,cells3+1),c,MPI_DOUBLE,rank_t,1_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(3),err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  call MPI_Irecv(IPfluct_padded(:,:,:,1),       c,MPI_DOUBLE,rank_b,1_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(4),err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  call MPI_Waitall(4,request,status,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  ! ToDo
#else
  if (any(status(MPI_ERROR,:) /= 0)) error stop 'MPI error'
#endif

 !--------------------------------------------------------------------------------------------------
 ! calculate nodal displacements
  nodeCoords = 0.0_pReal
  do j = 0,cells(2); do k = 0,cells3; do i = 0,cells(1)
    nodeCoords(1:3,i+1,j+1,k+1) = matmul(Favg,step*(real([i,j,k+cells3Offset],pReal)))
    averageFluct: do n = 1,8
      me = [i+neighbor(1,n),j+neighbor(2,n),k+neighbor(3,n)]
      nodeCoords(1:3,i+1,j+1,k+1) = nodeCoords(1:3,i+1,j+1,k+1) &
                                  + IPfluct_padded(1:3,modulo(me(1)-1,cells(1))+1,modulo(me(2)-1,cells(2))+1,me(3)+1)*0.125_pReal
    end do averageFluct
  end do; end do; end do

 !--------------------------------------------------------------------------------------------------
 ! calculate cell center displacements
  do k = 1,cells3; do j = 1,cells(2); do i = 1,cells(1)
    IPcoords(1:3,i,j,k) = vectorField_real(1:3,i,j,k) &
                        + matmul(Favg,step*(real([i,j,k+cells3Offset],pReal)-0.5_pReal))
  end do; end do; end do

  call discretization_setNodeCoords(reshape(NodeCoords,[3,(cells(1)+1)*(cells(2)+1)*(cells3+1)]))
  call discretization_setIPcoords  (reshape(IPcoords,  [3,cells(1)*cells(2)*cells3]))

end subroutine utilities_updateCoords


!---------------------------------------------------------------------------------------------------
!> @brief Write out the current reference stiffness for restart.
!---------------------------------------------------------------------------------------------------
subroutine utilities_saveReferenceStiffness()

  integer :: &
    fileUnit,ierr


  if (worldrank == 0) then
    print'(/,1x,a)', '... writing reference stiffness data required for restart to file .........'; flush(IO_STDOUT)
    open(newunit=fileUnit, file=getSolverJobName()//'.C_ref',&
         status='replace',access='stream',action='write',iostat=ierr)
    if (ierr /=0) call IO_error(100,ext_msg='could not open file '//getSolverJobName()//'.C_ref')
    write(fileUnit) C_ref*wgt
    close(fileUnit)
  end if

end subroutine utilities_saveReferenceStiffness


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of forward-backward transform.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  real(pReal), allocatable, dimension(:,:,:,:,:) :: tensorField_real_
  real(pReal), allocatable, dimension(:,:,:,:) :: vectorField_real_
  real(pReal), allocatable, dimension(:,:,:) :: scalarField_real_
  real(pReal), dimension(3,3) :: tensorSum
  real(pReal), dimension(3) :: vectorSum
  real(pReal) :: scalarSum
  real(pReal), dimension(3,3) :: r
  integer(MPI_INTEGER_KIND) :: err_MPI


  call random_number(tensorField_real)
  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  tensorField_real_ = tensorField_real
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)
  call MPI_Allreduce(sum(sum(sum(tensorField_real_,dim=5),dim=4),dim=3),tensorSum,9_MPI_INTEGER_KIND, &
                     MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  if (worldrank==0) then
    if (any(dNeq(tensorSum/tensorField_fourier(:,:,1,1,1)%re,1.0_pReal,1.0e-12_pReal))) &
      error stop 'mismatch avg tensorField FFT <-> real'
  end if
  call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  if (maxval(abs(tensorField_real_ - tensorField_real*wgt))>5.0e-15_pReal) &
    error stop 'mismatch tensorField FFT/invFFT <-> real'

  call random_number(vectorField_real)
  vectorField_real(1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  vectorField_real_ = vectorField_real
  call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)
  call MPI_Allreduce(sum(sum(sum(vectorField_real_,dim=4),dim=3),dim=2),vectorSum,3_MPI_INTEGER_KIND, &
                     MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  if (worldrank==0) then
    if (any(dNeq(vectorSum/vectorField_fourier(:,1,1,1)%re,1.0_pReal,1.0e-12_pReal))) &
      error stop 'mismatch avg vectorField FFT <-> real'
  end if
  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  vectorField_real(1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  if (maxval(abs(vectorField_real_ - vectorField_real*wgt))>5.0e-15_pReal) &
    error stop 'mismatch vectorField FFT/invFFT <-> real'

  call random_number(scalarField_real)
  scalarField_real(cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  scalarField_real_ = scalarField_real
  call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)
  call MPI_Allreduce(sum(sum(sum(scalarField_real_,dim=3),dim=2),dim=1),scalarSum,1_MPI_INTEGER_KIND, &
                     MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'
  if (worldrank==0) then
    if (dNeq(scalarSum/scalarField_fourier(1,1,1)%re,1.0_pReal,1.0e-12_pReal)) &
      error stop 'mismatch avg scalarField FFT <-> real'
  end if
  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  scalarField_real(cells(1)+1:cells1Red*2,:,:) = 0.0_pReal
  if (maxval(abs(scalarField_real_ - scalarField_real*wgt))>5.0e-15_pReal) &
    error stop 'mismatch scalarField FFT/invFFT <-> real'

  call random_number(r)
  call MPI_Bcast(r,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  if (err_MPI /= 0_MPI_INTEGER_KIND) error stop 'MPI error'

  scalarField_real_ = r(1,1)
  if (maxval(abs(utilities_scalarGradient(scalarField_real_)))>5.0e-9_pReal)   error stop 'non-zero grad(const)'

  vectorField_real_ = spread(spread(spread(r(1,:),2,cells(1)),3,cells(2)),4,cells3)
  if (maxval(abs(utilities_vectorDivergence(vectorField_real_)))>5.0e-9_pReal) error stop 'non-zero div(const)'

  tensorField_real_ = spread(spread(spread(r,3,cells(1)),4,cells(2)),5,cells3)
  if (utilities_divergenceRMS(tensorField_real_)>5.0e-14_pReal) error stop 'non-zero RMS div(const)'
  if (utilities_curlRMS(tensorField_real_)>5.0e-14_pReal)       error stop 'non-zero RMS curl(const)'

  if (cells(1) > 2 .and.  spectral_derivative_ID == DERIVATIVE_CONTINUOUS_ID) then
    scalarField_real_ = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
    vectorField_real_ = utilities_scalarGradient(scalarField_real_)/TAU*geomSize(1)
    scalarField_real_ = -spread(spread(planeSine  (cells(1)),2,cells(2)),3,cells3)
    if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-14_pReal) error stop 'grad cosine'
    scalarField_real_ = spread(spread(planeSine  (cells(1)),2,cells(2)),3,cells3)
    vectorField_real_ = utilities_scalarGradient(scalarField_real_)/TAU*geomSize(1)
    scalarField_real_ = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
    if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-14_pReal) error stop 'grad sine'

    vectorField_real_(2:3,:,:,:) = 0.0_pReal
    vectorField_real_(1,:,:,:) = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
    scalarField_real_ = utilities_vectorDivergence(vectorField_real_)/TAU*geomSize(1)
    vectorField_real_(1,:,:,:) =-spread(spread(planeSine(  cells(1)),2,cells(2)),3,cells3)
    if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-14_pReal) error stop 'div cosine'
    vectorField_real_(2:3,:,:,:) = 0.0_pReal
    vectorField_real_(1,:,:,:) = spread(spread(planeSine(  cells(1)),2,cells(2)),3,cells3)
    scalarField_real_ = utilities_vectorDivergence(vectorField_real_)/TAU*geomSize(1)
    vectorField_real_(1,:,:,:) = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
    if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-14_pReal) error stop 'div sine'
  end if

  contains

    function planeCosine(n)
      integer, intent(in) :: n
      real(pReal), dimension(n) :: planeCosine


      planeCosine = cos(real(math_range(n),pReal)/real(n,pReal)*TAU-TAU/real(n*2,pReal))

    end function planeCosine

    function planeSine(n)
      integer, intent(in) :: n
      real(pReal), dimension(n) :: planeSine


      planeSine = sin(real(math_range(n),pReal)/real(n,pReal)*TAU-TAU/real(n*2,pReal))

    end function planeSine

end subroutine selfTest

end module spectral_utilities
