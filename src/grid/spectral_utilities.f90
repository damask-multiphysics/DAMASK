!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Yi Hu, Max-Planck-Institut für Eisenforschung GmbH
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
  real(pREAL), protected, public               :: wgt                                               !< weighting factor 1/Nelems
  real(pREAL), protected, public, dimension(3) :: scaledGeomSize                                    !< scaled geometry size for calculation of divergence
  integer :: &
    cells1Red, &                                                                                    !< cells(1)/2+1
    cells2 = 0, &                                                                                   !< (local) cells in 2nd direction
    cells2Offset                                                                                    !< (local) cells offset in 2nd direction

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW

  real(C_DOUBLE),            dimension(:,:,:,:,:),     pointer     :: tensorField_real              !< tensor field in real space
  real(C_DOUBLE),            dimension(:,:,:,:),       pointer     :: vectorField_real              !< vector field in real space
  real(C_DOUBLE),            dimension(:,:,:),         pointer     :: scalarField_real              !< scalar field in real space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:,:),     pointer     :: tensorField_fourier           !< tensor field in Fourier space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:),       pointer     :: vectorField_fourier           !< vector field in Fourier space
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:),         pointer     :: scalarField_fourier           !< scalar field in Fourier space
  complex(pREAL),            dimension(:,:,:,:,:,:,:), allocatable :: Gamma_hat                     !< gamma operator (field) for spectral method
  complex(pREAL),            dimension(:,:,:,:,:,:,:), allocatable :: G_hat                         !< G operator (field) for Galerkin method
  complex(pREAL),            dimension(:,:,:,:),       allocatable :: xi1st                         !< wave vector field for first derivatives
  complex(pREAL),            dimension(:,:,:,:),       allocatable :: xi2nd                         !< wave vector field for second derivatives
  real(pREAL),               dimension(3,3,3,3)                    :: C_ref                         !< mechanic reference stiffness


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
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from spectral solver variants
    integer :: &
       iterationsNeeded  = 0
    logical :: &
       converged         = .true., &
       stagConverged     = .true.
  end type tSolutionState

  type :: tNumerics
    integer :: &
      divergence_correction                                                                         !< scale divergence/curl calculation
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
    utilities_G_Convolution, &
    utilities_GreenConvolution, &
    utilities_divergenceRMS, &
    utilities_curlRMS, &
    utilities_scalarGradient, &
    utilities_vectorDivergence, &
    utilities_updateCoords

contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate all neccessary fields and create plans for FFTW.
!--------------------------------------------------------------------------------------------------
subroutine spectral_utilities_init(active_Gamma, active_G, active_parabolic)

  logical, intent(in) :: active_Gamma, active_G, active_parabolic

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
  type(tDict) , pointer :: &
    num_solver, &
    num_grid, &
    num_grid_fft

  print'(/,1x,a)', '<<<+-  spectral_utilities init  -+>>>'

  print'(/,1x,a)', 'M. Diehl, Diploma Thesis TU München, 2010'
  print'(  1x,a)', 'https://doi.org/10.13140/2.1.3234.3840'//IO_EOL

  print'(  1x,a)', 'P. Eisenlohr et al., International Journal of Plasticity 46:37–53, 2013'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2012.09.012'//IO_EOL

  print'(  1x,a)', 'P. Shanthraj et al., International Journal of Plasticity 66:31–45, 2015'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2014.02.006'//IO_EOL

  print'(  1x,a)', 'P. Shanthraj et al., Handbook of Mechanics of Materials, 2019'
  print'(  1x,a)', 'https://doi.org/10.1007/978-981-10-6855-3_80'//IO_EOL

  print'(  1x,a)', 'M. Frigo and S.G. Johnson, Proceedings of the IEEE 93(2):216–231, 2005'
  print'(  1x,a)', 'https://doi.org/10.1109/jproc.2004.840301'

  num_solver      => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_grid        => num_solver%get_dict('grid',defaultVal=emptyDict)
  num_grid_fft    => num_grid%get_dict('FFT',defaultVal=emptyDict)

  call PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                                num_grid%get_asStr('PETSc_options',defaultVal=''),err_PETSc)
  CHKERRQ(err_PETSc)

  cells1Red = cells(1)/2 + 1
  wgt = real(product(cells),pREAL)**(-1)

  num%memory_efficient      = num_grid_fft%get_asBool('memory_efficient',     defaultVal=.true.)

!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and
! resolution-independent divergence
  select case (num_grid_fft%get_asStr('divergence_correction',defaultVal='grid+size'))
    case ('none')
      do j = 1, 3
       if (j /= minloc(geomSize,1) .and. j /= maxloc(geomSize,1)) &
         scaledGeomSize = geomSize/geomSize(j)
      end do
    case ('grid+size', 'size+grid')
      do j = 1, 3
       if (      j /= int(minloc(geomSize/real(cells,pREAL),1)) &
           .and. j /= int(maxloc(geomSize/real(cells,pREAL),1))) &
         scaledGeomSize = geomSize/geomSize(j)*real(cells(j),pREAL)
      end do
    case ('size')
      scaledGeomSize = geomSize
    case default
      call IO_error(301,ext_msg=trim(num_grid_fft%get_asStr('divergence_correction')))
  end select

  select case (num_grid_fft%get_asStr('derivative',defaultVal='continuous'))
    case ('continuous')
      spectral_derivative_ID = DERIVATIVE_CONTINUOUS_ID
    case ('central_difference')
      spectral_derivative_ID = DERIVATIVE_CENTRAL_DIFF_ID
    case ('FWBW_difference')
      spectral_derivative_ID = DERIVATIVE_FWBW_DIFF_ID
    case default
      call IO_error(892,ext_msg=trim(num_grid_fft%get_asStr('derivative')))
  end select


  select case(num_grid_fft%get_asStr('FFTW_plan_mode',defaultVal='FFTW_MEASURE'))
    case('FFTW_ESTIMATE')                                                                           ! ordered from slow execution (but fast plan creation) to fast execution
      FFTW_planner_flag = FFTW_ESTIMATE
    case('FFTW_MEASURE')
      FFTW_planner_flag = FFTW_MEASURE
    case('FFTW_PATIENT')
      FFTW_planner_flag = FFTW_PATIENT
    case('FFTW_EXHAUSTIVE')
      FFTW_planner_flag = FFTW_EXHAUSTIVE
    case default
      call IO_warning(47,'using default FFTW_MEASURE instead of "'//trim(num_grid_fft%get_asStr('FFTW_plan_mode'))//'"')
      FFTW_planner_flag = FFTW_MEASURE
  end select

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
  if (pREAL /= C_DOUBLE .or. kind(1) /= C_INT) error stop 'C and Fortran datatypes do not match'
  call fftw_set_timelimit(num_grid_fft%get_asReal('FFTW_timelimit',defaultVal=300.0_pREAL))

  print'(/,1x,a)', 'FFTW initialized'; flush(IO_STDOUT)

  cellsFFTW = int(cells,C_INTPTR_T)


!--------------------------------------------------------------------------------------------------
! set up FFTW data structures for tensor fields
! ToDo: FEM uses tensor field to calculate coordinates, this is not needed
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
! set up FFTW data structures for vector fields
! ToDo: FEM uses vector field to calculate coordinates, this is not needed
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
! set up FFTW data structures for scalar fields
  if (active_parabolic) then
    N = fftw_mpi_local_size_3d_transposed(cellsFFTW(3),cellsFFTW(2),int(cells1Red,C_INTPTR_T), &
                                          PETSC_COMM_WORLD,cells3FFTW,cells3_offset,cells2FFTW,cells2_offset)
    if (int(cells3FFTW) /= cells3) error stop 'domain decomposition mismatch (scalar, real space)'
    if (int(cells2FFTW) /= cells2) error stop 'domain decomposition mismatch (scalar, Fourier space)'

    scalarField = fftw_alloc_complex(N)
    call c_f_pointer(scalarField,scalarField_real, &
                     [int(cells1Red*2,C_INTPTR_T),cellsFFTW(2),cells3FFTW])
    call c_f_pointer(scalarField,scalarField_fourier, &
                     [int(cells1Red,  C_INTPTR_T),cellsFFTW(3),cells2FFTW])

    planScalarForth = fftw_mpi_plan_dft_r2c_3d(cellsFFTW(3),cellsFFTW(2),cellsFFTW(1), &
                                               scalarField_real,scalarField_fourier, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_OUT)
    if (.not. c_associated(planScalarForth)) error stop 'FFTW error r2c scalar'
    planScalarBack  = fftw_mpi_plan_dft_c2r_3d(cellsFFTW(3),cellsFFTW(2),cellsFFTW(1), &
                                               scalarField_fourier,scalarField_real, &
                                               PETSC_COMM_WORLD,FFTW_planner_flag+FFTW_MPI_TRANSPOSED_IN)
    if (.not. c_associated(planScalarBack))  error stop 'FFTW error c2r scalar'
  end if

!--------------------------------------------------------------------------------------------------
! discrete angular frequencies, ordered as in FFTW (wrap around)
  allocate (xi1st (3,cells1Red,cells(3),cells2),source = cmplx(0.0_pREAL,0.0_pREAL,pREAL))          ! frequencies for first derivatives, only half the size for first dimension
  allocate (xi2nd (3,cells1Red,cells(3),cells2),source = cmplx(0.0_pREAL,0.0_pREAL,pREAL))          ! frequencies for second derivatives, only half the size for first dimension

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
              xi1st(1:3,i,k,j-cells2Offset) = cmplx(0.0_pREAL,0.0_pREAL,pREAL)
            elsewhere
              xi1st(1:3,i,k,j-cells2Offset) = xi2nd(1:3,i,k,j-cells2Offset)
            endwhere
  end do; end do; end do

  if (active_Gamma) then
    if (num%memory_efficient) then                                                                  ! allocate just single fourth order tensor
      allocate (Gamma_hat(3,3,3,3,1,1,1), source = cmplx(0.0_pREAL,0.0_pREAL,pREAL))
    else                                                                                            ! precalculation of Gamma_hat field
      allocate (Gamma_hat(3,3,3,3,cells1Red,cells(3),cells2), source = cmplx(0.0_pREAL,0.0_pREAL,pREAL))
    end if
  end if

  if (active_G) G_hat = G_hat_init()

  call selfTest()

end subroutine spectral_utilities_init


!---------------------------------------------------------------------------------------------------
!> @brief Update reference stiffness and potentially precalculated gamma operator.
!> @details Set the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of an on-the-fly calculation, only the reference stiffness is updated.
!---------------------------------------------------------------------------------------------------
subroutine utilities_updateGamma(C)

  real(pREAL), intent(in), dimension(3,3,3,3) :: C                                                  !< input stiffness to store as reference stiffness

  complex(pREAL),              dimension(3,3) :: temp33_cmplx, xiDyad_cmplx
  real(pREAL),                 dimension(6,6) :: A, A_inv
  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err


  C_ref = C/wgt

  if (.not. num%memory_efficient) then
    Gamma_hat = cmplx(0.0_pREAL,0.0_pREAL,pREAL)                                                    ! for the singular point and any non invertible A
    !$OMP PARALLEL DO PRIVATE(l,m,n,o,temp33_cmplx,xiDyad_cmplx,A,A_inv,err)
    do j = cells2Offset+1, cells2Offset+cells2; do k = 1, cells(3); do i = 1, cells1Red
      if (any([i,j,k] /= 1)) then                                                                   ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
#ifndef __INTEL_COMPILER
        do concurrent(l = 1:3, m = 1:3)
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j-cells2Offset))*xi1st(m,i,k,j-cells2Offset)
        end do
        do concurrent(l = 1:3, m = 1:3)
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pREAL,pREAL)*xiDyad_cmplx)
        end do
#else
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j-cells2Offset))*xi1st(m,i,k,j-cells2Offset)
        forall(l = 1:3, m = 1:3) &
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pREAL,pREAL)*xiDyad_cmplx)
#endif
        A(1:3,1:3) = temp33_cmplx%re; A(4:6,4:6) =  temp33_cmplx%re
        A(1:3,4:6) = temp33_cmplx%im; A(4:6,1:3) = -temp33_cmplx%im
        if (abs(math_det33(A(1:3,1:3))) > 1.e-16_pREAL) then
          call math_invert(A_inv, err, A)
          temp33_cmplx = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pREAL)
#ifndef __INTEL_COMPILER
          do concurrent(l=1:3, m=1:3, n=1:3, o=1:3)
            Gamma_hat(l,m,n,o,i,k,j-cells2Offset) = temp33_cmplx(l,n) * xiDyad_cmplx(o,m)
          end do
#else
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            Gamma_hat(l,m,n,o,i,k,j-cells2Offset) = temp33_cmplx(l,n) * xiDyad_cmplx(o,m)
#endif
        end if
      end if
    end do; end do; end do
    !$OMP END PARALLEL DO
  end if

end subroutine utilities_updateGamma


!--------------------------------------------------------------------------------------------------
!> @brief Calculate gamma_hat * field_real (convolution).
!> @details The average value equals the given aim.
!--------------------------------------------------------------------------------------------------
function utilities_GammaConvolution(field, fieldAim) result(gammaField)

  real(pREAL), intent(in), dimension(3,3,cells(1),cells(2),cells3) :: field
  real(pREAL), intent(in), dimension(3,3) :: fieldAim                                               !< desired average value of the field after convolution
  real(pREAL),             dimension(3,3,cells(1),cells(2),cells3) :: gammaField

  complex(pREAL), dimension(3,3) :: temp33_cmplx, xiDyad_cmplx
  real(pREAL),    dimension(6,6) :: A, A_inv
  integer :: &
    i, j, k, &
    l, m, n, o
  logical :: err


  print'(/,1x,a)', '... doing gamma convolution ...............................................'
  flush(IO_STDOUT)

  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  memoryEfficient: if (num%memory_efficient) then
    !$OMP PARALLEL DO PRIVATE(l,m,n,o,temp33_cmplx,xiDyad_cmplx,A,A_inv,err,Gamma_hat)
    do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
      if (any([i,j+cells2Offset,k] /= 1)) then                                                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
#ifndef __INTEL_COMPILER
        do concurrent(l = 1:3, m = 1:3)
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j))*xi1st(m,i,k,j)
        end do
        do concurrent(l = 1:3, m = 1:3)
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pREAL,pREAL)*xiDyad_cmplx)
        end do
#else
        forall(l = 1:3, m = 1:3) &
          xiDyad_cmplx(l,m) = conjg(-xi1st(l,i,k,j))*xi1st(m,i,k,j)
        forall(l = 1:3, m = 1:3) &
          temp33_cmplx(l,m) = sum(cmplx(C_ref(l,1:3,m,1:3),0.0_pREAL,pREAL)*xiDyad_cmplx)
#endif
        A(1:3,1:3) = temp33_cmplx%re; A(4:6,4:6) =  temp33_cmplx%re
        A(1:3,4:6) = temp33_cmplx%im; A(4:6,1:3) = -temp33_cmplx%im
        if (abs(math_det33(A(1:3,1:3))) > 1.e-16_pREAL) then
          call math_invert(A_inv, err, A)
          temp33_cmplx = cmplx(A_inv(1:3,1:3),A_inv(1:3,4:6),pREAL)
#ifndef __INTEL_COMPILER
          do concurrent(l=1:3, m=1:3, n=1:3, o=1:3)
            Gamma_hat(l,m,n,o,1,1,1) = temp33_cmplx(l,n)*xiDyad_cmplx(o,m)
          end do
          do concurrent(l = 1:3, m = 1:3)
            temp33_cmplx(l,m) = sum(Gamma_hat(l,m,1:3,1:3,1,1,1)*tensorField_fourier(1:3,1:3,i,k,j))
          end do
#else
          forall(l=1:3, m=1:3, n=1:3, o=1:3) &
            Gamma_hat(l,m,n,o,1,1,1) = temp33_cmplx(l,n)*xiDyad_cmplx(o,m)
          forall(l = 1:3, m = 1:3) &
            temp33_cmplx(l,m) = sum(Gamma_hat(l,m,1:3,1:3,1,1,1)*tensorField_fourier(1:3,1:3,i,k,j))
#endif
          tensorField_fourier(1:3,1:3,i,k,j) = temp33_cmplx
        else
          tensorField_fourier(1:3,1:3,i,k,j) = cmplx(0.0_pREAL,0.0_pREAL,pREAL)
        end if
      end if
    end do; end do; end do
    !$OMP END PARALLEL DO
  else memoryEfficient
    !$OMP PARALLEL DO PRIVATE(l,m,temp33_cmplx)
    do j = 1, cells2;  do k = 1, cells(3);  do i = 1,cells1Red
#ifndef __INTEL_COMPILER
      do concurrent(l = 1:3, m = 1:3)
        temp33_cmplx(l,m) = sum(Gamma_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
      end do
#else
      forall(l = 1:3, m = 1:3) &
        temp33_cmplx(l,m) = sum(Gamma_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
#endif
      tensorField_fourier(1:3,1:3,i,k,j) = temp33_cmplx
    end do; end do; end do
    !$OMP END PARALLEL DO
  end if memoryEfficient

  if (cells3Offset == 0) tensorField_fourier(1:3,1:3,1,1,1) = cmplx(fieldAim,0.0_pREAL,pREAL)

  call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
  gammaField = tensorField_real(1:3,1:3,1:cells(1),1:cells(2),1:cells3)

end function utilities_GammaConvolution


!--------------------------------------------------------------------------------------------------
!> @brief Assemble G.
!--------------------------------------------------------------------------------------------------
function G_hat_init() result(G_hat_)

  complex(pREAL), dimension(:,:,:,:,:,:,:), allocatable :: G_hat_                                   !< G operator (field) for Galerkin method

  integer :: &
    i, j, k, &
    l, m, n, o
  complex(pREAL) :: xi_norm_2
  complex(pREAL), dimension(3,3), parameter :: delta = cmplx(math_I3,0.0_pREAL,pREAL)


  allocate(G_hat_(3,3,3,3,cells1Red,cells(3),cells2),source=cmplx(.0_pREAL,.0_pREAL,pREAL))

  !$OMP PARALLEL DO PRIVATE(l,m,n,o,xi_norm_2)
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    if (any([i,j+cells2Offset,k] /= 1)) then                                                        ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
      xi_norm_2 = cmplx(abs(dot_product(xi1st(:,i,k,j), xi1st(:,i,k,j))),0.0_pREAL,pREAL)
      if (xi_norm_2%re > 1.e-16_pREAL) then
#ifndef __INTEL_COMPILER
        do concurrent(l=1:3, m=1:3, n=1:3, o=1:3)
            G_hat_(l,m,n,o,i,k,j) = delta(l,n)*conjg(-xi1st(m,i,k,j))*xi1st(o,i,k,j)/xi_norm_2
        end do
#else
        forall(l=1:3, m=1:3, n=1:3, o=1:3)
            G_hat_(l,m,n,o,i,k,j) = delta(l,n)*conjg(-xi1st(m,i,k,j))*xi1st(o,i,k,j)/xi_norm_2
        end forall
#endif
      end if
    end if
  end do; end do; end do
  !$OMP END PARALLEL DO

end function G_hat_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate G * field_real (convolution).
!> @details G*field_real = Fourier_inv( G_hat : Fourier(field_real) )
!> @details Yi: make tensor field compatible
!> @details G_hat index according to S Lucarini et al. MSMSE 2021
!> @details fieldAim is for impose dP of stress bc in formResidual
!> @details fieldAim is not needed for stress bc in formJacobian GK_op
!--------------------------------------------------------------------------------------------------
function utilities_G_Convolution(field,stress_mask,fieldAim) result(G_Field)

  real(pREAL), intent(in), dimension(3,3,cells(1),cells(2),cells3) :: field
  logical,     intent(in), dimension(3,3) :: stress_mask                                            !< impose the mask component <=> G* in Lucarini
  real(pREAL), intent(in), dimension(3,3), optional :: fieldAim                                     !< desired average value of the field after convolution

  real(pREAL),             dimension(3,3,cells(1),cells(2),cells3) :: G_Field

  complex(pREAL), dimension(3,3) :: temp33_cmplx
  integer :: &
    i, j, k, &
    l, m


  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  !$OMP PARALLEL DO PRIVATE(l,m,temp33_cmplx)
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    if (any([i,j+cells2Offset,k] /= 1)) then                                                        ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
#ifndef __INTEL_COMPILER
      do concurrent(l=1:3, m=1:3)
        temp33_cmplx(l,m) = sum(G_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
      end do
#else
      forall(l=1:3, m=1:3)
        temp33_cmplx(l,m) = sum(G_hat(l,m,1:3,1:3,i,k,j)*tensorField_fourier(1:3,1:3,i,k,j))
      end forall
#endif
      tensorField_fourier(1:3,1:3,i,k,j) = temp33_cmplx
    end if
  end do; end do; end do
  !$OMP END PARALLEL DO

  ! Apply stress bounday conditions
  if (cells3Offset == 0) then
    if (present(fieldAim)) &
      tensorField_fourier(1:3,1:3,1,1,1) = cmplx(fieldAim/wgt,0.0_pREAL,pREAL)
    where (stress_mask) tensorField_fourier(1:3,1:3,1,1,1) = cmplx(0.0_pREAL,0.0_pREAL,pREAL)
  end if

  call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
  G_Field = tensorField_real(1:3,1:3,1:cells(1),1:cells(2),1:cells3)

end function utilities_G_Convolution


!--------------------------------------------------------------------------------------------------
!> @brief Convolution of Greens' operator for damage/thermal.
!--------------------------------------------------------------------------------------------------
function utilities_GreenConvolution(field, D_ref, mu_ref, Delta_t) result(greenField)

  real(pREAL), intent(in), dimension(cells(1),cells(2),cells3) :: field
  real(pREAL), dimension(3,3), intent(in) :: D_ref
  real(pREAL),                 intent(in) :: mu_ref, Delta_t
  real(pREAL), dimension(cells(1),cells(2),cells3) :: greenField

  complex(pREAL)                          :: GreenOp_hat
  integer                                 :: i, j, k


  scalarField_real(cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  scalarField_real(1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)

  !$OMP PARALLEL DO PRIVATE(GreenOp_hat)
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    GreenOp_hat = cmplx(wgt,0.0_pREAL,pREAL) &
                / (cmplx(mu_ref,0.0_pREAL,pREAL) + cmplx(Delta_t,0.0_pREAL,pREAL) &
                   * sum(conjg(xi1st(1:3,i,k,j))* matmul(cmplx(D_ref,0.0_pREAL,pREAL),xi1st(1:3,i,k,j))))
    scalarField_fourier(i,k,j) = scalarField_fourier(i,k,j)*GreenOp_hat
  end do; end do; end do
  !$OMP END PARALLEL DO

  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  greenField = scalarField_real(1:cells(1),1:cells(2),1:cells3)

end function utilities_GreenConvolution


!--------------------------------------------------------------------------------------------------
!> @brief Calculate root mean square of divergence.
!--------------------------------------------------------------------------------------------------
real(pREAL) function utilities_divergenceRMS(tensorField)

  real(pREAL), dimension(3,3,cells(1),cells(2),cells3), intent(in) :: tensorField

  integer :: i, j, k
  integer(MPI_INTEGER_KIND) :: err_MPI
  complex(pREAL), dimension(3) :: rescaledGeom


  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = tensorField
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pREAL,pREAL)

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
  utilities_divergenceRMS = 0.0_pREAL
  do j = 1, cells2; do k = 1, cells(3)
    do i = 2, cells1Red -1                                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
      utilities_divergenceRMS = utilities_divergenceRMS &
            + 2.0_pREAL*(sum (real(matmul(tensorField_fourier(1:3,1:3,i,k,j), &                     ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2, i.e. do not take square root and square again
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
  call parallelization_chkerr(err_MPI)
  utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                     ! RMS in real space calculated with Parsevals theorem from Fourier space
  if (cells(1) == 1) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pREAL                  ! counted twice in case of cells(1) == 1

end function utilities_divergenceRMS


!--------------------------------------------------------------------------------------------------
!> @brief Calculate root mean square of curl.
!--------------------------------------------------------------------------------------------------
real(pREAL) function utilities_curlRMS(tensorField)

  real(pREAL), dimension(3,3,cells(1),cells(2),cells3), intent(in) :: tensorField

  integer  ::  i, j, k, l
  integer(MPI_INTEGER_KIND) :: err_MPI
  complex(pREAL), dimension(3,3) :: curl_fourier
  complex(pREAL), dimension(3)   :: rescaledGeom


  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = tensorField
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

  rescaledGeom = cmplx(geomSize/scaledGeomSize,0.0_pREAL,pREAL)

!--------------------------------------------------------------------------------------------------
! calculating max curl criterion in Fourier space
  utilities_curlRMS = 0.0_pREAL

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
                        +2.0_pREAL*sum(curl_fourier%re**2+curl_fourier%im**2)                       ! Has somewhere a conj. complex counterpart. Therefore count it twice.
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
  call parallelization_chkerr(err_MPI)
  utilities_curlRMS = sqrt(utilities_curlRMS) * wgt                                                 ! RMS in real space calculated with Parsevals theorem from Fourier space
  if (cells(1) == 1) utilities_curlRMS = utilities_curlRMS * 0.5_pREAL                              ! counted twice in case of cells(1) == 1

end function utilities_curlRMS


!--------------------------------------------------------------------------------------------------
!> @brief Calculate gradient of scalar field.
!--------------------------------------------------------------------------------------------------
function utilities_scalarGradient(field) result(grad)

  real(pREAL), intent(in), dimension(  cells(1),cells(2),cells3) :: field
  real(pREAL),             dimension(3,cells(1),cells(2),cells3) :: grad

  integer :: i, j, k


  scalarField_real(cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
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

  real(pREAL), intent(in), dimension(3,cells(1),cells(2),cells3) :: field
  real(pREAL),             dimension(  cells(1),cells(2),cells3) :: div


  vectorField_real(1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  vectorField_real(1:3,1:cells(1),            1:cells(2),1:cells3) = field
  call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)
  scalarField_fourier(1:cells1Red,1:cells(3),1:cells2) = sum(vectorField_fourier(1:3,1:cells1Red,1:cells(3),1:cells2) &
                                                             *conjg(-xi1st),1)                      ! ToDo: use "xi1st" instead of "conjg(-xi1st)"?
  call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
  div = scalarField_real(1:cells(1),1:cells(2),1:cells3)*wgt

end function utilities_vectorDivergence


!--------------------------------------------------------------------------------------------------
!> @brief Calculate Filter for Fourier convolution.
!> @details this is the full operator to calculate derivatives, i.e. 2 \pi i k for the
! standard approach
!--------------------------------------------------------------------------------------------------
pure function utilities_getFreqDerivative(k_s)

  integer, intent(in),  dimension(3) :: k_s                                                         !< indices of frequency

  complex(pREAL),       dimension(3) :: utilities_getFreqDerivative


  select case (spectral_derivative_ID)
    case (DERIVATIVE_CONTINUOUS_ID)
      utilities_getFreqDerivative = cmplx(0.0_pREAL, TAU*real(k_s,pREAL)/geomSize,pREAL)

    case (DERIVATIVE_CENTRAL_DIFF_ID)
      utilities_getFreqDerivative = cmplx(0.0_pREAL, sin(TAU*real(k_s,pREAL)/real(cells,pREAL)), pREAL)/ &
                                    cmplx(2.0_pREAL*geomSize/real(cells,pREAL), 0.0_pREAL, pREAL)

    case (DERIVATIVE_FWBW_DIFF_ID)
      utilities_getFreqDerivative(1) = &
                               cmplx(cos(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)) - 1.0_pREAL, &
                                     sin(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)), pREAL)/ &
                               cmplx(4.0_pREAL*geomSize(1)/real(cells(1),pREAL), 0.0_pREAL, pREAL)
      utilities_getFreqDerivative(2) = &
                               cmplx(cos(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)) - 1.0_pREAL, &
                                     sin(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)), pREAL)/ &
                               cmplx(4.0_pREAL*geomSize(2)/real(cells(2),pREAL), 0.0_pREAL, pREAL)
      utilities_getFreqDerivative(3) = &
                               cmplx(cos(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(1),pREAL)/real(cells(1),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)) + 1.0_pREAL, &
                                     sin(TAU*real(k_s(2),pREAL)/real(cells(2),pREAL)), pREAL)* &
                               cmplx(cos(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)) - 1.0_pREAL, &
                                     sin(TAU*real(k_s(3),pREAL)/real(cells(3),pREAL)), pREAL)/ &
                               cmplx(4.0_pREAL*geomSize(3)/real(cells(3),pREAL), 0.0_pREAL, pREAL)
  end select

end function utilities_getFreqDerivative


!--------------------------------------------------------------------------------------------------
!> @brief Calculate coordinates in current configuration for given defgrad field
! using integration in Fourier space.
!--------------------------------------------------------------------------------------------------
subroutine utilities_updateCoords(F)

  real(pREAL),   dimension(3,3,cells(1),cells(2),cells3), intent(in) :: F

  real(pREAL),   dimension(3,  cells(1),cells(2),cells3)             :: x_p                         !< Point/cell center coordinates
  real(pREAL),   dimension(3,  cells(1),cells(2),0:cells3+1)         :: u_tilde_p_padded            !< Fluctuation of cell center displacement (padded along z for MPI)
  real(pREAL),   dimension(3,  cells(1)+1,cells(2)+1,cells3+1)       :: x_n                         !< Node coordinates
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
  real(pREAL),   dimension(3)   :: step
  real(pREAL),   dimension(3,3) :: Favg
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


  step = geomSize/real(cells, pREAL)

  tensorField_real(1:3,1:3,1:cells(1),            1:cells(2),1:cells3) = F
  tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,1:cells(2),1:cells3) = 0.0_pREAL
  call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)

 !--------------------------------------------------------------------------------------------------
 ! average F
  if (cells3Offset == 0) Favg = tensorField_fourier(1:3,1:3,1,1,1)%re*wgt
  call MPI_Bcast(Favg,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

 !--------------------------------------------------------------------------------------------------
 ! integration in Fourier space to get fluctuations of cell center displacements
  !$OMP PARALLEL DO
  do j = 1, cells2; do k = 1, cells(3); do i = 1, cells1Red
    if (any([i,j+cells2Offset,k] /= 1)) then
      vectorField_fourier(1:3,i,k,j) = matmul(tensorField_fourier(1:3,1:3,i,k,j),xi2nd(1:3,i,k,j)) &
                                     / sum(conjg(-xi2nd(1:3,i,k,j))*xi2nd(1:3,i,k,j))
    else
      vectorField_fourier(1:3,i,k,j) = cmplx(0.0,0.0,pREAL)
    end if
  end do; end do; end do
  !$OMP END PARALLEL DO

  call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
  u_tilde_p_padded(1:3,1:cells(1),1:cells(2),1:cells3) = vectorField_real(1:3,1:cells(1),1:cells(2),1:cells3) * wgt

 !--------------------------------------------------------------------------------------------------
 ! pad cell center fluctuations along z-direction (needed when running MPI simulation)
  c = product(shape(u_tilde_p_padded(:,:,:,1)))                                                     !< amount of data to transfer
  rank_t = modulo(worldrank+1_MPI_INTEGER_KIND,worldsize)
  rank_b = modulo(worldrank-1_MPI_INTEGER_KIND,worldsize)

  ! send bottom layer to process below
  call MPI_Isend(u_tilde_p_padded(:,:,:,1),       c,MPI_DOUBLE,rank_b,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(1),err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Irecv(u_tilde_p_padded(:,:,:,cells3+1),c,MPI_DOUBLE,rank_t,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(2),err_MPI)
  call parallelization_chkerr(err_MPI)

  ! send top layer to process above
  call MPI_Isend(u_tilde_p_padded(:,:,:,cells3)  ,c,MPI_DOUBLE,rank_t,1_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(3),err_MPI)
  call parallelization_chkerr(err_MPI)
  call MPI_Irecv(u_tilde_p_padded(:,:,:,0),       c,MPI_DOUBLE,rank_b,1_MPI_INTEGER_KIND,MPI_COMM_WORLD,request(4),err_MPI)
  call parallelization_chkerr(err_MPI)

  call MPI_Waitall(4,request,status,err_MPI)
  call parallelization_chkerr(err_MPI)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  ! ToDo
#else
  if (any(status(MPI_ERROR,:) /= 0)) error stop 'MPI error'
#endif

 !--------------------------------------------------------------------------------------------------
 ! calculate nodal positions
  x_n = 0.0_pREAL
  do j = 0,cells(2); do k = 0,cells3; do i = 0,cells(1)
    x_n(1:3,i+1,j+1,k+1) = matmul(Favg,step*(real([i,j,k+cells3Offset],pREAL)))
    averageFluct: do n = 1,8
      me = [i+neighbor(1,n),j+neighbor(2,n),k+neighbor(3,n)]
      x_n(1:3,i+1,j+1,k+1) = x_n(1:3,i+1,j+1,k+1) &
                           + u_tilde_p_padded(1:3,modulo(me(1)-1,cells(1))+1,modulo(me(2)-1,cells(2))+1,me(3))*0.125_pREAL
    end do averageFluct
  end do; end do; end do

 !--------------------------------------------------------------------------------------------------
 ! calculate cell center/point positions
  do k = 1,cells3; do j = 1,cells(2); do i = 1,cells(1)
    x_p(1:3,i,j,k) = u_tilde_p_padded(1:3,i,j,k) &
                   + matmul(Favg,step*(real([i,j,k+cells3Offset],pREAL)-0.5_pREAL))
  end do; end do; end do

  call discretization_setNodeCoords(reshape(x_n,[3,(cells(1)+1)*(cells(2)+1)*(cells3+1)]))
  call discretization_setIPcoords  (reshape(x_p,[3,cells(1)*cells(2)*cells3]))

end subroutine utilities_updateCoords


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of forward-backward transform.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  real(pREAL), allocatable, dimension(:,:,:,:,:) :: tensorField_real_
  real(pREAL), allocatable, dimension(:,:,:,:) :: vectorField_real_
  real(pREAL), allocatable, dimension(:,:,:) :: scalarField_real_
  real(pREAL), dimension(3,3) :: tensorSum
  real(pREAL), dimension(3) :: vectorSum
  real(pREAL) :: scalarSum
  real(pREAL), dimension(3,3) :: r
  integer(MPI_INTEGER_KIND) :: err_MPI

  call random_number(r)
  call MPI_Bcast(r,9_MPI_INTEGER_KIND,MPI_DOUBLE,0_MPI_INTEGER_KIND,MPI_COMM_WORLD,err_MPI)
  call parallelization_chkerr(err_MPI)

  if (associated(tensorField_real)) then
    call random_number(tensorField_real)
    tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    tensorField_real_ = tensorField_real
    call fftw_mpi_execute_dft_r2c(planTensorForth,tensorField_real,tensorField_fourier)
    call MPI_Allreduce(sum(sum(sum(tensorField_real_,dim=5),dim=4),dim=3),tensorSum,9_MPI_INTEGER_KIND, &
                       MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    if (worldrank==0) then
      if (any(dNeq(tensorSum/tensorField_fourier(:,:,1,1,1)%re,1.0_pREAL,1.0e-12_pREAL))) &
        error stop 'mismatch avg tensorField FFT <-> real'
    end if
    call fftw_mpi_execute_dft_c2r(planTensorBack,tensorField_fourier,tensorField_real)
    tensorField_real(1:3,1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    if (maxval(abs(tensorField_real_ - tensorField_real*wgt))>5.0e-15_pREAL) &
      error stop 'mismatch tensorField FFT/invFFT <-> real'

    tensorField_real_ = spread(spread(spread(r,3,cells(1)),4,cells(2)),5,cells3)
    if (utilities_divergenceRMS(tensorField_real_)>5.0e-14_pREAL) error stop 'non-zero RMS div(const)'
    if (utilities_curlRMS(tensorField_real_)>5.0e-14_pREAL)       error stop 'non-zero RMS curl(const)'
  end if

  if (associated(vectorField_real)) then
    call random_number(vectorField_real)
    vectorField_real(1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    vectorField_real_ = vectorField_real
    call fftw_mpi_execute_dft_r2c(planVectorForth,vectorField_real,vectorField_fourier)
    call MPI_Allreduce(sum(sum(sum(vectorField_real_,dim=4),dim=3),dim=2),vectorSum,3_MPI_INTEGER_KIND, &
                       MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    if (worldrank==0) then
      if (any(dNeq(vectorSum/vectorField_fourier(:,1,1,1)%re,1.0_pREAL,1.0e-12_pREAL))) &
        error stop 'mismatch avg vectorField FFT <-> real'
    end if
    call fftw_mpi_execute_dft_c2r(planVectorBack,vectorField_fourier,vectorField_real)
    vectorField_real(1:3,cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    if (maxval(abs(vectorField_real_ - vectorField_real*wgt))>5.0e-15_pREAL) &
      error stop 'mismatch vectorField FFT/invFFT <-> real'
  end if

  if (associated(scalarField_real)) then
    call random_number(scalarField_real)
    scalarField_real(cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    scalarField_real_ = scalarField_real
    call fftw_mpi_execute_dft_r2c(planScalarForth,scalarField_real,scalarField_fourier)
    call MPI_Allreduce(sum(sum(sum(scalarField_real_,dim=3),dim=2),dim=1),scalarSum,1_MPI_INTEGER_KIND, &
                       MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,err_MPI)
    call parallelization_chkerr(err_MPI)
    if (worldrank==0) then
      if (dNeq(scalarSum/scalarField_fourier(1,1,1)%re,1.0_pREAL,1.0e-12_pREAL)) &
        error stop 'mismatch avg scalarField FFT <-> real'
    end if
    call fftw_mpi_execute_dft_c2r(planScalarBack,scalarField_fourier,scalarField_real)
    scalarField_real(cells(1)+1:cells1Red*2,:,:) = 0.0_pREAL
    if (maxval(abs(scalarField_real_ - scalarField_real*wgt))>5.0e-15_pREAL) &
      error stop 'mismatch scalarField FFT/invFFT <-> real'
  end if

  if (associated(scalarField_real) .and. associated(vectorField_real)) then
    scalarField_real_ = r(1,1)
    if (maxval(abs(utilities_scalarGradient(scalarField_real_)))>5.0e-9_pREAL)   error stop 'non-zero grad(const)'

    vectorField_real_ = spread(spread(spread(r(1,:),2,cells(1)),3,cells(2)),4,cells3)
    if (maxval(abs(utilities_vectorDivergence(vectorField_real_)))>5.0e-9_pREAL) error stop 'non-zero div(const)'


    if (cells(1) > 2 .and.  spectral_derivative_ID == DERIVATIVE_CONTINUOUS_ID) then
      scalarField_real_ = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
      vectorField_real_ = utilities_scalarGradient(scalarField_real_)/TAU*geomSize(1)
      scalarField_real_ = -spread(spread(planeSine  (cells(1)),2,cells(2)),3,cells3)
      if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-12_pREAL) error stop 'grad cosine'
      scalarField_real_ = spread(spread(planeSine  (cells(1)),2,cells(2)),3,cells3)
      vectorField_real_ = utilities_scalarGradient(scalarField_real_)/TAU*geomSize(1)
      scalarField_real_ = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
      if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-12_pREAL) error stop 'grad sine'

      vectorField_real_(2:3,:,:,:) = 0.0_pREAL
      vectorField_real_(1,:,:,:) = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
      scalarField_real_ = utilities_vectorDivergence(vectorField_real_)/TAU*geomSize(1)
      vectorField_real_(1,:,:,:) =-spread(spread(planeSine(  cells(1)),2,cells(2)),3,cells3)
      if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-12_pREAL) error stop 'div cosine'
      vectorField_real_(2:3,:,:,:) = 0.0_pREAL
      vectorField_real_(1,:,:,:) = spread(spread(planeSine(  cells(1)),2,cells(2)),3,cells3)
      scalarField_real_ = utilities_vectorDivergence(vectorField_real_)/TAU*geomSize(1)
      vectorField_real_(1,:,:,:) = spread(spread(planeCosine(cells(1)),2,cells(2)),3,cells3)
      if (maxval(abs(vectorField_real_(1,:,:,:) - scalarField_real_))>5.0e-12_pREAL) error stop 'div sine'
    end if
  end if

  contains

    function planeCosine(n)
      integer, intent(in) :: n
      real(pREAL), dimension(n) :: planeCosine


      planeCosine = cos(real(math_range(n),pREAL)/real(n,pREAL)*TAU-TAU/real(n*2,pREAL))

    end function planeCosine

    function planeSine(n)
      integer, intent(in) :: n
      real(pREAL), dimension(n) :: planeSine


      planeSine = sin(real(math_range(n),pREAL)/real(n,pREAL)*TAU-TAU/real(n*2,pREAL))

    end function planeSine

end subroutine selfTest

end module spectral_utilities
