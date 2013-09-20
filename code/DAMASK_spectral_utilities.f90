! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
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

 implicit none
 private
#ifdef PETSc
#include <finclude/petscsys.h>
#endif
 logical,       public                                         :: cutBack =.false.                  !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill
!--------------------------------------------------------------------------------------------------
! grid related information information
 integer(pInt), public,  dimension(3) :: grid                                                       !< grid points as specified in geometry file
 real(pReal),   public                :: wgt                                                        !< weighting factor 1/Nelems
 real(pReal),   public,  dimension(3) :: geomSize                                                   !< size of geometry as specified in geometry file
 
!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 integer(pInt), public                                         :: grid1Red                          !< grid(1)/2
 real(pReal),   public,  dimension(:,:,:,:,:),     pointer     :: field_real                        !< real representation (some stress or deformation) of field_fourier
 complex(pReal),public,  dimension(:,:,:,:,:),     pointer     :: field_fourier                     !< field on which the Fourier transform operates
 real(pReal),   private, dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat                         !< gamma operator (field) for spectral method
 real(pReal),   private, dimension(:,:,:,:),       allocatable :: xi                                !< wave vector field for divergence and for gamma operator
 real(pReal),   private, dimension(3,3,3,3)                    :: C_ref                             !< reference stiffness
 real(pReal),   private, dimension(3)                          :: scaledGeomSize                    !< scaled geometry size for calculation of divergence (Basic, Basic PETSc)

!--------------------------------------------------------------------------------------------------
! debug fftw 
 complex(pReal),private, dimension(:,:,:), pointer :: scalarField_real, &                           !< scalar field real representation for debug of FFTW
                                                      scalarField_fourier                           !< scalar field complex representation for debug of FFTW
 
!--------------------------------------------------------------------------------------------------
! debug divergence
 real(pReal),   private, dimension(:,:,:,:), pointer     :: divReal                                 !< scalar field real representation for debugging divergence calculation
 complex(pReal),private, dimension(:,:,:,:), pointer     :: divFourier                              !< scalar field real representation for debugging divergence calculation

!--------------------------------------------------------------------------------------------------
! plans for FFTW
 type(C_PTR),   private :: &
   planForth, &                                                                                     !< FFTW plan P(x) to P(k)
   planBack, &                                                                                      !< FFTW plan F(k) to F(x)
   planDebugForth, &                                                                                !< FFTW plan for scalar field (proof that order of usual transform is correct)
   planDebugBack, &                                                                                 !< FFTW plan for scalar field inverse (proof that order of usual transform is correct)
   planDiv                                                                                          !< plan for FFTW in case of debugging divergence calculation

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
   logical       :: termIll           = .false.   
   integer(pInt) :: iterationsNeeded  = 0_pInt
 end type tSolutionState

 type, public :: tBoundaryCondition                                                                 !< set of parameters defining a boundary condition
   real(pReal), dimension(3,3) :: values      = 0.0_pReal
   real(pReal), dimension(3,3) :: maskFloat   = 0.0_pReal
   logical,     dimension(3,3) :: maskLogical = .false.
   character(len=64)           :: myType      = 'None'
 end type tBoundaryCondition
 
 public :: &
   utilities_init, &
   utilities_updateGamma, &
   utilities_FFTforward, &
   utilities_FFTbackward, &
   utilities_fourierConvolution, &
   utilities_inverseLaplace, &
   utilities_divergenceRMS, &
   utilities_curlRMS, &
   utilities_maskedCompliance, &
   utilities_constitutiveResponse, &
   utilities_calculateRate, &
   utilities_forwardField, &
   utilities_destroy
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
 use DAMASK_interface, only: &
  geometryFile
 use IO, only: &
   IO_error, &
   IO_warning, &
   IO_timeStamp, &
   IO_open_file
 use numerics, only: &                        
   DAMASK_NumThreadsInt, &
   fftw_planner_flag, &
   fftw_timelimit, &
   memory_efficient, &
   petsc_options, &
   divergence_correction
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
 use math                                                                                           ! must use the whole module for use of FFTW
 use mesh, only: &
  mesh_spectral_getSize, &
  mesh_spectral_getGrid

 implicit none
#ifdef PETSc
 external :: &
   PETScOptionsClear, &
   PETScOptionsInsertString, &
   MPI_Abort
 PetscErrorCode :: ierr
#endif  
 integer(pInt)               :: i, j, k
 integer(pInt), parameter    :: fileUnit = 228_pInt
 integer(pInt), dimension(3) :: k_s
 type(C_PTR) :: &
   tensorField, &                                                                                   !< field cotaining data for FFTW in real and fourier space (in place)
   scalarField_realC, &                                                                             !< field cotaining data for FFTW in real space when debugging FFTW (no in place)
   scalarField_fourierC, &                                                                          !< field cotaining data for FFTW in fourier space when debugging FFTW (no in place)
   div                                                                                              !< field cotaining data for FFTW in real and fourier space when debugging divergence (in place)
 write(6,'(/,a)')   ' <<<+-  DAMASK_spectral_utilities init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_SPECTRAL),debug_LEVELBASIC)         /= 0
 debugDivergence = iand(debug_level(debug_SPECTRAL),debug_SPECTRALDIVERGENCE) /= 0
 debugFFTW       = iand(debug_level(debug_SPECTRAL),debug_SPECTRALFFTW)       /= 0
 debugRotation   = iand(debug_level(debug_SPECTRAL),debug_SPECTRALROTATION)   /= 0
 debugPETSc      = iand(debug_level(debug_SPECTRAL),debug_SPECTRALPETSC)      /= 0
#ifdef PETSc
 if(debugPETSc) write(6,'(3(/,a),/)') &
                ' Initializing PETSc with debug options: ', &
                trim(PETScDebug), &
                ' add more using the PETSc_Options keyword in numerics.config '
 flush(6)
 call PetscOptionsClear(ierr); CHKERRQ(ierr)
 if(debugPETSc) call PetscOptionsInsertString(trim(PETSCDEBUG),ierr); CHKERRQ(ierr)
 call PetscOptionsInsertString(trim(petsc_options),ierr); CHKERRQ(ierr)
#else
 if(debugPETSc) call IO_warning(41_pInt, ext_msg='debug PETSc')
#endif

 call IO_open_file(fileUnit,geometryFile)                                                           ! parse info from geometry file...
 grid = mesh_spectral_getGrid(fileUnit)
 grid1Red = grid(1)/2_pInt + 1_pInt
 wgt = 1.0/real(product(grid),pReal)
 geomSize = mesh_spectral_getSize(fileUnit)
 close(fileUnit)

 write(6,'(a,3(i12  ))')  ' grid     a b c: ', grid
 write(6,'(a,3(es12.5))') ' size     x y z: ', geomSize

!--------------------------------------------------------------------------------------------------
! scale dimension to calculate either uncorrected, dimension-independent, or dimension- and reso-
! lution-independent divergence
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
! allocation
 allocate (xi(3,grid1Red,grid(2),grid(3)),source = 0.0_pReal)                                       ! frequencies, only half the size for first dimension
 tensorField = fftw_alloc_complex(int(grid1Red*grid(2)*grid(3)*9_pInt,C_SIZE_T))                    ! allocate aligned data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, field_real, [grid(1)+2_pInt-mod(grid(1),2_pInt),grid(2),grid(3),3,3])! place a pointer for a real representation on tensorField
 call c_f_pointer(tensorField, field_fourier,[grid1Red,                        grid(2),grid(3),3,3])! place a pointer for a complex representation on tensorField

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(0_pInt,ext_msg='Fortran to C')             ! check for correct precision in C
!$ if(DAMASK_NumThreadsInt > 0_pInt) then
!$   i = fftw_init_threads()                                                                        ! returns 0 in case of problem
!$   if (i == 0_pInt) call IO_error(error_ID = 809_pInt)
!$   call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
!$ endif
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans for the convolution
 planForth =  fftw_plan_many_dft_r2c(3,[grid(3),grid(2) ,grid(1)], 9, &                             ! dimensions,  logical length in each dimension in reversed order,  no. of transforms
                            field_real,[grid(3),grid(2) ,grid(1)+2_pInt-mod(grid(1),2_pInt)], &     ! input data,  physical length in each dimension in reversed order
                                     1, grid(3)*grid(2)*(grid(1)+2_pInt-mod(grid(1),2_pInt)), &     ! striding,    product of physical length in the 3 dimensions
                         field_fourier,[grid(3),grid(2) ,grid1Red], &                               ! output data, physical length in each dimension in reversed order
                                     1, grid(3)*grid(2)* grid1Red,       fftw_planner_flag)         ! striding,    product of physical length in the 3 dimensions,      planner precision

 planBack  =  fftw_plan_many_dft_c2r(3,[grid(3),grid(2) ,grid(1)], 9, &                             ! dimensions,  logical length in each dimension in reversed order,  no. of transforms
                         field_fourier,[grid(3),grid(2) ,grid1Red], &                               ! input data,  physical length in each dimension in reversed order
                                     1, grid(3)*grid(2)* grid1Red, &                                ! striding,    product of physical length in the 3 dimensions
                            field_real,[grid(3),grid(2) ,grid(1)+2_pInt-mod(grid(1),2_pInt)], &     ! output data, physical length in each dimension in reversed order
                                     1, grid(3)*grid(2)*(grid(1)+2_pInt-mod(grid(1),2_pInt)), &     ! striding,    product of physical length in the 3 dimensions
                                                                       fftw_planner_flag)           ! planner precision

!--------------------------------------------------------------------------------------------------
! depending on debug options, allocate more memory and create additional plans 
 if (debugDivergence) then
   div = fftw_alloc_complex(int(grid1Red*grid(2)*grid(3)*3_pInt,C_SIZE_T))
   call c_f_pointer(div,divReal,   [grid(1)+2_pInt-mod(grid(1),2_pInt),grid(2),grid(3),3])
   call c_f_pointer(div,divFourier,[grid1Red,                          grid(2),grid(3),3])
   planDiv  = fftw_plan_many_dft_c2r(3,[grid(3),grid(2) ,grid(1)],3,&
                            divFourier,[grid(3),grid(2) ,grid1Red],&
                                     1, grid(3)*grid(2)* grid1Red,&
                               divReal,[grid(3),grid(2) ,grid(1)+2_pInt-mod(grid(1),2_pInt)], &
                                     1, grid(3)*grid(2)*(grid(1)+2_pInt-mod(grid(1),2_pInt)), &
                                                                       fftw_planner_flag)
 endif

 if (debugFFTW) then
   scalarField_realC    = fftw_alloc_complex(int(product(grid),C_SIZE_T))                           ! allocate data for real representation (no in place transform)
   scalarField_fourierC = fftw_alloc_complex(int(product(grid),C_SIZE_T))                           ! allocate data for fourier representation (no in place transform)
   call c_f_pointer(scalarField_realC,    scalarField_real,    grid)                                ! place a pointer for a real representation
   call c_f_pointer(scalarField_fourierC, scalarField_fourier, grid)                                ! place a pointer for a fourier representation
   planDebugForth = fftw_plan_dft_3d(grid(3),grid(2),grid(1),&                                      ! reversed order (C style)
                                      scalarField_real,scalarField_fourier,-1,fftw_planner_flag)    ! input, output, forward FFT(-1), planner precision
   planDebugBack  = fftw_plan_dft_3d(grid(3),grid(2),grid(1),&                                      ! reversed order (C style)
                                      scalarField_fourier,scalarField_real,+1,fftw_planner_flag)    ! input, output, backward (1), planner precision
 endif 

 if (debugGeneral) write(6,'(/,a)') ' FFTW initialized'
 flush(6)
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 do k = 1_pInt, grid(3)
   k_s(3) = k - 1_pInt
   if(k > grid(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - grid(3)                                        ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
     do j = 1_pInt, grid(2)
       k_s(2) = j - 1_pInt
       if(j > grid(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - grid(2)                                    ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
         do i = 1_pInt, grid1Red
           k_s(1) = i - 1_pInt                                                                      ! symmetry, junst running from 0,1,...,N/2,N/2+1
           xi(1:3,i,j,k) = real(k_s, pReal)/scaledGeomSize                                          ! if divergence_correction is set, frequencies are calculated on unit length
 enddo; enddo; enddo
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(3,3,3,3,1,1,1), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(3,3,3,3,grid1Red ,grid(2),grid(3)), source = 0.0_pReal)     
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
   memory_efficient
 use math, only: &
   math_inv33

 implicit none
 real(pReal), intent(in), dimension(3,3,3,3) :: C                                                   !< input stiffness to store as reference stiffness
 logical    , intent(in)                     :: saveReference                                       !< save reference stiffness to file for restart
 real(pReal),                 dimension(3,3) :: temp33_Real, xiDyad
 real(pReal)                                 :: filter                                              !< weighting of current component
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o
  
 C_ref = C
 if (saveReference) then
   write(6,'(/,a)') ' writing reference stiffness to file'
   flush(6)
   call IO_write_jobRealFile(777,'C_ref',size(C_ref))
   write (777,rec=1) C_ref
   close(777)
 endif
 
 if(.not. memory_efficient) then                                                   
   do k = 1_pInt, grid(3); do j = 1_pInt, grid(2); do i = 1_pInt, grid1Red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       filter = utilities_getFilter(xi(1:3,i,j,k))                                                  ! weighting factor computed by getFilter function
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o, i,j,k) =  filter*temp33_Real(l,n)*xiDyad(m,o)
     endif  
   enddo; enddo; enddo
   gamma_hat(1:3,1:3,1:3,1:3, 1,1,1) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
 endif

end subroutine utilities_updateGamma


!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> @details Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTforward()
 use math

 implicit none
 integer(pInt)  :: row, column                                                                      ! if debug FFTW, compare 3D array field of row and column
 integer(pInt), dimension(2:3,2) :: Nyquist                                                         ! highest frequencies to be removed (1 if even, 2 if odd)
 real(pReal), dimension(2) :: myRand                                                                ! random numbers

!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
 if (debugFFTW) then
   call random_number(myRand)                                                                       ! two numbers: 0 <= x < 1
   row    = nint(myRand(1)*2_pReal + 1_pReal,pInt)
   column = nint(myRand(2)*2_pReal + 1_pReal,pInt)
   scalarField_real = cmplx(field_real(1:grid(1),1:grid(2),1:grid(3),row,column),0.0_pReal,pReal)   ! store the selected component 
 endif

!--------------------------------------------------------------------------------------------------
! doing the FFT
 call fftw_execute_dft_r2c(planForth,field_real,field_fourier)
  
!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
 if (debugFFTW) then
   call fftw_execute_dft(planDebugForth,scalarField_real,scalarField_fourier)
   where(abs(scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3))) > tiny(1.0_pReal))               ! avoid division by zero
     scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3)) = &  
                         (scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3))-&
                                field_fourier(1:grid1Red,1:grid(2),1:grid(3),row,column))/&
                          scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3))
   else where
     scalarField_real = cmplx(0.0,0.0,pReal)
   end where
   write(6,'(/,a,i1,1x,i1,a)') ' .. checking FT results of compontent ', row, column, ' ..'
   write(6,'(/,a,2(es11.4,1x))')  ' max FT relative error = ',&                                     ! print real and imaginary part seperately
     maxval(real (scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3)))),&
     maxval(aimag(scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3))))
   flush(6)
 endif

!--------------------------------------------------------------------------------------------------
! removing highest frequencies
 Nyquist(2,1:2) = [grid(2)/2_pInt + 1_pInt, grid(2)/2_pInt + 1_pInt + mod(grid(2),2_pInt)]
 Nyquist(3,1:2) = [grid(3)/2_pInt + 1_pInt, grid(3)/2_pInt + 1_pInt + mod(grid(3),2_pInt)]

 if(grid(1)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   field_fourier (grid1Red,  1:grid(2),                 1:grid(3),                 1:3,1:3) &
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 if(grid(2)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   field_fourier (1:grid1Red,Nyquist(2,1):Nyquist(2,2),1:grid(3),                 1:3,1:3) & 
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 if(grid(3)/=1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
   field_fourier (1:grid1Red,1:grid(2),                 Nyquist(3,1):Nyquist(3,2),1:3,1:3) &
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
end subroutine utilities_FFTforward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @details Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTbackward()
 use math                                                                                           !< must use the whole module for use of FFTW

 implicit none
 integer(pInt) :: row, column                                                                       !< if debug FFTW, compare 3D array field of row and column
 integer(pInt) ::  i, j, k, m, n
 real(pReal), dimension(2) :: myRand 

!--------------------------------------------------------------------------------------------------
! unpack FFT data for conj complex symmetric part. This data is not transformed when using c2r
 if (debugFFTW) then
   call random_number(myRand)                                                                       ! two numbers: 0 <= x < 1
   row    = nint(myRand(1)*2_pReal + 1_pReal,pInt)
   column = nint(myRand(2)*2_pReal + 1_pReal,pInt)
   scalarField_fourier(1:grid1Red,1:grid(2),1:grid(3)) &
                                         = field_fourier(1:grid1Red,1:grid(2),1:grid(3),row,column)
   do i = 0_pInt, grid(1)/2_pInt-2_pInt + mod(grid(1),2_pInt)
    m = 1_pInt
    do k = 1_pInt, grid(3)
      n = 1_pInt
      do j = 1_pInt, grid(2)
        scalarField_fourier(grid(1)-i,j,k) = conjg(scalarField_fourier(2+i,n,m))
        if(n == 1_pInt) n = grid(2) + 1_pInt
        n = n-1_pInt
     enddo
     if(m == 1_pInt) m = grid(3) + 1_pInt
     m = m -1_pInt
   enddo; enddo
 endif
 
!--------------------------------------------------------------------------------------------------
! doing the iFFT
 call fftw_execute_dft_c2r(planBack,field_fourier,field_real)                                       ! back transform of fluct deformation gradient

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
 if (debugFFTW) then
   call fftw_execute_dft(planDebugBack,scalarField_fourier,scalarField_real)
   where(abs(real(scalarField_real,pReal)) > tiny(1.0_pReal))                                       ! avoid division by zero
     scalarField_real = (scalarField_real &
                         - cmplx(field_real(1:grid(1),1:grid(2),1:grid(3),row,column), 0.0, pReal))/ &
                         scalarField_real
   else where
     scalarField_real = cmplx(0.0,0.0,pReal)
   end where
   write(6,'(/,a,i1,1x,i1,a)') ' ... checking iFT results of compontent ', row, column, ' ..'
   write(6,'(/,a,es11.4)')     ' max iFT relative error = ', maxval(real(scalarField_real,pReal))
   flush(6)
 endif
 
 field_real = field_real * wgt                                                                      ! normalize the result by number of elements

end subroutine utilities_FFTbackward


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution with inverse laplace kernel
!--------------------------------------------------------------------------------------------------
subroutine utilities_inverseLaplace()
 use math, only: &
   math_inv33, &
   PI

 implicit none  
 integer(pInt) :: i, j, k
 integer(pInt), dimension(3) :: k_s

 write(6,'(/,a)') ' ... doing inverse laplace .................................................'
 flush(6)
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation (mechanical equilibrium)
 do k = 1_pInt, grid(3)
  k_s(3) = k - 1_pInt
  if(k > grid(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - grid(3)                                        ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
  do j = 1_pInt, grid(2)
    k_s(2) = j - 1_pInt
    if(j > grid(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - grid(2)                                      ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
    do i = 1_pInt, grid1Red
      k_s(1) = i - 1_pInt  
      if (any(k_s /= 0_pInt)) field_fourier(i,j,k, 1:3,1:3) =  &
                              field_fourier(i,j,k, 1:3,1:3)/ &
                                        cmplx(-sum((2.0_pReal*PI*k_s/geomSize)*&
                                                   (2.0_pReal*PI*k_s/geomSize)),0.0_pReal,pReal)   ! symmetry, junst running from 0,1,...,N/2,N/2+1
enddo; enddo; enddo
field_fourier(1,1,1,1:3,1:3) = cmplx(0.0_pReal,0.0_pReal,pReal)

end subroutine utilities_inverseLaplace
 

!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real, ensuring that average value = fieldAim
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierConvolution(fieldAim)
 use numerics, only: &
   memory_efficient
 use math, only: &
   math_inv33

 implicit none  
 real(pReal), intent(in), dimension(3,3) :: fieldAim                                                !< desired average value of the field after convolution
 real(pReal),             dimension(3,3) :: xiDyad, temp33_Real
 real(pReal)                             :: filter                                                  !< weighting of current component
 complex(pReal),          dimension(3,3) :: temp33_complex
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o

 write(6,'(/,a)') ' ... doing convolution .....................................................'
 flush(6)
 
!--------------------------------------------------------------------------------------------------
! do the actual spectral method calculation (mechanical equilibrium)
 if(memory_efficient) then                                                                          ! memory saving version, on-the-fly calculation of gamma_hat
   do k = 1_pInt, grid(3); do j = 1_pInt, grid(2) ;do i = 1_pInt, grid1Red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,1:3,m,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       filter = utilities_getFilter(xi(1:3,i,j,k))                                                  ! weighting factor computed by getFilter function
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o, 1,1,1) =  filter*temp33_Real(l,n)*xiDyad(m,o)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, 1,1,1) * field_fourier(i,j,k,1:3,1:3))
       field_fourier(i,j,k,1:3,1:3) = temp33_Complex 
     endif             
   enddo; enddo; enddo
 else                                                                                               ! use precalculated gamma-operator
   do k = 1_pInt, grid(3);  do j = 1_pInt, grid(2);  do i = 1_pInt,grid1Red
     forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
       temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, i,j,k) * field_fourier(i,j,k,1:3,1:3))
     field_fourier(i,j,k, 1:3,1:3) = temp33_Complex
   enddo; enddo; enddo
 endif
 field_fourier(1,1,1,1:3,1:3) = cmplx(fieldAim*real(product(grid),pReal),0.0_pReal,pReal)           ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1  

end subroutine utilities_fourierConvolution
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS()
 use math                                                                                           !< must use the whole module for use of FFTW

 implicit none
 integer(pInt) :: i, j, k 
 real(pReal) :: &
   err_real_div_RMS, &                                                                              !< RMS of divergence in real space
   err_div_max, &                                                                                   !< maximum value of divergence in Fourier space
   err_real_div_max                                                                                 !< maximum value of divergence in real space
 complex(pReal), dimension(3) ::  temp3_complex

 write(6,'(/,a)') ' ... calculating divergence ................................................'
 flush(6)

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
 utilities_divergenceRMS = 0.0_pReal
 do k = 1_pInt, grid(3); do j = 1_pInt, grid(2)
   do i = 2_pInt, grid1Red -1_pInt                                                                  ! Has somewhere a conj. complex counterpart. Therefore count it twice.
     utilities_divergenceRMS = utilities_divergenceRMS &
           + 2.0_pReal*(sum (real(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),&               ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                           xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&                    ! --> sum squared L_2 norm of vector 
                       +sum(aimag(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),& 
                                                          xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
   enddo
   utilities_divergenceRMS = utilities_divergenceRMS &                                              ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart (if grid(1) /= 1)
              + sum( real(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                  xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
              + sum(aimag(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                  xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
              + sum( real(math_mul33x3_complex(field_fourier(grid1Red,j,k,1:3,1:3),&
                                                  xi(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal)&
              + sum(aimag(math_mul33x3_complex(field_fourier(grid1Red,j,k,1:3,1:3),&
                                                  xi(1:3,grid1Red,j,k))*TWOPIIMG)**2.0_pReal)
 enddo; enddo
 if(grid(1) == 1_pInt) utilities_divergenceRMS = utilities_divergenceRMS * 0.5_pReal                ! counted twice in case of grid(1) == 1
 utilities_divergenceRMS = sqrt(utilities_divergenceRMS) * wgt                                      ! RMS in real space calculated with Parsevals theorem from Fourier space

!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
 if (debugDivergence) then                                                                          ! calculate divergence again
   err_div_max = 0.0_pReal
   do k = 1_pInt, grid(3); do j = 1_pInt, grid(2); do i = 1_pInt, grid1Red
     temp3_Complex = math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3)*wgt,&                        ! weighting P_fourier
                                             xi(1:3,i,j,k))*TWOPIIMG
     err_div_max = max(err_div_max,sum(abs(temp3_Complex)**2.0_pReal))
     divFourier(i,j,k,1:3) = temp3_Complex                                                          ! need divergence NOT squared
   enddo; enddo; enddo
   
   call fftw_execute_dft_c2r(planDiv,divFourier,divReal)                                            ! already weighted

   err_real_div_RMS = sqrt(wgt*sum(divReal**2.0_pReal))                                             ! RMS in real space
   err_real_div_max = sqrt(maxval(sum(divReal**2.0_pReal,dim=4)))                                   ! max in real space                                       
   err_div_max      = sqrt(    err_div_max)                                                         ! max in Fourier space
   
   write(6,'(/,1x,a,es11.4)')      'error divergence  FT  RMS = ',utilities_divergenceRMS
   write(6,'(1x,a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
   write(6,'(1x,a,es11.4)')        'error divergence  FT  max = ',err_div_max
   write(6,'(1x,a,es11.4)')        'error divergence Real max = ',err_real_div_max
   flush(6)
 endif

end function utilities_divergenceRMS
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate max of curl of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_curlRMS()
 use math                                                                                           !< must use the whole module for use of FFTW

 implicit none
 integer(pInt)  ::  i, j, k, l 
 complex(pReal), dimension(3,3) ::  curl_fourier

 write(6,'(/,a)') ' ... calculating curl ......................................................'
 flush(6)
!--------------------------------------------------------------------------------------------------
! calculating max curl criterion in Fourier space
 utilities_curlRMS = 0.0_pReal
 
 do k = 1_pInt, grid(3); do j = 1_pInt, grid(2); 
 do i = 2_pInt, grid1Red - 1_pInt
   do l = 1_pInt, 3_pInt
     curl_fourier(l,1) = (+field_fourier(i,j,k,l,3)*xi(2,i,j,k)&
                          -field_fourier(i,j,k,l,2)*xi(3,i,j,k))*TWOPIIMG
     curl_fourier(l,2) = (+field_fourier(i,j,k,l,1)*xi(3,i,j,k)&
                          -field_fourier(i,j,k,l,3)*xi(1,i,j,k))*TWOPIIMG
     curl_fourier(l,3) = (+field_fourier(i,j,k,l,2)*xi(1,i,j,k)&
                          -field_fourier(i,j,k,l,1)*xi(2,i,j,k))*TWOPIIMG
   enddo
   utilities_curlRMS = utilities_curlRMS + &
                       2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)
 enddo 
 do l = 1_pInt, 3_pInt
   curl_fourier = (+field_fourier(1,j,k,l,3)*xi(2,1,j,k)&
                   -field_fourier(1,j,k,l,2)*xi(3,1,j,k))*TWOPIIMG
   curl_fourier = (+field_fourier(1,j,k,l,1)*xi(3,1,j,k)&
                   -field_fourier(1,j,k,l,3)*xi(1,1,j,k))*TWOPIIMG
   curl_fourier = (+field_fourier(1,j,k,l,2)*xi(1,1,j,k)&
                   -field_fourier(1,j,k,l,1)*xi(2,1,j,k))*TWOPIIMG
 enddo
 utilities_curlRMS = utilities_curlRMS + &
                     2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal)   
 do l = 1_pInt, 3_pInt  
   curl_fourier = (+field_fourier(grid1Red,j,k,l,3)*xi(2,grid1Red,j,k)&
                   -field_fourier(grid1Red,j,k,l,2)*xi(3,grid1Red,j,k))*TWOPIIMG
   curl_fourier = (+field_fourier(grid1Red,j,k,l,1)*xi(3,grid1Red,j,k)&
                   -field_fourier(grid1Red,j,k,l,3)*xi(1,grid1Red,j,k))*TWOPIIMG
   curl_fourier = (+field_fourier(grid1Red,j,k,l,2)*xi(1,grid1Red,j,k)&
                   -field_fourier(grid1Red,j,k,l,1)*xi(2,grid1Red,j,k))*TWOPIIMG
 enddo
 utilities_curlRMS = utilities_curlRMS + &
                     2.0_pReal*sum(real(curl_fourier)**2.0_pReal + aimag(curl_fourier)**2.0_pReal) 
 enddo; enddo
 utilities_curlRMS = sqrt(utilities_curlRMS) * wgt
 if(grid(1) == 1_pInt) utilities_curlRMS = utilities_curlRMS * 0.5_pReal                           ! counted twice in case of grid(1) == 1

end function utilities_curlRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculates mask compliance tensor used to adjust F to fullfill stress BC
!--------------------------------------------------------------------------------------------------
function utilities_maskedCompliance(rot_BC,mask_stress,C)
 use IO, only: &
   IO_error
 use math, only: &
   math_Plain3333to99, &
   math_plain99to3333, &
   math_rotate_forward3333, &
   math_rotate_forward33, &
   math_invert

 implicit none
 real(pReal),              dimension(3,3,3,3) :: utilities_maskedCompliance                        !< masked compliance
 real(pReal), intent(in) , dimension(3,3,3,3) :: C                                                !< current average stiffness
 real(pReal), intent(in) , dimension(3,3)     :: rot_BC                                           !< rotation of load frame
 logical,     intent(in),  dimension(3,3)     :: mask_stress                                     !< mask of stress BC
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

   if(debugGeneral) then 
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
   if(debugGeneral .or. errmatinv) then                                                             ! report
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
 if(debugGeneral) &                                                                                 ! report
   write(6,'(/,a,/,9(9(2x,f12.7,1x)/),/)',advance='no') ' Masked Compliance (load) * GPa =', &
                                                    transpose(temp99_Real*1.e9_pReal)
 flush(6)
 utilities_maskedCompliance = math_Plain99to3333(temp99_Real)

end function utilities_maskedCompliance 


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(F_lastInc,F,temperature,timeinc,&
                                          P,C_volAvg,C_minmaxAvg,P_av,forwardData,rotation_BC)
 use debug, only: &
   debug_reset, &
   debug_info
 use numerics, only: &
   usePingPong
 use math, only: &
   math_transpose33, &
   math_rotate_forward33, &
   math_det33
 use FEsolving, only: &
   restartWrite
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_COLLECT, &
   CPFEM_CALCRESULTS, &
   CPFEM_AGERESULTS, &
   CPFEM_BACKUPJACOBIAN, &
   CPFEM_RESTOREJACOBIAN
 use homogenization, only: &
   materialpoint_F0, &
   materialpoint_F, &
   materialpoint_Temperature, &
   materialpoint_P, &
   materialpoint_dPdF
 
 implicit none
 real(pReal), intent(inout)                                      :: temperature                     !< temperature (no field)
 real(pReal), intent(in), dimension(3,3,grid(1),grid(2),grid(3)) :: &
   F_lastInc, &                                                                                     !< target deformation gradient
   F                                                                                                !< previous deformation gradient
 real(pReal), intent(in)                                         :: timeinc                         !< loading time
 logical,     intent(in)                                         :: forwardData                     !< age results
 real(pReal), intent(in), dimension(3,3)                         :: rotation_BC                     !< rotation of load frame
 
 real(pReal),intent(out), dimension(3,3,3,3)                     :: C_volAvg, C_minmaxAvg           !< average stiffness
 real(pReal),intent(out), dimension(3,3)                         :: P_av                            !< average PK stress
 real(pReal),intent(out), dimension(3,3,grid(1),grid(2),grid(3)) :: P                               !< PK stress
 
 integer(pInt) :: &
   calcMode, &                                                                                      !< CPFEM mode for calculation
   collectMode, &                                                                                   !< CPFEM mode for collection
   j,k
 real(pReal), dimension(3,3,3,3) :: max_dPdF, min_dPdF
 real(pReal)   :: max_dPdF_norm, min_dPdF_norm, defgradDetMin, defgradDetMax, defgradDet

 write(6,'(/,a)') ' ... evaluating constitutive response ......................................'
 calcMode    = CPFEM_CALCRESULTS
 collectMode = CPFEM_COLLECT
 if (forwardData) then                                                                              ! aging results
   calcMode    = ior(calcMode,    CPFEM_AGERESULTS)             
   collectMode = ior(collectMode, CPFEM_BACKUPJACOBIAN)
 endif
 if (cutBack) then                                                                                  ! restore saved variables
  collectMode = ior(collectMode , CPFEM_RESTOREJACOBIAN)
  collectMode = iand(collectMode, not(CPFEM_BACKUPJACOBIAN))
  calcMode    = iand(calcMode,    not(CPFEM_AGERESULTS)) 
 endif

 call CPFEM_general(collectMode,usePingPong,F_lastInc(1:3,1:3,1,1,1),F(1:3,1:3,1,1,1), &            ! collect mode handles Jacobian backup / restoration
                   temperature,timeinc,1_pInt,1_pInt)
 
 materialpoint_F0 = reshape(F_lastInc, [3,3,1,product(grid)])
 materialpoint_F  = reshape(F,         [3,3,1,product(grid)])
 materialpoint_Temperature = temperature

 call debug_reset()

!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
 if(debugGeneral) then
   defgradDetMax = -huge(1.0_pReal)
   defgradDetMin = +huge(1.0_pReal)
   do j = 1_pInt, product(grid)
     defgradDet = math_det33(materialpoint_F(1:3,1:3,1,j))
     defgradDetMax = max(defgradDetMax,defgradDet)
     defgradDetMin = min(defgradDetMin,defgradDet) 
   end do
   write(6,'(a,1x,es11.4)') ' max determinant of deformation =', defgradDetMax
   write(6,'(a,1x,es11.4)') ' min determinant of deformation =', defgradDetMin
   flush(6)
 endif
  
 call CPFEM_general(calcMode,usePingPong,F_lastInc(1:3,1:3,1,1,1), F(1:3,1:3,1,1,1), &              ! first call calculates everything
                    temperature,timeinc,1_pInt,1_pInt)
 
 max_dPdF = 0.0_pReal
 max_dPdF_norm = 0.0_pReal
 min_dPdF = huge(1.0_pReal)
 min_dPdF_norm = huge(1.0_pReal)
 do k = 1_pInt, product(grid)
   if (max_dPdF_norm < sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)) then
     max_dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)
     max_dPdF_norm = sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)
   endif  
   if (min_dPdF_norm > sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)) then
     min_dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)
     min_dPdF_norm = sum(materialpoint_dPdF(1:3,1:3,1:3,1:3,1,k)**2.0_pReal)
   endif  
 end do

 P = reshape(materialpoint_P, [3,3,grid(1),grid(2),grid(3)])
 C_volAvg = sum(sum(materialpoint_dPdF,dim=6),dim=5) * wgt
 C_minmaxAvg = 0.5_pReal*(max_dPdF + min_dPdF)
 
 call debug_info()
 
 restartWrite = .false.                                                                             ! reset restartWrite status
 cutBack = .false.                                                                                  ! reset cutBack status
 
 P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt                                                    ! average of P 
 if (debugRotation) &
 write(6,'(/,a,/,3(3(2x,f12.4,1x)/))',advance='no') ' Piola--Kirchhoff stress (lab) / MPa =',&
                                                     math_transpose33(P_av)*1.e-6_pReal
 P_av = math_rotate_forward33(P_av,rotation_BC)
 write(6,'(/,a,/,3(3(2x,f12.4,1x)/))',advance='no') ' Piola--Kirchhoff stress / MPa =',&
                                                     math_transpose33(P_av)*1.e-6_pReal

end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief calculates forward rate, either guessing or just add delta/timeinc
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(avRate,timeinc_old,guess,field_lastInc,field)

 implicit none
 real(pReal), intent(in), dimension(3,3)                      :: avRate                             !< homogeneous addon
 real(pReal), intent(in) :: &
   timeinc_old                                                                                      !< timeinc of last step
 logical, intent(in) :: &
   guess                                                                                            !< guess along former trajectory
 real(pReal), intent(in), dimension(3,3,grid(1),grid(2),grid(3)) :: &
   field_lastInc, &                                                                                 !< data of previous step
   field                                                                                            !< data of current step
 real(pReal),             dimension(3,3,grid(1),grid(2),grid(3)) :: utilities_calculateRate
 
 if(guess) then
   utilities_calculateRate = (field-field_lastInc) / timeinc_old
 else
   utilities_calculateRate = spread(spread(spread(avRate,3,grid(1)),4,grid(2)),5,grid(3))
 endif

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, if aim is given, 
!> ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
pure function utilities_forwardField(timeinc,field_lastInc,rate,aim)

 implicit none
 real(pReal), intent(in) :: & 
   timeinc                                                                                          !< timeinc of current step
 real(pReal), intent(in),           dimension(3,3,grid(1),grid(2),grid(3)) :: &
   field_lastInc, &                                                                                 !< initial field
   rate                                                                                             !< rate by which to forward
 real(pReal), intent(in), optional, dimension(3,3) :: &
   aim                                                                                              !< average field value aim
 real(pReal),                       dimension(3,3,grid(1),grid(2),grid(3)) :: utilities_forwardField
 real(pReal),                       dimension(3,3)                      :: fieldDiff                !< <a + adot*t> - aim
 
 utilities_forwardField = field_lastInc + rate*timeinc
 if (present(aim)) then                                                                             !< correct to match average
   fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt - aim
   utilities_forwardField = utilities_forwardField - &
                          spread(spread(spread(fieldDiff,3,grid(1)),4,grid(2)),5,grid(3))
 endif

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief calculates filter for fourier convolution depending on type given in numerics.config
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_getFilter(k)
 use IO, only: &
   IO_error
 use numerics, only: &                        
   myfilter
 use math, only: &
   PI
  
 implicit none
 real(pReal),intent(in), dimension(3) :: k                                                          !< indices of frequency
  
 utilities_getFilter = 1.0_pReal

 select case (myfilter)
    case ('none')                                                                                   ! default, no weighting
    case ('cosine')                                                                                 ! cosine curve with 1 for avg and zero for highest freq
      utilities_getFilter = product(1.0_pReal + cos(PI*k*scaledGeomSize/grid))/8.0_pReal
    case ('gradient')                                                                               ! gradient, might need grid scaling as for cosine filter
      utilities_getFilter = 1.0_pReal/(1.0_pReal + &
                                       (k(1)*k(1) + k(2)*k(2) + k(3)*k(3)))
    case default
      call IO_error(error_ID = 892_pInt, ext_msg = trim(myfilter))
  end select 

end function utilities_getFilter


!--------------------------------------------------------------------------------------------------
!> @brief cleans up
!--------------------------------------------------------------------------------------------------
subroutine utilities_destroy()
 use math

 implicit none

 if (debugDivergence) call fftw_destroy_plan(planDiv)
 if (debugFFTW) call fftw_destroy_plan(planDebugForth)
 if (debugFFTW) call fftw_destroy_plan(planDebugBack)
 call fftw_destroy_plan(planForth)
 call fftw_destroy_plan(planBack)

end subroutine utilities_destroy


end module DAMASK_spectral_utilities
