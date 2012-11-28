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
 logical,       public                                         :: cutBack =.false.                  !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 real(pReal),   public,  dimension(:,:,:,:,:),     pointer     :: field_real                        !< real representation (some stress or deformation) of field_fourier
 complex(pReal),private, dimension(:,:,:,:,:),     pointer     :: field_fourier                     !< field on which the Fourier transform operates
 real(pReal),   private, dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat                         !< gamma operator (field) for spectral method
 real(pReal),   private, dimension(:,:,:,:),       allocatable :: xi                                !< wave vector field for divergence and for gamma operator
 real(pReal),   private, dimension(3,3,3,3)                    :: C_ref                             !< reference stiffness

!--------------------------------------------------------------------------------------------------
! debug fftw 
 complex(pReal),private, dimension(:,:,:), pointer :: scalarField_real, &                           !< scalar field real representation for debug of FFTW
                                                      scalarField_fourier                           !< scalar field complex representation for debug of FFTW
 
!--------------------------------------------------------------------------------------------------
! debug divergence
 real(pReal),   private, dimension(:,:,:,:), pointer     :: divergence_real                         !< scalar field real representation for debugging divergence calculation
 complex(pReal),private, dimension(:,:,:,:), pointer     :: divergence_fourier                      !< scalar field real representation for debugging divergence calculation
 real(pReal),   private, dimension(:,:,:,:), allocatable :: divergence_post                         !< data of divergence calculation using function from core modules (serves as a reference)

!--------------------------------------------------------------------------------------------------
! plans for FFTW
 type(C_PTR),   private                            :: plan_scalarField_forth, plan_scalarField_back !< plans for FFTW in case of debugging the Fourier transform
 type(C_PTR),   private                            :: plan_forward, plan_backward                   !< plans for FFTW
 type(C_PTR),   private                            :: plan_divergence                               !< plan for FFTW in case of debugging divergence calculation

!--------------------------------------------------------------------------------------------------
! variables controlling debugging
 logical, public :: &
   debugGeneral, &                                                                                  !< general debugging of spectral solver
   debugDivergence, &                                                                               !< debugging of divergence calculation (comparison to function used for post processing)
   debugRestart, &                                                                                  !< debbuging of restart features
   debugFFTW, &                                                                                     !< doing additional FFT on scalar field and compare to results of strided 3D FFT
   debugRotation                                                                                    !< also printing out results in lab frame

!--------------------------------------------------------------------------------------------------
! derived types
 type tSolutionState                                                                                 !< return type of solution from spectral solver variants
   logical       :: converged         = .true.   
   logical       :: regrid            = .false.   
   logical       :: termIll           = .false.   
   integer(pInt) :: iterationsNeeded  = 0_pInt
 end type tSolutionState

 type tBoundaryCondition                                                                             !< set of parameters defining a boundary condition
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
   utilities_divergenceRMS, &
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
 use IO, only: &
   IO_error
 use numerics, only: &                        
   DAMASK_NumThreadsInt, &
   fftw_planner_flag, &
   fftw_timelimit, &
   memory_efficient
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_levelBasic, &
   debug_spectralDivergence, &
   debug_spectralRestart, &
   debug_spectralFFTW, &
   debug_spectralRotation
 use mesh, only: &
   res, &
   res1_red, &
   virt_dim
 use math                                                                                           ! must use the whole module for use of FFTW
 
 implicit none
 integer(pInt)               :: i, j, k
 integer(pInt), dimension(3) :: k_s
 type(C_PTR) :: &
   tensorField, &                                                                                   !< field cotaining data for FFTW in real and fourier space (in place)
   scalarField_realC, &                                                                             !< field cotaining data for FFTW in real space when debugging FFTW (no in place)
   scalarField_fourierC, &                                                                          !< field cotaining data for FFTW in fourier space when debugging FFTW (no in place)
   divergence                                                                                       !< field cotaining data for FFTW in real and fourier space when debugging divergence (in place)
 
 write(6,'(/,a)') ' <<<+-  DAMASK_spectral_utilities init  -+>>>'
 write(6,'(a)')   ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)')   ''

!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_spectral),debug_levelBasic)         /= 0
 debugDivergence = iand(debug_level(debug_spectral),debug_spectralDivergence) /= 0
 debugRestart    = iand(debug_level(debug_spectral),debug_spectralRestart)    /= 0
 debugFFTW       = iand(debug_level(debug_spectral),debug_spectralFFTW)       /= 0
 debugRotation   = iand(debug_level(debug_spectral),debug_spectralRotation)   /= 0
 
!--------------------------------------------------------------------------------------------------
! allocation
 allocate (xi(3,res1_red,res(2),res(3)),source = 0.0_pReal)                                         ! frequencies, only half the size for first dimension
 tensorField = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))                      ! allocate aligned data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, field_real,      [ res(1)+2_pInt,res(2),res(3),3,3])                 ! place a pointer for a real representation on tensorField
 call c_f_pointer(tensorField, field_fourier,   [ res1_red,     res(2),res(3),3,3])                 ! place a pointer for a complex representation on tensorField

!--------------------------------------------------------------------------------------------------
! general initialization of FFTW (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=808_pInt)                         ! check for correct precision in C
!$ if(DAMASK_NumThreadsInt > 0_pInt) then
!$   i = fftw_init_threads()                                                                        ! returns 0 in case of problem
!$   if (i == 0_pInt) call IO_error(error_ID = 809_pInt)
!$   call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
!$ endif
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans for the convolution
 plan_forward =  fftw_plan_many_dft_r2c(3, [res(3),res(2) ,res(1)],        9,&                      ! dimensions,  logical length in each dimension in reversed order,  no. of transforms
                               field_real, [res(3),res(2) ,res(1)+2_pInt],&                         ! input data,  physical length in each dimension in reversed order
                                        1,  res(3)*res(2)*(res(1)+2_pInt),&                         ! striding,    product of physical length in the 3 dimensions
                            field_fourier, [res(3),res(2) ,res1_red],&                              ! output data, physical length in each dimension in reversed order
                                        1,  res(3)*res(2)* res1_red,       fftw_planner_flag)       ! striding,    product of physical length in the 3 dimensions,      planner precision

 plan_backward = fftw_plan_many_dft_c2r(3, [res(3),res(2) ,res(1)],        9,&                      ! dimensions,  logical length in each dimension in reversed order,  no. of transforms
                            field_fourier, [res(3),res(2) ,res1_red],&                              ! input data,  physical length in each dimension in reversed order
                                        1,  res(3)*res(2)* res1_red,&                               ! striding,    product of physical length in the 3 dimensions
                               field_real, [res(3),res(2) ,res(1)+2_pInt],&                         ! output data, physical length in each dimension in reversed order
                                        1,  res(3)*res(2)*(res(1)+2_pInt), fftw_planner_flag)       ! striding,    product of physical length in the 3 dimensions,      planner precision

!--------------------------------------------------------------------------------------------------
! depending on debug options, allocate more memory and create additional plans 
 if (debugDivergence) then
   divergence = fftw_alloc_complex(int(res1_red*res(2)*res(3)*3_pInt,C_SIZE_T))
   call c_f_pointer(divergence, divergence_real,    [ res(1)+2_pInt,res(2),res(3),3])
   call c_f_pointer(divergence, divergence_fourier, [ res1_red,     res(2),res(3),3])
   allocate (divergence_post(res(1),res(2),res(3),3),source = 0.0_pReal)
   plan_divergence = fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],3,&
                           divergence_fourier,[ res(3),res(2) ,res1_red],&
                                            1,  res(3)*res(2)* res1_red,&
                              divergence_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                            1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 endif

 if (debugFFTW) then
   scalarField_realC    = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))                    ! allocate data for real representation (no in place transform)
   scalarField_fourierC = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))                    ! allocate data for fourier representation (no in place transform)
   call c_f_pointer(scalarField_realC,    scalarField_real,    [res(1),res(2),res(3)])              ! place a pointer for a real representation
   call c_f_pointer(scalarField_fourierC, scalarField_fourier, [res(1),res(2),res(3)])              ! place a pointer for a fourier representation
   plan_scalarField_forth = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 ! reversed order (C style)
                                      scalarField_real,scalarField_fourier,-1,fftw_planner_flag)    ! input, output, forward FFT(-1), planner precision
   plan_scalarField_back  = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 ! reversed order (C style)
                                      scalarField_fourier,scalarField_real,+1,fftw_planner_flag)    ! input, output, backward (1), planner precision
 endif 

 if (debugGeneral) write(6,'(a)') 'FFTW initialized'
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 do k = 1_pInt, res(3)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)                                          ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2)                                      ! running from 0,1,...,N/2,N/2+1,-N/2,-N/2+1,...,-1
         do i = 1_pInt, res1_red
           k_s(1) = i - 1_pInt                                                                      ! symmetry, junst running from 0,1,...,N/2,N/2+1
           xi(1:3,i,j,k) = real(k_s, pReal)/virt_dim                                                ! if divergence_correction is set, frequencies are calculated on unit length
 enddo; enddo; enddo
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(3,3,3,3,1,1,1), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(3,3,3,3,res1_red ,res(2),res(3)), source =0.0_pReal)     
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
  IO_write_jobBinaryFile
 use numerics, only: &
   memory_efficient
 use math, only: &
   math_inv33
 use mesh, only: &
   res, &
   res1_red

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
   write(6,'(a)') 'writing reference stiffness to file'
   call IO_write_jobBinaryFile(777,'C_ref',size(C_ref))
   write (777,rec=1) C_ref
   close(777)
 endif
 
 if(.not. memory_efficient) then                                                   
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
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
!> @detailed Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTforward(row,column)
 use math
 use mesh, only : &
   virt_dim, &
   res, &
   res1_red

 implicit none
 integer(pInt), intent(in), optional :: row, column                                                 !< if debug FFTW, compare 3D array field of row and column
  
!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
  if (debugFFTW) then
    if (.not. present(row) .or. .not. present(column)) stop
    scalarField_real(1:res(1),1:res(2),1:res(3)) =&                                                 ! store the selected component
           cmplx(field_real(1:res(1),1:res(2),1:res(3),row,column),0.0_pReal,pReal)
  endif
  
!--------------------------------------------------------------------------------------------------
! call function to calculate divergence from math (for post processing) to check results
  if (debugDivergence) &
    divergence_post = math_divergenceFFT(virt_dim,field_real(1:res(1),1:res(2),1:res(3),1:3,1:3))   ! some elements are padded
  
!--------------------------------------------------------------------------------------------------
! doing the FFT
  call fftw_execute_dft_r2c(plan_forward,field_real,field_fourier)
  
!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
  if (debugFFTW) then
    call fftw_execute_dft(plan_scalarField_forth,scalarField_real,scalarField_fourier)
    write(6,'(a,i1,1x,i1)') 'checking FT results of compontent ', row, column
    write(6,'(a,2(es11.4,1x))')  'max FT relative error = ',&                                       ! print real and imaginary part seperately
      maxval( real((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-& 
                          field_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                    scalarField_fourier(1:res1_red,1:res(2),1:res(3)))), &
      maxval(aimag((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-&
                          field_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                    scalarField_fourier(1:res1_red,1:res(2),1:res(3))))
  endif

!--------------------------------------------------------------------------------------------------
! removing highest frequencies
 field_fourier  (  res1_red,1:res(2) ,             1:res(3)              ,1:3,1:3)&
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 field_fourier  (1:res1_red,  res(2)/2_pInt+1_pInt,1:res(3)              ,1:3,1:3)& 
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
 if(res(3)>1_pInt) &                                                                               ! do not delete the whole slice in case of 2D calculation
  field_fourier (1:res1_red,1:res(2),                res(3)/2_pInt+1_pInt,1:3,1:3)&
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
end subroutine utilities_FFTforward


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> @detailed Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine utilities_FFTbackward(row,column)
 use math                                                                                           !< must use the whole module for use of FFTW
 use mesh, only: &
   wgt, &
   res, &
   res1_red

 implicit none
 integer(pInt), intent(in), optional :: row, column                                                 !< if debug FFTW, compare 3D array field of row and column
 integer(pInt) ::  i, j, k, m, n
  
!--------------------------------------------------------------------------------------------------
! unpack FFT data for conj complex symmetric part. This data is not transformed when using c2r
 if (debugFFTW) then
   scalarField_fourier = field_fourier(1:res1_red,1:res(2),1:res(3),row,column)
   do i = 0_pInt, res(1)/2_pInt-2_pInt
    m = 1_pInt
    do k = 1_pInt, res(3)
      n = 1_pInt
      do j = 1_pInt, res(2)
        scalarField_fourier(res(1)-i,j,k) = conjg(scalarField_fourier(2+i,n,m))
        if(n == 1_pInt) n = res(2) + 1_pInt
        n = n-1_pInt
     enddo
     if(m == 1_pInt) m = res(3) + 1_pInt
     m = m -1_pInt
   enddo; enddo
 endif
 
!--------------------------------------------------------------------------------------------------
! doing the iFFT
 call fftw_execute_dft_c2r(plan_backward,field_fourier,field_real)                                  ! back transform of fluct deformation gradient

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
  if (debugFFTW) then
    write(6,'(a,i1,1x,i1)') 'checking iFT results of compontent ', row, column
    call fftw_execute_dft(plan_scalarField_back,scalarField_fourier,scalarField_real)
    write(6,'(a,es11.4)') 'max iFT relative error = ',&
        maxval((real(scalarField_real(1:res(1),1:res(2),1:res(3)))-&
                field_real(1:res(1),1:res(2),1:res(3),row,column))/&
                real(scalarField_real(1:res(1),1:res(2),1:res(3))))
  endif
  
  field_real = field_real * wgt                                                                     ! normalize the result by number of elements

!--------------------------------------------------------------------------------------------------
! calculate some additional output
 ! if(debugGeneral) then
 !   maxCorrectionSkew = 0.0_pReal
 !   maxCorrectionSym  = 0.0_pReal
 !   temp33_Real = 0.0_pReal
 !   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
 !     maxCorrectionSym  = max(maxCorrectionSym,&
 !                             maxval(math_symmetric33(field_real(i,j,k,1:3,1:3))))
 !     maxCorrectionSkew = max(maxCorrectionSkew,&
 !                             maxval(math_skew33(field_real(i,j,k,1:3,1:3))))
 !     temp33_Real = temp33_Real + field_real(i,j,k,1:3,1:3)
 !   enddo; enddo; enddo
 !   write(6,'(a,1x,es11.4)') 'max symmetric correction of deformation =',&
 !                                 maxCorrectionSym*wgt
 !   write(6,'(a,1x,es11.4)') 'max skew      correction of deformation =',&
 !                                 maxCorrectionSkew*wgt
 !   write(6,'(a,1x,es11.4)') 'max sym/skew of avg correction =         ',&
 !                                 maxval(math_symmetric33(temp33_real))/&
 !                                 maxval(math_skew33(temp33_real))
 ! endif


end subroutine utilities_FFTbackward


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real, ensuring that average value = fieldAim
!--------------------------------------------------------------------------------------------------
subroutine utilities_fourierConvolution(fieldAim)
 use numerics, only: &
   memory_efficient
 use math, only: &
   math_inv33
 use mesh, only: &
   mesh_NcpElems, &
   res, &
   res1_red

 implicit none  
 real(pReal), intent(in), dimension(3,3) :: fieldAim                                                !< desired average value of the field after convolution
 real(pReal),             dimension(3,3) :: xiDyad, temp33_Real
 real(pReal)                             :: filter                                                  !< weighting of current component
 complex(pReal),          dimension(3,3) :: temp33_complex
 integer(pInt) :: &
   i, j, k, &
   l, m, n, o

 write(6,'(/,a)') '... doing convolution .....................................................'
  
!--------------------------------------------------------------------------------------------------
! to the actual spectral method calculation (mechanical equilibrium)
 if(memory_efficient) then                                                                          ! memory saving version, on-the-fly calculation of gamma_hat
   do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res1_red
       if(any([i,j,k] /= 1_pInt)) then                                                              ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
         temp33_Real = math_inv33(temp33_Real)
         filter = utilities_getFilter(xi(1:3,i,j,k))                                                ! weighting factor computed by getFilter function
         forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
           gamma_hat(l,m,n,o, 1,1,1) =  filter*temp33_Real(l,n)*xiDyad(m,o)
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, 1,1,1) * field_fourier(i,j,k,1:3,1:3))
         field_fourier(i,j,k,1:3,1:3) = temp33_Complex 
     endif             
   enddo; enddo; enddo
 else                                                                                               ! use precalculated gamma-operator
   do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
     forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt) &
       temp33_Complex(m,n) = sum(gamma_hat(m,n,1:3,1:3, i,j,k) * field_fourier(i,j,k,1:3,1:3))
     field_fourier(i,j,k, 1:3,1:3) = temp33_Complex
   enddo; enddo; enddo
 endif
 field_fourier(1,1,1,1:3,1:3) = cmplx(fieldAim*real(mesh_NcpElems,pReal),0.0_pReal,pReal)           ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1  

end subroutine utilities_fourierConvolution
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_divergenceRMS()
 use math                                                                                           !< must use the whole module for use of FFTW
 use mesh, only: &
   wgt, &
   res, &
   res1_red

 implicit none
 integer(pInt) :: i, j, k 
 real(pReal) :: &
   err_div_RMS, &                                                                                   !< RMS of divergence in Fourier space
   err_real_div_RMS, &                                                                              !< RMS of divergence in real space
   err_post_div_RMS, &                                                                              !< RMS of divergence in Fourier space, calculated using function for post processing
   err_div_max, &                                                                                   !< maximum value of divergence in Fourier space
   err_real_div_max                                                                                 !< maximum value of divergence in real space
 complex(pReal), dimension(3) ::  temp3_complex

 write(6,'(/,a)') '... calculating divergence ................................................'

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
 utilities_divergenceRMS = 0.0_pReal
 do k = 1_pInt, res(3); do j = 1_pInt, res(2)
   do i = 2_pInt, res1_red -1_pInt                                                                 ! Has somewhere a conj. complex counterpart. Therefore count it twice.
     utilities_divergenceRMS = utilities_divergenceRMS &
           + 2.0_pReal*(sum (real(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),&              ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                           xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&                   ! --> sum squared L_2 norm of vector 
                       +sum(aimag(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),& 
                                                          xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
   enddo
   utilities_divergenceRMS = utilities_divergenceRMS &                                             ! these two layers (DC and Nyquist) do not have a conjugate complex counterpart
                 + sum( real(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                     xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                 + sum(aimag(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                     xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                 + sum( real(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                     xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)&
                 + sum(aimag(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                     xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)
 enddo; enddo

 utilities_divergenceRMS = sqrt(utilities_divergenceRMS) *wgt                                      ! RMS in real space calculated with Parsevals theorem from Fourier space
 
!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
 if (debugDivergence) then                                                                         ! calculate divergence again
   err_div_max = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     temp3_Complex = math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3)*wgt,&                       ! weighting P_fourier
                                             xi(1:3,i,j,k))*TWOPIIMG
     err_div_max = max(err_div_max,sum(abs(temp3_Complex)**2.0_pReal))
     divergence_fourier(i,j,k,1:3) = temp3_Complex                                                 ! need divergence NOT squared
   enddo; enddo; enddo
   
   call fftw_execute_dft_c2r(plan_divergence,divergence_fourier,divergence_real)                   ! already weighted

   err_real_div_RMS = sqrt(wgt*sum(divergence_real**2.0_pReal))                                    ! RMS in real space
   err_post_div_RMS = sqrt(wgt*sum(divergence_post**2.0_pReal))                                    ! RMS in real space from funtion in math.f90
   err_real_div_max = sqrt(maxval(sum(divergence_real**2.0_pReal,dim=4)))                          ! max in real space                                       
   err_div_max      = sqrt(    err_div_max)                                                        ! max in Fourier space
   
   write(6,'(1x,a,es11.4)')        'error divergence  FT  RMS = ',err_div_RMS
   write(6,'(1x,a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
   write(6,'(1x,a,es11.4)')        'error divergence post RMS = ',err_post_div_RMS
   write(6,'(1x,a,es11.4)')        'error divergence  FT  max = ',err_div_max
   write(6,'(1x,a,es11.4)')        'error divergence Real max = ',err_real_div_max
 endif

end function utilities_divergenceRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculates mask compliance tensor
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
   if(debugGeneral) &
     write(6,'(a,/,9(9(2x,f12.7,1x)/))',advance='no') 'Stiffness C rotated  / GPa =',&
                                                  transpose(temp99_Real)/1.e9_pReal
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
   if(errmatinv) call IO_error(error_ID=400_pInt)
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
! check if inversion was successfull
   sTimesC = matmul(c_reduced,s_reduced)
   do m=1_pInt, size_reduced
     do n=1_pInt, size_reduced
       if(m==n .and. abs(sTimesC(m,n)) > (1.0_pReal + 10.0e-12_pReal)) errmatinv = .true.           ! diagonal elements of S*C should be 1
       if(m/=n .and. abs(sTimesC(m,n)) > (0.0_pReal + 10.0e-12_pReal)) errmatinv = .true.           ! off diagonal elements of S*C should be 0
     enddo
   enddo
   if(debugGeneral .or. errmatinv) then                                                             ! report
     write(formatString, '(I16.16)') size_reduced
     formatString = '(a,/,'//trim(formatString)//'('//trim(formatString)//'(2x,es9.2,1x)/))'
     write(6,trim(formatString),advance='no') 'C * S', transpose(matmul(c_reduced,s_reduced))
     write(6,trim(formatString),advance='no') 'S', transpose(s_reduced)
   endif
   if(errmatinv) call IO_error(error_ID=400_pInt)
   deallocate(c_reduced)
   deallocate(s_reduced)
   deallocate(sTimesC)
 else
   temp99_real = 0.0_pReal
 endif
 if(debugGeneral) &                                                                                 ! report
   write(6,'(a,/,9(9(2x,f12.7,1x)/))',advance='no') 'Masked Compliance * GPa =', &
                                                  transpose(temp99_Real*1.e9_pReal)
 utilities_maskedCompliance = math_Plain99to3333(temp99_Real)

end function utilities_maskedCompliance 


!--------------------------------------------------------------------------------------------------
!> @brief calculates constitutive response
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(F_lastInc,F,temperature,timeinc,&
                                          P,C,P_av,forwardData,rotation_BC)
 use debug, only: &
   debug_reset, &
   debug_info
 use math, only: &
   math_transpose33, &
   math_rotate_forward33
 use FEsolving, only: &
   restartWrite
 use mesh, only: &
   res, &
   wgt
 use CPFEM, only: &
   CPFEM_general
 
 implicit none
 real(pReal), intent(inout)                                      :: temperature                     !< temperature (no field)
 real(pReal), intent(in),    dimension(3,3,res(1),res(2),res(3)) :: &
   F_lastInc, &                                                                                     !< target deformation gradient
   F                                                                                                !< previous deformation gradient
 real(pReal), intent(in)                                         :: timeinc                         !< loading time
 logical,     intent(in)                                         :: forwardData                     !< age results
 real(pReal), intent(in),    dimension(3,3)                      :: rotation_BC                     !< rotation of load frame
 
 real(pReal),intent(out),    dimension(3,3,3,3)                  :: C                               !< average stiffness
 real(pReal),intent(out),    dimension(3,3)                      :: P_av                            !< average PK stress
 real(pReal),intent(out),    dimension(3,3,res(1),res(2),res(3)) :: P                               !< PK stress
 
 integer(pInt) :: &
   i, j, k, &
   ielem, &
   calcMode, &                                                                                      !< CPFEM mode for calculation
   collectMode                                                                                      !< CPFEM mode for collection
 real(pReal), dimension(3,3,3,3) :: dPdF                                                            !< d P / d F
 real(pReal), dimension(6)       :: sigma                                                           !< cauchy stress in mandel notation
 real(pReal), dimension(6,6)     :: dsde                                                            !< d sigma / d Epsilon

 write(6,'(/,a,/)') '... evaluating constitutive response ......................................'
 if (forwardData) then                                                                              ! aging results
   calcMode = 1_pInt
   collectMode = 4_pInt
  else                                                                                              ! normal calculation
   calcMode = 2_pInt
   collectMode = 3_pInt
 endif
 if (cutBack) then                                                                                  ! restore saved variables
  calcMode = 2_pInt
  collectMode = 5_pInt
 endif
!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
  ! if(debugGeneral) then
    ! defgradDetMax = -huge(1.0_pReal)
    ! defgradDetMin = +huge(1.0_pReal)
    ! do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
      ! defgradDet = math_det33(F(i,j,k,1:3,1:3))
      ! defgradDetMax = max(defgradDetMax,defgradDet)
      ! defgradDetMin = min(defgradDetMin,defgradDet) 
    ! enddo; enddo; enddo

    ! write(6,'(a,1x,es11.4)') 'max determinant of deformation =', defgradDetMax
    ! write(6,'(a,1x,es11.4)') 'min determinant of deformation =', defgradDetMin
  ! endif
 if (DebugGeneral) write(6,*) 'collect mode: ', collectMode,' calc mode: ', calcMode
 flush(6)
 
 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt
   call CPFEM_general(collectMode,&                                                                 ! collect cycle
                       F_lastInc(1:3,1:3,i,j,k),F(1:3,1:3,i,j,k), &
                       temperature,timeinc,ielem,1_pInt,sigma,dsde,P(1:3,1:3,i,j,k),dPdF)
   collectMode = 3_pInt
 enddo; enddo; enddo

 P = 0.0_pReal                                                                                      ! needed because of the padding for FFTW
 C = 0.0_pReal
 ielem = 0_pInt 
 call debug_reset()
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt
   call CPFEM_general(calcMode,&                                                                    ! first element in first iteration retains CPFEM_mode 1, 
                      F_lastInc(1:3,1:3,i,j,k), F(1:3,1:3,i,j,k), &          ! others get 2 (saves winding forward effort)
                      temperature,timeinc,ielem,1_pInt,sigma,dsde,P(1:3,1:3,i,j,k),dPdF)
   calcMode = 2_pInt
   C = C + dPdF
 enddo; enddo; enddo
 C = C * wgt
 call debug_info()
 
 restartWrite = .false.                                                                             ! reset restartWrite status
 cutBack = .false.                                                                                  ! reset cutBack status
 
 P_av = sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt                                                    ! average of P 
 if (debugRotation) &
 write(6,'(a,/,3(3(2x,f12.7,1x)/))',advance='no') 'Piola-Kirchhoff stress (lab) / MPa =',&
                                                     math_transpose33(P_av)/1.e6_pReal
 P_av = math_rotate_forward33(P_av,rotation_BC)
 write(6,'(a,/,3(3(2x,f12.7,1x)/))',advance='no') 'Piola-Kirchhoff stress / MPa =',&
                                                     math_transpose33(P_av)/1.e6_pReal
end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief calculates forward rate, either guessing or just add delta/timeinc
!--------------------------------------------------------------------------------------------------
pure function utilities_calculateRate(avRate,timeinc_old,guess,field_lastInc,field)
 use mesh, only: &
   res
 
 implicit none
 real(pReal), intent(in), dimension(3,3)                      :: avRate                             !< homogeneous addon
 real(pReal), intent(in) :: &
   timeinc_old                                                                                      !< timeinc of last step
 logical, intent(in) :: &
   guess                                                                                            !< guess along former trajectory
 real(pReal), intent(in), dimension(3,3,res(1),res(2),res(3)) :: &
   field_lastInc, &                                                                                 !< data of previous step
   field                                                                                            !< data of current step
 real(pReal),             dimension(3,3,res(1),res(2),res(3)) :: utilities_calculateRate
 
 if(guess) then
   utilities_calculateRate = (field-field_lastInc) / timeinc_old
 else
   utilities_calculateRate = spread(spread(spread(avRate,3,res(1)),4,res(2)),5,res(3))
 endif

end function utilities_calculateRate


!--------------------------------------------------------------------------------------------------
!> @brief forwards a field with a pointwise given rate, ensures that the average matches the aim
!--------------------------------------------------------------------------------------------------
pure function utilities_forwardField(timeinc,aim,field_lastInc,rate)
 use mesh, only: &
   res, &
   wgt

 implicit none
 real(pReal), intent(in)                                      :: timeinc                            !< timeinc of current step
 real(pReal), intent(in), dimension(3,3)                      :: aim                                !< average field value aim
 real(pReal), intent(in), dimension(3,3,res(1),res(2),res(3)) :: &
   field_lastInc,&                                                                                  !< initial field
   rate                                                                                             !< rate by which to forward
 real(pReal),             dimension(3,3,res(1),res(2),res(3)) :: utilities_forwardField
 real(pReal),             dimension(3,3)                      :: fieldDiff                          !< <a + adot*t> - aim
 
 utilities_forwardField = field_lastInc + rate*timeinc
 fieldDiff = sum(sum(sum(utilities_forwardField,dim=5),dim=4),dim=3)*wgt - aim
 utilities_forwardField = utilities_forwardField - &
                          spread(spread(spread(fieldDiff,3,res(1)),4,res(2)),5,res(3))

end function utilities_forwardField


!--------------------------------------------------------------------------------------------------
!> @brief calculates filter for fourier convolution depending on type given in numerics.config
!--------------------------------------------------------------------------------------------------
real(pReal) function utilities_getFilter(k)
 use IO, only: &
   IO_error
 use numerics, only: &                        
   myfilter
 use mesh, only: &
   res
 use math, only: &
   PI
  
 implicit none
 real(pReal),intent(in), dimension(3) :: k                                                          !< indices of frequency
  
 select case (myfilter)
    case ('none')
      utilities_getFilter = 1.0_pReal
    case ('cosine')                                                                                 !< cosine curve with 1 for avg and zero for highest freq
      utilities_getFilter = (1.0_pReal + cos(PI*k(3)/res(3))) &
                           *(1.0_pReal + cos(PI*k(2)/res(2))) &
                           *(1.0_pReal + cos(PI*k(1)/res(1)))/8.0_pReal
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
 if (debugDivergence) call fftw_destroy_plan(plan_divergence)

 if (debugFFTW) call fftw_destroy_plan(plan_scalarField_forth)
 if (debugFFTW) call fftw_destroy_plan(plan_scalarField_back)
  
 call fftw_destroy_plan(plan_forward)
 call fftw_destroy_plan(plan_backward)

end subroutine utilities_destroy


end module DAMASK_spectral_utilities
