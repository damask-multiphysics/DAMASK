!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the different spectral solver variants
!--------------------------------------------------------------------------------------------------
module DAMASK_spectral_Utilities
 
 use prec, only: &
   pReal, &
   pInt
 use mesh,  only : &
   res, &
   res1_red, &
   geomdim, &
   mesh_NcpElems, &
   wgt
   
 use math
 
 use IO, only: &
   IO_error
 
 implicit none

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 type(C_PTR),   private                                        :: plan_forward, plan_backward                                                        ! plans for fftw
 real(pReal),   private, dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat                                   ! gamma operator (field) for spectral method
 real(pReal),   private, dimension(:,:,:,:),       allocatable :: xi                                          ! wave vector field for divergence and for gamma operator
 complex(pReal),private, dimension(:,:,:,:,:),     pointer     :: field_fourier
 real(pReal),   private, dimension(3,3,3,3)                    :: C_ref
 
 real(pReal),   public,  dimension(:,:,:,:,:),     pointer     :: field_real

!--------------------------------------------------------------------------------------------------
! debug fftw 
 type(C_PTR),   private  :: plan_scalarField_forth, plan_scalarField_back
 complex(pReal),private,  dimension(:,:,:), pointer :: scalarField_real
 complex(pReal),private,  dimension(:,:,:), pointer :: scalarField_fourier
 
!--------------------------------------------------------------------------------------------------
! debug divergence
 type(C_PTR), private :: plan_divergence
 real(pReal),    private, dimension(:,:,:,:), pointer :: divergence_real
 complex(pReal), private, dimension(:,:,:,:), pointer :: divergence_fourier
 real(pReal),    dimension(:,:,:,:), allocatable :: divergence_post

!--------------------------------------------------------------------------------------------------
!variables controlling debugging
 logical,public :: debugGeneral, debugDivergence, debugRestart, debugFFTW

!--------------------------------------------------------------------------------------------------
! derived types
 type solutionState 
   logical ::        converged       = .false.   
   logical ::        regrid          = .false.   
   logical ::        term_ill        = .false.   
 end type solutionState

 type boundaryCondition
   real(pReal), dimension(3,3) :: values = 0.0_pReal
   real(pReal), dimension(3,3) :: maskFloat = 0.0_pReal
   logical,     dimension(3,3) :: maskLogical = .false.
   character(len=64)           :: myType = 'None'
 end type boundaryCondition

contains 

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, sets debug flags, create plans for fftw
!> @details Sets the debug levels for general, divergence, restart and fftw from the biwise coding 
!> provided by the debug module to logicals.
!> Allocates all fields used by FFTW and create the corresponding plans depending on the debug
!> level chosen.
!> Initializes FFTW.
!--------------------------------------------------------------------------------------------------
subroutine Utilities_init()

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
   debug_spectralFFTW
 
 use mesh, only : &
   virt_dim
 
 implicit none
 integer(pInt) :: i, j, k
 integer(pInt), dimension(3) :: k_s
 !$ integer(pInt) :: ierr
 type(C_PTR) :: tensorField                                                                         ! field in real and fourier space
 type(C_PTR) :: scalarField_realC, scalarField_fourierC
 type(C_PTR) :: divergence
 

 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectralSolver Utilities init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ''

!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_spectral),debug_levelBasic)         /= 0
 debugDivergence = iand(debug_level(debug_spectral),debug_spectralDivergence) /= 0
 debugRestart    = iand(debug_level(debug_spectral),debug_spectralRestart)    /= 0
 debugFFTW       = iand(debug_level(debug_spectral),debug_spectralFFTW)       /= 0
 
!--------------------------------------------------------------------------------------------------
! allocation
 allocate (xi         (3,res1_red,res(2),res(3)),      source = 0.0_pReal)                           ! start out isothermally
 tensorField = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))                       ! allocate continous data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, field_real,      [ res(1)+2_pInt,res(2),res(3),3,3])                  ! place a pointer for a real representation on tensorField
 call c_f_pointer(tensorField, field_fourier,   [ res1_red,     res(2),res(3),3,3])                  ! place a pointer for a complex representation on tensorField

!--------------------------------------------------------------------------------------------------
! general initialization of fftw (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=808_pInt)                         ! check for correct precision in C
!$ if(DAMASK_NumThreadsInt > 0_pInt) then
!$   ierr = fftw_init_threads()
!$   if (ierr == 0_pInt) call IO_error(error_ID = 809_pInt)
!$   call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
!$ endif
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans
 plan_forward =    fftw_plan_many_dft_r2c(3,[ res(3),res(2) ,res(1)],9,&                             ! dimensions , length in each dimension in reversed order
                                 field_real,[ res(3),res(2) ,res(1)+2_pInt],&                        ! input data , physical length in each dimension in reversed order
                                         1,  res(3)*res(2)*(res(1)+2_pInt),&                        ! striding   , product of physical lenght in the 3 dimensions
                                field_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 plan_backward   =fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],9,&
                             field_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,&
                                field_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                         1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)

!--------------------------------------------------------------------------------------------------
! depending on (debug) options, allocate more memory and create additional plans 
 if (debugDivergence) then
   divergence = fftw_alloc_complex(int(res1_red*res(2)*res(3)*3_pInt,C_SIZE_T))
   call c_f_pointer(divergence, divergence_real,    [ res(1)+2_pInt,res(2),res(3),3])
   call c_f_pointer(divergence, divergence_fourier, [ res1_red,     res(2),res(3),3])
   allocate (divergence_post(res(1),res(2),res(3),3));  divergence_post = 0.0_pReal
   plan_divergence = fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],3,&
                           divergence_fourier,[ res(3),res(2) ,res1_red],&
                                            1,  res(3)*res(2)* res1_red,&
                              divergence_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                            1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 endif

 if (debugFFTW) then
   scalarField_realC    = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))                    ! do not do an inplace transform  
   scalarField_fourierC = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))
   call c_f_pointer(scalarField_realC,    scalarField_real,    [res(1),res(2),res(3)])
   call c_f_pointer(scalarField_fourierC, scalarField_fourier, [res(1),res(2),res(3)])
   plan_scalarField_forth = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_real,scalarField_fourier,-1,fftw_planner_flag)
   plan_scalarField_back  = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_fourier,scalarField_real,+1,fftw_planner_flag)
 endif 

 if (debugGeneral) write(6,'(a)') 'FFTW initialized'
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 do k = 1_pInt, res(3)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1_pInt, res1_red
           k_s(1) = i - 1_pInt
           xi(1:3,i,j,k) = real(k_s, pReal)/virt_dim
 enddo; enddo; enddo
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(3,3,3,3,1,1,1), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(3,3,3,3,res1_red ,res(2),res(3)), source =0.0_pReal)                                                  ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif

end subroutine Utilities_init

!--------------------------------------------------------------------------------------------------
!> @brief updates references stiffness and potentially precalculated gamma operator
!> @details Sets the current reference stiffness to the stiffness given as an argument.
!> If the gamma operator is precalculated, it is calculated with this stiffness.
!> In case of a on-the-fly calculation, only the reference stiffness is updated.
!> The gamma operator is filtered depening on the filter selected in numerics
!--------------------------------------------------------------------------------------------------
subroutine Utilities_updateGamma(C)
   
 use numerics, only: &
   memory_efficient
  
 implicit none

 real(pReal), dimension(3,3,3,3), intent(in) :: C
 real(pReal), dimension(3,3) :: temp33_Real, xiDyad
 real(pReal) :: filter
 integer(pInt) :: i, j, k,   l, m, n, o
  
 C_ref = C
 if(.not. memory_efficient) then                                                   
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       filter = Utilities_getFilter(xi(1:3,i,j,k))
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
         gamma_hat(l,m,n,o, i,j,k) =  filter*temp33_Real(l,n)*xiDyad(m,o)
     endif  
   enddo; enddo; enddo
   gamma_hat(1:3,1:3,1:3,1:3, 1,1,1) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif
end subroutine Utilities_updateGamma

!--------------------------------------------------------------------------------------------------
!> @brief forward FFT of data in field_real to field_fourier with highest freqs. removed
!> Does an unweighted FFT transform from real to complex.
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!--------------------------------------------------------------------------------------------------
subroutine Utilities_forwardFFT(row,column)
 use mesh, only : &
   virt_dim
  
  implicit none
  integer(pInt), intent(in), optional :: row, column
  
!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
  if (debugFFTW) then
    if (.not. present(row) .or. .not. present(column)) stop
    scalarField_real(1:res(1),1:res(2),1:res(3)) =&                                          ! store the selected component
           cmplx(field_real(1:res(1),1:res(2),1:res(3),row,column),0.0_pReal,pReal)
  endif
  
!--------------------------------------------------------------------------------------------------
! call function to calculate divergence from math (for post processing) to check results
  if (debugDivergence) &
    call divergence_fft(res,virt_dim,3_pInt,field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),divergence_post)
  
!--------------------------------------------------------------------------------------------------
! doing the FT
  call fftw_execute_dft_r2c(plan_forward,field_real,field_fourier)
  
!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
  if (debugFFTW) then
    call fftw_execute_dft(plan_scalarField_forth,scalarField_real,scalarField_fourier)
    write(6,'(a,i1,1x,i1)') 'checking FT results of compontent ', row, column
    write(6,'(a,2(es11.4,1x))')  'max FT relative error = ',&
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
 if(res(3)>1_pInt) &
  field_fourier (1:res1_red,1:res(2),                res(3)/2_pInt+1_pInt,1:3,1:3)&
                                                     = cmplx(0.0_pReal,0.0_pReal,pReal)
end subroutine Utilities_forwardFFT


!--------------------------------------------------------------------------------------------------
!> @brief backward FFT of data in field_fourier to field_real
!> Does an inverse FFT transform from complex to real
!> In case of debugging the FFT, also one component of the tensor (specified by row and column)
!> is independetly transformed complex to complex and compared to the whole tensor transform
!> results is weighted by number of points stored in wgt
!--------------------------------------------------------------------------------------------------
subroutine Utilities_backwardFFT(row,column)

  implicit none
  
  integer(pInt), intent(in), optional :: row, column
  integer(pInt) ::  i, j, k, m, n
  
!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
  if (debugFFTW) then
    scalarField_fourier = field_fourier(1:res1_red,1:res(2),1:res(3),row,column)
    do i = 0_pInt, res(1)/2_pInt-2_pInt                                                      ! unpack fft data for conj complex symmetric part
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


  call fftw_execute_dft_c2r(plan_backward,field_fourier,field_real)                      ! back transform of fluct deformation gradient

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
  field_real = field_real * wgt

end subroutine Utilities_backwardFFT


!--------------------------------------------------------------------------------------------------
!> @brief doing convolution gamma_hat * field_real with average value given by fieldAim
!--------------------------------------------------------------------------------------------------
subroutine Utilities_fourierConvolution(fieldAim)

 use numerics, only: &
   memory_efficient
    
 implicit none  
  
 real(pReal), dimension(3,3), intent(in) :: fieldAim
 real(pReal), dimension(3,3) :: xiDyad, temp33_Real
 real(pReal) :: filter
 integer(pInt) :: i, j, k,    l, m, n, o
 complex(pReal), dimension(3,3) ::  temp33_complex

 write(6,'(a)') ''
 write(6,'(a)') '... doing convolution .................'
 write(6,'(a)') ''
  
!--------------------------------------------------------------------------------------------------
! to the actual spectral method calculation (mechanical equilibrium)
 if(memory_efficient) then                                                                  ! memory saving version, on-the-fly calculation of gamma_hat
   do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res1_red
       if(any([i,j,k] /= 1_pInt)) then                                                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
         temp33_Real = math_inv33(temp33_Real)
         filter = Utilities_getFilter(xi(1:3,i,j,k))
         forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, o=1_pInt:3_pInt)&
           gamma_hat(l,m,n,o, 1,1,1) =  filter*temp33_Real(l,n)*xiDyad(m,o)
         forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
           temp33_Complex(l,m) = sum(gamma_hat(l,m,1:3,1:3, 1,1,1) * field_fourier(i,j,k,1:3,1:3))
         field_fourier(i,j,k,1:3,1:3) = temp33_Complex 
     endif             
   enddo; enddo; enddo
 else                                                                                       ! use precalculated gamma-operator
   do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
     forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt) &
       temp33_Complex(m,n) = sum(gamma_hat(m,n,1:3,1:3, i,j,k) * field_fourier(i,j,k,1:3,1:3))
     field_fourier(i,j,k, 1:3,1:3) = temp33_Complex
   enddo; enddo; enddo
 endif

 field_fourier(1,1,1,1:3,1:3) = cmplx(fieldAim*real(mesh_NcpElems,pReal),0.0_pReal,pReal)             ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1  

end subroutine Utilities_fourierConvolution
 

!--------------------------------------------------------------------------------------------------
!> @brief calculate root mean square of divergence of field_fourier
!--------------------------------------------------------------------------------------------------
real(pReal) function Utilities_divergenceRMS()
 
  integer(pInt) :: i, j, k 
  real(pReal) :: err_div_RMS, err_real_div_RMS, err_post_div_RMS,&
                 err_div_max, err_real_div_max
  complex(pReal), dimension(3) ::  temp3_complex

  write(6,'(a)') ''
  write(6,'(a)') '... calculating divergence .................'

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
  Utilities_divergenceRMS = 0.0_pReal
  do k = 1_pInt, res(3); do j = 1_pInt, res(2)
    do i = 2_pInt, res1_red -1_pInt                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
      Utilities_divergenceRMS = Utilities_divergenceRMS &
            + 2.0_pReal*(sum (real(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),&           ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                            xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&            ! --> sum squared L_2 norm of vector 
                        +sum(aimag(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),& 
                                                           xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
    enddo
    Utilities_divergenceRMS = Utilities_divergenceRMS &                                                              ! Those two layers (DC and Nyquist) do not have a conjugate complex counterpart
                  + sum( real(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                      xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                  + sum(aimag(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                      xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                  + sum( real(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                      xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)&
                  + sum(aimag(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                      xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)
  enddo; enddo

  Utilities_divergenceRMS = sqrt(Utilities_divergenceRMS) *wgt                                                        ! RMS in real space calculated with Parsevals theorem from Fourier space
  
!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
  if (debugDivergence) then                                                                  ! calculate divergence again
    err_div_max = 0.0_pReal
    do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
      temp3_Complex = math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3)*wgt,&                    ! weighting P_fourier
                                              xi(1:3,i,j,k))*TWOPIIMG
      err_div_max = max(err_div_max,sum(abs(temp3_Complex)**2.0_pReal))
      divergence_fourier(i,j,k,1:3) = temp3_Complex                                          ! need divergence NOT squared
    enddo; enddo; enddo
    
    call fftw_execute_dft_c2r(plan_divergence,divergence_fourier,divergence_real)            ! already weighted

    err_real_div_RMS =  sqrt(wgt*sum(divergence_real**2.0_pReal))                            ! RMS in real space
    err_post_div_RMS =  sqrt(wgt*sum(divergence_post**2.0_pReal))                            ! RMS in real space
    err_real_div_max =  sqrt(maxval(sum(divergence_real**2.0_pReal,dim=4)))                  ! max in real space                                       
    err_div_max      = sqrt(    err_div_max)                                                 ! max in Fourier space
    
    write(6,'(1x,a,es11.4)')        'error divergence  FT  RMS = ',err_div_RMS
    write(6,'(1x,a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
    write(6,'(1x,a,es11.4)')        'error divergence post RMS = ',err_post_div_RMS
    write(6,'(1x,a,es11.4)')        'error divergence  FT  max = ',err_div_max
    write(6,'(1x,a,es11.4)')        'error divergence Real max = ',err_real_div_max
  endif

end function Utilities_divergenceRMS


!--------------------------------------------------------------------------------------------------
!> @brief calculates mask compliance
!--------------------------------------------------------------------------------------------------
function Utilities_maskedCompliance(rot_BC,mask_stressVector,C)

 real(pReal), dimension(3,3,3,3) :: Utilities_maskedCompliance
 real(pReal), dimension(3,3,3,3), intent(in) :: C
 integer(pInt) :: i, j, k, m, n 
 real(pReal), dimension(3,3), intent(in) :: rot_BC
 logical, dimension(9), intent(in) :: mask_stressVector
 real(pReal), dimension(3,3,3,3) :: C_lastInc
 real(pReal), dimension(9,9) :: temp99_Real   
 integer(pInt) :: size_reduced = 0_pInt 
 real(pReal), dimension(:,:), allocatable ::  s_reduced, c_reduced                                  ! reduced compliance and stiffness (only for stress BC)
 logical :: errmatinv
 
 size_reduced = count(mask_stressVector)
 if(size_reduced > 0_pInt )then
   allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)

   C_lastInc = math_rotate_forward3333(C,rot_BC)                                                    ! calculate stiffness from former inc
   temp99_Real = math_Plain3333to99(C_lastInc)
   k = 0_pInt                                                                                       ! build reduced stiffness
   do n = 1_pInt,9_pInt
     if(mask_stressVector(n)) then
       k = k + 1_pInt
       j = 0_pInt
       do m = 1_pInt,9_pInt
         if(mask_stressVector(m)) then
           j = j + 1_pInt
           c_reduced(k,j) = temp99_Real(n,m)
   endif; enddo; endif; enddo
   call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)                               ! invert reduced stiffness
   if(errmatinv) call IO_error(error_ID=400_pInt)
   temp99_Real = 0.0_pReal                                                                          ! build full compliance
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
   deallocate(c_reduced)
   deallocate(s_reduced)
 else
   temp99_real = 0.0_pReal
 endif

 Utilities_maskedCompliance = math_Plain99to3333(temp99_Real)

end function Utilities_maskedCompliance 

subroutine Utilities_constitutiveResponse(coordinates,F_lastInc,F,temperature,timeinc,&
                                P,C,P_av,ForwardData,rotation_BC)
  use debug, only: &
    debug_reset, &
    debug_info 
  use CPFEM, only: &
    CPFEM_general
  use FEsolving, only: restartWrite
  
  implicit none
  
  real(pReal), dimension(res(1),res(2),res(3)) :: temperature
  real(pReal), dimension(res(1),res(2),res(3),3) :: coordinates
  
  real(pReal), dimension(3,3,res(1),res(2),res(3)) :: F,F_lastInc, P
  real(pReal) :: timeinc
  logical :: ForwardData
  integer(pInt) :: i, j, k, ielem
  integer(pInt) :: CPFEM_mode
  real(pReal), dimension(3,3,3,3) :: dPdF, C
  real(pReal), dimension(6)       :: sigma                                                                ! cauchy stress
  real(pReal), dimension(6,6)     :: dsde
  real(pReal), dimension(3,3)     :: P_av, rotation_BC
  
  write(6,'(a)') ''
  write(6,'(a)') '... evaluating constitutive response .................'
  write(6,'(a)') ''
  
  if (ForwardData) then
    CPFEM_mode = 1_pInt
   else
    CPFEM_mode = 2_pInt
  endif

  ielem = 0_pInt
  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    ielem = ielem + 1_pInt
    call CPFEM_general(3_pInt,&                                                              ! collect cycle
                        coordinates(i,j,k,1:3), F_lastInc(1:3,1:3,i,j,k),F(1:3,1:3,i,j,k), &
                        temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde,P(1:3,1:3,i,j,k),dPdF)
  enddo; enddo; enddo

  P = 0.0_pReal                                                                         ! needed because of the padding for FFTW
  C = 0.0_pReal
  ielem = 0_pInt 
  call debug_reset()
  do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
    ielem = ielem + 1_pInt
    call CPFEM_general(CPFEM_mode,&                                                          ! first element in first iteration retains CPFEM_mode 1, 
                       coordinates(i,j,k,1:3),F_lastInc(1:3,1:3,i,j,k), F(1:3,1:3,i,j,k), &  ! others get 2 (saves winding forward effort)
                       temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde,P(1:3,1:3,i,j,k),dPdF)
    CPFEM_mode = 2_pInt
    C = C + dPdF
  enddo; enddo; enddo
  call debug_info()
  
  P_av = math_rotate_forward33(sum(sum(sum(P,dim=5),dim=4),dim=3) * wgt,rotation_BC)               !average of P rotated
  restartWrite = .false.

  
  write (6,'(a,/,3(3(2x,f12.7,1x)/))',advance='no') ' Piola-Kirchhoff stress / MPa =',&
                                                      math_transpose33(P_av)/1.e6_pReal
  
  C = C * wgt
end subroutine Utilities_constitutiveResponse


subroutine Utilities_forwardField(delta_aim,timeinc,timeinc_old,guessmode,field_lastInc,field)
 
 real(pReal), intent(in), dimension(3,3) :: delta_aim

 real(pReal), intent(in) :: timeinc, timeinc_old, guessmode
 real(pReal), intent(inout), dimension(3,3,res(1),res(2),res(3)) :: field_lastInc,field
 
 if (guessmode == 1.0_pReal) then
   field = field + (field-field_lastInc) * timeinc/timeinc_old
   field_lastInc = (field + field_lastInc * timeinc/timeinc_old) /(1.0_pReal + timeinc/timeinc_old)
 else
   field_lastInc = field
   field = field + spread(spread(spread(delta_aim,3,res(1)),4,res(2)),5,res(3))
 endif

end subroutine Utilities_forwardField
 
real(pReal) function Utilities_getFilter(k)

  use numerics, only: &                        
   myfilter
  
  implicit none
  
  real(pReal), dimension(3),intent(in) :: k
  
  select case (myfilter)
       
    case ('none')
      Utilities_getFilter = 1.0_pReal
      
    case ('cosine')
      Utilities_getFilter = 0.125_pReal*(1.0_pReal + cos(pi*k(3)*geomdim(3)/(res(3)/2_pInt + 1_pInt))) &
                           *(1.0_pReal + cos(pi*k(2)*geomdim(2)/(res(2)/2_pInt + 1_pInt))) &
                           *(1.0_pReal + cos(pi*k(1)*geomdim(1)/(res(1)/2_pInt + 1_pInt)))
    
  end select 

end function Utilities_getFilter

subroutine Utilities_destroy()

  implicit none
  
  if (debugDivergence) call fftw_destroy_plan(plan_divergence)
  
  if (debugFFTW) then
    call fftw_destroy_plan(plan_scalarField_forth)
    call fftw_destroy_plan(plan_scalarField_back)
  endif 
  
  call fftw_destroy_plan(plan_forward)
  call fftw_destroy_plan(plan_backward)

end subroutine Utilities_destroy
       
end module DAMASK_spectral_Utilities
