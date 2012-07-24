! Copyright 2012 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
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
!##################################################################################################
!* $Id$
!##################################################################################################
! Material subroutine for BVP solution using spectral method
!
! Run 'DAMASK_spectral.exe --help' to get usage hints
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts,
!            D.D. Tjahjanto,
!            C. Kords,
!            M. Diehl,
!            R. Lebensohn
!
! MPI fuer Eisenforschung, Duesseldorf

module DAMASK_spectral_Utilities
 
 use prec, only: &
   pReal, &
   pInt

 use math
 
 use IO, only: &
   IO_error
 
 implicit none
 
 type solutionState                                 ! mask of stress boundary conditions
   logical ::        converged       = .false.   
   logical ::        regrid          = .false.   
   logical ::        term_ill        = .false.   
 end type solutionState

 character(len=5) :: solverType, parameter = 'basic'
!--------------------------------------------------------------------------------------------------
! common pointwise data
 real(pReal),    dimension(:,:,:,:,:), allocatable ::  F, F_lastInc, P
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature


!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 type(C_PTR) ::  plan_forward, plan_backward                                                        ! plans for fftw
 real(pReal),    dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                   ! gamma operator (field) for spectral method
 real(pReal),    dimension(:,:,:,:),       allocatable ::  xi                                          ! wave vector field for divergence and for gamma operator
 real(pReal),    dimension(:,:,:,:,:),     pointer :: field_real
 complex(pReal), dimension(:,:,:,:,:),     pointer :: field_fourier

!--------------------------------------------------------------------------------------------------
! debug fftw 
 type(C_PTR) :: plan_scalarField_forth, plan_scalarField_back
 complex(pReal), dimension(:,:,:), pointer :: scalarField_real
 complex(pReal), dimension(:,:,:), pointer :: scalarField_fourier
 
!--------------------------------------------------------------------------------------------------
! debug divergence
 type(C_PTR) :: plan_divergence
 real(pReal),    dimension(:,:,:,:), pointer :: divergence_real
 complex(pReal), dimension(:,:,:,:), pointer :: divergence_fourier
 real(pReal),    dimension(:,:,:,:), allocatable :: divergence_post

!--------------------------------------------------------------------------------------------------
!variables controlling debugging
 logical :: debugGeneral, debugDivergence, debugRestart, debugFFTW

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) :: &
   F_aim = math_I3, &
   F_aim_lastInc = math_I3, &
   P_av

 real(pReal), dimension(3,3,3,3) :: &
   C_ref = 0.0_pReal, &
   C = 0.0_pReal
   
 real(pReal), dimension(3) ::   geomdim = 0.0_pReal, virt_dim = 0.0_pReal             ! physical dimension of volume element per direction
 integer(pInt), dimension(3) :: res = 1_pInt 
 real(pReal) ::  wgt
 integer(pInt) :: res1_red, Npoints

contains 

subroutine Utilities_init(F,P,F_...)

 use DAMASK_interface, only: &
   getSolverJobName
   
 use mesh,  only : &
   mesh_spectral_getResolution, &
   mesh_spectral_getDimension
 
 use numerics, only: &
   divergence_correction, &                             
   DAMASK_NumThreadsInt, &
   fftw_planner_flag, &
   fftw_timelimit
   
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_levelBasic, &
   debug_spectralDivergence, &
   debug_spectralRestart, &
   debug_spectralFFTW

 use FEsolving, only: &
   restartInc
 use numerics, only: &
   memory_efficient
 
 use CPFEM, only: &
   CPFEM_general
 
 use IO, only: &
   IO_read_JobBinaryFile, &
   IO_write_JobBinaryFile
 
 implicit none

 real(pReal), dimension(3,3) :: temp33_Real, xiDyad
 integer(pInt) :: i, j, k, l, m, n, q, ierr
 integer(pInt), dimension(3) :: k_s
  
 type(C_PTR) :: tensorField                                                                         ! field in real and fourier space
 type(C_PTR) :: scalarField_realC, scalarField_fourierC
 type(C_PTR) :: divergence
 

 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectralSolver init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ''

!--------------------------------------------------------------------------------------------------
! set debugging parameters
 debugGeneral    = iand(debug_level(debug_spectral),debug_levelBasic)         /= 0
 debugDivergence = iand(debug_level(debug_spectral),debug_spectralDivergence) /= 0
 debugRestart    = iand(debug_level(debug_spectral),debug_spectralRestart)    /= 0
 debugFFTW       = iand(debug_level(debug_spectral),debug_spectralFFTW)       /= 0
 
!##################################################################################################
! initialization 
!##################################################################################################
 res     =   mesh_spectral_getResolution()
 geomdim = mesh_spectral_getDimension()
 res1_red = res(1)/2_pInt + 1_pInt
 Npoints = res(1)*res(2)*res(3)
 wgt = 1.0/real(Npoints,pReal)

 allocate (F          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (F_lastInc  (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (P          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (xi         (3,res1_red,res(2),res(3)),      source = 0.0_pReal)
 allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
 allocate (temperature(  res(1),  res(2),res(3)),      source = 0.0_pReal)                           ! start out isothermally
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
! init fields                 
 if (restartInc == 1_pInt) then                                                                     ! no deformation (no restart)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     F(i,j,k,1:3,1:3) = math_I3
     F_lastInc(i,j,k,1:3,1:3) = math_I3
     coordinates(i,j,k,1:3) = geomdim/real(res,pReal)*real([i,j,k],pReal) &
                            - geomdim/real(2_pInt*res,pReal)
   enddo; enddo; enddo
 elseif (restartInc > 1_pInt) then                                                                  ! using old values from file                                                      
   if (debugRestart) write(6,'(a,i6,a)') 'Reading values of increment ',&
                                             restartInc - 1_pInt,' from file' 
   call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                trim(getSolverJobName()),size(F))
   read (777,rec=1) F
   close (777)
   call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad_lastInc',&
                                                trim(getSolverJobName()),size(F_lastInc))
   read (777,rec=1) F_lastInc
   close (777)
   call IO_read_jobBinaryFile(777,'F_aim',trim(getSolverJobName()),size(F_aim))
   read (777,rec=1) F_aim
   close (777)
   call IO_read_jobBinaryFile(777,'F_aim_lastInc',trim(getSolverJobName()),size(F_aim_lastInc))
   read (777,rec=1) F_aim_lastInc
   close (777)

   coordinates = 0.0 ! change it later!!!
 endif
 call constitutiveResponse(.FALSE.,0.0_pReal)
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
 if (divergence_correction) then
   do i = 1_pInt, 3_pInt
    if (i /= minloc(geomdim,1) .and. i /= maxloc(geomdim,1)) virt_dim = geomdim/geomdim(i)
   enddo
 else
   virt_dim = geomdim
 endif

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

!--------------------------------------------------------------------------------------------------
! calculate the gamma operator
 if (restartInc == 1_pInt) then
   C_ref = C 
   call IO_write_jobBinaryFile(777,'C_ref',size(C_ref))
   write (777,rec=1) C_ref
   close(777)
 elseif (restartInc > 1_pInt) then
   call IO_read_jobBinaryFile(777,'C_ref',trim(getSolverJobName()),size(C_ref))
   read (777,rec=1) C_ref
   close (777)
 endif
 
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(res1_red ,res(2),res(3),3,3,3,3), source =0.0_pReal)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
       forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
         temp33_Real(l,m) = sum(C_ref(l,m,1:3,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, q=1_pInt:3_pInt)&
         gamma_hat(i,j,k, l,m,n,q) =  temp33_Real(l,n)*xiDyad(m,q)
     endif  
   enddo; enddo; enddo
   gamma_hat(1,1,1, 1:3,1:3,1:3,1:3) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif
end subroutine Utilities_init
 
real(pReal) function convolution(calcDivergence, field_aim,)
 use numerics, only: &
   memory_efficient, &
   err_div_tol
 real(pReal),    dimension(3,3) ::                         xiDyad                                      ! product of wave vectors
 real(pReal) :: err_div = 0.0_pReal
 real(pReal),    dimension(3,3) ::          temp33_Real
 integer(pInt) :: i, j, k, l, m, n, q

!--------------------------------------------------------------------------------------------------
!variables for additional output due to general debugging
 real(pReal) :: maxCorrectionSym, maxCorrectionSkew
   logical :: calcDivergence
   real(pReal), dimension(3,3) :: field_avg, field_aim
  integer(pInt) :: row, column
    real(pReal) :: field_av_L2, err_div_RMS, err_real_div_RMS, err_post_div_RMS,&
                err_div_max, err_real_div_max
                     complex(pReal), dimension(3) ::  temp3_complex
     complex(pReal), dimension(3,3) ::  temp33_complex
  
  !--------------------------------------------------------------------------------------------------
! actual spectral method         
         write(6,'(a)') ''
         write(6,'(a)') '... doing convolution .................'
!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
         if (debugFFTW) then
           row =  3 !  (mod(totalIncsCounter+iter-2_pInt,9_pInt))/3_pInt + 1_pInt                      ! go through the elements of the tensors, controlled by totalIncsCounter and iter, starting at 1
           column = 3 !(mod(totalIncsCounter+iter-2_pInt,3_pInt))        + 1_pInt
           scalarField_real(1:res(1),1:res(2),1:res(3)) =&                                          ! store the selected component
                  cmplx(field_real(1:res(1),1:res(2),1:res(3),row,column),0.0_pReal,pReal)
         endif

!--------------------------------------------------------------------------------------------------
! call function to calculate divergence from math (for post processing) to check results
         if (debugDivergence) &
           call divergence_fft(res,virt_dim,3_pInt,&
                               field_real(1:res(1),1:res(2),1:res(3),1:3,1:3),divergence_post)      ! padding
              
!--------------------------------------------------------------------------------------------------
! doing the FT because it simplifies calculation of average stress in real space also
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


!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
         if( calcDivergence) then
         field_avg = real(field_fourier(1,1,1,1:3,1:3),pReal)*wgt

         field_av_L2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(field_avg,&                    ! L_2 norm of average stress (http://mathworld.wolfram.com/SpectralNorm.html)
                                                     math_transpose33(field_avg)))))
         err_div_RMS = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2)
           do i = 2_pInt, res1_red -1_pInt                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
             err_div_RMS = err_div_RMS &
                   + 2.0_pReal*(sum (real(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),&           ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                   xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&            ! --> sum squared L_2 norm of vector 
                               +sum(aimag(math_mul33x3_complex(field_fourier(i,j,k,1:3,1:3),& 
                                                                  xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
           enddo
           err_div_RMS = err_div_RMS &                                                              ! Those two layers (DC and Nyquist) do not have a conjugate complex counterpart
                         + sum( real(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(field_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum( real(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(field_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)
         enddo; enddo

         err_div_RMS = sqrt(err_div_RMS)*wgt                                                        ! RMS in real space calculated with Parsevals theorem from Fourier space
         err_div = err_div_RMS/field_av_L2                                                        ! criterion to stop iterations


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

           err_real_div_RMS = 0.0_pReal
           err_post_div_RMS = 0.0_pReal
           err_real_div_max = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             err_real_div_RMS = err_real_div_RMS +   sum(divergence_real(i,j,k,1:3)**2.0_pReal)     ! avg of squared L_2 norm of div(stress) in real space
             err_post_div_RMS = err_post_div_RMS +   sum(divergence_post(i,j,k,1:3)**2.0_pReal)     ! avg of squared L_2 norm of div(stress) in real space
             err_real_div_max = max(err_real_div_max,sum(divergence_real(i,j,k,1:3)**2.0_pReal))    ! max of squared L_2 norm of div(stress) in real space
           enddo; enddo; enddo

           err_real_div_RMS = sqrt(wgt*err_real_div_RMS)                                            ! RMS in real space
           err_post_div_RMS = sqrt(wgt*err_post_div_RMS)                                            ! RMS in real space
           err_real_div_max = sqrt(    err_real_div_max)                                            ! max in real space
           err_div_max      = sqrt(    err_div_max)                                                 ! max in Fourier space
           
           write(6,'(a,es11.4)')        'error divergence  FT  RMS = ',err_div_RMS
           write(6,'(a,es11.4)')        'error divergence Real RMS = ',err_real_div_RMS
           write(6,'(a,es11.4)')        'error divergence post RMS = ',err_post_div_RMS
           write(6,'(a,es11.4)')        'error divergence  FT  max = ',err_div_max
           write(6,'(a,es11.4)')        'error divergence Real max = ',err_real_div_max
         endif
         write(6,'(a,f6.2,a,es11.4,a)') 'error divergence = ', err_div/err_div_tol,&
                                                            ' (',err_div,' N/m³)'
         end if
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
                 forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, q=1_pInt:3_pInt)&
                   gamma_hat(1,1,1, l,m,n,q) =  temp33_Real(l,n)*xiDyad(m,q)
                 forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                   temp33_Complex(l,m) = sum(gamma_hat(1,1,1, l,m, 1:3,1:3) *&
                                                                field_fourier(i,j,k,1:3,1:3))
                 field_fourier(i,j,k,1:3,1:3) = temp33_Complex 
             endif             
           enddo; enddo; enddo
   
         else                                                                                       ! use precalculated gamma-operator
           
           do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
             forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt) &
               temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n, 1:3,1:3) *&
                                                                field_fourier(i,j,k,1:3,1:3))
             field_fourier(i,j,k, 1:3,1:3) = temp33_Complex
           enddo; enddo; enddo

         endif
         field_fourier(1,1,1,1:3,1:3) = cmplx(field_aim,0.0_pReal,pReal)             ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
         if (debugFFTW) then
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
              scalarField_fourier(i,j,k) = field_fourier(i,j,k,row,column)
           enddo; enddo; enddo
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
!--------------------------------------------------------------------------------------------------
! doing the inverse FT
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

!--------------------------------------------------------------------------------------------------
! calculate some additional output
         if(debugGeneral) then
           maxCorrectionSkew = 0.0_pReal
           maxCorrectionSym  = 0.0_pReal
           temp33_Real = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             maxCorrectionSym  = max(maxCorrectionSym,&
                                     maxval(math_symmetric33(field_real(i,j,k,1:3,1:3))))
             maxCorrectionSkew = max(maxCorrectionSkew,&
                                     maxval(math_skew33(field_real(i,j,k,1:3,1:3))))
             temp33_Real = temp33_Real + field_real(i,j,k,1:3,1:3)
           enddo; enddo; enddo
           write(6,'(a,1x,es11.4)') 'max symmetric correction of deformation =',&
                                         maxCorrectionSym*wgt
           write(6,'(a,1x,es11.4)') 'max skew      correction of deformation =',&
                                         maxCorrectionSkew*wgt
           write(6,'(a,1x,es11.4)') 'max sym/skew of avg correction =         ',&
                                         maxval(math_symmetric33(temp33_real))/&
                                         maxval(math_skew33(temp33_real))
         endif
         field_real = field_real * wgt
         convolution = err_div/err_div_tol
end function convolution


function S_lastInc(rot_BC,mask_stressVector1)
 real(pReal), dimension(3,3,3,3) :: S_lastInc
 integer(pInt) :: i, j, k, m,n 
 real(pReal), dimension(3,3), intent(in) :: rot_BC
 logical, dimension(9), intent(in) :: mask_stressVector1
 real(pReal), dimension(3,3,3,3) :: C_lastInc
 real(pReal), dimension(9,9) ::           temp99_Real   
 integer(pInt) :: size_reduced = 0_pInt 
 real(pReal), dimension(:,:), allocatable ::  s_reduced, c_reduced                                  ! reduced compliance and stiffness (only for stress BC)
 logical :: errmatinv
 size_reduced = count(mask_stressVector1)
 if (allocated(c_reduced))        deallocate(c_reduced)
 if (allocated(c_reduced))        deallocate(c_reduced)
 allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
 allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)



 C_lastInc = math_rotate_forward3333(C,rot_BC)                               ! calculate stiffness from former inc
 temp99_Real = math_Plain3333to99(C_lastInc)
 k = 0_pInt                                                                                 ! build reduced stiffness
 do n = 1_pInt,9_pInt
   if(mask_stressVector1(n)) then
     k = k + 1_pInt
     j = 0_pInt
     do m = 1_pInt,9_pInt
       if(mask_stressVector1(m)) then
         j = j + 1_pInt
         c_reduced(k,j) = temp99_Real(n,m)
 endif; enddo; endif; enddo
 call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)                         ! invert reduced stiffness
 if(errmatinv) call IO_error(error_ID=400_pInt)
 temp99_Real = 0.0_pReal                                                                    ! build full compliance
  k = 0_pInt
  do n = 1_pInt,9_pInt
    if(mask_stressVector1(n)) then
      k = k + 1_pInt
      j = 0_pInt
      do m = 1_pInt,9_pInt
        if(mask_stressVector1(m)) then
          j = j + 1_pInt
          temp99_Real(n,m) = s_reduced(k,j)
 endif; enddo; endif; enddo
 S_lastInc = (math_Plain99to3333(temp99_Real))
 if (allocated(c_reduced))        deallocate(c_reduced)
 if (allocated(c_reduced))        deallocate(c_reduced)

end function S_lastInc


!--------------------------------------------------------------------------------------------------
! calculate reduced compliance

real(pReal) function BCcorrection(mask_stressVector,P_BC,S_lastInc)
use numerics, only: err_stress_tolrel, err_stress_tolabs
 logical, dimension(9) :: mask_stressVector
 real(pReal) ::  err_stress, err_stress_tol     
 real(pReal), dimension(3,3), parameter ::  ones = 1.0_pReal, zeroes = 0.0_pReal   
       real(pReal), dimension(3,3,3,3) :: S_lastInc
       real(pReal), dimension(3,3) :: &
       P_BC , &
   mask_stress, &
   mask_defgrad
        mask_stress = merge(ones,zeroes,reshape(mask_stressVector,[3,3]))                                   
 mask_defgrad  = merge(zeroes,ones,reshape(mask_stressVector,[3,3]))
 
!--------------------------------------------------------------------------------------------------
! stress BC handling
                                                    ! calculate stress BC if applied
           err_stress = maxval(abs(mask_stress * (P_av - P_BC)))                     ! maximum deviaton (tensor norm not applicable)
           err_stress_tol = min(maxval(abs(P_av)) * err_stress_tolrel,err_stress_tolabs)            ! don't use any tensor norm for the relative criterion because the comparison should be coherent
           write(6,'(a)') '' 
           write(6,'(a)') '... correcting deformation gradient to fulfill BCs ...............'
           write(6,'(a,f6.2,a,es11.4,a)') 'error stress = ', err_stress/err_stress_tol, &
                                                                            ' (',err_stress,' Pa)'  
           F_aim = F_aim - math_mul3333xx33(S_lastInc, ((P_av - P_BC)))              ! residual on given stress components
           write(6,'(a,1x,es11.4)')'determinant of new deformation = ',math_det33(F_aim)
          BCcorrection = err_stress/err_stress_tol
end function BCcorrection

subroutine constitutiveResponse(F,P,ForwardData,timeinc)
 use debug, only: &
   debug_reset, &
   debug_info 
    use CPFEM, only: &
   CPFEM_general
use FEsolving, only: restartWrite
  real(pReal) :: timeinc
  logical :: ForwardData
  integer(pInt) :: i, j, k, ielem
  integer(pInt) :: CPFEM_mode
   real(pReal), dimension(3,3,3,3) :: dPdF
 real(pReal), dimension(6)       :: sigma                                                                ! cauchy stress
 real(pReal), dimension(6,6)     :: dsde
 real(pReal), dimension(3,3)     :: P_av_lab, rotation_BC
  if (ForwardData) then
    CPFEM_mode = 1_pInt
   else
    CPFEM_mode = 2_pInt
  endif
   write(6,'(a)') ''
   write(6,'(a)') '... update stress field P(F) .....................................'
   ielem = 0_pInt
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     ielem = ielem + 1_pInt
     call CPFEM_general(3_pInt,&                                                              ! collect cycle
                        coordinates(i,j,k,1:3), F_lastInc(i,j,k,1:3,1:3),F(i,j,k,1:3,1:3), &
                        temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde,P(i,j,k,1:3,1:3),dPdF)
   enddo; enddo; enddo

   P = 0.0_pReal                                                                         ! needed because of the padding for FFTW
   C = 0.0_pReal
   P_av_lab = 0.0_pReal
   ielem = 0_pInt 
   call debug_reset()
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     ielem = ielem + 1_pInt
     call CPFEM_general(CPFEM_mode,&                                                          ! first element in first iteration retains CPFEM_mode 1, 
                        coordinates(i,j,k,1:3),F_lastInc(i,j,k,1:3,1:3), F(i,j,k,1:3,1:3), &  ! others get 2 (saves winding forward effort)
                        temperature(i,j,k),timeinc,ielem,1_pInt,sigma,dsde,P(i,j,k,1:3,1:3),dPdF)
     CPFEM_mode = 2_pInt
     C = C + dPdF
     P_av_lab = P_av_lab + P(i,j,k,1:3,1:3)
  enddo; enddo; enddo
  call debug_info()
  restartWrite = .false.
  P_av_lab = P_av_lab * wgt
  P_av = math_rotate_forward33(P_av_lab,rotation_BC)
  write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') 'Piola-Kirchhoff stress / MPa =',&
                                                          math_transpose33(P_av)/1.e6_pReal
  C = C * wgt
end subroutine constitutiveResponse
       
end module DAMASK_spectral_Utilities
