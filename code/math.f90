!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mathematical library, including random number generation and tensor represenations
!--------------------------------------------------------------------------------------------------
module math
 use, intrinsic :: iso_c_binding
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 real(pReal),    parameter, public :: PI = 3.14159265358979323846264338327950288419716939937510_pReal !< ratio of a circle's circumference to its diameter
 real(pReal),    parameter, public :: INDEG = 180.0_pReal/PI                                        !< conversion from radian into degree
 real(pReal),    parameter, public :: INRAD = PI/180.0_pReal                                        !< conversion from degree into radian
 complex(pReal), parameter, public :: TWOPIIMG = (0.0_pReal,2.0_pReal)* PI                          !< Re(0.0), Im(2xPi)

 real(pReal), dimension(3,3), parameter, public :: &
   MATH_I3 = reshape([&
     1.0_pReal,0.0_pReal,0.0_pReal, &
     0.0_pReal,1.0_pReal,0.0_pReal, &
     0.0_pReal,0.0_pReal,1.0_pReal  &
     ],[3,3])                                                                                       !< 3x3 Identity

 integer(pInt), dimension (2,6), parameter, private :: &
   mapMandel = reshape([&
     1_pInt,1_pInt, &
     2_pInt,2_pInt, &
     3_pInt,3_pInt, &
     1_pInt,2_pInt, &
     2_pInt,3_pInt, &
     1_pInt,3_pInt  &
     ],[2,6])                                                                                       !< arrangement in Mandel notation

 real(pReal), dimension(6), parameter, private :: &
   nrmMandel = [&
     1.0_pReal,                1.0_pReal,                1.0_pReal,&
     1.414213562373095_pReal,  1.414213562373095_pReal,  1.414213562373095_pReal ]                  !< weighting for Mandel notation (forward)

 real(pReal), dimension(6), parameter , public :: &
   invnrmMandel = [&
     1.0_pReal,                1.0_pReal,                1.0_pReal,&
     0.7071067811865476_pReal, 0.7071067811865476_pReal, 0.7071067811865476_pReal ]                 !< weighting for Mandel notation (backward)

 integer(pInt), dimension (2,6), parameter, private :: &
   mapVoigt = reshape([&
     1_pInt,1_pInt, &
     2_pInt,2_pInt, &
     3_pInt,3_pInt, &
     2_pInt,3_pInt, &
     1_pInt,3_pInt, &
     1_pInt,2_pInt  &
     ],[2,6])                                                                                       !< arrangement in Voigt notation

 real(pReal), dimension(6), parameter, private :: &
   nrmVoigt    = 1.0_pReal, &                                                                       !< weighting for Voigt notation (forward)
   invnrmVoigt = 1.0_pReal                                                                          !< weighting for Voigt notation (backward)

 integer(pInt), dimension (2,9), parameter, private :: &
   mapPlain = reshape([&
     1_pInt,1_pInt, &
     1_pInt,2_pInt, &
     1_pInt,3_pInt, &
     2_pInt,1_pInt, &
     2_pInt,2_pInt, &
     2_pInt,3_pInt, &
     3_pInt,1_pInt, &
     3_pInt,2_pInt, &
     3_pInt,3_pInt  &
     ],[2,9])                                                                                       !< arrangement in Plain notation

#ifdef Spectral
 include 'fftw3.f03'
#endif 

 public :: &
   math_init, &
   math_qsort, &
   math_range, &
   math_identity2nd, &
   math_identity4th, &
   math_civita, &
   math_delta, &
   math_crossproduct, &
   math_tensorproduct33, &
   math_mul3x3, &
   math_mul6x6, &
   math_mul33xx33, &
   math_mul3333xx33, &
   math_mul3333xx3333, &
   math_mul33x33, &
   math_mul66x66, &
   math_mul99x99, &
   math_mul33x3, &
   math_mul33x3_complex, &
   math_mul66x6 , &
   math_exp33 , &
   math_transpose33, &
   math_inv33, &
   math_invert33, &
   math_invSym3333, &
   math_invert, &
   math_symmetric33, &
   math_symmetric66, &
   math_skew33, &
   math_spherical33, &
   math_deviatoric33, &
   math_equivStrain33, &
   math_equivStress33, &
   math_trace33, &
   math_j3_33, &
   math_det33, &
   math_Plain33to9, &
   math_Plain9to33, &
   math_Mandel33to6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_Plain99to3333, &
   math_Mandel66toPlain66, &
   math_Plain66toMandel66, &
   math_Mandel3333to66, &
   math_Mandel66to3333, &
   math_Voigt66to3333, &
   math_qRand, &
   math_qMul, &
   math_qDot, &
   math_qConj, &
   math_qInv, &
   math_qRot, &
   math_RtoEuler, &
   math_RtoQ, &
   math_EulerToR, &
   math_EulerToQ, &
   math_EulerAxisAngleToR, &
   math_axisAngleToR, &
   math_EulerAxisAngleToQ, &
   math_axisAngleToQ, &
   math_qToRodrig, &
   math_qToEuler, &
   math_qToEulerAxisAngle, &
   math_qToAxisAngle, &
   math_qToR, &
   math_EulerMisorientation, &
   math_sampleRandomOri, &
   math_sampleGaussOri, &
   math_sampleFiberOri, &
   math_sampleGaussVar, &
   math_symmetricEulers, &
   math_spectralDecompositionSym33, &
   math_spectralDecompositionSym, &
   math_rotationalPart33, &
   math_invariants33, &
   math_eigenvaluesSym33, &
   math_factorial, &
   math_binomial, &
   math_multinomial, &
   math_volTetrahedron, &
   math_areaTriangle, &
   math_rotate_forward33, &
   math_rotate_backward33, &
   math_rotate_forward3333
#ifdef Spectral
 public :: &
   fftw_set_timelimit, &
   fftw_plan_dft_3d, &
   fftw_plan_many_dft_r2c, &
   fftw_plan_many_dft_c2r, &
   fftw_plan_with_nthreads, &
   fftw_init_threads, &
   fftw_alloc_complex, &
   fftw_execute_dft, &
   fftw_execute_dft_r2c, &
   fftw_execute_dft_c2r, &
   fftw_destroy_plan, &
   math_tensorAvg
#endif     
 private :: &
   math_partition, &
   halton, &
   halton_memory, &
   halton_ndim_set, &
   halton_seed_set, &
   i_to_halton, &
   prime
 external :: &
   dsyev, &
   dgetrf, &
   dgetri

contains

!--------------------------------------------------------------------------------------------------
!> @brief initialization of random seed generator
!--------------------------------------------------------------------------------------------------
subroutine math_init

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec,     only: tol_math_check
 use numerics, only: &
   worldrank, &
   fixedSeed
 use IO,       only: IO_error, IO_timeStamp

 implicit none
 integer(pInt) :: i
 real(pReal), dimension(3,3) :: R,R2
 real(pReal), dimension(3) ::   Eulers,v
 real(pReal), dimension(4) ::   q,q2,axisangle,randTest
! the following variables are system dependend and shound NOT be pInt
 integer :: randSize                                                                                ! gfortran requires a variable length to compile
 integer, dimension(:), allocatable :: randInit                                                     ! if recalculations of former randomness (with given seed) is necessary
                                                                                                    ! comment the first random_seed call out, set randSize to 1, and use ifort
 character(len=64) :: error_msg

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  math init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 call random_seed(size=randSize)
 if (allocated(randInit)) deallocate(randInit)
 allocate(randInit(randSize))
 if (fixedSeed > 0_pInt) then
   randInit(1:randSize) = int(fixedSeed)                                                            ! fixedSeed is of type pInt, randInit not
   call random_seed(put=randInit)
 else
   call random_seed()
   call random_seed(get = randInit)
   randInit(2:randSize) = randInit(1)
   call random_seed(put = randInit)
 endif

 do i = 1_pInt, 4_pInt
   call random_number(randTest(i))
 enddo

 mainProcess2: if (worldrank == 0) then 
   write(6,*) 'size  of random seed:    ', randSize
   do i =1, randSize
     write(6,*) 'value of random seed:    ', i, randInit(i)
   enddo
   write(6,'(a,4(/,26x,f17.14),/)') ' start of random sequence: ', randTest
 endif mainProcess2 

 call random_seed(put = randInit)

 call halton_seed_set(int(randInit(1), pInt))
 call halton_ndim_set(3_pInt)

 ! --- check rotation dictionary ---

 q = math_qRand()          ! random quaternion
 
 ! +++ q -> a -> q  +++
 axisangle = math_qToAxisAngle(q)
 q2 = math_axisAngleToQ(axisangle(1:3),axisangle(4))
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ q -> R -> q  +++
 R = math_qToR(q)
 q2 = math_RtoQ(R)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(402_pInt,ext_msg=error_msg)
 endif

 ! +++ q -> euler -> q  +++
 Eulers = math_qToEuler(q)
 q2 = math_EulerToQ(Eulers)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(403_pInt,ext_msg=error_msg)
 endif

 ! +++ R -> euler -> R  +++
 Eulers = math_RtoEuler(R)
 R2 = math_EulerToR(Eulers)
 if ( any(abs( R-R2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation ',maxval(abs( R-R2))
   call IO_error(404_pInt,ext_msg=error_msg)
 endif

 ! +++ check rotation sense of q and R +++
 q = math_qRand()          ! random quaternion
 call halton(3_pInt,v)     ! random vector
 R = math_qToR(q)
 if (any(abs(math_mul33x3(R,v) - math_qRot(q,v)) > tol_math_check)) then
   write(6,'(a,4(f8.3,1x))') 'q',q
   call IO_error(409_pInt)
 endif

end subroutine math_init
 

!--------------------------------------------------------------------------------------------------
!> @brief Quicksort algorithm for two-dimensional integer arrays
! Sorting is done with respect to array(1,:)
! and keeps array(2:N,:) linked to it.
!--------------------------------------------------------------------------------------------------
recursive subroutine math_qsort(a, istart, iend)

 implicit none
 integer(pInt), dimension(:,:), intent(inout) :: a
 integer(pInt), intent(in) :: istart,iend
 integer(pInt) :: ipivot

 if (istart < iend) then
   ipivot = math_partition(a,istart, iend)
   call math_qsort(a, istart, ipivot-1_pInt)
   call math_qsort(a, ipivot+1_pInt, iend)
 endif

end subroutine math_qsort


!--------------------------------------------------------------------------------------------------
!> @brief Partitioning required for quicksort
!--------------------------------------------------------------------------------------------------
integer(pInt) function math_partition(a, istart, iend)

 implicit none
 integer(pInt), dimension(:,:), intent(inout) :: a
 integer(pInt), intent(in) :: istart,iend
 integer(pInt) :: d,i,j,k,x,tmp

 d = int(size(a,1_pInt), pInt) ! number of linked data
! set the starting and ending points, and the pivot point

 i = istart

 j = iend
 x = a(1,istart)
 do
! find the first element on the right side less than or equal to the pivot point
   do j = j, istart, -1_pInt
     if (a(1,j) <= x) exit
   enddo
! find the first element on the left side greater than the pivot point
   do i = i, iend
     if (a(1,i) > x) exit
   enddo
   if (i < j) then ! if the indexes do not cross, exchange values
     do k = 1_pInt,d
      tmp = a(k,i)
      a(k,i) = a(k,j)
      a(k,j) = tmp
     enddo
   else           ! if they do cross, exchange left value with pivot and return with the partition index
     do k = 1_pInt,d
      tmp = a(k,istart)
      a(k,istart) = a(k,j)
      a(k,j) = tmp
     enddo
     math_partition = j
     return
   endif
 enddo

end function math_partition


!--------------------------------------------------------------------------------------------------
!> @brief range of integers starting at one
!--------------------------------------------------------------------------------------------------
pure function math_range(N)

 implicit none
 integer(pInt), intent(in) :: N                                                                     !< length of range
 integer(pInt) :: i
 integer(pInt), dimension(N) :: math_range

 math_range = [(i,i=1_pInt,N)]

end function math_range


!--------------------------------------------------------------------------------------------------
!> @brief second rank identity tensor of specified dimension
!--------------------------------------------------------------------------------------------------
pure function math_identity2nd(dimen)

 implicit none
 integer(pInt), intent(in) :: dimen                                                                 !< tensor dimension
 integer(pInt) :: i
 real(pReal), dimension(dimen,dimen) :: math_identity2nd

 math_identity2nd = 0.0_pReal
 forall (i=1_pInt:dimen) math_identity2nd(i,i) = 1.0_pReal

end function math_identity2nd

!--------------------------------------------------------------------------------------------------
!> @brief symmetric fourth rank identity tensor of specified dimension
!  from http://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)#Derivative_of_a_second-order_tensor_with_respect_to_itself
!--------------------------------------------------------------------------------------------------
pure function math_identity4th(dimen)

 implicit none
 integer(pInt), intent(in) :: dimen                                                                 !< tensor dimension
 integer(pInt) :: i,j,k,l
 real(pReal), dimension(dimen,dimen,dimen,dimen) ::  math_identity4th

 forall (i=1_pInt:dimen,j=1_pInt:dimen,k=1_pInt:dimen,l=1_pInt:dimen) math_identity4th(i,j,k,l) = &
        0.5_pReal*(math_I3(i,k)*math_I3(j,l)+math_I3(i,l)*math_I3(j,k))

end function math_identity4th


!--------------------------------------------------------------------------------------------------
!> @brief permutation tensor e_ijk used for computing cross product of two tensors
! e_ijk =  1 if even permutation of ijk
! e_ijk = -1 if odd permutation of ijk
! e_ijk =  0 otherwise
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_civita(i,j,k)

 implicit none
 integer(pInt), intent(in) :: i,j,k

 math_civita = 0.0_pReal
 if (((i == 1_pInt).and.(j == 2_pInt).and.(k == 3_pInt)) .or. &
     ((i == 2_pInt).and.(j == 3_pInt).and.(k == 1_pInt)) .or. &
     ((i == 3_pInt).and.(j == 1_pInt).and.(k == 2_pInt))) math_civita = 1.0_pReal
 if (((i == 1_pInt).and.(j == 3_pInt).and.(k == 2_pInt)) .or. &
     ((i == 2_pInt).and.(j == 1_pInt).and.(k == 3_pInt)) .or. &
     ((i == 3_pInt).and.(j == 2_pInt).and.(k == 1_pInt))) math_civita = -1.0_pReal

end function math_civita


!--------------------------------------------------------------------------------------------------
!> @brief kronecker delta function d_ij
! d_ij = 1 if i = j
! d_ij = 0 otherwise
! inspired by http://fortraninacworld.blogspot.de/2012/12/ternary-operator.html
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_delta(i,j)

 implicit none
 integer(pInt), intent (in) :: i,j

 math_delta = merge(0.0_pReal, 1.0_pReal, i /= j)

end function math_delta


!--------------------------------------------------------------------------------------------------
!> @brief cross product a x b
!--------------------------------------------------------------------------------------------------
pure function math_crossproduct(A,B)

 implicit none
 real(pReal), dimension(3), intent(in) ::  A,B
 real(pReal), dimension(3) ::  math_crossproduct

 math_crossproduct = [ A(2)*B(3) -A(3)*B(2), &
                       A(3)*B(1) -A(1)*B(3), &
                       A(1)*B(2) -A(2)*B(1) ]

end function math_crossproduct


!--------------------------------------------------------------------------------------------------
!> @brief tensor product a \otimes b
!--------------------------------------------------------------------------------------------------
pure function math_tensorproduct33(A,B)

 implicit none
 real(pReal), dimension(3,3) ::  math_tensorproduct33
 real(pReal), dimension(3), intent(in) ::  A,B
 integer(pInt) :: i,j

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_tensorproduct33(i,j) = A(i)*B(j)

end function math_tensorproduct33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3x3 = 1
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_mul3x3(A,B)

 implicit none
 real(pReal), dimension(3), intent(in) ::  A,B

 math_mul3x3 = sum(A*B)

end function math_mul3x3


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 6x6 = 1
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_mul6x6(A,B)

 implicit none
 real(pReal), dimension(6), intent(in) ::  A,B

 math_mul6x6 = sum(A*B)

end function math_mul6x6


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 33xx33 = 1 (double contraction --> ij * ij)
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_mul33xx33(A,B)

 implicit none
 real(pReal), dimension(3,3), intent(in) ::  A,B
 integer(pInt) :: i,j
 real(pReal), dimension(3,3) ::              C

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) C(i,j) = A(i,j) * B(i,j)
 math_mul33xx33 = sum(C)

end function math_mul33xx33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3333x33 = 33 (double contraction --> ijkl *kl = ij)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx33(A,B)

 implicit none
 real(pReal), dimension(3,3) :: math_mul3333xx33
 real(pReal), dimension(3,3,3,3), intent(in) ::  A
 real(pReal), dimension(3,3), intent(in) ::  B
 integer(pInt) :: i,j 

 forall(i = 1_pInt:3_pInt,j = 1_pInt:3_pInt) &
   math_mul3333xx33(i,j) = sum(A(i,j,1:3,1:3)*B(1:3,1:3))
   
end function math_mul3333xx33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3333x3333 = 3333 (ijkl *klmn = ijmn)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx3333(A,B)

 implicit none
 integer(pInt) :: i,j,k,l
 real(pReal), dimension(3,3,3,3), intent(in) ::  A
 real(pReal), dimension(3,3,3,3), intent(in) ::  B
 real(pReal), dimension(3,3,3,3) :: math_mul3333xx3333

 forall(i = 1_pInt:3_pInt,j = 1_pInt:3_pInt, k = 1_pInt:3_pInt, l= 1_pInt:3_pInt) &
   math_mul3333xx3333(i,j,k,l) = sum(A(i,j,1:3,1:3)*B(1:3,1:3,k,l))

end function math_mul3333xx3333


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 33x33 = 33
!--------------------------------------------------------------------------------------------------
pure function math_mul33x33(A,B)

 implicit none
 real(pReal), dimension(3,3) ::  math_mul33x33
 real(pReal), dimension(3,3), intent(in) ::  A,B
 integer(pInt) :: i,j

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) &
   math_mul33x33(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j)

end function math_mul33x33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 66x66 = 66
!--------------------------------------------------------------------------------------------------
pure function math_mul66x66(A,B)

 implicit none
 real(pReal), dimension(6,6) ::  math_mul66x66
 real(pReal), dimension(6,6), intent(in) ::  A,B
 integer(pInt) :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) math_mul66x66(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) + &
   A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)

end function math_mul66x66


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 99x99 = 99
!--------------------------------------------------------------------------------------------------
pure function math_mul99x99(A,B)

 implicit none
 real(pReal), dimension(9,9) ::  math_mul99x99
 real(pReal), dimension(9,9), intent(in) ::  A,B
 integer(pInt)  i,j

 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_mul99x99(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) + &
   A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j) + &
   A(i,7)*B(7,j) + A(i,8)*B(8,j) + A(i,9)*B(9,j)

end function math_mul99x99


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 33x3 = 3
!--------------------------------------------------------------------------------------------------
pure function math_mul33x3(A,B)

 implicit none
 real(pReal), dimension(3) ::  math_mul33x3
 real(pReal), dimension(3,3), intent(in) ::  A
 real(pReal), dimension(3),   intent(in) ::  B
 integer(pInt) :: i 

 forall (i=1_pInt:3_pInt) math_mul33x3(i) = sum(A(i,1:3)*B)

end function math_mul33x3


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication complex(33) x real(3) = complex(3)
!--------------------------------------------------------------------------------------------------
pure function math_mul33x3_complex(A,B)

 implicit none
 complex(pReal), dimension(3) ::  math_mul33x3_complex
 complex(pReal), dimension(3,3), intent(in) ::  A
 real(pReal),    dimension(3),   intent(in) ::  B
 integer(pInt) :: i 

 forall (i=1_pInt:3_pInt) math_mul33x3_complex(i) = sum(A(i,1:3)*cmplx(B,0.0_pReal,pReal))

end function math_mul33x3_complex


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 66x6 = 6
!--------------------------------------------------------------------------------------------------
pure function math_mul66x6(A,B)

 implicit none
 real(pReal), dimension(6) ::  math_mul66x6
 real(pReal), dimension(6,6), intent(in) ::  A
 real(pReal), dimension(6),   intent(in) ::  B
 integer(pInt) :: i

 forall (i=1_pInt:6_pInt) math_mul66x6(i) = &
   A(i,1)*B(1) + A(i,2)*B(2) + A(i,3)*B(3) + &
   A(i,4)*B(4) + A(i,5)*B(5) + A(i,6)*B(6)

end function math_mul66x6


!--------------------------------------------------------------------------------------------------
!> @brief 3x3 matrix exponential up to series approximation order n (default 5)
!--------------------------------------------------------------------------------------------------
pure function math_exp33(A,n)

 implicit none
 integer(pInt) :: i,order
 integer(pInt), intent(in), optional :: n
 real(pReal), dimension(3,3), intent(in) ::  A
 real(pReal), dimension(3,3) ::  B,math_exp33
 real(pReal) :: invfac

 order = merge(n,5_pInt,present(n))

 B = math_I3                                                                                        ! init
 invfac = 1.0_pReal                                                                                 ! 0!
 math_exp33 = B                                                                                     ! A^0 = eye2

 do i = 1_pInt,n
   invfac = invfac/real(i)                                                                          ! invfac = 1/i!
   B = math_mul33x33(B,A)
   math_exp33 = math_exp33 + invfac*B                                                               ! exp = SUM (A^i)/i!
 enddo

end function math_exp33


!--------------------------------------------------------------------------------------------------
!> @brief transposition of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_transpose33(A)

 implicit none 
 real(pReal),dimension(3,3) :: math_transpose33
 real(pReal),dimension(3,3),intent(in)  :: A
 integer(pInt) :: i,j

 forall(i=1_pInt:3_pInt, j=1_pInt:3_pInt) math_transpose33(i,j) = A(j,i)

end function math_transpose33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 33 matrix (function)
!   direct Cramer inversion of matrix A.
!   returns all zeroes if not possible, i.e. if det close to zero
!--------------------------------------------------------------------------------------------------
pure function math_inv33(A)

 implicit none
 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal) :: DetA
 real(pReal),dimension(3,3) :: math_inv33

 math_inv33(1,1) =  A(2,2) * A(3,3) - A(2,3) * A(3,2)
 math_inv33(2,1) = -A(2,1) * A(3,3) + A(2,3) * A(3,1)
 math_inv33(3,1) =  A(2,1) * A(3,2) - A(2,2) * A(3,1)

 DetA = A(1,1) * math_inv33(1,1) + A(1,2) * math_inv33(2,1) + A(1,3) * math_inv33(3,1)

 if (abs(DetA) > tiny(DetA)) then                                             ! use a real threshold here
   math_inv33(1,2) = -A(1,2) * A(3,3) + A(1,3) * A(3,2)
   math_inv33(2,2) =  A(1,1) * A(3,3) - A(1,3) * A(3,1)
   math_inv33(3,2) = -A(1,1) * A(3,2) + A(1,2) * A(3,1)

   math_inv33(1,3) =  A(1,2) * A(2,3) - A(1,3) * A(2,2)
   math_inv33(2,3) = -A(1,1) * A(2,3) + A(1,3) * A(2,1)
   math_inv33(3,3) =  A(1,1) * A(2,2) - A(1,2) * A(2,1)
   
   math_inv33 = math_inv33/DetA
 else
   math_inv33 = 0.0_pReal
 endif

end function math_inv33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 33 matrix (subroutine)
!   direct Cramer inversion of matrix A.
!   also returns determinant
!   returns error if not possible, i.e. if det close to zero
!--------------------------------------------------------------------------------------------------
pure subroutine math_invert33(A, InvA, DetA, error)

 implicit none
 logical, intent(out) :: error
 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal),dimension(3,3),intent(out) :: InvA
 real(pReal), intent(out) :: DetA

 InvA(1,1) =  A(2,2) * A(3,3) - A(2,3) * A(3,2)
 InvA(2,1) = -A(2,1) * A(3,3) + A(2,3) * A(3,1)
 InvA(3,1) =  A(2,1) * A(3,2) - A(2,2) * A(3,1)

 DetA = A(1,1) * InvA(1,1) + A(1,2) * InvA(2,1) + A(1,3) * InvA(3,1)

 if (abs(DetA) <= tiny(DetA)) then
   error = .true.
 else
   InvA(1,2) = -A(1,2) * A(3,3) + A(1,3) * A(3,2)
   InvA(2,2) =  A(1,1) * A(3,3) - A(1,3) * A(3,1)
   InvA(3,2) = -A(1,1) * A(3,2) + A(1,2) * A(3,1)

   InvA(1,3) =  A(1,2) * A(2,3) - A(1,3) * A(2,2)
   InvA(2,3) = -A(1,1) * A(2,3) + A(1,3) * A(2,1)
   InvA(3,3) =  A(1,1) * A(2,2) - A(1,2) * A(2,1)

   InvA = InvA/DetA
   error = .false.
 endif

end subroutine math_invert33


!--------------------------------------------------------------------------------------------------
!> @brief Inversion of symmetriced 3x3x3x3 tensor.
!--------------------------------------------------------------------------------------------------
function math_invSym3333(A)
 use IO, only: &
   IO_error

 implicit none
 real(pReal),dimension(3,3,3,3)            :: math_invSym3333

 real(pReal),dimension(3,3,3,3),intent(in) :: A

 integer(pInt) :: ierr
 integer(pInt), dimension(6)   :: ipiv6
 real(pReal),   dimension(6,6) :: temp66_Real
 real(pReal),   dimension(6)   :: work6

 temp66_real = math_Mandel3333to66(A)
#if(FLOAT==8)
 call dgetrf(6,6,temp66_real,6,ipiv6,ierr)
 call dgetri(6,temp66_real,6,ipiv6,work6,6,ierr)
#elif(FLOAT==4)
 call sgetrf(6,6,temp66_real,6,ipiv6,ierr)
 call sgetri(6,temp66_real,6,ipiv6,work6,6,ierr)
#endif
 if (ierr == 0_pInt) then
   math_invSym3333 = math_Mandel66to3333(temp66_real)
 else
   call IO_error(400_pInt, ext_msg = 'math_invSym3333')
 endif

end function math_invSym3333


!--------------------------------------------------------------------------------------------------
!> @brief invert matrix of arbitrary dimension
!--------------------------------------------------------------------------------------------------
subroutine math_invert(myDim,A, InvA, error)

 implicit none
 integer(pInt), intent(in)  :: myDim
 real(pReal), dimension(myDim,myDim), intent(in)  :: A


 integer(pInt) :: ierr
 integer(pInt), dimension(myDim)       :: ipiv
 real(pReal),   dimension(myDim)       :: work
 
 real(pReal), dimension(myDim,myDim), intent(out) :: invA
 logical, intent(out) :: error
 
 invA = A 
#if(FLOAT==8)
 call dgetrf(myDim,myDim,invA,myDim,ipiv,ierr)
 call dgetri(myDim,InvA,myDim,ipiv,work,myDim,ierr)
#elif(FLOAT==4)
 call sgetrf(myDim,myDim,invA,myDim,ipiv,ierr)
 call sgetri(myDim,InvA,myDim,ipiv,work,myDim,ierr)
#endif
 error = merge(.true.,.false., ierr /= 0_pInt)                                                      ! http://fortraninacworld.blogspot.de/2012/12/ternary-operator.html
 
end subroutine math_invert


!--------------------------------------------------------------------------------------------------
!> @brief symmetrize a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_symmetric33(m)

 implicit none
 real(pReal), dimension(3,3) :: math_symmetric33
 real(pReal), dimension(3,3), intent(in) :: m
 
 math_symmetric33 = 0.5_pReal * (m + transpose(m))

end function math_symmetric33


!--------------------------------------------------------------------------------------------------
!> @brief symmetrize a 66 matrix
!--------------------------------------------------------------------------------------------------
pure function math_symmetric66(m)

 implicit none
 real(pReal), dimension(6,6) :: math_symmetric66
 real(pReal), dimension(6,6), intent(in) :: m

 math_symmetric66 = 0.5_pReal * (m + transpose(m))

end function math_symmetric66


!--------------------------------------------------------------------------------------------------
!> @brief skew part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_skew33(m)

 implicit none
 real(pReal), dimension(3,3) :: math_skew33
 real(pReal), dimension(3,3), intent(in) :: m

 math_skew33 = m - math_symmetric33(m)

end function math_skew33

!--------------------------------------------------------------------------------------------------
!> @brief hydrostatic part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_spherical33(m)

 implicit none
 real(pReal), dimension(3,3) :: math_spherical33
 real(pReal), dimension(3,3), intent(in) :: m

 math_spherical33 = math_I3 * math_trace33(m)/3.0_pReal

end function math_spherical33


!--------------------------------------------------------------------------------------------------
!> @brief deviatoric part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_deviatoric33(m)

 implicit none
 real(pReal), dimension(3,3) :: math_deviatoric33
 real(pReal), dimension(3,3), intent(in) :: m

 math_deviatoric33 = m - math_spherical33(m)

end function math_deviatoric33


!--------------------------------------------------------------------------------------------------
!> @brief equivalent scalar quantity of a full symmetric strain tensor
!--------------------------------------------------------------------------------------------------
pure function math_equivStrain33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal), dimension(3) :: e,s
 real(pReal) :: math_equivStrain33
 real(pReal), parameter :: TWOTHIRD = 2.0_pReal/3.0_pReal

 e = [2.0_pReal*m(1,1)-m(2,2)-m(3,3), &
      2.0_pReal*m(2,2)-m(3,3)-m(1,1), &
      2.0_pReal*m(3,3)-m(1,1)-m(2,2)]/3.0_pReal
 s = [m(1,2),m(2,3),m(1,3)]*2.0_pReal

 math_equivStrain33 = TWOTHIRD*(1.50_pReal*(sum(e**2.0_pReal)) + &
                                0.75_pReal*(sum(s**2.0_pReal)))**(0.5_pReal)

end function math_equivStrain33


!--------------------------------------------------------------------------------------------------
!> @brief von Mises equivalent of a full symmetric stress tensor
!--------------------------------------------------------------------------------------------------
pure function math_equivStress33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: math_equivStress33

 math_equivStress33 =( ( (m(1,1)-m(2,2))**2.0_pReal + &
                         (m(2,2)-m(3,3))**2.0_pReal + &
                         (m(3,3)-m(1,1))**2.0_pReal + &
                         6.0_pReal*( m(1,2)**2.0_pReal + &
                                     m(2,3)**2.0_pReal + &
                                     m(1,3)**2.0_pReal &
                                   ) &
                       )**0.5_pReal &
                     )/sqrt(2.0_pReal)

end function math_equivStress33


!--------------------------------------------------------------------------------------------------
!> @brief trace of a 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_trace33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m

 math_trace33 = m(1,1) + m(2,2) + m(3,3)

end function math_trace33


!--------------------------------------------------------------------------------------------------
!> @brief invariant 3 of a 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_j3_33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m

 math_j3_33 = sum(m**3.0_pReal)/3.0_pReal

end function math_j3_33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_det33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m
 
 math_det33 = m(1,1)* (m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
            - m(1,2)* (m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
            + m(1,3)* (m(2,1)*m(3,2)-m(2,2)*m(3,1))

end function math_det33


!--------------------------------------------------------------------------------------------------
!> @brief convert 33 matrix into vector 9
!--------------------------------------------------------------------------------------------------
pure function math_Plain33to9(m33)

 implicit none
 real(pReal), dimension(9) :: math_Plain33to9
 real(pReal), dimension(3,3), intent(in) :: m33
 integer(pInt) :: i

 forall (i=1_pInt:9_pInt) math_Plain33to9(i) = m33(mapPlain(1,i),mapPlain(2,i))

end function math_Plain33to9


!--------------------------------------------------------------------------------------------------
!> @brief convert Plain 9 back to 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_Plain9to33(v9)

 implicit none
 real(pReal), dimension(3,3) :: math_Plain9to33
 real(pReal), dimension(9), intent(in) :: v9
 integer(pInt) :: i

 forall (i=1_pInt:9_pInt) math_Plain9to33(mapPlain(1,i),mapPlain(2,i)) = v9(i)

end function math_Plain9to33


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 33 matrix into Mandel vector 6
!--------------------------------------------------------------------------------------------------
pure function math_Mandel33to6(m33)

 implicit none
 real(pReal), dimension(6) :: math_Mandel33to6
 real(pReal), dimension(3,3), intent(in) :: m33
 
 integer(pInt) :: i

 forall (i=1_pInt:6_pInt) math_Mandel33to6(i) = nrmMandel(i)*m33(mapMandel(1,i),mapMandel(2,i))

end function math_Mandel33to6


!--------------------------------------------------------------------------------------------------
!> @brief convert Mandel 6 back to symmetric 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_Mandel6to33(v6)

 implicit none
 real(pReal), dimension(6), intent(in) :: v6
 real(pReal), dimension(3,3) :: math_Mandel6to33
 integer(pInt) :: i

 forall (i=1_pInt:6_pInt)
  math_Mandel6to33(mapMandel(1,i),mapMandel(2,i)) = invnrmMandel(i)*v6(i)
  math_Mandel6to33(mapMandel(2,i),mapMandel(1,i)) = invnrmMandel(i)*v6(i)
 end forall

end function math_Mandel6to33


!--------------------------------------------------------------------------------------------------
!> @brief convert 3333 tensor into plain matrix 99
!--------------------------------------------------------------------------------------------------
pure function math_Plain3333to99(m3333)

 implicit none
 real(pReal), dimension(3,3,3,3), intent(in) :: m3333
 real(pReal), dimension(9,9) :: math_Plain3333to99
 integer(pInt) :: i,j

 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_Plain3333to99(i,j) = &
   m3333(mapPlain(1,i),mapPlain(2,i),mapPlain(1,j),mapPlain(2,j))

end function math_Plain3333to99

!--------------------------------------------------------------------------------------------------
!> @brief plain matrix 99 into 3333 tensor
!--------------------------------------------------------------------------------------------------
pure function math_Plain99to3333(m99)

 implicit none
 real(pReal), dimension(9,9), intent(in) :: m99
 real(pReal), dimension(3,3,3,3) :: math_Plain99to3333
 integer(pInt) :: i,j

 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_Plain99to3333(mapPlain(1,i),mapPlain(2,i),&
     mapPlain(1,j),mapPlain(2,j)) = m99(i,j)

end function math_Plain99to3333


!--------------------------------------------------------------------------------------------------
!> @brief convert Mandel matrix 66 into Plain matrix 66
!--------------------------------------------------------------------------------------------------
pure function math_Mandel66toPlain66(m66)

 implicit none
 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(6,6) :: math_Mandel66toPlain66
 integer(pInt) :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) &
   math_Mandel66toPlain66(i,j) = invnrmMandel(i) * invnrmMandel(j) * m66(i,j)

end function math_Mandel66toPlain66


!--------------------------------------------------------------------------------------------------
!> @brief convert Plain matrix 66 into Mandel matrix 66
!--------------------------------------------------------------------------------------------------
pure function math_Plain66toMandel66(m66)

 implicit none
 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(6,6) :: math_Plain66toMandel66
 integer(pInt)  :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) &
   math_Plain66toMandel66(i,j) = nrmMandel(i) * nrmMandel(j) * m66(i,j)

end function math_Plain66toMandel66


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 3333 tensor into Mandel matrix 66
!--------------------------------------------------------------------------------------------------
pure function math_Mandel3333to66(m3333)

 implicit none

 real(pReal), dimension(3,3,3,3), intent(in) :: m3333
 real(pReal), dimension(6,6) :: math_Mandel3333to66
 integer(pInt) :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) math_Mandel3333to66(i,j) = &
   nrmMandel(i)*nrmMandel(j)*m3333(mapMandel(1,i),mapMandel(2,i),mapMandel(1,j),mapMandel(2,j))

end function math_Mandel3333to66


!--------------------------------------------------------------------------------------------------
!> @brief convert Mandel matrix 66 back to symmetric 3333 tensor
!--------------------------------------------------------------------------------------------------
pure function math_Mandel66to3333(m66)

 implicit none
 real(pReal), dimension(3,3,3,3) :: math_Mandel66to3333
 real(pReal), dimension(6,6), intent(in) :: m66
 integer(pInt) :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt)
   math_Mandel66to3333(mapMandel(1,i),mapMandel(2,i),mapMandel(1,j),mapMandel(2,j)) = &
                              invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(2,i),mapMandel(1,i),mapMandel(1,j),mapMandel(2,j)) = &
                              invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(1,i),mapMandel(2,i),mapMandel(2,j),mapMandel(1,j)) = &
                              invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(2,i),mapMandel(1,i),mapMandel(2,j),mapMandel(1,j)) = &
                              invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
 end forall

end function math_Mandel66to3333


!--------------------------------------------------------------------------------------------------
!> @brief convert Voigt matrix 66 back to symmetric 3333 tensor
!--------------------------------------------------------------------------------------------------
pure function math_Voigt66to3333(m66)

 implicit none
 real(pReal), dimension(3,3,3,3) :: math_Voigt66to3333
 real(pReal), dimension(6,6), intent(in) :: m66
 integer(pInt) :: i,j

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt)
   math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(1,j),mapVoigt(2,j)) =  &
                                                            invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(1,j),mapVoigt(2,j)) = &
                                                            invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(2,j),mapVoigt(1,j)) = &
                                                            invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(2,j),mapVoigt(1,j)) = & 
                                                            invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
 end forall

end function math_Voigt66to3333


!--------------------------------------------------------------------------------------------------
!> @brief random quaternion
!--------------------------------------------------------------------------------------------------
function math_qRand()

 implicit none
 real(pReal), dimension(4) :: math_qRand
 real(pReal), dimension(3) :: rnd

 call halton(3_pInt,rnd)
 math_qRand(1) = cos(2.0_pReal*PI*rnd(1))*sqrt(rnd(3))
 math_qRand(2) = sin(2.0_pReal*PI*rnd(2))*sqrt(1.0_pReal-rnd(3))
 math_qRand(3) = cos(2.0_pReal*PI*rnd(2))*sqrt(1.0_pReal-rnd(3))
 math_qRand(4) = sin(2.0_pReal*PI*rnd(1))*sqrt(rnd(3))

end function math_qRand


!--------------------------------------------------------------------------------------------------
!> @brief quaternion multiplication q1xq2 = q12
!--------------------------------------------------------------------------------------------------
pure function math_qMul(A,B)

 implicit none
 real(pReal), dimension(4) ::  math_qMul
 real(pReal), dimension(4), intent(in) ::  A, B

 math_qMul = [ A(1)*B(1) - A(2)*B(2) - A(3)*B(3) - A(4)*B(4), &
               A(1)*B(2) + A(2)*B(1) + A(3)*B(4) - A(4)*B(3), &
               A(1)*B(3) - A(2)*B(4) + A(3)*B(1) + A(4)*B(2), &
               A(1)*B(4) + A(2)*B(3) - A(3)*B(2) + A(4)*B(1) ]

end function math_qMul


!--------------------------------------------------------------------------------------------------
!> @brief quaternion dotproduct
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_qDot(A,B)

 implicit none
 real(pReal), dimension(4), intent(in) :: A, B

 math_qDot = sum(A*B)

end function math_qDot


!--------------------------------------------------------------------------------------------------
!> @brief quaternion conjugation
!--------------------------------------------------------------------------------------------------
pure function math_qConj(Q)

 implicit none
 real(pReal), dimension(4) ::  math_qConj
 real(pReal), dimension(4), intent(in) ::  Q

 math_qConj = [Q(1), -Q(2:4)]

end function math_qConj


!--------------------------------------------------------------------------------------------------
!> @brief quaternion norm
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_qNorm(Q)

 implicit none
 real(pReal), dimension(4), intent(in) ::  Q

 math_qNorm = norm2(Q)

end function math_qNorm


!--------------------------------------------------------------------------------------------------
!> @brief quaternion inversion
!--------------------------------------------------------------------------------------------------
pure function math_qInv(Q)

 implicit none
 real(pReal), dimension(4), intent(in) ::  Q
 real(pReal), dimension(4) ::  math_qInv
 real(pReal) :: squareNorm

 math_qInv = 0.0_pReal

 squareNorm = math_qDot(Q,Q)
 if (abs(squareNorm) > tiny(squareNorm)) &
   math_qInv = math_qConj(Q) / squareNorm

end function math_qInv


!--------------------------------------------------------------------------------------------------
!> @brief action of a quaternion on a vector (rotate vector v with Q)
!--------------------------------------------------------------------------------------------------
pure function math_qRot(Q,v)

 implicit none
 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3), intent(in) :: v
 real(pReal), dimension(3) :: math_qRot
 real(pReal), dimension(4,4) :: T
 integer(pInt) :: i, j

 do i = 1_pInt,4_pInt
   do j = 1_pInt,i
     T(i,j) = Q(i) * Q(j)
   enddo
 enddo

 math_qRot = [-v(1)*(T(3,3)+T(4,4)) + v(2)*(T(3,2)-T(4,1)) + v(3)*(T(4,2)+T(3,1)), &
               v(1)*(T(3,2)+T(4,1)) - v(2)*(T(2,2)+T(4,4)) + v(3)*(T(4,3)-T(2,1)), &
               v(1)*(T(4,2)-T(3,1)) + v(2)*(T(4,3)+T(2,1)) - v(3)*(T(2,2)+T(3,3))]

 math_qRot = 2.0_pReal * math_qRot + v

end function math_qRot


!--------------------------------------------------------------------------------------------------
!> @brief Euler angles (in radians) from rotation matrix
!> @details rotation matrix is meant to represent a PASSIVE rotation, 
!> composed of INTRINSIC rotations around the axes of the 
!> rotating reference frame 
!> (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_RtoEuler(R)

 implicit none
 real(pReal), dimension (3,3), intent(in) :: R
 real(pReal), dimension(3) :: math_RtoEuler
 real(pReal) :: sqhkl, squvw, sqhk

 sqhkl=sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3)+R(3,3)*R(3,3))
 squvw=sqrt(R(1,1)*R(1,1)+R(2,1)*R(2,1)+R(3,1)*R(3,1))
 sqhk =sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3))

! calculate PHI
 math_RtoEuler(2) = acos(math_limit(R(3,3)/sqhkl,-1.0_pReal, 1.0_pReal))

 if((math_RtoEuler(2) < 1.0e-8_pReal) .or. (pi-math_RtoEuler(2) < 1.0e-8_pReal)) then
   math_RtoEuler(3) = 0.0_pReal
   math_RtoEuler(1) = acos(math_limit(R(1,1)/squvw, -1.0_pReal, 1.0_pReal))
   if(R(2,1) > 0.0_pReal) math_RtoEuler(1) = 2.0_pReal*pi-math_RtoEuler(1)
 else
   math_RtoEuler(3) = acos(math_limit(R(2,3)/sqhk, -1.0_pReal, 1.0_pReal))
   if(R(1,3) < 0.0) math_RtoEuler(3) = 2.0_pReal*pi-math_RtoEuler(3)
   math_RtoEuler(1) = acos(math_limit(-R(3,2)/sin(math_RtoEuler(2)), -1.0_pReal, 1.0_pReal))
   if(R(3,1) < 0.0) math_RtoEuler(1) = 2.0_pReal*pi-math_RtoEuler(1)
 end if

end function math_RtoEuler


!--------------------------------------------------------------------------------------------------
!> @brief converts a rotation matrix into a quaternion (w+ix+jy+kz)
!> @details math adopted from http://arxiv.org/pdf/math/0701759v1.pdf
!--------------------------------------------------------------------------------------------------
pure function math_RtoQ(R)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: R
 real(pReal), dimension(4)   :: absQ, math_RtoQ
 real(pReal) :: max_absQ
 integer, dimension(1) :: largest !no pInt, maxloc returns integer default

 math_RtoQ = 0.0_pReal

 absQ  = [+ R(1,1) + R(2,2) + R(3,3), &
          + R(1,1) - R(2,2) - R(3,3), &
          - R(1,1) + R(2,2) - R(3,3), &
          - R(1,1) - R(2,2) + R(3,3)] + 1.0_pReal

 largest = maxloc(absQ)

 largestComponent: select case(largest(1))
   case (1) largestComponent
      !1----------------------------------
      math_RtoQ(2) = R(3,2) - R(2,3)
      math_RtoQ(3) = R(1,3) - R(3,1)
      math_RtoQ(4) = R(2,1) - R(1,2)
                                 
   case (2) largestComponent
      math_RtoQ(1) = R(3,2) - R(2,3)
      !2----------------------------------
      math_RtoQ(3) = R(2,1) + R(1,2)
      math_RtoQ(4) = R(1,3) + R(3,1)
                                 
   case (3) largestComponent
      math_RtoQ(1) = R(1,3) - R(3,1)
      math_RtoQ(2) = R(2,1) + R(1,2)
      !3----------------------------------
      math_RtoQ(4) = R(3,2) + R(2,3)
                                 
   case (4) largestComponent
      math_RtoQ(1) = R(2,1) - R(1,2)
      math_RtoQ(2) = R(1,3) + R(3,1)
      math_RtoQ(3) = R(2,3) + R(3,2)
      !4----------------------------------
 end select  largestComponent

 max_absQ = 0.5_pReal * sqrt(absQ(largest(1)))
 math_RtoQ = math_RtoQ * 0.25_pReal / max_absQ
 math_RtoQ(largest(1)) = max_absQ

end function math_RtoQ


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from Euler angles (in radians)
!> @details rotation matrix is meant to represent a PASSIVE rotation, 
!> @details composed of INTRINSIC rotations around the axes of the 
!> @details rotating reference frame 
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_EulerToR(Euler)

 implicit none
 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_EulerToR
 real(pReal) c1, c, c2, s1, s, s2

 C1 = cos(Euler(1))
 C  = cos(Euler(2))
 C2 = cos(Euler(3))
 S1 = sin(Euler(1))
 S  = sin(Euler(2))
 S2 = sin(Euler(3))

 math_EulerToR(1,1)=C1*C2-S1*S2*C
 math_EulerToR(1,2)=-C1*S2-S1*C2*C
 math_EulerToR(1,3)=S1*S
 math_EulerToR(2,1)=S1*C2+C1*S2*C
 math_EulerToR(2,2)=-S1*S2+C1*C2*C
 math_EulerToR(2,3)=-C1*S
 math_EulerToR(3,1)=S2*S
 math_EulerToR(3,2)=C2*S
 math_EulerToR(3,3)=C
 
 math_EulerToR = transpose(math_EulerToR)  ! convert to passive rotation

end function math_EulerToR


!--------------------------------------------------------------------------------------------------
!> @brief quaternion (w+ix+jy+kz) from 3-1-3 Euler angles (in radians)
!> @details quaternion is meant to represent a PASSIVE rotation, 
!> @details composed of INTRINSIC rotations around the axes of the 
!> @details rotating reference frame 
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_EulerToQ(eulerangles)

 implicit none
 real(pReal), dimension(3), intent(in) :: eulerangles
 real(pReal), dimension(4) :: math_EulerToQ
 real(pReal), dimension(3) :: halfangles
 real(pReal) :: c, s

 halfangles = 0.5_pReal * eulerangles

 c = cos(halfangles(2))
 s = sin(halfangles(2))

 math_EulerToQ= [cos(halfangles(1)+halfangles(3)) * c, &
                 cos(halfangles(1)-halfangles(3)) * s, &
                 sin(halfangles(1)-halfangles(3)) * s, &
                 sin(halfangles(1)+halfangles(3)) * c ]
 math_EulerToQ = math_qConj(math_EulerToQ)                 ! convert to passive rotation

end function math_EulerToQ


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from axis and angle (in radians)
!> @details rotation matrix is meant to represent a ACTIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from http://mathworld.wolfram.com/RodriguesRotationFormula.html
!--------------------------------------------------------------------------------------------------
pure function math_axisAngleToR(axis,omega)

 implicit none
 real(pReal), dimension(3,3) :: math_axisAngleToR
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: axisNrm
 real(pReal) :: norm,s,c,c1

 norm = sqrt(math_mul3x3(axis,axis))
 if (norm > 1.0e-8_pReal) then                             ! non-zero rotation
   axisNrm = axis/norm                                     ! normalize axis to be sure

   s = sin(omega)
   c = cos(omega)
   c1 = 1.0_pReal - c

   math_axisAngleToR(1,1) =  c + c1*axisNrm(1)**2.0_pReal
   math_axisAngleToR(1,2) = -s*axisNrm(3) + c1*axisNrm(1)*axisNrm(2)
   math_axisAngleToR(1,3) =  s*axisNrm(2) + c1*axisNrm(1)*axisNrm(3)
                             
   math_axisAngleToR(2,1) =  s*axisNrm(3) + c1*axisNrm(2)*axisNrm(1)
   math_axisAngleToR(2,2) =  c + c1*axisNrm(2)**2.0_pReal
   math_axisAngleToR(2,3) = -s*axisNrm(1) + c1*axisNrm(2)*axisNrm(3)
                             
   math_axisAngleToR(3,1) = -s*axisNrm(2) + c1*axisNrm(3)*axisNrm(1)
   math_axisAngleToR(3,2) =  s*axisNrm(1) + c1*axisNrm(3)*axisNrm(2)
   math_axisAngleToR(3,3) =  c + c1*axisNrm(3)**2.0_pReal
 else
   math_axisAngleToR = math_I3
 endif
 
end function math_axisAngleToR


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from axis and angle (in radians)
!> @details rotation matrix is meant to represent a PASSIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_EulerAxisAngleToR(axis,omega)

 implicit none
 real(pReal), dimension(3,3) :: math_EulerAxisAngleToR
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega

 math_EulerAxisAngleToR = transpose(math_axisAngleToR(axis,omega))  ! convert to passive rotation

end function math_EulerAxisAngleToR


!--------------------------------------------------------------------------------------------------
!> @brief quaternion (w+ix+jy+kz) from Euler axis and angle (in radians)
!> @details quaternion is meant to represent a PASSIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from 
!> @details http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29#Rodrigues_parameters
!--------------------------------------------------------------------------------------------------
pure function math_EulerAxisAngleToQ(axis,omega)

 implicit none
 real(pReal), dimension(4) :: math_EulerAxisAngleToQ
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega

 math_EulerAxisAngleToQ = math_qConj(math_axisAngleToQ(axis,omega))                                   ! convert to passive rotation

end function math_EulerAxisAngleToQ


!--------------------------------------------------------------------------------------------------
!> @brief quaternion (w+ix+jy+kz) from axis and angle (in radians)
!> @details quaternion is meant to represent an ACTIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from 
!> @details http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29#Rodrigues_parameters
!--------------------------------------------------------------------------------------------------
pure function math_axisAngleToQ(axis,omega)

 implicit none
 real(pReal), dimension(4) :: math_axisAngleToQ
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: axisNrm
 real(pReal) :: norm

 norm = sqrt(math_mul3x3(axis,axis))
 rotation: if (norm > 1.0e-8_pReal) then
   axisNrm = axis/norm                                                                                ! normalize axis to be sure
   math_axisAngleToQ = [cos(0.5_pReal*omega), sin(0.5_pReal*omega) * axisNrm(1:3)]
 else rotation
   math_axisAngleToQ = [1.0_pReal,0.0_pReal,0.0_pReal,0.0_pReal]
 endif rotation

end function math_axisAngleToQ


!--------------------------------------------------------------------------------------------------
!> @brief orientation matrix from quaternion (w+ix+jy+kz)
!> @details taken from http://arxiv.org/pdf/math/0701759v1.pdf
!> @details see also http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
!--------------------------------------------------------------------------------------------------
pure function math_qToR(q)

 implicit none
 real(pReal), dimension(4), intent(in) :: q
 real(pReal), dimension(3,3) :: math_qToR, T,S
 integer(pInt) :: i, j

 forall (i = 1_pInt:3_pInt, j = 1_pInt:3_pInt) &
   T(i,j) = q(i+1_pInt) * q(j+1_pInt)

 S = reshape( [0.0_pReal,     -q(4),      q(3), &
                    q(4), 0.0_pReal,     -q(2), &
                   -q(3),      q(2), 0.0_pReal],[3,3])                                             ! notation is transposed

 math_qToR = (2.0_pReal * q(1)*q(1) - 1.0_pReal) * math_I3 &
           + 2.0_pReal * T - 2.0_pReal * q(1) * S

end function math_qToR


!--------------------------------------------------------------------------------------------------
!> @brief 3-1-3 Euler angles (in radians) from quaternion (w+ix+jy+kz)
!> @details quaternion is meant to represent a PASSIVE rotation, 
!> @details composed of INTRINSIC rotations around the axes of the 
!> @details rotating reference frame 
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_qToEuler(qPassive)

 implicit none
 real(pReal), dimension(4), intent(in) :: qPassive
 real(pReal), dimension(4) :: q
 real(pReal), dimension(3) :: math_qToEuler

 q = math_qConj(qPassive)    ! convert to active rotation, since formulas are defined for active rotations

 math_qToEuler(2) = acos(1.0_pReal-2.0_pReal*(q(2)*q(2)+q(3)*q(3)))

 if (abs(math_qToEuler(2)) < 1.0e-6_pReal) then
   math_qToEuler(1) = sign(2.0_pReal*acos(math_limit(q(1),-1.0_pReal, 1.0_pReal)),q(4))
   math_qToEuler(3) = 0.0_pReal
 else
   math_qToEuler(1) = atan2(q(1)*q(3)+q(2)*q(4), q(1)*q(2)-q(3)*q(4))
   math_qToEuler(3) = atan2(-q(1)*q(3)+q(2)*q(4), q(1)*q(2)+q(3)*q(4))
 endif

 math_qToEuler = merge(math_qToEuler + [2.0_pReal*PI, PI, 2.0_pReal*PI], &                           ! ensure correct range
                       math_qToEuler,                                    math_qToEuler<0.0_pReal)

end function math_qToEuler


!--------------------------------------------------------------------------------------------------
!> @brief axis-angle (x, y, z, ang in radians) from quaternion (w+ix+jy+kz)
!> @details quaternion is meant to represent an ACTIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from 
!> @details http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29#Rodrigues_parameters
!--------------------------------------------------------------------------------------------------
pure function math_qToAxisAngle(Q)

 implicit none
 real(pReal), dimension(4), intent(in) :: Q
 real(pReal) :: halfAngle, sinHalfAngle
 real(pReal), dimension(4) :: math_qToAxisAngle

 halfAngle = acos(max(-1.0_pReal, min(1.0_pReal, Q(1))))            ! limit to [-1,1] --> 0 to 180 deg
 sinHalfAngle = sin(halfAngle)

 if (sinHalfAngle <= 1.0e-4_pReal) then                              ! very small rotation angle?
   math_qToAxisAngle = 0.0_pReal
 else
   math_qToAxisAngle= [ Q(2:4)/sinHalfAngle, halfAngle*2.0_pReal] 
 endif

end function math_qToAxisAngle


!--------------------------------------------------------------------------------------------------
!> @brief Euler axis-angle (x, y, z, ang in radians) from quaternion (w+ix+jy+kz)
!> @details quaternion is meant to represent a PASSIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!--------------------------------------------------------------------------------------------------
pure function math_qToEulerAxisAngle(qPassive)

 implicit none
 real(pReal), dimension(4), intent(in) :: qPassive
 real(pReal), dimension(4) :: q
 real(pReal), dimension(4) :: math_qToEulerAxisAngle

 q = math_qConj(qPassive)                        ! convert to active rotation
 math_qToEulerAxisAngle = math_qToAxisAngle(q)
 
end function math_qToEulerAxisAngle


!--------------------------------------------------------------------------------------------------
!> @brief Rodrigues vector (x, y, z) from unit quaternion (w+ix+jy+kz)
!--------------------------------------------------------------------------------------------------
pure function math_qToRodrig(Q)
 use prec, only: &
   DAMASK_NaN, &
   tol_math_check

 implicit none
 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3) :: math_qToRodrig

 math_qToRodrig = merge(Q(2:4)/Q(1),DAMASK_NaN,abs(Q(1)) > tol_math_check)                          ! NaN for 180 deg since Rodrig is unbound

end function math_qToRodrig


!--------------------------------------------------------------------------------------------------
!> @brief misorientation angle between two sets of Euler angles
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_EulerMisorientation(EulerA,EulerB)

 implicit none
 real(pReal), dimension(3), intent(in) :: EulerA,EulerB
 real(pReal), dimension(3,3) :: r
 real(pReal) :: tr

 r = math_mul33x33(math_EulerToR(EulerB),transpose(math_EulerToR(EulerA)))

 tr = (math_trace33(r)-1.0_pReal)*0.4999999_pReal
 math_EulerMisorientation = abs(0.5_pReal*PI-asin(tr))

end function math_EulerMisorientation


!--------------------------------------------------------------------------------------------------
!> @brief draw a random sample from Euler space
!--------------------------------------------------------------------------------------------------
function math_sampleRandomOri()

 implicit none
 real(pReal), dimension(3) :: math_sampleRandomOri, rnd

 call halton(3_pInt,rnd)
 math_sampleRandomOri  = [rnd(1)*2.0_pReal*PI, &
                          acos(2.0_pReal*rnd(2)-1.0_pReal), &
                          rnd(3)*2.0_pReal*PI]

end function math_sampleRandomOri


!--------------------------------------------------------------------------------------------------
!> @brief draw a random sample from Gauss component with noise (in radians) half-width
!--------------------------------------------------------------------------------------------------
function math_sampleGaussOri(center,noise)
 use prec, only: &
   tol_math_check

 implicit none
 real(pReal), intent(in) :: noise
 real(pReal), dimension(3), intent(in) :: center 
 real(pReal) :: cosScatter,scatter
 real(pReal), dimension(3) :: math_sampleGaussOri, disturb
 real(pReal), dimension(3), parameter :: ORIGIN = [0.0_pReal,0.0_pReal,0.0_pReal]
 real(pReal), dimension(5) :: rnd
 integer(pInt) :: i

 if (abs(noise) < tol_math_check) then
   math_sampleGaussOri = center
   return
 endif

! Helming uses different distribution with Bessel functions
! therefore the gauss scatter width has to be scaled differently
 scatter = 0.95_pReal * noise
 cosScatter = cos(scatter)

 do
   call halton(5_pInt,rnd)
   forall (i=1_pInt:3_pInt) rnd(i) = 2.0_pReal*rnd(i)-1.0_pReal  ! expand 1:3 to range [-1,+1]
   disturb  = [ scatter * rnd(1), &                                                       ! phi1
                sign(1.0_pReal,rnd(2))*acos(cosScatter+(1.0_pReal-cosScatter)*rnd(4)), &  ! Phi
                scatter * rnd(2)]                                                         ! phi2
   if (rnd(5) <= exp(-1.0_pReal*(math_EulerMisorientation(ORIGIN,disturb)/scatter)**2_pReal)) exit
 enddo

 math_sampleGaussOri = math_RtoEuler(math_mul33x33(math_EulerToR(disturb),math_EulerToR(center)))

end function math_sampleGaussOri


!--------------------------------------------------------------------------------------------------
!> @brief draw a random sample from Fiber component with noise (in radians)
!--------------------------------------------------------------------------------------------------
function math_sampleFiberOri(alpha,beta,noise)
 use prec, only: &
   tol_math_check

 implicit none
 real(pReal), dimension(3) :: math_sampleFiberOri, fiberInC,fiberInS,axis
 real(pReal), dimension(2), intent(in) :: alpha,beta
 real(pReal), dimension(6) :: rnd
 real(pReal), dimension(3,3) :: oRot,fRot,pRot
 real(pReal) :: noise, scatter, cos2Scatter, angle
 integer(pInt), dimension(2,3), parameter :: ROTMAP = reshape([2_pInt,3_pInt,&
                                                               3_pInt,1_pInt,&
                                                               1_pInt,2_pInt],[2,3])
 integer(pInt) :: i

! Helming uses different distribution with Bessel functions
! therefore the gauss scatter width has to be scaled differently
 scatter = 0.95_pReal * noise
 cos2Scatter = cos(2.0_pReal*scatter)

! fiber axis in crystal coordinate system
 fiberInC = [ sin(alpha(1))*cos(alpha(2)) , &
              sin(alpha(1))*sin(alpha(2)), &
              cos(alpha(1))]
! fiber axis in sample coordinate system
 fiberInS = [ sin(beta(1))*cos(beta(2)), &
              sin(beta(1))*sin(beta(2)), &
              cos(beta(1))]

! ---# rotation matrix from sample to crystal system #---
 angle = -acos(dot_product(fiberInC,fiberInS))
 if(abs(angle) > tol_math_check) then
!   rotation axis between sample and crystal system (cross product)
   forall(i=1_pInt:3_pInt) axis(i) = fiberInC(ROTMAP(1,i))*fiberInS(ROTMAP(2,i))-fiberInC(ROTMAP(2,i))*fiberInS(ROTMAP(1,i))
   oRot = math_EulerAxisAngleToR(math_crossproduct(fiberInC,fiberInS),angle)
 else
   oRot = math_I3
 end if

! ---# rotation matrix about fiber axis (random angle) #---
 do
   call halton(6_pInt,rnd)
   fRot = math_EulerAxisAngleToR(fiberInS,rnd(1)*2.0_pReal*pi)

! ---# rotation about random axis perpend to fiber #---
! random axis pependicular to fiber axis
   axis(1:2) = rnd(2:3)
   if (abs(fiberInS(3)) > tol_math_check) then
     axis(3)=-(axis(1)*fiberInS(1)+axis(2)*fiberInS(2))/fiberInS(3)
   else if(abs(fiberInS(2)) > tol_math_check) then
     axis(3)=axis(2)
     axis(2)=-(axis(1)*fiberInS(1)+axis(3)*fiberInS(3))/fiberInS(2)
   else if(abs(fiberInS(1)) > tol_math_check) then
     axis(3)=axis(1)
     axis(1)=-(axis(2)*fiberInS(2)+axis(3)*fiberInS(3))/fiberInS(1)
   end if

! scattered rotation angle
   if (noise > 0.0_pReal) then
     angle = acos(cos2Scatter+(1.0_pReal-cos2Scatter)*rnd(4))
     if (rnd(5) <= exp(-1.0_pReal*(angle/scatter)**2.0_pReal)) exit
   else
     angle = 0.0_pReal
     exit
   end if
 enddo
 if (rnd(6) <= 0.5) angle = -angle
 
 pRot = math_EulerAxisAngleToR(axis,angle)

! ---# apply the three rotations #---
 math_sampleFiberOri = math_RtoEuler(math_mul33x33(pRot,math_mul33x33(fRot,oRot)))

end function math_sampleFiberOri


!--------------------------------------------------------------------------------------------------
!> @brief draw a random sample from Gauss variable
!--------------------------------------------------------------------------------------------------
real(pReal) function math_sampleGaussVar(meanvalue, stddev, width)
 use prec, only: &
   tol_math_check

 implicit none
 real(pReal), intent(in) ::            meanvalue, &      ! meanvalue of gauss distribution
                                       stddev            ! standard deviation of gauss distribution
 real(pReal), intent(in), optional ::  width             ! width of considered values as multiples of standard deviation
 real(pReal), dimension(2) ::          rnd               ! random numbers
 real(pReal) ::                        scatter, &        ! normalized scatter around meanvalue
                                       myWidth

 if (abs(stddev) < tol_math_check) then
   math_sampleGaussVar = meanvalue
   return
 endif

 myWidth = merge(width,3.0_pReal,present(width))                                                    ! use +-3*sigma as default value for scatter if not given

 do
   call halton(2_pInt, rnd)
   scatter = myWidth * (2.0_pReal * rnd(1) - 1.0_pReal)
   if (rnd(2) <= exp(-0.5_pReal * scatter ** 2.0_pReal)) exit                                       ! test if scattered value is drawn
 enddo

 math_sampleGaussVar = scatter * stddev

end function math_sampleGaussVar


!--------------------------------------------------------------------------------------------------
!> @brief symmetrically equivalent Euler angles for given sample symmetry 1:triclinic, 2:monoclinic, 4:orthotropic
!--------------------------------------------------------------------------------------------------
pure function math_symmetricEulers(sym,Euler)

 implicit none
 integer(pInt), intent(in) :: sym
 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_symmetricEulers
 integer(pInt) :: i,j

 math_symmetricEulers(1,1) = PI+Euler(1)
 math_symmetricEulers(2,1) = Euler(2)
 math_symmetricEulers(3,1) = Euler(3)

 math_symmetricEulers(1,2) = PI-Euler(1)
 math_symmetricEulers(2,2) = PI-Euler(2)
 math_symmetricEulers(3,2) = PI+Euler(3)

 math_symmetricEulers(1,3) = 2.0_pReal*PI-Euler(1)
 math_symmetricEulers(2,3) = PI-Euler(2)
 math_symmetricEulers(3,3) = PI+Euler(3)

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_symmetricEulers(j,i) = modulo(math_symmetricEulers(j,i),2.0_pReal*pi)

 select case (sym)
   case (4_pInt) ! all done

   case (2_pInt)  ! return only first
     math_symmetricEulers(1:3,2:3) = 0.0_pReal

   case default         ! return blank
     math_symmetricEulers = 0.0_pReal
 end select

end function math_symmetricEulers


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric matrix m
!--------------------------------------------------------------------------------------------------
subroutine math_spectralDecompositionSym(m,values,vectors,error)

 implicit none
 real(pReal), dimension(:,:),                  intent(in)  :: m
 real(pReal), dimension(size(m,1)),            intent(out) :: values
 real(pReal), dimension(size(m,1),size(m,1)),  intent(out) :: vectors
 logical, intent(out) :: error

 integer(pInt) :: info
 real(pReal), dimension((64+2)*size(m,1)) :: work                                                    ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f

 vectors = M                                                                                         ! copy matrix to input (doubles as output) array
#if(FLOAT==8)
 call dsyev('V','U',size(m,1),vectors,size(m,1),values,work,(64+2)*size(m,1),info)
#elif(FLOAT==4)
 call ssyev('V','U',size(m,1),vectors,size(m,1),values,work,(64+2)*size(m,1),info)
#endif
 error = (info == 0_pInt)

end subroutine math_spectralDecompositionSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric 3x3 matrix m using an analytical expression
!> and the general LAPACK powered version as fallback
!> @author Joachim Kopp, Max–Planck–Institut für Kernphysik, Heidelberg (Copyright (C) 2006)
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details See http://arxiv.org/abs/physics/0610206 (DSYEVH3)
!--------------------------------------------------------------------------------------------------
subroutine math_spectralDecompositionSym33(m,values,vectors)
 
 implicit none
 real(pReal), dimension(3,3),intent(in)  :: m
 real(pReal), dimension(3),  intent(out) :: values
 real(pReal), dimension(3,3),intent(out) :: vectors
 real(pReal) :: T, U, norm, threshold
 logical :: error

 values = math_eigenvaluesSym33(m)

 vectors(1:3,2) = [ m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2), &
                    m(1, 3) * m(1, 2) - m(2, 3) * m(1, 1), &
                    m(1, 2)**2_pInt]

 T = maxval(abs(values))
 U = MAX(T, T**2_pInt)
 threshold = sqrt(5.0e-14_pReal * U**2_pInt)

! Calculate first eigenvector by the formula v[0] = (m - lambda[0]).e1 x (m - lambda[0]).e2
 vectors(1:3,1) = [ vectors(1,2) + m(1, 3) * values(1), &
                    vectors(2,2) + m(2, 3) * values(1), &
                    (m(1,1) - values(1)) * (m(2,2) - values(1)) - vectors(3,2)]
 norm = norm2(vectors(1:3, 1))

 fallback1: if(norm < threshold) then
   call math_spectralDecompositionSym(m,values,vectors,error)
   return
 endif fallback1

 vectors(1:3,1) = vectors(1:3, 1) / norm 
 
! Calculate second eigenvector by the formula v[1] = (m - lambda[1]).e1 x (m - lambda[1]).e2
 vectors(1:3,2) = [ vectors(1,2) + m(1, 3) * values(2), &
                    vectors(2,2) + m(2, 3) * values(2), &
                    (m(1,1) - values(2)) * (m(2,2) - values(2)) - vectors(3,2)]
 norm = norm2(vectors(1:3, 2))

 fallback2: if(norm < threshold) then
   call math_spectralDecompositionSym(m,values,vectors,error)
   return
 endif fallback2
 vectors(1:3,2) = vectors(1:3, 2) / norm

! Calculate third eigenvector according to  v[2] = v[0] x v[1]
 vectors(1:3,3) = math_crossproduct(vectors(1:3,1),vectors(1:3,2))

end subroutine math_spectralDecompositionSym33


!--------------------------------------------------------------------------------------------------
!> @brief rotational part from polar decomposition of tensor m
!--------------------------------------------------------------------------------------------------
function math_rotationalPart33(m)
 use IO, only: &
   IO_warning

 implicit none
 real(pReal), intent(in), dimension(3,3) :: m
 real(pReal), dimension(3,3) :: math_rotationalPart33
 real(pReal), dimension(3,3) :: U, mTm , Uinv, EB
 real(pReal), dimension(3) :: EV

 mTm = math_mul33x33(math_transpose33(m),m)
 call math_spectralDecompositionSym33(mTm,EV,EB)

 U = sqrt(EV(1)) * math_tensorproduct33(EB(1:3,1),EB(1:3,1)) &
   + sqrt(EV(2)) * math_tensorproduct33(EB(1:3,2),EB(1:3,2)) &
   + sqrt(EV(3)) * math_tensorproduct33(EB(1:3,3),EB(1:3,3))

 Uinv = math_inv33(U)
 if (all(abs(Uinv) <= tiny(Uinv))) then                                                             ! math_inv33 returns zero when failed, avoid floating point equality comparison
   math_rotationalPart33 = math_I3
   call IO_warning(650_pInt)
 else
   math_rotationalPart33 = math_mul33x33(m,Uinv)
 endif

end function math_rotationalPart33


!--------------------------------------------------------------------------------------------------
!> @brief Eigenvalues of symmetric matrix m
! will return NaN  on error
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym(m)
 use prec, only: &
   DAMASK_NaN
 
 implicit none
 real(pReal), dimension(:,:),                  intent(in)  :: m
 real(pReal), dimension(size(m,1))                         :: math_eigenvaluesSym
 real(pReal), dimension(size(m,1),size(m,1))               :: vectors

 integer(pInt) :: info
 real(pReal), dimension((64+2)*size(m,1)) :: work                                                    ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f

 vectors = m                                                                                         ! copy matrix to input (doubles as output) array
#if(FLOAT==8)
 call dsyev('N','U',size(m,1),vectors,size(m,1),math_eigenvaluesSym,work,(64+2)*size(m,1),info)
#elif(FLOAT==4)
 call ssyev('N','U',size(m,1),vectors,size(m,1),math_eigenvaluesSym,work,(64+2)*size(m,1),info)
#endif
 if (info /= 0_pInt) math_eigenvaluesSym = DAMASK_NaN

end function math_eigenvaluesSym

!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues of symmetric 3x3 matrix m using an analytical expression
!> @author Joachim Kopp, Max–Planck–Institut für Kernphysik, Heidelberg (Copyright (C) 2006)
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details See http://arxiv.org/abs/physics/0610206 (DSYEVC3)
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in)  :: m
 real(pReal), dimension(3)                :: math_eigenvaluesSym33

 real(pReal), dimension(0:2) :: c
 real(pReal) :: p, q, phi

 c(2) = - math_trace33(m)
 c(1) = m(1,1)*m(2,2)    + m(1,1)*m(3,3)    + m(2,2)*m(3,3) &
      -(m(1,2)**2 + m(2,3)**2 + m(1,3)**2)
 c(0) = m(1,1)*m(2,3)**2 + m(2,2)*m(1,3)**2 + m(3,3)*m(1,2)**2 &
      -(m(1,1)*m(2,2)*m(3,3)  + 2.0_pReal * m(1,3)*m(1,2)*m(2,3))

 p = c(2)**2 -3.0_pReal*c(1)
 q = -13.5_pReal*c(0) -c(2)**3 + 4.5_pReal*c(2)*c(1)

 phi = atan2(  sqrt(abs(27.0_pReal * (0.25_pReal*c(1)**2 *(p-c(1)) +c(0)*(q + 6.75_pReal*c(0)))  )), q)&
       *(1.0_pReal/3.0_pReal)

 math_eigenvaluesSym33 = ([ cos(phi)*2.0_pReal, &
                           -cos(phi)-sqrt(1.0_pReal/3.0_pReal)*sin(phi), &
                           -cos(phi)+sqrt(1.0_pReal/3.0_pReal)*sin(phi) &
                          ] *  sqrt(p) -c(2))/3.0_pReal

                         

end function math_eigenvaluesSym33

!--------------------------------------------------------------------------------------------------
!> @brief invariants of matrix m
!--------------------------------------------------------------------------------------------------
pure function math_invariants33(m)

 implicit none
 real(pReal), dimension(3,3) , intent(in) :: m
 real(pReal), dimension(3) :: math_invariants33

 math_invariants33(1) = math_trace33(m)
 math_invariants33(2) = math_invariants33(1)**2.0_pReal/2.0_pReal &
                      -(m(1,1)**2.0_pReal+m(2,2)**2.0_pReal+m(3,3)**2.0_pReal)* 0.5_pReal &
                      - m(1,2)*m(2,1) -m(1,3)*m(3,1) -m(2,3)*m(3,2)
 math_invariants33(3) = math_det33(m)

end function math_invariants33


!--------------------------------------------------------------------------------------------------
!> @brief computes the next element in the Halton sequence.
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine halton(ndim, r)
 
 implicit none
 integer(pInt), intent(in)                   :: ndim                                                  !< dimension of the element
 real(pReal),   intent(out), dimension(ndim) :: r                                                     !< next element of the current Halton sequence
 integer(pInt), dimension(ndim) :: base
 integer(pInt) :: seed
 integer(pInt), dimension(1) :: value_halton

 call halton_memory ('GET', 'SEED', 1_pInt, value_halton)
 seed = value_halton(1)

 call halton_memory ('GET', 'BASE', ndim, base)

 call i_to_halton (seed, base, ndim, r)

 value_halton(1) = 1_pInt
 call halton_memory ('INC', 'SEED', 1_pInt, value_halton)

end subroutine halton


!--------------------------------------------------------------------------------------------------
!> @brief sets or returns quantities associated with the Halton sequence.
!> @details If action_halton is 'SET' and action_halton is 'BASE', then NDIM is input, and
!> @details is the number of entries in value_halton to be put into BASE.
!> @details If action_halton is 'SET', then on input, value_halton contains values to be assigned
!> @details to the internal variable.
!> @details If action_halton is 'GET', then on output, value_halton contains the values of
!> @details the specified internal variable.
!> @details If action_halton is 'INC', then on input, value_halton contains the increment to
!> @details be added to the specified internal variable.
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine halton_memory (action_halton, name_halton, ndim, value_halton)

 implicit none
 character(len = *), intent(in) :: & 
   action_halton, &                                                                                 !< desired action: GET the value of a particular quantity, SET the value of a particular quantity, INC the value of a particular quantity (only for SEED)
   name_halton                                                                                      !< name of the quantity: BASE: Halton base(s), NDIM: spatial dimension, SEED: current Halton seed
 integer(pInt), dimension(*), intent(inout) :: value_halton
 integer(pInt), allocatable, save, dimension(:) :: base
 logical, save :: first_call = .true.
 integer(pInt), intent(in) :: ndim                                                                  !< dimension of the quantity
 integer(pInt):: i
 integer(pInt), save :: ndim_save = 0_pInt, seed = 1_pInt

 if (first_call) then
   ndim_save = 1_pInt
   allocate(base(ndim_save))
   base(1) = 2_pInt
   first_call = .false.
 endif
 
!--------------------------------------------------------------------------------------------------
! Set
 if(action_halton(1:1) == 'S' .or. action_halton(1:1) == 's') then

   if(name_halton(1:1) == 'B' .or. name_halton(1:1) == 'b') then

     if(ndim_save /= ndim) then
       deallocate(base)
       ndim_save = ndim
       allocate(base(ndim_save))
     endif

     base(1:ndim) = value_halton(1:ndim)

   elseif(name_halton(1:1) == 'N' .or. name_halton(1:1) == 'n') then

     if(ndim_save /= value_halton(1)) then
       deallocate(base)
       ndim_save = value_halton(1)
       allocate(base(ndim_save))
       do i = 1_pInt, ndim_save
         base(i) = prime (i)
       enddo
     else
       ndim_save = value_halton(1)
     endif
   elseif(name_halton(1:1) == 'S' .or. name_halton(1:1) == 's') then
     seed = value_halton(1)
 endif
 
!--------------------------------------------------------------------------------------------------
! Get
 elseif(action_halton(1:1) == 'G' .or. action_halton(1:1) == 'g') then
   if(name_halton(1:1) == 'B' .or. name_halton(1:1) == 'b') then
     if(ndim /= ndim_save) then
  deallocate(base)
  ndim_save = ndim
  allocate(base(ndim_save))
  do i = 1_pInt, ndim_save
    base(i) = prime(i)
  enddo
     endif
     value_halton(1:ndim_save) = base(1:ndim_save)
   elseif(name_halton(1:1) == 'N' .or. name_halton(1:1) == 'n') then
     value_halton(1) = ndim_save
   elseif(name_halton(1:1) == 'S' .or. name_halton(1:1) == 's') then
     value_halton(1) = seed
   endif
   
!--------------------------------------------------------------------------------------------------
!   Increment
 elseif(action_halton(1:1) == 'I' .or. action_halton(1:1) == 'i') then
   if(name_halton(1:1) == 'S' .or. name_halton(1:1) == 's') then
     seed = seed + value_halton(1)
   end if
 endif

end subroutine halton_memory


!--------------------------------------------------------------------------------------------------
!> @brief sets the dimension for a Halton sequence
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine halton_ndim_set (ndim)

 implicit none
 integer(pInt), intent(in) :: ndim                                                                  !< dimension of the Halton vectors
 integer(pInt) :: value_halton(1)

 value_halton(1) = ndim
 call halton_memory ('SET', 'NDIM', 1_pInt, value_halton)

end subroutine halton_ndim_set


!--------------------------------------------------------------------------------------------------
!> @brief  sets the seed for the Halton sequence.
!> @details Calling HALTON repeatedly returns the elements of the Halton sequence in order, 
!> @details starting with element number 1.
!> @details An internal counter, called SEED, keeps track of the next element to return. Each time 
!> @details is computed, and then SEED is incremented by 1.
!> @details To restart the Halton sequence, it is only necessary to reset SEED to 1. It might also 
!> @details be desirable to reset SEED to some other value. This routine allows the user to specify 
!> @details any value of SEED.
!> @details The default value of SEED is 1, which restarts the Halton sequence.
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine halton_seed_set(seed)
 implicit none

 integer(pInt), parameter :: NDIM = 1_pInt
 integer(pInt), intent(in) :: seed                                                                  !< seed for the Halton sequence.
 integer(pInt) :: value_halton(ndim)

 value_halton(1) = seed
 call halton_memory ('SET', 'SEED', NDIM, value_halton)

end subroutine halton_seed_set


!--------------------------------------------------------------------------------------------------
!> @brief computes an element of a Halton sequence.
!> @details Only the absolute value of SEED is considered. SEED = 0 is allowed, and returns R = 0.
!> @details Halton Bases should be distinct prime numbers. This routine only checks that each base 
!> @details is greater than 1.
!> @details Reference:
!> @details J.H. Halton: On the efficiency of certain quasi-random sequences of points in evaluating 
!> @details multi-dimensional integrals, Numerische Mathematik, Volume 2, pages 84-90, 1960.
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine i_to_halton (seed, base, ndim, r)
 use IO, only: &
   IO_error
 
 implicit none
 integer(pInt), intent(in) :: ndim                                                                  !< dimension of the sequence
 integer(pInt), intent(in), dimension(ndim) :: base                                                 !< Halton bases
 real(pReal), dimension(ndim) ::  base_inv
 integer(pInt), dimension(ndim) :: digit
 real(pReal), dimension(ndim), intent(out) ::r                                                      !< the SEED-th element of the Halton sequence for the given bases
 integer(pInt) , intent(in):: seed                                                                  !< index of the desired element
 integer(pInt), dimension(ndim) :: seed2

 seed2(1:ndim) = abs(seed)

 r(1:ndim) = 0.0_pReal

 if (any (base(1:ndim) <= 1_pInt)) call IO_error(error_ID=405_pInt)

 base_inv(1:ndim) = 1.0_pReal / real (base(1:ndim), pReal)

 do while ( any ( seed2(1:ndim) /= 0_pInt) )
   digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim))
   r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), pReal) * base_inv(1:ndim)
   base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), pReal)
   seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
 enddo

end subroutine i_to_halton


!--------------------------------------------------------------------------------------------------
!> @brief returns any of the first 1500 prime numbers.
!> @details n <= 0 returns 1500, the index of the largest prime (12553) available.
!> @details n = 0 is legal, returning PRIME = 1.
!> @details Reference:
!> @details Milton Abramowitz and Irene Stegun: Handbook of Mathematical Functions,
!> @details US Department of Commerce, 1964, pages 870-873.
!> @details Daniel Zwillinger: CRC Standard Mathematical Tables and Formulae,
!> @details 30th Edition, CRC Press, 1996, pages 95-98.
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
integer(pInt) function prime(n)
 use IO, only: &
   IO_error

 implicit none
 integer(pInt), intent(in) :: n                                                                     !< index of the desired prime number
 integer(pInt), parameter :: PRIME_MAX = 1500_pInt
 integer(pInt), save :: icall = 0_pInt
 integer(pInt), save, dimension(PRIME_MAX) :: npvec

 if (icall == 0_pInt) then
   icall = 1_pInt

   npvec = [&
          2_pInt,    3_pInt,    5_pInt,    7_pInt,   11_pInt,   13_pInt,   17_pInt,   19_pInt,   23_pInt,   29_pInt, &
         31_pInt,   37_pInt,   41_pInt,   43_pInt,   47_pInt,   53_pInt,   59_pInt,   61_pInt,   67_pInt,   71_pInt, &
         73_pInt,   79_pInt,   83_pInt,   89_pInt,   97_pInt,  101_pInt,  103_pInt,  107_pInt,  109_pInt,  113_pInt, &
        127_pInt,  131_pInt,  137_pInt,  139_pInt,  149_pInt,  151_pInt,  157_pInt,  163_pInt,  167_pInt,  173_pInt, &
        179_pInt,  181_pInt,  191_pInt,  193_pInt,  197_pInt,  199_pInt,  211_pInt,  223_pInt,  227_pInt,  229_pInt, &
        233_pInt,  239_pInt,  241_pInt,  251_pInt,  257_pInt,  263_pInt,  269_pInt,  271_pInt,  277_pInt,  281_pInt, &
        283_pInt,  293_pInt,  307_pInt,  311_pInt,  313_pInt,  317_pInt,  331_pInt,  337_pInt,  347_pInt,  349_pInt, &
        353_pInt,  359_pInt,  367_pInt,  373_pInt,  379_pInt,  383_pInt,  389_pInt,  397_pInt,  401_pInt,  409_pInt, &
        419_pInt,  421_pInt,  431_pInt,  433_pInt,  439_pInt,  443_pInt,  449_pInt,  457_pInt,  461_pInt,  463_pInt, &
        467_pInt,  479_pInt,  487_pInt,  491_pInt,  499_pInt,  503_pInt,  509_pInt,  521_pInt,  523_pInt,  541_pInt, &
   ! 101:200
         547_pInt,  557_pInt,  563_pInt,  569_pInt,  571_pInt,  577_pInt,  587_pInt,  593_pInt,  599_pInt,  601_pInt, &
         607_pInt,  613_pInt,  617_pInt,  619_pInt,  631_pInt,  641_pInt,  643_pInt,  647_pInt,  653_pInt,  659_pInt, &
         661_pInt,  673_pInt,  677_pInt,  683_pInt,  691_pInt,  701_pInt,  709_pInt,  719_pInt,  727_pInt,  733_pInt, &
         739_pInt,  743_pInt,  751_pInt,  757_pInt,  761_pInt,  769_pInt,  773_pInt,  787_pInt,  797_pInt,  809_pInt, &
         811_pInt,  821_pInt,  823_pInt,  827_pInt,  829_pInt,  839_pInt,  853_pInt,  857_pInt,  859_pInt,  863_pInt, &
         877_pInt,  881_pInt,  883_pInt,  887_pInt,  907_pInt,  911_pInt,  919_pInt,  929_pInt,  937_pInt,  941_pInt, &
         947_pInt,  953_pInt,  967_pInt,  971_pInt,  977_pInt,  983_pInt,  991_pInt,  997_pInt, 1009_pInt, 1013_pInt, &
        1019_pInt, 1021_pInt, 1031_pInt, 1033_pInt, 1039_pInt, 1049_pInt, 1051_pInt, 1061_pInt, 1063_pInt, 1069_pInt, &
        1087_pInt, 1091_pInt, 1093_pInt, 1097_pInt, 1103_pInt, 1109_pInt, 1117_pInt, 1123_pInt, 1129_pInt, 1151_pInt, &
        1153_pInt, 1163_pInt, 1171_pInt, 1181_pInt, 1187_pInt, 1193_pInt, 1201_pInt, 1213_pInt, 1217_pInt, 1223_pInt, &
   ! 201:300
        1229_pInt, 1231_pInt, 1237_pInt, 1249_pInt, 1259_pInt, 1277_pInt, 1279_pInt, 1283_pInt, 1289_pInt, 1291_pInt, &
        1297_pInt, 1301_pInt, 1303_pInt, 1307_pInt, 1319_pInt, 1321_pInt, 1327_pInt, 1361_pInt, 1367_pInt, 1373_pInt, &
        1381_pInt, 1399_pInt, 1409_pInt, 1423_pInt, 1427_pInt, 1429_pInt, 1433_pInt, 1439_pInt, 1447_pInt, 1451_pInt, &
        1453_pInt, 1459_pInt, 1471_pInt, 1481_pInt, 1483_pInt, 1487_pInt, 1489_pInt, 1493_pInt, 1499_pInt, 1511_pInt, &
        1523_pInt, 1531_pInt, 1543_pInt, 1549_pInt, 1553_pInt, 1559_pInt, 1567_pInt, 1571_pInt, 1579_pInt, 1583_pInt, &
        1597_pInt, 1601_pInt, 1607_pInt, 1609_pInt, 1613_pInt, 1619_pInt, 1621_pInt, 1627_pInt, 1637_pInt, 1657_pInt, &
        1663_pInt, 1667_pInt, 1669_pInt, 1693_pInt, 1697_pInt, 1699_pInt, 1709_pInt, 1721_pInt, 1723_pInt, 1733_pInt, &
        1741_pInt, 1747_pInt, 1753_pInt, 1759_pInt, 1777_pInt, 1783_pInt, 1787_pInt, 1789_pInt, 1801_pInt, 1811_pInt, &
        1823_pInt, 1831_pInt, 1847_pInt, 1861_pInt, 1867_pInt, 1871_pInt, 1873_pInt, 1877_pInt, 1879_pInt, 1889_pInt, &
        1901_pInt, 1907_pInt, 1913_pInt, 1931_pInt, 1933_pInt, 1949_pInt, 1951_pInt, 1973_pInt, 1979_pInt, 1987_pInt, &
   ! 301:400
        1993_pInt, 1997_pInt, 1999_pInt, 2003_pInt, 2011_pInt, 2017_pInt, 2027_pInt, 2029_pInt, 2039_pInt, 2053_pInt, &
        2063_pInt, 2069_pInt, 2081_pInt, 2083_pInt, 2087_pInt, 2089_pInt, 2099_pInt, 2111_pInt, 2113_pInt, 2129_pInt, &
        2131_pInt, 2137_pInt, 2141_pInt, 2143_pInt, 2153_pInt, 2161_pInt, 2179_pInt, 2203_pInt, 2207_pInt, 2213_pInt, &
        2221_pInt, 2237_pInt, 2239_pInt, 2243_pInt, 2251_pInt, 2267_pInt, 2269_pInt, 2273_pInt, 2281_pInt, 2287_pInt, &
        2293_pInt, 2297_pInt, 2309_pInt, 2311_pInt, 2333_pInt, 2339_pInt, 2341_pInt, 2347_pInt, 2351_pInt, 2357_pInt, &
        2371_pInt, 2377_pInt, 2381_pInt, 2383_pInt, 2389_pInt, 2393_pInt, 2399_pInt, 2411_pInt, 2417_pInt, 2423_pInt, &
        2437_pInt, 2441_pInt, 2447_pInt, 2459_pInt, 2467_pInt, 2473_pInt, 2477_pInt, 2503_pInt, 2521_pInt, 2531_pInt, &
        2539_pInt, 2543_pInt, 2549_pInt, 2551_pInt, 2557_pInt, 2579_pInt, 2591_pInt, 2593_pInt, 2609_pInt, 2617_pInt, &
        2621_pInt, 2633_pInt, 2647_pInt, 2657_pInt, 2659_pInt, 2663_pInt, 2671_pInt, 2677_pInt, 2683_pInt, 2687_pInt, &
        2689_pInt, 2693_pInt, 2699_pInt, 2707_pInt, 2711_pInt, 2713_pInt, 2719_pInt, 2729_pInt, 2731_pInt, 2741_pInt, &
   ! 401:500
        2749_pInt, 2753_pInt, 2767_pInt, 2777_pInt, 2789_pInt, 2791_pInt, 2797_pInt, 2801_pInt, 2803_pInt, 2819_pInt, &
        2833_pInt, 2837_pInt, 2843_pInt, 2851_pInt, 2857_pInt, 2861_pInt, 2879_pInt, 2887_pInt, 2897_pInt, 2903_pInt, &
        2909_pInt, 2917_pInt, 2927_pInt, 2939_pInt, 2953_pInt, 2957_pInt, 2963_pInt, 2969_pInt, 2971_pInt, 2999_pInt, &
        3001_pInt, 3011_pInt, 3019_pInt, 3023_pInt, 3037_pInt, 3041_pInt, 3049_pInt, 3061_pInt, 3067_pInt, 3079_pInt, &
        3083_pInt, 3089_pInt, 3109_pInt, 3119_pInt, 3121_pInt, 3137_pInt, 3163_pInt, 3167_pInt, 3169_pInt, 3181_pInt, &
        3187_pInt, 3191_pInt, 3203_pInt, 3209_pInt, 3217_pInt, 3221_pInt, 3229_pInt, 3251_pInt, 3253_pInt, 3257_pInt, &
        3259_pInt, 3271_pInt, 3299_pInt, 3301_pInt, 3307_pInt, 3313_pInt, 3319_pInt, 3323_pInt, 3329_pInt, 3331_pInt, &
        3343_pInt, 3347_pInt, 3359_pInt, 3361_pInt, 3371_pInt, 3373_pInt, 3389_pInt, 3391_pInt, 3407_pInt, 3413_pInt, &
        3433_pInt, 3449_pInt, 3457_pInt, 3461_pInt, 3463_pInt, 3467_pInt, 3469_pInt, 3491_pInt, 3499_pInt, 3511_pInt, &
        3517_pInt, 3527_pInt, 3529_pInt, 3533_pInt, 3539_pInt, 3541_pInt, 3547_pInt, 3557_pInt, 3559_pInt, 3571_pInt, &
   ! 501:600
        3581_pInt, 3583_pInt, 3593_pInt, 3607_pInt, 3613_pInt, 3617_pInt, 3623_pInt, 3631_pInt, 3637_pInt, 3643_pInt, &
        3659_pInt, 3671_pInt, 3673_pInt, 3677_pInt, 3691_pInt, 3697_pInt, 3701_pInt, 3709_pInt, 3719_pInt, 3727_pInt, &
        3733_pInt, 3739_pInt, 3761_pInt, 3767_pInt, 3769_pInt, 3779_pInt, 3793_pInt, 3797_pInt, 3803_pInt, 3821_pInt, &
        3823_pInt, 3833_pInt, 3847_pInt, 3851_pInt, 3853_pInt, 3863_pInt, 3877_pInt, 3881_pInt, 3889_pInt, 3907_pInt, &
        3911_pInt, 3917_pInt, 3919_pInt, 3923_pInt, 3929_pInt, 3931_pInt, 3943_pInt, 3947_pInt, 3967_pInt, 3989_pInt, &
        4001_pInt, 4003_pInt, 4007_pInt, 4013_pInt, 4019_pInt, 4021_pInt, 4027_pInt, 4049_pInt, 4051_pInt, 4057_pInt, &
        4073_pInt, 4079_pInt, 4091_pInt, 4093_pInt, 4099_pInt, 4111_pInt, 4127_pInt, 4129_pInt, 4133_pInt, 4139_pInt, &
        4153_pInt, 4157_pInt, 4159_pInt, 4177_pInt, 4201_pInt, 4211_pInt, 4217_pInt, 4219_pInt, 4229_pInt, 4231_pInt, &
        4241_pInt, 4243_pInt, 4253_pInt, 4259_pInt, 4261_pInt, 4271_pInt, 4273_pInt, 4283_pInt, 4289_pInt, 4297_pInt, &
        4327_pInt, 4337_pInt, 4339_pInt, 4349_pInt, 4357_pInt, 4363_pInt, 4373_pInt, 4391_pInt, 4397_pInt, 4409_pInt, &
   ! 601:700
        4421_pInt, 4423_pInt, 4441_pInt, 4447_pInt, 4451_pInt, 4457_pInt, 4463_pInt, 4481_pInt, 4483_pInt, 4493_pInt, &
        4507_pInt, 4513_pInt, 4517_pInt, 4519_pInt, 4523_pInt, 4547_pInt, 4549_pInt, 4561_pInt, 4567_pInt, 4583_pInt, &
        4591_pInt, 4597_pInt, 4603_pInt, 4621_pInt, 4637_pInt, 4639_pInt, 4643_pInt, 4649_pInt, 4651_pInt, 4657_pInt, &
        4663_pInt, 4673_pInt, 4679_pInt, 4691_pInt, 4703_pInt, 4721_pInt, 4723_pInt, 4729_pInt, 4733_pInt, 4751_pInt, &
        4759_pInt, 4783_pInt, 4787_pInt, 4789_pInt, 4793_pInt, 4799_pInt, 4801_pInt, 4813_pInt, 4817_pInt, 4831_pInt, &
        4861_pInt, 4871_pInt, 4877_pInt, 4889_pInt, 4903_pInt, 4909_pInt, 4919_pInt, 4931_pInt, 4933_pInt, 4937_pInt, &
        4943_pInt, 4951_pInt, 4957_pInt, 4967_pInt, 4969_pInt, 4973_pInt, 4987_pInt, 4993_pInt, 4999_pInt, 5003_pInt, &
        5009_pInt, 5011_pInt, 5021_pInt, 5023_pInt, 5039_pInt, 5051_pInt, 5059_pInt, 5077_pInt, 5081_pInt, 5087_pInt, &
        5099_pInt, 5101_pInt, 5107_pInt, 5113_pInt, 5119_pInt, 5147_pInt, 5153_pInt, 5167_pInt, 5171_pInt, 5179_pInt, &
        5189_pInt, 5197_pInt, 5209_pInt, 5227_pInt, 5231_pInt, 5233_pInt, 5237_pInt, 5261_pInt, 5273_pInt, 5279_pInt, &
   ! 701:800
        5281_pInt, 5297_pInt, 5303_pInt, 5309_pInt, 5323_pInt, 5333_pInt, 5347_pInt, 5351_pInt, 5381_pInt, 5387_pInt, &
        5393_pInt, 5399_pInt, 5407_pInt, 5413_pInt, 5417_pInt, 5419_pInt, 5431_pInt, 5437_pInt, 5441_pInt, 5443_pInt, &
        5449_pInt, 5471_pInt, 5477_pInt, 5479_pInt, 5483_pInt, 5501_pInt, 5503_pInt, 5507_pInt, 5519_pInt, 5521_pInt, &
        5527_pInt, 5531_pInt, 5557_pInt, 5563_pInt, 5569_pInt, 5573_pInt, 5581_pInt, 5591_pInt, 5623_pInt, 5639_pInt, &
        5641_pInt, 5647_pInt, 5651_pInt, 5653_pInt, 5657_pInt, 5659_pInt, 5669_pInt, 5683_pInt, 5689_pInt, 5693_pInt, &
        5701_pInt, 5711_pInt, 5717_pInt, 5737_pInt, 5741_pInt, 5743_pInt, 5749_pInt, 5779_pInt, 5783_pInt, 5791_pInt, &
        5801_pInt, 5807_pInt, 5813_pInt, 5821_pInt, 5827_pInt, 5839_pInt, 5843_pInt, 5849_pInt, 5851_pInt, 5857_pInt, &
        5861_pInt, 5867_pInt, 5869_pInt, 5879_pInt, 5881_pInt, 5897_pInt, 5903_pInt, 5923_pInt, 5927_pInt, 5939_pInt, &
        5953_pInt, 5981_pInt, 5987_pInt, 6007_pInt, 6011_pInt, 6029_pInt, 6037_pInt, 6043_pInt, 6047_pInt, 6053_pInt, &
        6067_pInt, 6073_pInt, 6079_pInt, 6089_pInt, 6091_pInt, 6101_pInt, 6113_pInt, 6121_pInt, 6131_pInt, 6133_pInt, &
   ! 801:900
        6143_pInt, 6151_pInt, 6163_pInt, 6173_pInt, 6197_pInt, 6199_pInt, 6203_pInt, 6211_pInt, 6217_pInt, 6221_pInt, &
        6229_pInt, 6247_pInt, 6257_pInt, 6263_pInt, 6269_pInt, 6271_pInt, 6277_pInt, 6287_pInt, 6299_pInt, 6301_pInt, &
        6311_pInt, 6317_pInt, 6323_pInt, 6329_pInt, 6337_pInt, 6343_pInt, 6353_pInt, 6359_pInt, 6361_pInt, 6367_pInt, &
        6373_pInt, 6379_pInt, 6389_pInt, 6397_pInt, 6421_pInt, 6427_pInt, 6449_pInt, 6451_pInt, 6469_pInt, 6473_pInt, &
        6481_pInt, 6491_pInt, 6521_pInt, 6529_pInt, 6547_pInt, 6551_pInt, 6553_pInt, 6563_pInt, 6569_pInt, 6571_pInt, &
        6577_pInt, 6581_pInt, 6599_pInt, 6607_pInt, 6619_pInt, 6637_pInt, 6653_pInt, 6659_pInt, 6661_pInt, 6673_pInt, &
        6679_pInt, 6689_pInt, 6691_pInt, 6701_pInt, 6703_pInt, 6709_pInt, 6719_pInt, 6733_pInt, 6737_pInt, 6761_pInt, &
        6763_pInt, 6779_pInt, 6781_pInt, 6791_pInt, 6793_pInt, 6803_pInt, 6823_pInt, 6827_pInt, 6829_pInt, 6833_pInt, &
        6841_pInt, 6857_pInt, 6863_pInt, 6869_pInt, 6871_pInt, 6883_pInt, 6899_pInt, 6907_pInt, 6911_pInt, 6917_pInt, &
        6947_pInt, 6949_pInt, 6959_pInt, 6961_pInt, 6967_pInt, 6971_pInt, 6977_pInt, 6983_pInt, 6991_pInt, 6997_pInt, &
   ! 901:1000
        7001_pInt, 7013_pInt, 7019_pInt, 7027_pInt, 7039_pInt, 7043_pInt, 7057_pInt, 7069_pInt, 7079_pInt, 7103_pInt, &
        7109_pInt, 7121_pInt, 7127_pInt, 7129_pInt, 7151_pInt, 7159_pInt, 7177_pInt, 7187_pInt, 7193_pInt, 7207_pInt, &
        7211_pInt, 7213_pInt, 7219_pInt, 7229_pInt, 7237_pInt, 7243_pInt, 7247_pInt, 7253_pInt, 7283_pInt, 7297_pInt, &
        7307_pInt, 7309_pInt, 7321_pInt, 7331_pInt, 7333_pInt, 7349_pInt, 7351_pInt, 7369_pInt, 7393_pInt, 7411_pInt, &
        7417_pInt, 7433_pInt, 7451_pInt, 7457_pInt, 7459_pInt, 7477_pInt, 7481_pInt, 7487_pInt, 7489_pInt, 7499_pInt, &
        7507_pInt, 7517_pInt, 7523_pInt, 7529_pInt, 7537_pInt, 7541_pInt, 7547_pInt, 7549_pInt, 7559_pInt, 7561_pInt, &
        7573_pInt, 7577_pInt, 7583_pInt, 7589_pInt, 7591_pInt, 7603_pInt, 7607_pInt, 7621_pInt, 7639_pInt, 7643_pInt, &
        7649_pInt, 7669_pInt, 7673_pInt, 7681_pInt, 7687_pInt, 7691_pInt, 7699_pInt, 7703_pInt, 7717_pInt, 7723_pInt, &
        7727_pInt, 7741_pInt, 7753_pInt, 7757_pInt, 7759_pInt, 7789_pInt, 7793_pInt, 7817_pInt, 7823_pInt, 7829_pInt, &
        7841_pInt, 7853_pInt, 7867_pInt, 7873_pInt, 7877_pInt, 7879_pInt, 7883_pInt, 7901_pInt, 7907_pInt, 7919_pInt, &
   ! 1001:1100
        7927_pInt, 7933_pInt, 7937_pInt, 7949_pInt, 7951_pInt, 7963_pInt, 7993_pInt, 8009_pInt, 8011_pInt, 8017_pInt, &
        8039_pInt, 8053_pInt, 8059_pInt, 8069_pInt, 8081_pInt, 8087_pInt, 8089_pInt, 8093_pInt, 8101_pInt, 8111_pInt, &
        8117_pInt, 8123_pInt, 8147_pInt, 8161_pInt, 8167_pInt, 8171_pInt, 8179_pInt, 8191_pInt, 8209_pInt, 8219_pInt, &
        8221_pInt, 8231_pInt, 8233_pInt, 8237_pInt, 8243_pInt, 8263_pInt, 8269_pInt, 8273_pInt, 8287_pInt, 8291_pInt, &
        8293_pInt, 8297_pInt, 8311_pInt, 8317_pInt, 8329_pInt, 8353_pInt, 8363_pInt, 8369_pInt, 8377_pInt, 8387_pInt, &
        8389_pInt, 8419_pInt, 8423_pInt, 8429_pInt, 8431_pInt, 8443_pInt, 8447_pInt, 8461_pInt, 8467_pInt, 8501_pInt, &
        8513_pInt, 8521_pInt, 8527_pInt, 8537_pInt, 8539_pInt, 8543_pInt, 8563_pInt, 8573_pInt, 8581_pInt, 8597_pInt, &
        8599_pInt, 8609_pInt, 8623_pInt, 8627_pInt, 8629_pInt, 8641_pInt, 8647_pInt, 8663_pInt, 8669_pInt, 8677_pInt, &
        8681_pInt, 8689_pInt, 8693_pInt, 8699_pInt, 8707_pInt, 8713_pInt, 8719_pInt, 8731_pInt, 8737_pInt, 8741_pInt, &
        8747_pInt, 8753_pInt, 8761_pInt, 8779_pInt, 8783_pInt, 8803_pInt, 8807_pInt, 8819_pInt, 8821_pInt, 8831_pInt, &
   ! 1101:1200
        8837_pInt, 8839_pInt, 8849_pInt, 8861_pInt, 8863_pInt, 8867_pInt, 8887_pInt, 8893_pInt, 8923_pInt, 8929_pInt, &
        8933_pInt, 8941_pInt, 8951_pInt, 8963_pInt, 8969_pInt, 8971_pInt, 8999_pInt, 9001_pInt, 9007_pInt, 9011_pInt, &
        9013_pInt, 9029_pInt, 9041_pInt, 9043_pInt, 9049_pInt, 9059_pInt, 9067_pInt, 9091_pInt, 9103_pInt, 9109_pInt, &
        9127_pInt, 9133_pInt, 9137_pInt, 9151_pInt, 9157_pInt, 9161_pInt, 9173_pInt, 9181_pInt, 9187_pInt, 9199_pInt, &
        9203_pInt, 9209_pInt, 9221_pInt, 9227_pInt, 9239_pInt, 9241_pInt, 9257_pInt, 9277_pInt, 9281_pInt, 9283_pInt, &
        9293_pInt, 9311_pInt, 9319_pInt, 9323_pInt, 9337_pInt, 9341_pInt, 9343_pInt, 9349_pInt, 9371_pInt, 9377_pInt, &
        9391_pInt, 9397_pInt, 9403_pInt, 9413_pInt, 9419_pInt, 9421_pInt, 9431_pInt, 9433_pInt, 9437_pInt, 9439_pInt, &
        9461_pInt, 9463_pInt, 9467_pInt, 9473_pInt, 9479_pInt, 9491_pInt, 9497_pInt, 9511_pInt, 9521_pInt, 9533_pInt, &
        9539_pInt, 9547_pInt, 9551_pInt, 9587_pInt, 9601_pInt, 9613_pInt, 9619_pInt, 9623_pInt, 9629_pInt, 9631_pInt, &
        9643_pInt, 9649_pInt, 9661_pInt, 9677_pInt, 9679_pInt, 9689_pInt, 9697_pInt, 9719_pInt, 9721_pInt, 9733_pInt, &
   ! 1201:1300
         9739_pInt, 9743_pInt, 9749_pInt, 9767_pInt, 9769_pInt, 9781_pInt, 9787_pInt, 9791_pInt, 9803_pInt, 9811_pInt, &
         9817_pInt, 9829_pInt, 9833_pInt, 9839_pInt, 9851_pInt, 9857_pInt, 9859_pInt, 9871_pInt, 9883_pInt, 9887_pInt, &
         9901_pInt, 9907_pInt, 9923_pInt, 9929_pInt, 9931_pInt, 9941_pInt, 9949_pInt, 9967_pInt, 9973_pInt,10007_pInt, &
        10009_pInt,10037_pInt,10039_pInt,10061_pInt,10067_pInt,10069_pInt,10079_pInt,10091_pInt,10093_pInt,10099_pInt, &
        10103_pInt,10111_pInt,10133_pInt,10139_pInt,10141_pInt,10151_pInt,10159_pInt,10163_pInt,10169_pInt,10177_pInt, &
        10181_pInt,10193_pInt,10211_pInt,10223_pInt,10243_pInt,10247_pInt,10253_pInt,10259_pInt,10267_pInt,10271_pInt, &
        10273_pInt,10289_pInt,10301_pInt,10303_pInt,10313_pInt,10321_pInt,10331_pInt,10333_pInt,10337_pInt,10343_pInt, &
        10357_pInt,10369_pInt,10391_pInt,10399_pInt,10427_pInt,10429_pInt,10433_pInt,10453_pInt,10457_pInt,10459_pInt, &
        10463_pInt,10477_pInt,10487_pInt,10499_pInt,10501_pInt,10513_pInt,10529_pInt,10531_pInt,10559_pInt,10567_pInt, &
        10589_pInt,10597_pInt,10601_pInt,10607_pInt,10613_pInt,10627_pInt,10631_pInt,10639_pInt,10651_pInt,10657_pInt, &
   ! 1301:1400
        10663_pInt,10667_pInt,10687_pInt,10691_pInt,10709_pInt,10711_pInt,10723_pInt,10729_pInt,10733_pInt,10739_pInt, &
        10753_pInt,10771_pInt,10781_pInt,10789_pInt,10799_pInt,10831_pInt,10837_pInt,10847_pInt,10853_pInt,10859_pInt, &
        10861_pInt,10867_pInt,10883_pInt,10889_pInt,10891_pInt,10903_pInt,10909_pInt,19037_pInt,10939_pInt,10949_pInt, &
        10957_pInt,10973_pInt,10979_pInt,10987_pInt,10993_pInt,11003_pInt,11027_pInt,11047_pInt,11057_pInt,11059_pInt, &
        11069_pInt,11071_pInt,11083_pInt,11087_pInt,11093_pInt,11113_pInt,11117_pInt,11119_pInt,11131_pInt,11149_pInt, &
        11159_pInt,11161_pInt,11171_pInt,11173_pInt,11177_pInt,11197_pInt,11213_pInt,11239_pInt,11243_pInt,11251_pInt, &
        11257_pInt,11261_pInt,11273_pInt,11279_pInt,11287_pInt,11299_pInt,11311_pInt,11317_pInt,11321_pInt,11329_pInt, &
        11351_pInt,11353_pInt,11369_pInt,11383_pInt,11393_pInt,11399_pInt,11411_pInt,11423_pInt,11437_pInt,11443_pInt, &
        11447_pInt,11467_pInt,11471_pInt,11483_pInt,11489_pInt,11491_pInt,11497_pInt,11503_pInt,11519_pInt,11527_pInt, &
        11549_pInt,11551_pInt,11579_pInt,11587_pInt,11593_pInt,11597_pInt,11617_pInt,11621_pInt,11633_pInt,11657_pInt, &
   ! 1401:1500
        11677_pInt,11681_pInt,11689_pInt,11699_pInt,11701_pInt,11717_pInt,11719_pInt,11731_pInt,11743_pInt,11777_pInt, &
        11779_pInt,11783_pInt,11789_pInt,11801_pInt,11807_pInt,11813_pInt,11821_pInt,11827_pInt,11831_pInt,11833_pInt, &
        11839_pInt,11863_pInt,11867_pInt,11887_pInt,11897_pInt,11903_pInt,11909_pInt,11923_pInt,11927_pInt,11933_pInt, &
        11939_pInt,11941_pInt,11953_pInt,11959_pInt,11969_pInt,11971_pInt,11981_pInt,11987_pInt,12007_pInt,12011_pInt, &
        12037_pInt,12041_pInt,12043_pInt,12049_pInt,12071_pInt,12073_pInt,12097_pInt,12101_pInt,12107_pInt,12109_pInt, &
        12113_pInt,12119_pInt,12143_pInt,12149_pInt,12157_pInt,12161_pInt,12163_pInt,12197_pInt,12203_pInt,12211_pInt, &
        12227_pInt,12239_pInt,12241_pInt,12251_pInt,12253_pInt,12263_pInt,12269_pInt,12277_pInt,12281_pInt,12289_pInt, &
        12301_pInt,12323_pInt,12329_pInt,12343_pInt,12347_pInt,12373_pInt,12377_pInt,12379_pInt,12391_pInt,12401_pInt, &
        12409_pInt,12413_pInt,12421_pInt,12433_pInt,12437_pInt,12451_pInt,12457_pInt,12473_pInt,12479_pInt,12487_pInt, &
        12491_pInt,12497_pInt,12503_pInt,12511_pInt,12517_pInt,12527_pInt,12539_pInt,12541_pInt,12547_pInt,12553_pInt]
 endif

 if(n < 0_pInt) then
   prime = PRIME_MAX
 else if (n == 0_pInt) then
   prime = 1_pInt
 else if (n <= PRIME_MAX) then
   prime = npvec(n)
 else
   prime = -1_pInt
   call IO_error(error_ID=406_pInt)
 end if

end function prime


!--------------------------------------------------------------------------------------------------
!> @brief factorial
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function math_factorial(n)

 implicit none
 integer(pInt), intent(in) :: n
 integer(pInt) :: i
 
 math_factorial = product([(i, i=1,n)])

end function math_factorial


!--------------------------------------------------------------------------------------------------
!> @brief binomial coefficient
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function math_binomial(n,k)

 implicit none
 integer(pInt), intent(in) :: n, k
 integer(pInt) :: i, j
 
 j = min(k,n-k)
 math_binomial = product([(i, i=n, n-j+1, -1)])/math_factorial(j)

end function math_binomial


!--------------------------------------------------------------------------------------------------
!> @brief multinomial coefficient
!--------------------------------------------------------------------------------------------------
integer(pInt) pure function math_multinomial(alpha)

 implicit none
 integer(pInt), intent(in), dimension(:) :: alpha
 integer(pInt) :: i
 
 math_multinomial = 1_pInt
 do i = 1, size(alpha)
   math_multinomial = math_multinomial*math_binomial(sum(alpha(1:i)),alpha(i))
 enddo  

end function math_multinomial


!--------------------------------------------------------------------------------------------------
!> @brief volume of tetrahedron given by four vertices
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_volTetrahedron(v1,v2,v3,v4)

 implicit none
 real(pReal), dimension (3), intent(in) :: v1,v2,v3,v4
 real(pReal), dimension (3,3) :: m

 m(1:3,1) = v1-v2
 m(1:3,2) = v2-v3
 m(1:3,3) = v3-v4

 math_volTetrahedron = math_det33(m)/6.0_pReal

end function math_volTetrahedron


!--------------------------------------------------------------------------------------------------
!> @brief area of triangle given by three vertices
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_areaTriangle(v1,v2,v3)

 implicit none
 real(pReal), dimension (3), intent(in) :: v1,v2,v3

 math_areaTriangle = 0.5_pReal * norm2(math_crossproduct(v1-v2,v1-v3))

end function math_areaTriangle


!--------------------------------------------------------------------------------------------------
!> @brief rotate 33 tensor forward
!--------------------------------------------------------------------------------------------------
pure function math_rotate_forward33(tensor,rot_tensor)

 implicit none

 real(pReal), dimension(3,3) ::  math_rotate_forward33
 real(pReal), dimension(3,3), intent(in) :: tensor, rot_tensor

 math_rotate_forward33 = math_mul33x33(rot_tensor,&
                         math_mul33x33(tensor,math_transpose33(rot_tensor)))

end function math_rotate_forward33


!--------------------------------------------------------------------------------------------------
!> @brief rotate 33 tensor backward
!--------------------------------------------------------------------------------------------------
pure function math_rotate_backward33(tensor,rot_tensor)

 implicit none
 real(pReal), dimension(3,3) ::  math_rotate_backward33
 real(pReal), dimension(3,3), intent(in) :: tensor, rot_tensor

 math_rotate_backward33 = math_mul33x33(math_transpose33(rot_tensor),&
                           math_mul33x33(tensor,rot_tensor))

end function math_rotate_backward33


!--------------------------------------------------------------------------------------------------
!> @brief rotate 3333 tensor C'_ijkl=g_im*g_jn*g_ko*g_lp*C_mnop
!--------------------------------------------------------------------------------------------------
pure function math_rotate_forward3333(tensor,rot_tensor)

 implicit none
 real(pReal), dimension(3,3,3,3) ::  math_rotate_forward3333
 real(pReal), dimension(3,3), intent(in) :: rot_tensor
 real(pReal), dimension(3,3,3,3), intent(in) :: tensor
 integer(pInt) :: i,j,k,l,m,n,o,p

 math_rotate_forward3333= 0.0_pReal

 do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt; do k = 1_pInt,3_pInt; do l = 1_pInt,3_pInt
   do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt; do p = 1_pInt,3_pInt
     math_rotate_forward3333(i,j,k,l) = math_rotate_forward3333(i,j,k,l) &
                                      + rot_tensor(m,i) * rot_tensor(n,j) &
                                      * rot_tensor(o,k) * rot_tensor(p,l) * tensor(m,n,o,p)
 enddo; enddo; enddo; enddo; enddo; enddo; enddo; enddo

end function math_rotate_forward3333


!--------------------------------------------------------------------------------------------------
!> @brief calculate average of tensor field
!--------------------------------------------------------------------------------------------------
function math_tensorAvg(field)
 
 implicit none
 real(pReal), dimension(3,3) :: math_tensorAvg
 real(pReal), intent(in), dimension(:,:,:,:,:) :: field
 real(pReal) :: wgt
 
 wgt = 1.0_pReal/real(size(field,3)*size(field,4)*size(field,5), pReal) 
 math_tensorAvg = sum(sum(sum(field,dim=5),dim=4),dim=3)*wgt

end function math_tensorAvg


!--------------------------------------------------------------------------------------------------
!> @brief limits a scalar value to a certain range (either one or two sided)
! Will return NaN if left > right
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_limit(a, left, right)
 use prec, only: &
   DAMASK_NaN

 implicit none
 real(pReal), intent(in) :: a
 real(pReal), intent(in), optional :: left, right

   
 math_limit = min ( &
                   max (merge(left, -huge(a), present(left)), a), &
                        merge(right, huge(a), present(right)) &
                  )

 if (present(left) .and. present(right)) math_limit = merge (DAMASK_NaN,math_limit, left>right)

end function math_limit
 
end module math
