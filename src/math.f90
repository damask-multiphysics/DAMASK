!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mathematical library, including random number generation and tensor representations
!--------------------------------------------------------------------------------------------------
module math
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 real(pReal),    parameter, public :: PI = 3.141592653589793_pReal                                  !< ratio of a circle's circumference to its diameter
 real(pReal),    parameter, public :: INDEG = 180.0_pReal/PI                                        !< conversion from radian into degree
 real(pReal),    parameter, public :: INRAD = PI/180.0_pReal                                        !< conversion from degree into radian
 complex(pReal), parameter, public :: TWOPIIMG = (0.0_pReal,2.0_pReal)*(PI,0.0_pReal)               !< Re(0.0), Im(2xPi)

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
   math_eigenvectorBasisSym33, &
   math_eigenvectorBasisSym33_log, &
   math_eigenvectorBasisSym, &
   math_eigenValuesVectorsSym33, &
   math_eigenValuesVectorsSym, &
   math_rotationalPart33, &
   math_invariantsSym33, &
   math_eigenvaluesSym33, &
   math_factorial, &
   math_binomial, &
   math_multinomial, &
   math_volTetrahedron, &
   math_areaTriangle, &
   math_rotate_forward33, &
   math_rotate_backward33, &
   math_rotate_forward3333, &
   math_limit
 private :: &
   halton, &
   halton_memory, &
   halton_ndim_set, &
   halton_seed_set

contains

!--------------------------------------------------------------------------------------------------
!> @brief initialization of random seed generator
!--------------------------------------------------------------------------------------------------
subroutine math_init

#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use numerics, only: randomSeed
 use IO,       only: IO_timeStamp

 implicit none
 integer(pInt) :: i
 real(pReal), dimension(4) ::   randTest
! the following variables are system dependend and shound NOT be pInt
 integer :: randSize                                                                                ! gfortran requires a variable length to compile
 integer, dimension(:), allocatable :: randInit                                                     ! if recalculations of former randomness (with given seed) is necessary
                                                                                                    ! comment the first random_seed call out, set randSize to 1, and use ifort
 write(6,'(/,a)')   ' <<<+-  math init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 call random_seed(size=randSize)
 if (allocated(randInit)) deallocate(randInit)
 allocate(randInit(randSize))
 if (randomSeed > 0_pInt) then
   randInit(1:randSize) = int(randomSeed)                                                           ! randomSeed is of type pInt, randInit not
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

 write(6,'(a,I2)') ' size  of random seed:    ', randSize
 do i = 1_pInt,randSize
   write(6,'(a,I2,I14)') ' value of random seed:    ', i, randInit(i)
 enddo
 write(6,'(a,4(/,26x,f17.14),/)') ' start of random sequence: ', randTest

 call random_seed(put = randInit)

 call halton_seed_set(int(randInit(1), pInt))
 call halton_ndim_set(3_pInt)

 call math_check()

end subroutine math_init
 
!--------------------------------------------------------------------------------------------------
!> @brief check correctness of (some) math functions
!--------------------------------------------------------------------------------------------------
subroutine math_check

 use prec,     only: tol_math_check
 use IO,       only: IO_error

 implicit none
 character(len=64) :: error_msg

 real(pReal), dimension(3,3) :: R,R2
 real(pReal), dimension(3) ::   Eulers,v
 real(pReal), dimension(4) ::   q,q2,axisangle

 ! --- check rotation dictionary ---

 q = math_qRand()          ! random quaternion
 
 ! +++ q -> a -> q  +++
 axisangle = math_qToAxisAngle(q)
 q2 = math_axisAngleToQ(axisangle(1:3),axisangle(4))
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) &
          'quat -> axisAngle -> quat maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ q -> R -> q  +++
 R = math_qToR(q)
 q2 = math_RtoQ(R)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) &
          'quat -> R -> quat maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ q -> euler -> q  +++
 Eulers = math_qToEuler(q)
 q2 = math_EulerToQ(Eulers)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) &
          'quat -> euler -> quat maximum deviation ',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ R -> euler -> R  +++
 Eulers = math_RtoEuler(R)
 R2 = math_EulerToR(Eulers)
 if ( any(abs( R-R2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) &
          'R -> euler -> R maximum deviation ',maxval(abs( R-R2))
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ check rotation sense of q and R +++
 call halton(3_pInt,v)     ! random vector
 R = math_qToR(q)
 if (any(abs(math_mul33x3(R,v) - math_qRot(q,v)) > tol_math_check)) then
   write (error_msg, '(a)' ) 'R(q)*v has different sense than q*v'
   call IO_error(401_pInt,ext_msg=error_msg)
 endif

 ! +++ check vector expansion +++
 if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,3.0_pReal,3.0_pReal,3.0_pReal] - &
             math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1_pInt,2_pInt,3_pInt,0_pInt])) > tol_math_check)) then
   write (error_msg, '(a)' ) 'math_expand [1,2,3] by [1,2,3,0] => [1,2,2,3,3,3]'
   call IO_error(401_pInt,ext_msg=error_msg)
 endif
 if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal] - &
             math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1_pInt,2_pInt])) > tol_math_check)) then
   write (error_msg, '(a)' ) 'math_expand [1,2,3] by [1,2] => [1,2,2]'
   call IO_error(401_pInt,ext_msg=error_msg)
 endif
 if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,1.0_pReal,1.0_pReal,1.0_pReal] - &
             math_expand([1.0_pReal,2.0_pReal],[1_pInt,2_pInt,3_pInt])) > tol_math_check)) then
   write (error_msg, '(a)' ) 'math_expand [1,2] by [1,2,3] => [1,2,2,1,1,1]'
   call IO_error(401_pInt,ext_msg=error_msg)
 endif
 
end subroutine math_check
 

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
   ipivot = qsort_partition(a,istart, iend)
   call math_qsort(a, istart, ipivot-1_pInt)
   call math_qsort(a, ipivot+1_pInt, iend)
 endif

!--------------------------------------------------------------------------------------------------
 contains

 !-------------------------------------------------------------------------------------------------
 !> @brief Partitioning required for quicksort
 !-------------------------------------------------------------------------------------------------
 integer(pInt) function qsort_partition(a, istart, iend)
 
   implicit none
   integer(pInt), dimension(:,:), intent(inout) :: a
   integer(pInt), intent(in) :: istart,iend
   integer(pInt) :: i,j,k,tmp
  
   do
  ! find the first element on the right side less than or equal to the pivot point
     do j = iend, istart, -1_pInt
       if (a(1,j) <= a(1,istart)) exit
     enddo
  ! find the first element on the left side greater than the pivot point
     do i = istart, iend
       if (a(1,i) > a(1,istart)) exit
     enddo
     if (i < j) then ! if the indexes do not cross, exchange values
       do k = 1_pInt, int(size(a,1_pInt), pInt)
         tmp = a(k,i)
         a(k,i) = a(k,j)
         a(k,j) = tmp
       enddo
     else           ! if they do cross, exchange left value with pivot and return with the partition index
       do k = 1_pInt, int(size(a,1_pInt), pInt)
         tmp = a(k,istart)
         a(k,istart) = a(k,j)
         a(k,j) = tmp
       enddo
       qsort_partition = j
       return
     endif
   enddo
 
 end function qsort_partition

end subroutine math_qsort


!--------------------------------------------------------------------------------------------------
!> @brief vector expansion
!> @details takes a set of numbers (a,b,c,...) and corresponding multiples (x,y,z,...)
!> to return a vector of x times a, y times b, z times c, ... 
!--------------------------------------------------------------------------------------------------
pure function math_expand(what,how)

 implicit none
 real(pReal),   dimension(:), intent(in) :: what
 integer(pInt), dimension(:), intent(in) :: how
 real(pReal), dimension(sum(how)) ::  math_expand
 integer(pInt) :: i

 do i = 1_pInt, size(how)
   math_expand(sum(how(1:i-1))+1:sum(how(1:i))) = what(mod(i-1_pInt,size(what))+1_pInt)
 enddo

end function math_expand


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
!> @brief tensor product A \otimes B of arbitrary sized vectors A and B
!--------------------------------------------------------------------------------------------------
pure function math_tensorproduct(A,B)

 implicit none
 real(pReal), dimension(:), intent(in) ::  A,B
 real(pReal), dimension(size(A,1),size(B,1)) ::  math_tensorproduct
 integer(pInt) :: i,j

 forall (i=1_pInt:size(A,1),j=1_pInt:size(B,1)) math_tensorproduct(i,j) = A(i)*B(j)

end function math_tensorproduct


!--------------------------------------------------------------------------------------------------
!> @brief tensor product A \otimes B of leght-3 vectors A and B
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
 integer(pInt) :: i
 integer(pInt), intent(in), optional :: n
 real(pReal), dimension(3,3), intent(in) :: A
 real(pReal), dimension(3,3) :: B, math_exp33
 real(pReal) :: invFac

 B = math_I3                                                                                        ! init
 invFac = 1.0_pReal                                                                                 ! 0!
 math_exp33 = B                                                                                     ! A^0 = eye2

 do i = 1_pInt, merge(n,5_pInt,present(n))
   invFac = invFac/real(i,pReal)                                                                    ! invfac = 1/i!
   B = math_mul33x33(B,A)
   math_exp33 = math_exp33 + invFac*B                                                               ! exp = SUM (A^i)/i!
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
 use prec, only: &
   dNeq0

 implicit none
 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal) :: DetA
 real(pReal),dimension(3,3) :: math_inv33

 math_inv33(1,1) =  A(2,2) * A(3,3) - A(2,3) * A(3,2)
 math_inv33(2,1) = -A(2,1) * A(3,3) + A(2,3) * A(3,1)
 math_inv33(3,1) =  A(2,1) * A(3,2) - A(2,2) * A(3,1)

 DetA = A(1,1) * math_inv33(1,1) + A(1,2) * math_inv33(2,1) + A(1,3) * math_inv33(3,1)

 if (dNeq0(DetA)) then
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
 use prec, only: &
   dEq0

 implicit none
 logical, intent(out) :: error
 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal),dimension(3,3),intent(out) :: InvA
 real(pReal), intent(out) :: DetA

 InvA(1,1) =  A(2,2) * A(3,3) - A(2,3) * A(3,2)
 InvA(2,1) = -A(2,1) * A(3,3) + A(2,3) * A(3,1)
 InvA(3,1) =  A(2,1) * A(3,2) - A(2,2) * A(3,1)

 DetA = A(1,1) * InvA(1,1) + A(1,2) * InvA(2,1) + A(1,3) * InvA(3,1)

 if (dEq0(DetA)) then
   InvA = 0.0_pReal
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
 external :: &
  dgetrf, &
  dgetri

 temp66_real = math_Mandel3333to66(A)
 call dgetrf(6,6,temp66_real,6,ipiv6,ierr)
 call dgetri(6,temp66_real,6,ipiv6,work6,6,ierr)
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
 external :: &
  dgetrf, &
  dgetri
 
 invA = A 
 call dgetrf(myDim,myDim,invA,myDim,ipiv,ierr)
 call dgetri(myDim,InvA,myDim,ipiv,work,myDim,ierr)
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
!> @brief determinant of a symmetric 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_detSym33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m
 
  math_detSym33 = -(m(1,1)*m(2,3)**2_pInt + m(2,2)*m(1,3)**2_pInt + m(3,3)*m(1,2)**2_pInt) &
                  + m(1,1)*m(2,2)*m(3,3)  + 2.0_pReal * m(1,2)*m(1,3)*m(2,3)

end function  math_detSym33


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
! http://math.stackexchange.com/questions/131336/uniform-random-quaternion-in-a-restricted-angle-range
! K. Shoemake. Uniform random rotations. In D. Kirk, editor, Graphics Gems III, pages 124-132. 
! Academic, New York, 1992.
!--------------------------------------------------------------------------------------------------
function math_qRand()

 implicit none
 real(pReal), dimension(4) :: math_qRand
 real(pReal), dimension(3) :: rnd

 call halton(3_pInt,rnd)
 math_qRand = [cos(2.0_pReal*PI*rnd(1))*sqrt(rnd(3)), &
               sin(2.0_pReal*PI*rnd(2))*sqrt(1.0_pReal-rnd(3)), &
               cos(2.0_pReal*PI*rnd(2))*sqrt(1.0_pReal-rnd(3)), &
               sin(2.0_pReal*PI*rnd(1))*sqrt(rnd(3))]

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
 use prec, only: &
   dNeq0

 implicit none
 real(pReal), dimension(4), intent(in) ::  Q
 real(pReal), dimension(4) ::  math_qInv
 real(pReal) :: squareNorm

 math_qInv = 0.0_pReal

 squareNorm = math_qDot(Q,Q)
 if (dNeq0(squareNorm)) math_qInv = math_qConj(Q) / squareNorm

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
!> @brief rotation matrix from Bunge-Euler (3-1-3) angles (in radians)
!> @details rotation matrix is meant to represent a PASSIVE rotation, composed of INTRINSIC
!> @details rotations around the axes of the details rotating reference frame.
!> @details similar to eu2om from "D Rowenhorst et al. Consistent representations of and conversions
!> @details between 3D rotations, Model. Simul. Mater. Sci. Eng. 23-8 (2015)", but R is transposed
!--------------------------------------------------------------------------------------------------
pure function math_EulerToR(Euler)

 implicit none
 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_EulerToR
 real(pReal) :: c1, C, c2, s1, S, s2

 c1 = cos(Euler(1))
 C  = cos(Euler(2))
 c2 = cos(Euler(3))
 s1 = sin(Euler(1))
 S  = sin(Euler(2))
 s2 = sin(Euler(3))

 math_EulerToR(1,1) =  c1*c2 -s1*C*s2
 math_EulerToR(1,2) = -c1*s2 -s1*C*c2
 math_EulerToR(1,3) =  s1*S

 math_EulerToR(2,1) =  s1*c2 +c1*C*s2
 math_EulerToR(2,2) = -s1*s2 +c1*C*c2
 math_EulerToR(2,3) = -c1*S

 math_EulerToR(3,1) =  S*s2
 math_EulerToR(3,2) =  S*c2
 math_EulerToR(3,3) =  C
 
 math_EulerToR = transpose(math_EulerToR)  ! convert to passive rotation

end function math_EulerToR


!--------------------------------------------------------------------------------------------------
!> @brief quaternion (w+ix+jy+kz) from Bunge-Euler (3-1-3) angles (in radians)
!> @details rotation matrix is meant to represent a PASSIVE rotation, composed of INTRINSIC
!> @details rotations around the axes of the details rotating reference frame.
!> @details similar to eu2qu from "D Rowenhorst et al. Consistent representations of and
!> @details conversions between 3D rotations, Model. Simul. Mater. Sci. Eng. 23-8 (2015)", but
!> @details Q is conjucated and Q is not reversed for Q(0) < 0.
!--------------------------------------------------------------------------------------------------
pure function math_EulerToQ(eulerangles)

 implicit none
 real(pReal), dimension(3), intent(in) :: eulerangles
 real(pReal), dimension(4) :: math_EulerToQ
 real(pReal) :: c, s, sigma, delta

 c = cos(0.5_pReal * eulerangles(2))
 s = sin(0.5_pReal * eulerangles(2))
 sigma = 0.5_pReal * (eulerangles(1)+eulerangles(3))
 delta = 0.5_pReal * (eulerangles(1)-eulerangles(3))
 
 math_EulerToQ= [c * cos(sigma), &
                 s * cos(delta), &
                 s * sin(delta), &
                 c * sin(sigma) ]
 math_EulerToQ = math_qConj(math_EulerToQ)                                                          ! convert to passive rotation

end function math_EulerToQ


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from axis and angle (in radians)
!> @details rotation matrix is meant to represent a ACTIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from http://mathworld.wolfram.com/RodriguesRotationFormula.html
!> @details equivalent to eu2om (P=-1) from "D Rowenhorst et al. Consistent representations of and
!> @details conversions between 3D rotations, Model. Simul. Mater. Sci. Eng. 23-8 (2015)"
!--------------------------------------------------------------------------------------------------
pure function math_axisAngleToR(axis,omega)

 implicit none
 real(pReal), dimension(3,3) :: math_axisAngleToR
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: n
 real(pReal) :: norm,s,c,c1

 norm = norm2(axis)
 wellDefined: if (norm > 1.0e-8_pReal) then
   n = axis/norm                                           ! normalize axis to be sure

   s = sin(omega)
   c = cos(omega)
   c1 = 1.0_pReal - c

   math_axisAngleToR(1,1) =  c + c1*n(1)**2.0_pReal
   math_axisAngleToR(1,2) =  c1*n(1)*n(2) - s*n(3)
   math_axisAngleToR(1,3) =  c1*n(1)*n(3) + s*n(2)
                             
   math_axisAngleToR(2,1) =  c1*n(1)*n(2) + s*n(3)
   math_axisAngleToR(2,2) =  c + c1*n(2)**2.0_pReal
   math_axisAngleToR(2,3) =  c1*n(2)*n(3) - s*n(1)
                             
   math_axisAngleToR(3,1) =  c1*n(1)*n(3) - s*n(2)
   math_axisAngleToR(3,2) =  c1*n(2)*n(3) + s*n(1)
   math_axisAngleToR(3,3) =  c + c1*n(3)**2.0_pReal
 else wellDefined
   math_axisAngleToR = math_I3
 endif wellDefined
 
end function math_axisAngleToR


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from axis and angle (in radians)
!> @details rotation matrix is meant to represent a PASSIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details eq-uivalent to eu2qu (P=+1) from "D Rowenhorst et al. Consistent representations of and
!> @details conversions between 3D rotations, Model. Simul. Mater. Sci. Eng. 23-8 (2015)"
!--------------------------------------------------------------------------------------------------
pure function math_EulerAxisAngleToR(axis,omega)

 implicit none
 real(pReal), dimension(3,3) :: math_EulerAxisAngleToR
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega

 math_EulerAxisAngleToR = transpose(math_axisAngleToR(axis,omega))                                  ! convert to passive rotation

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

 math_EulerAxisAngleToQ = math_qConj(math_axisAngleToQ(axis,omega))                                 ! convert to passive rotation

end function math_EulerAxisAngleToQ


!--------------------------------------------------------------------------------------------------
!> @brief quaternion (w+ix+jy+kz) from axis and angle (in radians)
!> @details quaternion is meant to represent an ACTIVE rotation
!> @details (see http://en.wikipedia.org/wiki/Euler_angles for definitions)
!> @details formula for active rotation taken from
!> @details http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29#Rodrigues_parameters
!> @details equivalent to eu2qu (P=+1) from "D Rowenhorst et al. Consistent representations of and
!> @details conversions between 3D rotations, Model. Simul. Mater. Sci. Eng. 23-8 (2015)"
!--------------------------------------------------------------------------------------------------
pure function math_axisAngleToQ(axis,omega)

 implicit none
 real(pReal), dimension(4) :: math_axisAngleToQ
 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: axisNrm
 real(pReal) :: norm

 norm = norm2(axis)
 wellDefined: if (norm > 1.0e-8_pReal) then
   axisNrm = axis/norm                                                                              ! normalize axis to be sure
   math_axisAngleToQ = [cos(0.5_pReal*omega), sin(0.5_pReal*omega) * axisNrm(1:3)]
 else wellDefined
   math_axisAngleToQ = [1.0_pReal,0.0_pReal,0.0_pReal,0.0_pReal]
 endif wellDefined

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
                   -q(3),      q(2), 0.0_pReal],[3,3])                                              ! notation is transposed

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

 math_qToEuler(2) = acos(1.0_pReal-2.0_pReal*(q(2)**2+q(3)**2))

 if (abs(math_qToEuler(2)) < 1.0e-6_pReal) then
   math_qToEuler(1) = sign(2.0_pReal*acos(math_limit(q(1),-1.0_pReal, 1.0_pReal)),q(4))
   math_qToEuler(3) = 0.0_pReal
 else
   math_qToEuler(1) = atan2(+q(1)*q(3)+q(2)*q(4), q(1)*q(2)-q(3)*q(4))
   math_qToEuler(3) = atan2(-q(1)*q(3)+q(2)*q(4), q(1)*q(2)+q(3)*q(4))
 endif

 math_qToEuler = merge(math_qToEuler + [2.0_pReal*PI, PI, 2.0_pReal*PI], &                          ! ensure correct range
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

 halfAngle = acos(math_limit(Q(1),-1.0_pReal,1.0_pReal))
 sinHalfAngle = sin(halfAngle)

 smallRotation: if (sinHalfAngle <= 1.0e-4_pReal) then
   math_qToAxisAngle = 0.0_pReal
 else smallRotation
   math_qToAxisAngle= [ Q(2:4)/sinHalfAngle, halfAngle*2.0_pReal] 
 endif smallRotation

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
 use, intrinsic :: &
   IEEE_arithmetic
 use prec, only: &
   tol_math_check

 implicit none
 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3) :: math_qToRodrig

 math_qToRodrig = merge(Q(2:4)/Q(1),IEEE_value(1.0_pReal,IEEE_quiet_NaN),abs(Q(1)) > tol_math_check)! NaN for 180 deg since Rodrig is unbound

end function math_qToRodrig


!--------------------------------------------------------------------------------------------------
!> @brief misorientation angle between two sets of Euler angles
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_EulerMisorientation(EulerA,EulerB)

 implicit none
 real(pReal), dimension(3), intent(in) :: EulerA,EulerB
 real(pReal) :: cosTheta

 cosTheta = (math_trace33(math_mul33x33(math_EulerToR(EulerB), &
                              transpose(math_EulerToR(EulerA)))) - 1.0_pReal) * 0.5_pReal

 math_EulerMisorientation = acos(math_limit(cosTheta,-1.0_pReal,1.0_pReal))

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
 real(pReal), dimension(3), parameter :: ORIGIN = 0.0_pReal
 real(pReal), dimension(5) :: rnd

 noScatter: if (abs(noise) < tol_math_check) then
   math_sampleGaussOri = center
 else noScatter
  ! Helming uses different distribution with Bessel functions
  ! therefore the gauss scatter width has to be scaled differently
   scatter = 0.95_pReal * noise
   cosScatter = cos(scatter)

   do
     call halton(5_pInt,rnd)
     rnd(1:3) = 2.0_pReal*rnd(1:3)-1.0_pReal                                                        ! expand 1:3 to range [-1,+1]
     disturb  = [ scatter * rnd(1), &                                                               ! phi1
                  sign(1.0_pReal,rnd(2))*acos(cosScatter+(1.0_pReal-cosScatter)*rnd(4)), &          ! Phi
                  scatter * rnd(3)]                                                                 ! phi2
     if (rnd(5) <= exp(-1.0_pReal*(math_EulerMisorientation(ORIGIN,disturb)/scatter)**2_pReal)) exit
   enddo

   math_sampleGaussOri = math_RtoEuler(math_mul33x33(math_EulerToR(disturb),math_EulerToR(center)))
 endif noScatter

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
!> @brief symmetrically equivalent Euler angles for given sample symmetry 
!> @detail 1 (equivalent to != 2 and !=4):triclinic, 2:monoclinic, 4:orthotropic
!--------------------------------------------------------------------------------------------------
pure function math_symmetricEulers(sym,Euler)

 implicit none
 integer(pInt), intent(in) :: sym                                                                   !< symmetry Class
 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_symmetricEulers

 math_symmetricEulers = transpose(reshape([PI+Euler(1), PI-Euler(1), 2.0_pReal*PI-Euler(1), &
                                              Euler(2), PI-Euler(2), PI          -Euler(2), &
                                              Euler(3), PI+Euler(3), PI          +Euler(3)],[3,3])) ! transpose is needed to have symbolic notation instead of column-major

 math_symmetricEulers = modulo(math_symmetricEulers,2.0_pReal*pi)

 select case (sym)
   case (4_pInt)                                                                                    ! orthotropic: all done

   case (2_pInt)                                                                                    ! monoclinic: return only first
     math_symmetricEulers(1:3,2:3) = 0.0_pReal

   case default                                                                                     ! triclinic: return blank
     math_symmetricEulers = 0.0_pReal
 end select

end function math_symmetricEulers


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric matrix m
!--------------------------------------------------------------------------------------------------
subroutine math_eigenValuesVectorsSym(m,values,vectors,error)

 implicit none
 real(pReal), dimension(:,:),                  intent(in)  :: m
 real(pReal), dimension(size(m,1)),            intent(out) :: values
 real(pReal), dimension(size(m,1),size(m,1)),  intent(out) :: vectors
 logical, intent(out) :: error
 integer(pInt) :: info
 real(pReal), dimension((64+2)*size(m,1)) :: work                                                    ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f
 external :: &
  dsyev

 vectors = m                                                                                         ! copy matrix to input (doubles as output) array
 call dsyev('V','U',size(m,1),vectors,size(m,1),values,work,(64+2)*size(m,1),info)
 error = (info == 0_pInt)

end subroutine math_eigenValuesVectorsSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric 33 matrix m using an analytical expression
!> and the general LAPACK powered version for arbritrary sized matrices as fallback
!> @author Joachim Kopp, Max–Planck–Institut für Kernphysik, Heidelberg (Copyright (C) 2006)
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details See http://arxiv.org/abs/physics/0610206 (DSYEVH3)
!--------------------------------------------------------------------------------------------------
subroutine math_eigenValuesVectorsSym33(m,values,vectors)
 
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
 U = max(T, T**2_pInt)
 threshold = sqrt(5.68e-14_pReal * U**2_pInt)

! Calculate first eigenvector by the formula v[0] = (m - lambda[0]).e1 x (m - lambda[0]).e2
 vectors(1:3,1) = [ vectors(1,2) + m(1, 3) * values(1), &
                    vectors(2,2) + m(2, 3) * values(1), &
                    (m(1,1) - values(1)) * (m(2,2) - values(1)) - vectors(3,2)]
 norm = norm2(vectors(1:3, 1))

 fallback1: if(norm < threshold) then
   call math_eigenValuesVectorsSym(m,values,vectors,error)
   return
 endif fallback1

 vectors(1:3,1) = vectors(1:3, 1) / norm 
 
! Calculate second eigenvector by the formula v[1] = (m - lambda[1]).e1 x (m - lambda[1]).e2
 vectors(1:3,2) = [ vectors(1,2) + m(1, 3) * values(2), &
                    vectors(2,2) + m(2, 3) * values(2), &
                    (m(1,1) - values(2)) * (m(2,2) - values(2)) - vectors(3,2)]
 norm = norm2(vectors(1:3, 2))

 fallback2: if(norm < threshold) then
   call math_eigenValuesVectorsSym(m,values,vectors,error)
   return
 endif fallback2
 vectors(1:3,2) = vectors(1:3, 2) / norm

! Calculate third eigenvector according to  v[2] = v[0] x v[1]
 vectors(1:3,3) = math_crossproduct(vectors(1:3,1),vectors(1:3,2))

end subroutine math_eigenValuesVectorsSym33


!--------------------------------------------------------------------------------------------------
!> @brief eigenvector basis of symmetric matrix m
!--------------------------------------------------------------------------------------------------
function math_eigenvectorBasisSym(m)

 implicit none
 real(pReal), dimension(:,:),      intent(in)  :: m
 real(pReal), dimension(size(m,1))             :: values
 real(pReal), dimension(size(m,1),size(m,1))   :: vectors
 real(pReal), dimension(size(m,1),size(m,1))   :: math_eigenvectorBasisSym
 logical :: error
 integer(pInt) :: i

 math_eigenvectorBasisSym = 0.0_pReal
 call math_eigenValuesVectorsSym(m,values,vectors,error)
 if(error) return
 
 do i=1_pInt, size(m,1)
   math_eigenvectorBasisSym = math_eigenvectorBasisSym &
                            + sqrt(values(i)) * math_tensorproduct(vectors(:,i),vectors(:,i))
 enddo

end function math_eigenvectorBasisSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvector basis of symmetric 33 matrix m
!--------------------------------------------------------------------------------------------------
function math_eigenvectorBasisSym33(m)

 implicit none
 real(pReal), dimension(3,3)              :: math_eigenvectorBasisSym33
 real(pReal), dimension(3)                :: invariants, values
 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: P, Q, rho, phi
 real(pReal), parameter :: TOL=1.e-14_pReal
 real(pReal), dimension(3,3,3) :: N, EB

 invariants = math_invariantsSym33(m)
 EB = 0.0_pReal

 P = invariants(2)-invariants(1)**2.0_pReal/3.0_pReal
 Q = -2.0_pReal/27.0_pReal*invariants(1)**3.0_pReal+product(invariants(1:2))/3.0_pReal-invariants(3)

 threeSimilarEigenvalues: if(all(abs([P,Q]) < TOL)) then
   values = invariants(1)/3.0_pReal
!   this is not really correct, but at least the basis is correct
   EB(1,1,1)=1.0_pReal
   EB(2,2,2)=1.0_pReal
   EB(3,3,3)=1.0_pReal
 else threeSimilarEigenvalues
   rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   phi=acos(math_limit(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
   values = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                             [cos(phi/3.0_pReal), &
                              cos((phi+2.0_pReal*PI)/3.0_pReal), &
                              cos((phi+4.0_pReal*PI)/3.0_pReal) &
                             ] + invariants(1)/3.0_pReal
   N(1:3,1:3,1) = m-values(1)*math_I3
   N(1:3,1:3,2) = m-values(2)*math_I3
   N(1:3,1:3,3) = m-values(3)*math_I3
   twoSimilarEigenvalues: if(abs(values(1)-values(2)) < TOL) then 
     EB(1:3,1:3,3)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                               ((values(3)-values(1))*(values(3)-values(2)))
     EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,3)
   elseif(abs(values(2)-values(3)) < TOL) then twoSimilarEigenvalues
     EB(1:3,1:3,1)=math_mul33x33(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                               ((values(1)-values(2))*(values(1)-values(3)))
     EB(1:3,1:3,2)=math_I3-EB(1:3,1:3,1)
   elseif(abs(values(3)-values(1)) < TOL) then twoSimilarEigenvalues 
     EB(1:3,1:3,2)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                               ((values(2)-values(1))*(values(2)-values(3)))
     EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,2)
   else twoSimilarEigenvalues
     EB(1:3,1:3,1)=math_mul33x33(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                               ((values(1)-values(2))*(values(1)-values(3)))
     EB(1:3,1:3,2)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                               ((values(2)-values(1))*(values(2)-values(3)))
     EB(1:3,1:3,3)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                               ((values(3)-values(1))*(values(3)-values(2)))
   endif twoSimilarEigenvalues
 endif threeSimilarEigenvalues

 math_eigenvectorBasisSym33 = sqrt(values(1)) * EB(1:3,1:3,1) &
                            + sqrt(values(2)) * EB(1:3,1:3,2) &
                            + sqrt(values(3)) * EB(1:3,1:3,3)

end function math_eigenvectorBasisSym33


!--------------------------------------------------------------------------------------------------
!> @brief logarithm eigenvector basis of symmetric 33 matrix m
!--------------------------------------------------------------------------------------------------
function math_eigenvectorBasisSym33_log(m)

 implicit none
 real(pReal), dimension(3,3)              :: math_eigenvectorBasisSym33_log
 real(pReal), dimension(3)                :: invariants, values
 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: P, Q, rho, phi
 real(pReal), parameter :: TOL=1.e-14_pReal
 real(pReal), dimension(3,3,3) :: N, EB

 invariants = math_invariantsSym33(m)
 EB = 0.0_pReal

 P = invariants(2)-invariants(1)**2.0_pReal/3.0_pReal
 Q = -2.0_pReal/27.0_pReal*invariants(1)**3.0_pReal+product(invariants(1:2))/3.0_pReal-invariants(3)

 threeSimilarEigenvalues: if(all(abs([P,Q]) < TOL)) then
   values = invariants(1)/3.0_pReal
!   this is not really correct, but at least the basis is correct
   EB(1,1,1)=1.0_pReal
   EB(2,2,2)=1.0_pReal
   EB(3,3,3)=1.0_pReal
 else threeSimilarEigenvalues
   rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   phi=acos(math_limit(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
   values = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                             [cos(phi/3.0_pReal), &
                              cos((phi+2.0_pReal*PI)/3.0_pReal), &
                              cos((phi+4.0_pReal*PI)/3.0_pReal) &
                             ] + invariants(1)/3.0_pReal
   N(1:3,1:3,1) = m-values(1)*math_I3
   N(1:3,1:3,2) = m-values(2)*math_I3
   N(1:3,1:3,3) = m-values(3)*math_I3
   twoSimilarEigenvalues: if(abs(values(1)-values(2)) < TOL) then 
     EB(1:3,1:3,3)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                               ((values(3)-values(1))*(values(3)-values(2)))
     EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,3)
   elseif(abs(values(2)-values(3)) < TOL) then twoSimilarEigenvalues
     EB(1:3,1:3,1)=math_mul33x33(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                               ((values(1)-values(2))*(values(1)-values(3)))
     EB(1:3,1:3,2)=math_I3-EB(1:3,1:3,1)
   elseif(abs(values(3)-values(1)) < TOL) then twoSimilarEigenvalues 
     EB(1:3,1:3,2)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                               ((values(2)-values(1))*(values(2)-values(3)))
     EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,2)
   else twoSimilarEigenvalues
     EB(1:3,1:3,1)=math_mul33x33(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                               ((values(1)-values(2))*(values(1)-values(3)))
     EB(1:3,1:3,2)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                               ((values(2)-values(1))*(values(2)-values(3)))
     EB(1:3,1:3,3)=math_mul33x33(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                               ((values(3)-values(1))*(values(3)-values(2)))
   endif twoSimilarEigenvalues
 endif threeSimilarEigenvalues

 math_eigenvectorBasisSym33_log = log(sqrt(values(1))) * EB(1:3,1:3,1) &
                            + log(sqrt(values(2))) * EB(1:3,1:3,2) &
                            + log(sqrt(values(3))) * EB(1:3,1:3,3)

end function math_eigenvectorBasisSym33_log

!--------------------------------------------------------------------------------------------------
!> @brief rotational part from polar decomposition of 33 tensor m
!--------------------------------------------------------------------------------------------------
function math_rotationalPart33(m)
 use prec, only: &
   dEq0
 use IO, only: &
   IO_warning

 implicit none
 real(pReal), intent(in), dimension(3,3) :: m
 real(pReal), dimension(3,3) :: math_rotationalPart33
 real(pReal), dimension(3,3) :: U , Uinv

 U = math_eigenvectorBasisSym33(math_mul33x33(transpose(m),m))
 Uinv = math_inv33(U)

 inversionFailed: if (all(dEq0(Uinv))) then
   math_rotationalPart33 = math_I3
   call IO_warning(650_pInt)
 else inversionFailed
   math_rotationalPart33 = math_mul33x33(m,Uinv)
 endif inversionFailed

end function math_rotationalPart33


!--------------------------------------------------------------------------------------------------
!> @brief Eigenvalues of symmetric matrix m
! will return NaN on error
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym(m)
 use, intrinsic :: &
   IEEE_arithmetic
 
 implicit none
 real(pReal), dimension(:,:),                  intent(in)  :: m
 real(pReal), dimension(size(m,1))                         :: math_eigenvaluesSym
 real(pReal), dimension(size(m,1),size(m,1))               :: vectors
 integer(pInt) :: info
 real(pReal), dimension((64+2)*size(m,1)) :: work                                                    ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f
 external :: &
  dsyev

 vectors = m                                                                                         ! copy matrix to input (doubles as output) array
 call dsyev('N','U',size(m,1),vectors,size(m,1),math_eigenvaluesSym,work,(64+2)*size(m,1),info)
 if (info /= 0_pInt) math_eigenvaluesSym = IEEE_value(1.0_pReal,IEEE_quiet_NaN)

end function math_eigenvaluesSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues of symmetric 33 matrix m using an analytical expression
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details similar to http://arxiv.org/abs/physics/0610206 (DSYEVC3)
!> but apparently more stable solution and has general LAPACK powered version for arbritrary sized
!> matrices as fallback
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym33(m)

 implicit none
 real(pReal), intent(in), dimension(3,3) :: m
 real(pReal), dimension(3) :: math_eigenvaluesSym33,invariants
 real(pReal) :: P, Q, rho, phi
 real(pReal), parameter :: TOL=1.e-14_pReal

 invariants = math_invariantsSym33(m)                                                               ! invariants are coefficients in characteristic polynomial apart for the sign of c0 and c2 in http://arxiv.org/abs/physics/0610206

 P = invariants(2)-invariants(1)**2.0_pReal/3.0_pReal                                               ! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)
 Q = -2.0_pReal/27.0_pReal*invariants(1)**3.0_pReal+product(invariants(1:2))/3.0_pReal-invariants(3)! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)

 if(all(abs([P,Q]) < TOL)) then
   math_eigenvaluesSym33 = math_eigenvaluesSym(m)
 else
   rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   phi=acos(math_limit(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
   math_eigenvaluesSym33 = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                         [cos(phi/3.0_pReal), &
                          cos((phi+2.0_pReal*PI)/3.0_pReal), &
                          cos((phi+4.0_pReal*PI)/3.0_pReal) &
                         ] + invariants(1)/3.0_pReal
 endif

end function math_eigenvaluesSym33


!--------------------------------------------------------------------------------------------------
!> @brief invariants of symmetrix 33 matrix m
!--------------------------------------------------------------------------------------------------
pure function math_invariantsSym33(m)

 implicit none
 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal), dimension(3) :: math_invariantsSym33

 math_invariantsSym33(1) = math_trace33(m)
 math_invariantsSym33(2) = m(1,1)*m(2,2) + m(1,1)*m(3,3) + m(2,2)*m(3,3) &
                         -(m(1,2)**2     + m(1,3)**2     + m(2,3)**2)
 math_invariantsSym33(3) = math_detSym33(m)

end function math_invariantsSym33


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

!--------------------------------------------------------------------------------------------------
 contains
 
 !-------------------------------------------------------------------------------------------------
 !> @brief computes an element of a Halton sequence.
 !> @details Only the absolute value of SEED is considered. SEED = 0 is allowed, and returns R = 0.
 !> @details Halton Bases should be distinct prime numbers. This routine only checks that each base 
 !> @details is greater than 1.
 !> @details Reference:
 !> @details J.H. Halton: On the efficiency of certain quasi-random sequences of points in evaluating 
 !> @details multi-dimensional integrals, Numerische Mathematik, Volume 2, pages 84-90, 1960.
 !> @author John Burkardt
 !-------------------------------------------------------------------------------------------------
 subroutine i_to_halton (seed, base, ndim, r)
   use IO, only: &
     IO_error
   
   implicit none
   integer(pInt), intent(in) :: &
     ndim, &                                                                                        !< dimension of the sequence
     seed                                                                                           !< index of the desired element
   integer(pInt), intent(in),  dimension(ndim) :: base                                              !< Halton bases
   real(pReal),   intent(out), dimension(ndim) :: r                                                 !< the SEED-th element of the Halton sequence for the given bases
  
   real(pReal),                dimension(ndim) ::  base_inv
   integer(pInt),              dimension(ndim) :: &
     digit, &
     seed2
  
   seed2 = abs(seed)
   r = 0.0_pReal
  
   if (any (base(1:ndim) <= 1_pInt)) call IO_error(error_ID=405_pInt)
  
   base_inv(1:ndim) = 1.0_pReal / real (base(1:ndim), pReal)
  
   do while ( any ( seed2(1:ndim) /= 0_pInt) )
     digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim))
     r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), pReal) * base_inv(1:ndim)
     base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), pReal)
     seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
   enddo
  
 end subroutine i_to_halton


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
  use IO, only: &
    IO_lc

 implicit none
 character(len = *), intent(in) :: & 
   action_halton, &                                                                                 !< desired action: GET the value of a particular quantity, SET the value of a particular quantity, INC the value of a particular quantity (only for SEED)
   name_halton                                                                                      !< name of the quantity: BASE: Halton base(s), NDIM: spatial dimension, SEED: current Halton seed
 integer(pInt), dimension(*), intent(inout) :: value_halton
 integer(pInt), allocatable, save, dimension(:) :: base
 logical, save :: first_call = .true.
 integer(pInt), intent(in) :: ndim                                                                  !< dimension of the quantity
 integer(pInt), save :: ndim_save = 0_pInt, seed = 1_pInt
 integer(pInt) :: i

 if (first_call) then
   ndim_save = 1_pInt
   allocate(base(ndim_save))
   base(1) = 2_pInt
   first_call = .false.
 endif
 
!--------------------------------------------------------------------------------------------------
! Set
 actionHalton: if(IO_lc(action_halton(1:1)) == 's') then

   nameSet: if(IO_lc(name_halton(1:1)) == 'b') then
     if(ndim_save /= ndim) ndim_save = ndim
     base = value_halton(1:ndim)
   elseif(IO_lc(name_halton(1:1)) == 'n') then nameSet
     if(ndim_save /= value_halton(1)) then
       ndim_save = value_halton(1)
       base = [(prime(i),i=1_pInt,ndim_save)]
     else
       ndim_save = value_halton(1)
     endif
   elseif(IO_lc(name_halton(1:1)) == 's') then nameSet
     seed = value_halton(1)
   endif nameSet
 
!--------------------------------------------------------------------------------------------------
! Get
 elseif(IO_lc(action_halton(1:1)) == 'g') then actionHalton
   nameGet: if(IO_lc(name_halton(1:1)) == 'b') then
     if(ndim /= ndim_save) then
       ndim_save = ndim
       base = [(prime(i),i=1_pInt,ndim_save)]
     endif
     value_halton(1:ndim_save) = base(1:ndim_save)
   elseif(IO_lc(name_halton(1:1)) == 'n') then nameGet
     value_halton(1) = ndim_save
   elseif(IO_lc(name_halton(1:1)) == 's') then nameGet
     value_halton(1) = seed
   endif nameGet
   
!--------------------------------------------------------------------------------------------------
!   Increment
 elseif(IO_lc(action_halton(1:1)) == 'i') then actionHalton
   if(IO_lc(name_halton(1:1)) == 's') seed = seed + value_halton(1)
 endif actionHalton

!--------------------------------------------------------------------------------------------------
 contains
 
 !--------------------------------------------------------------------------------------------------
 !> @brief returns any of the first 1500 prime numbers.
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
   integer(pInt), intent(in) :: n                                                                   !< index of the desired prime number
   integer(pInt), dimension(0:1500), parameter :: &
     npvec = int([&
             1, &
             2,     3,     5,     7,    11,    13,    17,    19,    23,    29, &
            31,    37,    41,    43,    47,    53,    59,    61,    67,    71, &
            73,    79,    83,    89,    97,   101,   103,   107,   109,   113, &
           127,   131,   137,   139,   149,   151,   157,   163,   167,   173, &
           179,   181,   191,   193,   197,   199,   211,   223,   227,   229, &
           233,   239,   241,   251,   257,   263,   269,   271,   277,   281, &
           283,   293,   307,   311,   313,   317,   331,   337,   347,   349, &
           353,   359,   367,   373,   379,   383,   389,   397,   401,   409, &
           419,   421,   431,   433,   439,   443,   449,   457,   461,   463, &
           467,   479,   487,   491,   499,   503,   509,   521,   523,   541, &
     ! 101:200
           547,   557,   563,   569,   571,   577,   587,   593,   599,   601, &
           607,   613,   617,   619,   631,   641,   643,   647,   653,   659, &
           661,   673,   677,   683,   691,   701,   709,   719,   727,   733, &
           739,   743,   751,   757,   761,   769,   773,   787,   797,   809, &
           811,   821,   823,   827,   829,   839,   853,   857,   859,   863, &
           877,   881,   883,   887,   907,   911,   919,   929,   937,   941, &
           947,   953,   967,   971,   977,   983,   991,   997,  1009,  1013, &
          1019,  1021,  1031,  1033,  1039,  1049,  1051,  1061,  1063,  1069, &
          1087,  1091,  1093,  1097,  1103,  1109,  1117,  1123,  1129,  1151, &
          1153,  1163,  1171,  1181,  1187,  1193,  1201,  1213,  1217,  1223, &
     ! 201:300
          1229,  1231,  1237,  1249,  1259,  1277,  1279,  1283,  1289,  1291, &
          1297,  1301,  1303,  1307,  1319,  1321,  1327,  1361,  1367,  1373, &
          1381,  1399,  1409,  1423,  1427,  1429,  1433,  1439,  1447,  1451, &
          1453,  1459,  1471,  1481,  1483,  1487,  1489,  1493,  1499,  1511, &
          1523,  1531,  1543,  1549,  1553,  1559,  1567,  1571,  1579,  1583, &
          1597,  1601,  1607,  1609,  1613,  1619,  1621,  1627,  1637,  1657, &
          1663,  1667,  1669,  1693,  1697,  1699,  1709,  1721,  1723,  1733, &
          1741,  1747,  1753,  1759,  1777,  1783,  1787,  1789,  1801,  1811, &
          1823,  1831,  1847,  1861,  1867,  1871,  1873,  1877,  1879,  1889, &
          1901,  1907,  1913,  1931,  1933,  1949,  1951,  1973,  1979,  1987, &
     ! 301:400
          1993,  1997,  1999,  2003,  2011,  2017,  2027,  2029,  2039,  2053, &
          2063,  2069,  2081,  2083,  2087,  2089,  2099,  2111,  2113,  2129, &
          2131,  2137,  2141,  2143,  2153,  2161,  2179,  2203,  2207,  2213, &
          2221,  2237,  2239,  2243,  2251,  2267,  2269,  2273,  2281,  2287, &
          2293,  2297,  2309,  2311,  2333,  2339,  2341,  2347,  2351,  2357, &
          2371,  2377,  2381,  2383,  2389,  2393,  2399,  2411,  2417,  2423, &
          2437,  2441,  2447,  2459,  2467,  2473,  2477,  2503,  2521,  2531, &
          2539,  2543,  2549,  2551,  2557,  2579,  2591,  2593,  2609,  2617, &
          2621,  2633,  2647,  2657,  2659,  2663,  2671,  2677,  2683,  2687, &
          2689,  2693,  2699,  2707,  2711,  2713,  2719,  2729,  2731,  2741, &
     ! 401:500
          2749,  2753,  2767,  2777,  2789,  2791,  2797,  2801,  2803,  2819, &
          2833,  2837,  2843,  2851,  2857,  2861,  2879,  2887,  2897,  2903, &
          2909,  2917,  2927,  2939,  2953,  2957,  2963,  2969,  2971,  2999, &
          3001,  3011,  3019,  3023,  3037,  3041,  3049,  3061,  3067,  3079, &
          3083,  3089,  3109,  3119,  3121,  3137,  3163,  3167,  3169,  3181, &
          3187,  3191,  3203,  3209,  3217,  3221,  3229,  3251,  3253,  3257, &
          3259,  3271,  3299,  3301,  3307,  3313,  3319,  3323,  3329,  3331, &
          3343,  3347,  3359,  3361,  3371,  3373,  3389,  3391,  3407,  3413, &
          3433,  3449,  3457,  3461,  3463,  3467,  3469,  3491,  3499,  3511, &
          3517,  3527,  3529,  3533,  3539,  3541,  3547,  3557,  3559,  3571, &
     ! 501:600
          3581,  3583,  3593,  3607,  3613,  3617,  3623,  3631,  3637,  3643, &
          3659,  3671,  3673,  3677,  3691,  3697,  3701,  3709,  3719,  3727, &
          3733,  3739,  3761,  3767,  3769,  3779,  3793,  3797,  3803,  3821, &
          3823,  3833,  3847,  3851,  3853,  3863,  3877,  3881,  3889,  3907, &
          3911,  3917,  3919,  3923,  3929,  3931,  3943,  3947,  3967,  3989, &
          4001,  4003,  4007,  4013,  4019,  4021,  4027,  4049,  4051,  4057, &
          4073,  4079,  4091,  4093,  4099,  4111,  4127,  4129,  4133,  4139, &
          4153,  4157,  4159,  4177,  4201,  4211,  4217,  4219,  4229,  4231, &
          4241,  4243,  4253,  4259,  4261,  4271,  4273,  4283,  4289,  4297, &
          4327,  4337,  4339,  4349,  4357,  4363,  4373,  4391,  4397,  4409, &
     ! 601:700
          4421,  4423,  4441,  4447,  4451,  4457,  4463,  4481,  4483,  4493, &
          4507,  4513,  4517,  4519,  4523,  4547,  4549,  4561,  4567,  4583, &
          4591,  4597,  4603,  4621,  4637,  4639,  4643,  4649,  4651,  4657, &
          4663,  4673,  4679,  4691,  4703,  4721,  4723,  4729,  4733,  4751, &
          4759,  4783,  4787,  4789,  4793,  4799,  4801,  4813,  4817,  4831, &
          4861,  4871,  4877,  4889,  4903,  4909,  4919,  4931,  4933,  4937, &
          4943,  4951,  4957,  4967,  4969,  4973,  4987,  4993,  4999,  5003, &
          5009,  5011,  5021,  5023,  5039,  5051,  5059,  5077,  5081,  5087, &
          5099,  5101,  5107,  5113,  5119,  5147,  5153,  5167,  5171,  5179, &
          5189,  5197,  5209,  5227,  5231,  5233,  5237,  5261,  5273,  5279, &
     ! 701:800
          5281,  5297,  5303,  5309,  5323,  5333,  5347,  5351,  5381,  5387, &
          5393,  5399,  5407,  5413,  5417,  5419,  5431,  5437,  5441,  5443, &
          5449,  5471,  5477,  5479,  5483,  5501,  5503,  5507,  5519,  5521, &
          5527,  5531,  5557,  5563,  5569,  5573,  5581,  5591,  5623,  5639, &
          5641,  5647,  5651,  5653,  5657,  5659,  5669,  5683,  5689,  5693, &
          5701,  5711,  5717,  5737,  5741,  5743,  5749,  5779,  5783,  5791, &
          5801,  5807,  5813,  5821,  5827,  5839,  5843,  5849,  5851,  5857, &
          5861,  5867,  5869,  5879,  5881,  5897,  5903,  5923,  5927,  5939, &
          5953,  5981,  5987,  6007,  6011,  6029,  6037,  6043,  6047,  6053, &
          6067,  6073,  6079,  6089,  6091,  6101,  6113,  6121,  6131,  6133, &
     ! 801:900
          6143,  6151,  6163,  6173,  6197,  6199,  6203,  6211,  6217,  6221, &
          6229,  6247,  6257,  6263,  6269,  6271,  6277,  6287,  6299,  6301, &
          6311,  6317,  6323,  6329,  6337,  6343,  6353,  6359,  6361,  6367, &
          6373,  6379,  6389,  6397,  6421,  6427,  6449,  6451,  6469,  6473, &
          6481,  6491,  6521,  6529,  6547,  6551,  6553,  6563,  6569,  6571, &
          6577,  6581,  6599,  6607,  6619,  6637,  6653,  6659,  6661,  6673, &
          6679,  6689,  6691,  6701,  6703,  6709,  6719,  6733,  6737,  6761, &
          6763,  6779,  6781,  6791,  6793,  6803,  6823,  6827,  6829,  6833, &
          6841,  6857,  6863,  6869,  6871,  6883,  6899,  6907,  6911,  6917, &
          6947,  6949,  6959,  6961,  6967,  6971,  6977,  6983,  6991,  6997, &
     ! 901:1000 
          7001,  7013,  7019,  7027,  7039,  7043,  7057,  7069,  7079,  7103, &
          7109,  7121,  7127,  7129,  7151,  7159,  7177,  7187,  7193,  7207, &
          7211,  7213,  7219,  7229,  7237,  7243,  7247,  7253,  7283,  7297, &
          7307,  7309,  7321,  7331,  7333,  7349,  7351,  7369,  7393,  7411, &
          7417,  7433,  7451,  7457,  7459,  7477,  7481,  7487,  7489,  7499, &
          7507,  7517,  7523,  7529,  7537,  7541,  7547,  7549,  7559,  7561, &
          7573,  7577,  7583,  7589,  7591,  7603,  7607,  7621,  7639,  7643, &
          7649,  7669,  7673,  7681,  7687,  7691,  7699,  7703,  7717,  7723, &
          7727,  7741,  7753,  7757,  7759,  7789,  7793,  7817,  7823,  7829, &
          7841,  7853,  7867,  7873,  7877,  7879,  7883,  7901,  7907,  7919, &
     ! 1001:1100
          7927,  7933,  7937,  7949,  7951,  7963,  7993,  8009,  8011,  8017, &
          8039,  8053,  8059,  8069,  8081,  8087,  8089,  8093,  8101,  8111, &
          8117,  8123,  8147,  8161,  8167,  8171,  8179,  8191,  8209,  8219, &
          8221,  8231,  8233,  8237,  8243,  8263,  8269,  8273,  8287,  8291, &
          8293,  8297,  8311,  8317,  8329,  8353,  8363,  8369,  8377,  8387, &
          8389,  8419,  8423,  8429,  8431,  8443,  8447,  8461,  8467,  8501, &
          8513,  8521,  8527,  8537,  8539,  8543,  8563,  8573,  8581,  8597, &
          8599,  8609,  8623,  8627,  8629,  8641,  8647,  8663,  8669,  8677, &
          8681,  8689,  8693,  8699,  8707,  8713,  8719,  8731,  8737,  8741, &
          8747,  8753,  8761,  8779,  8783,  8803,  8807,  8819,  8821,  8831, &
     ! 1101:1200
          8837,  8839,  8849,  8861,  8863,  8867,  8887,  8893,  8923,  8929, &
          8933,  8941,  8951,  8963,  8969,  8971,  8999,  9001,  9007,  9011, &
          9013,  9029,  9041,  9043,  9049,  9059,  9067,  9091,  9103,  9109, &
          9127,  9133,  9137,  9151,  9157,  9161,  9173,  9181,  9187,  9199, &
          9203,  9209,  9221,  9227,  9239,  9241,  9257,  9277,  9281,  9283, &
          9293,  9311,  9319,  9323,  9337,  9341,  9343,  9349,  9371,  9377, &
          9391,  9397,  9403,  9413,  9419,  9421,  9431,  9433,  9437,  9439, &
          9461,  9463,  9467,  9473,  9479,  9491,  9497,  9511,  9521,  9533, &
          9539,  9547,  9551,  9587,  9601,  9613,  9619,  9623,  9629,  9631, &
          9643,  9649,  9661,  9677,  9679,  9689,  9697,  9719,  9721,  9733, &
     ! 1201:1300
          9739,  9743,  9749,  9767,  9769,  9781,  9787,  9791,  9803,  9811, &
          9817,  9829,  9833,  9839,  9851,  9857,  9859,  9871,  9883,  9887, &
          9901,  9907,  9923,  9929,  9931,  9941,  9949,  9967,  9973, 10007, &
         10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099, &
         10103, 10111, 10133, 10139, 10141, 10151, 10159, 10163, 10169, 10177, &
         10181, 10193, 10211, 10223, 10243, 10247, 10253, 10259, 10267, 10271, &
         10273, 10289, 10301, 10303, 10313, 10321, 10331, 10333, 10337, 10343, &
         10357, 10369, 10391, 10399, 10427, 10429, 10433, 10453, 10457, 10459, &
         10463, 10477, 10487, 10499, 10501, 10513, 10529, 10531, 10559, 10567, &
         10589, 10597, 10601, 10607, 10613, 10627, 10631, 10639, 10651, 10657, &
     ! 1301:1400
         10663, 10667, 10687, 10691, 10709, 10711, 10723, 10729, 10733, 10739, &
         10753, 10771, 10781, 10789, 10799, 10831, 10837, 10847, 10853, 10859, &
         10861, 10867, 10883, 10889, 10891, 10903, 10909, 19037, 10939, 10949, &
         10957, 10973, 10979, 10987, 10993, 11003, 11027, 11047, 11057, 11059, &
         11069, 11071, 11083, 11087, 11093, 11113, 11117, 11119, 11131, 11149, &
         11159, 11161, 11171, 11173, 11177, 11197, 11213, 11239, 11243, 11251, &
         11257, 11261, 11273, 11279, 11287, 11299, 11311, 11317, 11321, 11329, &
         11351, 11353, 11369, 11383, 11393, 11399, 11411, 11423, 11437, 11443, &
         11447, 11467, 11471, 11483, 11489, 11491, 11497, 11503, 11519, 11527, &
         11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657, &
     ! 1401:1500
         11677, 11681, 11689, 11699, 11701, 11717, 11719, 11731, 11743, 11777, &
         11779, 11783, 11789, 11801, 11807, 11813, 11821, 11827, 11831, 11833, &
         11839, 11863, 11867, 11887, 11897, 11903, 11909, 11923, 11927, 11933, &
         11939, 11941, 11953, 11959, 11969, 11971, 11981, 11987, 12007, 12011, &
         12037, 12041, 12043, 12049, 12071, 12073, 12097, 12101, 12107, 12109, &
         12113, 12119, 12143, 12149, 12157, 12161, 12163, 12197, 12203, 12211, &
         12227, 12239, 12241, 12251, 12253, 12263, 12269, 12277, 12281, 12289, &
         12301, 12323, 12329, 12343, 12347, 12373, 12377, 12379, 12391, 12401, &
         12409, 12413, 12421, 12433, 12437, 12451, 12457, 12473, 12479, 12487, &
         12491, 12497, 12503, 12511, 12517, 12527, 12539, 12541, 12547, 12553],pInt)
  
   if (n < size(npvec)) then
     prime = npvec(n)
   else
     call IO_error(error_ID=406_pInt)
   end if
  
 end function prime

end subroutine halton_memory


!--------------------------------------------------------------------------------------------------
!> @brief sets the dimension for a Halton sequence
!> @author John Burkardt
!--------------------------------------------------------------------------------------------------
subroutine halton_ndim_set(ndim)

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
!> @brief limits a scalar value to a certain range (either one or two sided)
! Will return NaN if left > right
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_limit(a, left, right)
 use, intrinsic :: &
   IEEE_arithmetic

 implicit none
 real(pReal), intent(in) :: a
 real(pReal), intent(in), optional :: left, right

   
 math_limit = min ( &
                   max (merge(left, -huge(a), present(left)), a), &
                        merge(right, huge(a), present(right)) &
                  )

 if (present(left) .and. present(right)) &
   math_limit = merge (IEEE_value(1.0_pReal,IEEE_quiet_NaN),math_limit, left>right)

end function math_limit


!--------------------------------------------------------------------------------------------------
!> @brief Modified Bessel I function of order 0
!> @author John Burkardt
!> @details original version available on https://people.sc.fsu.edu/~jburkardt/f_src/toms715/toms715.html
!--------------------------------------------------------------------------------------------------
real(pReal) function bessel_i0 (x)
 use, intrinsic :: IEEE_ARITHMETIC

 implicit none
 real(pReal), intent(in) :: x
 integer(pInt) :: i
 real(pReal) :: sump_p, sump_q, xAbs, xx
 real(pReal), parameter, dimension(15) :: p_small = real( &
   [-5.2487866627945699800e-18, -1.5982226675653184646e-14, -2.6843448573468483278e-11, &
    -3.0517226450451067446e-08, -2.5172644670688975051e-05, -1.5453977791786851041e-02, &
    -7.0935347449210549190e+00, -2.4125195876041896775e+03, -5.9545626019847898221e+05, &
    -1.0313066708737980747e+08, -1.1912746104985237192e+10, -8.4925101247114157499e+11, &
    -3.2940087627407749166e+13, -5.5050369673018427753e+14, -2.2335582639474375249e+15], pReal)
 real(pReal), parameter, dimension(5) :: q_small = real( &
   [-3.7277560179962773046e+03,  6.5158506418655165707e+06, -6.5626560740833869295e+09, &
     3.7604188704092954661e+12, -9.7087946179594019126e+14], pReal)
 real(pReal), parameter, dimension(8) :: p_large = real( &
   [-3.9843750000000000000e-01,  2.9205384596336793945e+00, -2.4708469169133954315e+00, &
     4.7914889422856814203e-01, -3.7384991926068969150e-03, -2.6801520353328635310e-03, &
     9.9168777670983678974e-05, -2.1877128189032726730e-06], pReal)
 real(pReal), parameter, dimension(7) :: q_large = real( &
   [-3.1446690275135491500e+01,  8.5539563258012929600e+01, -6.0228002066743340583e+01, &
     1.3982595353892851542e+01, -1.1151759188741312645e+00,  3.2547697594819615062e-02, &
    -5.5194330231005480228e-04], pReal)


 xAbs = abs(x)

 argRange: if (xAbs < 5.55e-17_pReal) then
   bessel_i0 = 1.0_pReal
 else if (xAbs < 15.0_pReal) then argRange
   xx = xAbs**2.0_pReal
   sump_p = p_small(1)
   do i = 2, 15
     sump_p = sump_p * xx + p_small(i)
   end do
   xx = xx - 225.0_pReal
   sump_q = ((((xx+q_small(1))*xx+q_small(2))*xx+q_small(3))*xx+q_small(4))*xx+q_small(5)
   bessel_i0 = sump_p / sump_q
 else if (xAbs <= 713.986_pReal) then argRange
   xx = 1.0_pReal / xAbs - 2.0_pReal/30.0_pReal
   sump_p = ((((((p_large(1)*xx+p_large(2))*xx+p_large(3))*xx+p_large(4))*xx+ &
                 p_large(5))*xx+p_large(6))*xx+p_large(7))*xx+p_large(8)
   sump_q = ((((((xx+q_large(1))*xx+q_large(2))*xx+q_large(3))*xx+ &
                  q_large(4))*xx+q_large(5))*xx+q_large(6))*xx+q_large(7)
   bessel_i0 = sump_p / sump_q

   avoidOverflow: if (xAbs > 698.986_pReal) then
     bessel_i0 = ((bessel_i0*exp(xAbs-40.0_pReal)-p_large(1)*exp(xAbs-40.0_pReal))/sqrt(xAbs))*exp(40.0)
   else avoidOverflow
     bessel_i0 = ((bessel_i0*exp(xAbs)-p_large(1)*exp(xAbs))/sqrt(xAbs))
   endif avoidOverflow

 else argRange
   bessel_i0 =  IEEE_value(bessel_i0,IEEE_positive_inf)
 end if argRange

end function bessel_i0


!--------------------------------------------------------------------------------------------------
!> @brief Modified Bessel I function of order 1
!> @author John Burkardt
!> @details original version available on https://people.sc.fsu.edu/~jburkardt/f_src/toms715/toms715.html
!--------------------------------------------------------------------------------------------------
real(pReal) function bessel_i1 (x)
 use, intrinsic :: IEEE_ARITHMETIC

 implicit none
 real(pReal), intent(in) :: x
 integer(pInt) :: i
 real(pReal) :: sump_p, sump_q, xAbs, xx
 real(pReal), dimension(15), parameter ::  p_small = real( &
   [-1.9705291802535139930e-19, -6.5245515583151902910e-16, -1.1928788903603238754e-12, &
    -1.4831904935994647675e-09, -1.3466829827635152875e-06, -9.1746443287817501309e-04, &
    -4.7207090827310162436e-01, -1.8225946631657315931e+02, -5.1894091982308017540e+04, &
    -1.0588550724769347106e+07, -1.4828267606612366099e+09, -1.3357437682275493024e+11, &
    -6.9876779648010090070e+12, -1.7732037840791591320e+14, -1.4577180278143463643e+15], pReal)
 real(pReal), dimension(5), parameter ::  q_small = real( &
   [-4.0076864679904189921e+03,  7.4810580356655069138e+06, -8.0059518998619764991e+09, &
     4.8544714258273622913e+12, -1.3218168307321442305e+15], pReal)
 real(pReal), dimension(8), parameter ::  p_large = real( &
   [-6.0437159056137600000e-02,  4.5748122901933459000e-01, -4.2843766903304806403e-01, &
     9.7356000150886612134e-02, -3.2457723974465568321e-03, -3.6395264712121795296e-04, &
     1.6258661867440836395e-05, -3.6347578404608223492e-07], pReal)
 real(pReal), dimension(6), parameter ::  q_large = real( &
   [-3.8806586721556593450e+00,  3.2593714889036996297e+00, -8.5017476463217924408e-01, &
     7.4212010813186530069e-02, -2.2835624489492512649e-03,  3.7510433111922824643e-05], pReal)
 real(pReal), parameter :: pbar = 3.98437500e-01


 xAbs = abs(x)

 argRange: if (xAbs < 5.55e-17_pReal) then
   bessel_i1 = 0.5_pReal * xAbs
 else if (xAbs < 15.0_pReal) then argRange
   xx = xAbs**2.0_pReal
   sump_p = p_small(1)
   do i = 2, 15
     sump_p = sump_p * xx + p_small(i)
   end do
   xx = xx - 225.0_pReal
   sump_q = ((((xx+q_small(1))*xx+q_small(2))*xx+q_small(3))*xx+q_small(4)) * xx + q_small(5)
   bessel_i1 = (sump_p / sump_q) * xAbs
 else if (xAbs <= 713.986_pReal) then argRange
   xx = 1.0_pReal / xAbs - 2.0_pReal/30.0_pReal
   sump_p = ((((((p_large(1)*xx+p_large(2))*xx+p_large(3))*xx+p_large(4))*xx+&
                 p_large(5))*xx+p_large(6))*xx+p_large(7))*xx+p_large(8)
   sump_q = (((((xx+q_large(1))*xx+q_large(2))*xx+q_large(3))*xx+ q_large(4))*xx+q_large(5))*xx+q_large(6)
   bessel_i1 = sump_p / sump_q

   avoidOverflow: if (xAbs > 698.986_pReal) then
     bessel_i1 = ((bessel_i1 * exp(xAbs-40.0_pReal) + pbar * exp(xAbs-40.0_pReal)) / sqrt(xAbs)) * exp(40.0_pReal)
   else avoidOverflow
     bessel_i1 = ((bessel_i1 * exp(xAbs) + pbar * exp(xAbs)) / sqrt(xAbs))
   endif avoidOverflow

 else argRange
   bessel_i1 =  IEEE_value(bessel_i1,IEEE_positive_inf)
 end if argRange

 if (x < 0.0_pReal) bessel_i1 = -bessel_i1

end function bessel_i1
 
end module math
