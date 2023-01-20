!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mathematical library, including random number generation and tensor representations
!--------------------------------------------------------------------------------------------------
module math
  use prec
  use IO
  use config
  use YAML_types
  use parallelization
  use LAPACK_interface

#ifdef PETSC
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif
#endif

  implicit none(type,external)
  public
#if __INTEL_COMPILER >= 1900
  ! do not make use of associated entities available to other modules
  private :: &
    IO, &
    config
#endif

  real(pReal), parameter :: &
    PI = acos(-1.0_pReal), &                                                                        !< ratio of a circle's circumference to its diameter
    TAU = 2.0_pReal*PI, &                                                                           !< ratio of a circle's circumference to its radius
    INDEG = 360.0_pReal/TAU, &                                                                      !< conversion from radian to degree
    INRAD = TAU/360.0_pReal                                                                         !< conversion from degree to radian

  real(pReal), dimension(3,3), parameter :: &
    math_I3 = reshape([&
      1.0_pReal,0.0_pReal,0.0_pReal, &
      0.0_pReal,1.0_pReal,0.0_pReal, &
      0.0_pReal,0.0_pReal,1.0_pReal  &
      ],shape(math_I3))                                                                             !< 3x3 Identity

  real(pReal), dimension(*), parameter, private :: &
    NRMMANDEL = [1.0_pReal, 1.0_pReal,1.0_pReal, sqrt(2.0_pReal), sqrt(2.0_pReal), sqrt(2.0_pReal)] !< forward weighting for Mandel notation

  real(pReal), dimension(*), parameter, private :: &
    INVNRMMANDEL = 1.0_pReal/NRMMANDEL                                                              !< backward weighting for Mandel notation

  integer, dimension (2,6), parameter, private :: &
    MAPNYE = reshape([&
      1,1, &
      2,2, &
      3,3, &
      1,2, &
      2,3, &
      1,3  &
      ],shape(MAPNYE))                                                                              !< arrangement in Nye notation.

  integer, dimension (2,6), parameter, private :: &
    MAPVOIGT = reshape([&
      1,1, &
      2,2, &
      3,3, &
      2,3, &
      1,3, &
      1,2  &
      ],shape(MAPVOIGT))                                                                            !< arrangement in Voigt notation

  integer, dimension (2,9), parameter, private :: &
    MAPPLAIN = reshape([&
      1,1, &
      1,2, &
      1,3, &
      2,1, &
      2,2, &
      2,3, &
      3,1, &
      3,2, &
      3,3  &
      ],shape(MAPPLAIN))                                                                            !< arrangement in Plain notation

!---------------------------------------------------------------------------------------------------
 private :: &
   selfTest

contains

!--------------------------------------------------------------------------------------------------
!> @brief initialization of random seed generator and internal checks
!--------------------------------------------------------------------------------------------------
subroutine math_init()

  real(pReal), dimension(4) :: randTest
  integer :: randSize
  integer, dimension(:), allocatable :: seed
  type(tDict), pointer :: &
    num_generic


  print'(/,1x,a)', '<<<+-  math init  -+>>>'; flush(IO_STDOUT)

  num_generic => config_numerics%get_dict('generic',defaultVal=emptyDict)

  call random_seed(size=randSize)
  allocate(seed(randSize))

  if (num_generic%contains('random_seed')) then
    seed = num_generic%get_as1dInt('random_seed',requiredSize=randSize)
  else
    call random_seed()
    call random_seed(get = seed)
  end if

  call random_seed(put = seed + worldrank*42_MPI_INTEGER_KIND)
  call random_number(randTest)

  print'(/,a,i2)',              ' size  of random seed:     ', randSize
  print*,                       'value of random seed:     ', seed
  print'(  a,4(/,26x,f17.14))', ' start of random sequence: ', randTest

  call selfTest()

end subroutine math_init


!--------------------------------------------------------------------------------------------------
!> @brief Sorting of two-dimensional integer arrays
!> @details Based on quicksort.
!  Sorting is done with respect to array(sortDim,:) and keeps array(/=sortDim,:) linked to it.
!  Default: sortDim=1
!--------------------------------------------------------------------------------------------------
pure recursive subroutine math_sort(a, istart, iend, sortDim)

  integer, dimension(:,:), intent(inout) :: a
  integer, intent(in),optional :: istart,iend, sortDim

  integer :: ipivot,s,e,d


  if (present(istart)) then
    s = istart
  else
    s = lbound(a,2)
  end if

  if (present(iend)) then
    e = iend
  else
    e = ubound(a,2)
  end if

  if (present(sortDim)) then
    d = sortDim
  else
    d = 1
  end if

  if (s < e) then
    call qsort_partition(a,ipivot, s,e, d)
    call math_sort(a, s, ipivot-1, d)
    call math_sort(a, ipivot+1, e, d)
  end if


  contains

  !-------------------------------------------------------------------------------------------------
  !> @brief Partitioning required for quicksort
  !-------------------------------------------------------------------------------------------------
  pure subroutine qsort_partition(a,p, istart, iend, sort)

    integer, dimension(:,:), intent(inout) :: a
    integer,                 intent(out)   :: p                                                     ! Pivot element
    integer,                 intent(in)    :: istart,iend,sort

    integer, dimension(size(a,1)) :: tmp
    integer :: i,j


    do
      ! find the first element on the right side less than or equal to the pivot point
      do j = iend, istart, -1
        if (a(sort,j) <= a(sort,istart)) exit
      end do
      ! find the first element on the left side greater than the pivot point
      do i = istart, iend
        if (a(sort,i) > a(sort,istart)) exit
      end do
      cross: if (i >= j) then ! exchange left value with pivot and return with the partition index
        tmp         = a(:,istart)
        a(:,istart) = a(:,j)
        a(:,j)      = tmp
        p           = j
        return
      else cross              ! exchange values
        tmp    = a(:,i)
        a(:,i) = a(:,j)
        a(:,j) = tmp
      end if cross
    end do

  end subroutine qsort_partition

end subroutine math_sort


!--------------------------------------------------------------------------------------------------
!> @brief vector expansion
!> @details takes a set of numbers (a,b,c,...) and corresponding multiples (x,y,z,...)
!> to return a vector of x times a, y times b, z times c, ...
!--------------------------------------------------------------------------------------------------
pure function math_expand(what,how)

  real(pReal),   dimension(:), intent(in) :: what
  integer,       dimension(:), intent(in) :: how
  real(pReal), dimension(sum(how)) ::  math_expand

  integer :: i


  if (sum(how) == 0) return

  do i = 1, size(how)
    math_expand(sum(how(1:i-1))+1:sum(how(1:i))) = what(mod(i-1,size(what))+1)
  end do

end function math_expand


!--------------------------------------------------------------------------------------------------
!> @brief range of integers starting at one
!--------------------------------------------------------------------------------------------------
pure function math_range(N)

  integer, intent(in) :: N                                                                          !< length of range
  integer, dimension(N) :: math_range

  integer :: i


  math_range = [(i,i=1,N)]

end function math_range


!--------------------------------------------------------------------------------------------------
!> @brief Rank two identity tensor of specified dimension.
!--------------------------------------------------------------------------------------------------
pure function math_eye(d)

  integer, intent(in) :: d                                                                          !< tensor dimension
  real(pReal), dimension(d,d) :: math_eye

  integer :: i


  math_eye = 0.0_pReal
  do i=1,d
    math_eye(i,i) = 1.0_pReal
  end do

end function math_eye


!--------------------------------------------------------------------------------------------------
!> @brief Symmetric rank four identity tensor.
!  from http://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)#Derivative_of_a_second-order_tensor_with_respect_to_itself
!--------------------------------------------------------------------------------------------------
pure function math_identity4th()

  real(pReal), dimension(3,3,3,3) :: math_identity4th

  integer :: i,j,k,l


#ifndef __INTEL_COMPILER
  do concurrent(i=1:3, j=1:3, k=1:3, l=1:3)
    math_identity4th(i,j,k,l) = 0.5_pReal*(math_I3(i,k)*math_I3(j,l)+math_I3(i,l)*math_I3(j,k))
  end do
#else
  forall(i=1:3, j=1:3, k=1:3, l=1:3) &
    math_identity4th(i,j,k,l) = 0.5_pReal*(math_I3(i,k)*math_I3(j,l)+math_I3(i,l)*math_I3(j,k))
#endif

end function math_identity4th


!--------------------------------------------------------------------------------------------------
!> @brief permutation tensor e_ijk
! e_ijk =  1 if even permutation of ijk
! e_ijk = -1 if odd permutation of ijk
! e_ijk =  0 otherwise
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_LeviCivita(i,j,k)

  integer, intent(in) :: i,j,k

  integer :: o


  if     (any([(all(cshift([i,j,k],o) == [1,2,3]),o=0,2)])) then
    math_LeviCivita = +1.0_pReal
  elseif (any([(all(cshift([i,j,k],o) == [3,2,1]),o=0,2)])) then
    math_LeviCivita = -1.0_pReal
  else
    math_LeviCivita =  0.0_pReal
  end if

end function math_LeviCivita


!--------------------------------------------------------------------------------------------------
!> @brief kronecker delta function d_ij
! d_ij = 1 if i = j
! d_ij = 0 otherwise
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_delta(i,j)

  integer, intent (in) :: i,j


  math_delta = merge(0.0_pReal, 1.0_pReal, i /= j)

end function math_delta


!--------------------------------------------------------------------------------------------------
!> @brief cross product a x b
!--------------------------------------------------------------------------------------------------
pure function math_cross(A,B)

  real(pReal), dimension(3), intent(in) ::  A,B
  real(pReal), dimension(3) :: math_cross


  math_cross = [ A(2)*B(3) -A(3)*B(2), &
                 A(3)*B(1) -A(1)*B(3), &
                 A(1)*B(2) -A(2)*B(1) ]

end function math_cross


!--------------------------------------------------------------------------------------------------
!> @brief outer product of arbitrary sized vectors (A ⊗ B / i,j)
!--------------------------------------------------------------------------------------------------
pure function math_outer(A,B)

  real(pReal), dimension(:), intent(in) ::  A,B
  real(pReal), dimension(size(A,1),size(B,1)) ::  math_outer

  integer :: i,j


#ifndef __INTEL_COMPILER
  do concurrent(i=1:size(A,1), j=1:size(B,1))
    math_outer(i,j) = A(i)*B(j)
  end do
#else
  forall(i=1:size(A,1), j=1:size(B,1)) math_outer(i,j) = A(i)*B(j)
#endif

end function math_outer


!--------------------------------------------------------------------------------------------------
!> @brief inner product of arbitrary sized vectors (A · B / i,i)
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_inner(A,B)

  real(pReal), dimension(:),         intent(in) :: A
  real(pReal), dimension(size(A,1)), intent(in) :: B


  math_inner = sum(A*B)

end function math_inner


!--------------------------------------------------------------------------------------------------
!> @brief double contraction of 3x3 matrices (A : B / ij,ij)
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_tensordot(A,B)

  real(pReal), dimension(3,3), intent(in) :: A,B


  math_tensordot = sum(A*B)

end function math_tensordot


!--------------------------------------------------------------------------------------------------
!> @brief matrix double contraction 3333x33 = 33 (ijkl,kl)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx33(A,B)

  real(pReal), dimension(3,3,3,3), intent(in) :: A
  real(pReal), dimension(3,3),     intent(in) :: B
  real(pReal), dimension(3,3) :: math_mul3333xx33

  integer :: i,j


#ifndef __INTEL_COMPILER
  do concurrent(i=1:3, j=1:3)
    math_mul3333xx33(i,j) = sum(A(i,j,1:3,1:3)*B(1:3,1:3))
  end do
#else
  forall (i=1:3, j=1:3) math_mul3333xx33(i,j) = sum(A(i,j,1:3,1:3)*B(1:3,1:3))
#endif

end function math_mul3333xx33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3333x3333 = 3333 (ijkl,klmn)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx3333(A,B)

  real(pReal), dimension(3,3,3,3), intent(in) :: A
  real(pReal), dimension(3,3,3,3), intent(in) :: B
  real(pReal), dimension(3,3,3,3) :: math_mul3333xx3333

  integer :: i,j,k,l


#ifndef __INTEL_COMPILER
  do concurrent(i=1:3, j=1:3, k=1:3, l=1:3)
    math_mul3333xx3333(i,j,k,l) = sum(A(i,j,1:3,1:3)*B(1:3,1:3,k,l))
  end do
#else
  forall(i=1:3, j=1:3, k=1:3, l=1:3) math_mul3333xx3333(i,j,k,l) = sum(A(i,j,1:3,1:3)*B(1:3,1:3,k,l))
#endif

end function math_mul3333xx3333


!--------------------------------------------------------------------------------------------------
!> @brief 3x3 matrix exponential up to series approximation order n (default 5)
!--------------------------------------------------------------------------------------------------
pure function math_exp33(A,n)

  real(pReal), dimension(3,3), intent(in) :: A
  integer,                     intent(in), optional :: n
  real(pReal), dimension(3,3) :: B, math_exp33

  real(pReal) :: invFac
  integer     :: n_,i


  if (present(n)) then
    n_ = n
  else
    n_ = 5
  end if

  invFac     = 1.0_pReal                                                                            ! 0!
  B          = math_I3
  math_exp33 = math_I3                                                                              ! A^0 = I

  do i = 1, n_
    invFac = invFac/real(i,pReal)                                                                   ! invfac = 1/(i!)
    B = matmul(B,A)
    math_exp33 = math_exp33 + invFac*B                                                              ! exp = SUM (A^i)/(i!)
  end do

end function math_exp33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 3x3 matrix (function)
!> @details Direct Cramer inversion of matrix A. Returns all zeroes if not possible, i.e.
! if determinant is close to zero
!--------------------------------------------------------------------------------------------------
pure function math_inv33(A)

  real(pReal), dimension(3,3), intent(in) :: A
  real(pReal), dimension(3,3) :: math_inv33

  real(pReal) :: DetA
  logical     :: error


  call math_invert33(math_inv33,DetA,error,A)
  if (error) math_inv33 = 0.0_pReal

end function math_inv33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 3x3 matrix (subroutine)
!> @details Direct Cramer inversion of matrix A. Also returns determinant
!  Returns an error if not possible, i.e. if determinant is close to zero
!--------------------------------------------------------------------------------------------------
pure subroutine math_invert33(InvA,DetA,error, A)

  real(pReal), dimension(3,3), intent(out) :: InvA
  real(pReal),                 intent(out), optional :: DetA
  logical,                     intent(out) :: error
  real(pReal), dimension(3,3), intent(in)  :: A

  real(pReal) :: Det


  InvA(1,1) =  A(2,2) * A(3,3) - A(2,3) * A(3,2)
  InvA(2,1) = -A(2,1) * A(3,3) + A(2,3) * A(3,1)
  InvA(3,1) =  A(2,1) * A(3,2) - A(2,2) * A(3,1)

  Det = A(1,1) * InvA(1,1) + A(1,2) * InvA(2,1) + A(1,3) * InvA(3,1)

  if (dEq0(Det)) then
    InvA = 0.0_pReal
    if (present(DetA)) DetA = 0.0_pReal
    error = .true.
  else
    InvA(1,2) = -A(1,2) * A(3,3) + A(1,3) * A(3,2)
    InvA(2,2) =  A(1,1) * A(3,3) - A(1,3) * A(3,1)
    InvA(3,2) = -A(1,1) * A(3,2) + A(1,2) * A(3,1)

    InvA(1,3) =  A(1,2) * A(2,3) - A(1,3) * A(2,2)
    InvA(2,3) = -A(1,1) * A(2,3) + A(1,3) * A(2,1)
    InvA(3,3) =  A(1,1) * A(2,2) - A(1,2) * A(2,1)

    InvA = InvA/Det
    if (present(DetA)) DetA = Det
    error = .false.
  end if

end subroutine math_invert33


!--------------------------------------------------------------------------------------------------
!> @brief Inversion of symmetriced 3x3x3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_invSym3333(A)

  real(pReal),dimension(3,3,3,3) :: math_invSym3333

  real(pReal),dimension(3,3,3,3),intent(in) :: A

  integer,     dimension(6)   :: ipiv6
  real(pReal), dimension(6,6) :: temp66
  real(pReal), dimension(6*6) :: work
  integer                     :: ierr_i, ierr_f


  temp66 = math_sym3333to66(A)
  call dgetrf(6,6,temp66,6,ipiv6,ierr_i)
  call dgetri(6,temp66,6,ipiv6,work,size(work,1),ierr_f)
  if (ierr_i /= 0 .or. ierr_f /= 0) then
    error stop 'matrix inversion error'
  else
    math_invSym3333 = math_66toSym3333(temp66)
  end if

end function math_invSym3333


!--------------------------------------------------------------------------------------------------
!> @brief invert quadratic matrix of arbitrary dimension
!--------------------------------------------------------------------------------------------------
pure subroutine math_invert(InvA, error, A)

  real(pReal), dimension(:,:),                 intent(in)  :: A
  real(pReal), dimension(size(A,1),size(A,1)), intent(out) :: invA
  logical,                                     intent(out) :: error

  integer,     dimension(size(A,1))    :: ipiv
  real(pReal), dimension(size(A,1)**2) :: work
  integer                              :: ierr


  invA = A
  call dgetrf(size(A,1),size(A,1),invA,size(A,1),ipiv,ierr)
  error = (ierr /= 0)
  call dgetri(size(A,1),InvA,size(A,1),ipiv,work,size(work,1),ierr)
  error = error .or. (ierr /= 0)

end subroutine math_invert


!--------------------------------------------------------------------------------------------------
!> @brief symmetrize a 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_symmetric33(m)

  real(pReal), dimension(3,3) :: math_symmetric33
  real(pReal), dimension(3,3), intent(in) :: m


  math_symmetric33 = 0.5_pReal * (m + transpose(m))

end function math_symmetric33


!--------------------------------------------------------------------------------------------------
!> @brief skew part of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_skew33(m)

  real(pReal), dimension(3,3) :: math_skew33
  real(pReal), dimension(3,3), intent(in) :: m


  math_skew33 = m - math_symmetric33(m)

end function math_skew33


!--------------------------------------------------------------------------------------------------
!> @brief hydrostatic part of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_spherical33(m)

  real(pReal), dimension(3,3) :: math_spherical33
  real(pReal), dimension(3,3), intent(in) :: m


  math_spherical33 = math_I3 * math_trace33(m)/3.0_pReal

end function math_spherical33


!--------------------------------------------------------------------------------------------------
!> @brief deviatoric part of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_deviatoric33(m)

  real(pReal), dimension(3,3) :: math_deviatoric33
  real(pReal), dimension(3,3), intent(in) :: m


  math_deviatoric33 = m - math_spherical33(m)

end function math_deviatoric33


!--------------------------------------------------------------------------------------------------
!> @brief trace of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_trace33(m)

  real(pReal), dimension(3,3), intent(in) :: m


  math_trace33 = m(1,1) + m(2,2) + m(3,3)

end function math_trace33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_det33(m)

  real(pReal), dimension(3,3), intent(in) :: m


  math_det33 = m(1,1)* (m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
             - m(1,2)* (m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
             + m(1,3)* (m(2,1)*m(3,2)-m(2,2)*m(3,1))

end function math_det33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a symmetric 3x3 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_detSym33(m)

  real(pReal), dimension(3,3), intent(in) :: m


  math_detSym33 = -(m(1,1)*m(2,3)**2 + m(2,2)*m(1,3)**2 + m(3,3)*m(1,2)**2) &
                  + m(1,1)*m(2,2)*m(3,3) + 2.0_pReal * m(1,2)*m(1,3)*m(2,3)

end function  math_detSym33


!--------------------------------------------------------------------------------------------------
!> @brief convert 3x3 matrix into vector 9
!--------------------------------------------------------------------------------------------------
pure function math_33to9(m33)

  real(pReal), dimension(9)               :: math_33to9
  real(pReal), dimension(3,3), intent(in) :: m33

  integer :: i


  math_33to9 = [(m33(MAPPLAIN(1,i),MAPPLAIN(2,i)),i=1,9)]

end function math_33to9


!--------------------------------------------------------------------------------------------------
!> @brief convert 9 vector into 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_9to33(v9)

  real(pReal), dimension(3,3)           :: math_9to33
  real(pReal), dimension(9), intent(in) :: v9

  integer :: i


  do i = 1, 9
    math_9to33(MAPPLAIN(1,i),MAPPLAIN(2,i)) = v9(i)
  end do

end function math_9to33


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 3x3 matrix into 6 vector
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_sym33to6(m33,weighted)

  real(pReal), dimension(6)               :: math_sym33to6
  real(pReal), dimension(3,3), intent(in) :: m33                                                    !< symmetric 3x3 matrix (no internal check)
  logical,     optional,       intent(in) :: weighted                                               !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i


  if (present(weighted)) then
    w = merge(NRMMANDEL,1.0_pReal,weighted)
  else
    w = NRMMANDEL
  end if

  math_sym33to6 = [(w(i)*m33(MAPNYE(1,i),MAPNYE(2,i)),i=1,6)]

end function math_sym33to6


!--------------------------------------------------------------------------------------------------
!> @brief convert 6 vector into symmetric 3x3 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_6toSym33(v6,weighted)

  real(pReal), dimension(3,3)           :: math_6toSym33
  real(pReal), dimension(6), intent(in) :: v6                                                       !< 6 vector
  logical,     optional,     intent(in) :: weighted                                                 !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i


  if (present(weighted)) then
    w = merge(INVNRMMANDEL,1.0_pReal,weighted)
  else
    w = INVNRMMANDEL
  end if

  do i=1,6
    math_6toSym33(MAPNYE(1,i),MAPNYE(2,i)) = w(i)*v6(i)
    math_6toSym33(MAPNYE(2,i),MAPNYE(1,i)) = w(i)*v6(i)
  end do

end function math_6toSym33


!--------------------------------------------------------------------------------------------------
!> @brief convert 3x3x3x3 matrix into 9x9 matrix
!--------------------------------------------------------------------------------------------------
pure function math_3333to99(m3333)

  real(pReal), dimension(9,9)                 :: math_3333to99
  real(pReal), dimension(3,3,3,3), intent(in) :: m3333

  integer :: i,j


#ifndef __INTEL_COMPILER
  do concurrent(i=1:9, j=1:9)
    math_3333to99(i,j) = m3333(MAPPLAIN(1,i),MAPPLAIN(2,i),MAPPLAIN(1,j),MAPPLAIN(2,j))
  end do
#else
  forall(i=1:9, j=1:9) math_3333to99(i,j) = m3333(MAPPLAIN(1,i),MAPPLAIN(2,i),MAPPLAIN(1,j),MAPPLAIN(2,j))
#endif

end function math_3333to99


!--------------------------------------------------------------------------------------------------
!> @brief convert 9x9 matrix into 3x3x3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_99to3333(m99)

  real(pReal), dimension(3,3,3,3)         :: math_99to3333
  real(pReal), dimension(9,9), intent(in) :: m99

  integer :: i,j


#ifndef __INTEL_COMPILER
  do concurrent(i=1:9, j=1:9)
    math_99to3333(MAPPLAIN(1,i),MAPPLAIN(2,i),MAPPLAIN(1,j),MAPPLAIN(2,j)) = m99(i,j)
  end do
#else
  forall(i=1:9, j=1:9) math_99to3333(MAPPLAIN(1,i),MAPPLAIN(2,i),MAPPLAIN(1,j),MAPPLAIN(2,j)) = m99(i,j)
#endif

end function math_99to3333


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 3x3x3x3 matrix into 6x6 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only rearranges order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_sym3333to66(m3333,weighted)

  real(pReal), dimension(6,6)                 :: math_sym3333to66
  real(pReal), dimension(3,3,3,3), intent(in) :: m3333                                              !< symmetric 3x3x3x3 matrix (no internal check)
  logical,     optional,           intent(in) :: weighted                                           !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i,j


  if (present(weighted)) then
    w = merge(NRMMANDEL,1.0_pReal,weighted)
  else
    w = NRMMANDEL
  end if

#ifndef __INTEL_COMPILER
  do concurrent(i=1:6, j=1:6)
    math_sym3333to66(i,j) = w(i)*w(j)*m3333(MAPNYE(1,i),MAPNYE(2,i),MAPNYE(1,j),MAPNYE(2,j))
  end do
#else
  forall(i=1:6, j=1:6) math_sym3333to66(i,j) = w(i)*w(j)*m3333(MAPNYE(1,i),MAPNYE(2,i),MAPNYE(1,j),MAPNYE(2,j))
#endif

end function math_sym3333to66


!--------------------------------------------------------------------------------------------------
!> @brief convert 6x6 matrix into symmetric 3x3x3x3 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only rearranges order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_66toSym3333(m66,weighted)

  real(pReal), dimension(3,3,3,3)            :: math_66toSym3333
  real(pReal), dimension(6,6),    intent(in) :: m66                                                 !< 6x6 matrix
  logical,     optional,          intent(in) :: weighted                                            !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i,j


  if (present(weighted)) then
    w = merge(INVNRMMANDEL,1.0_pReal,weighted)
  else
    w = INVNRMMANDEL
  end if

  do i=1,6; do j=1,6
    math_66toSym3333(MAPNYE(1,i),MAPNYE(2,i),MAPNYE(1,j),MAPNYE(2,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(MAPNYE(2,i),MAPNYE(1,i),MAPNYE(1,j),MAPNYE(2,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(MAPNYE(1,i),MAPNYE(2,i),MAPNYE(2,j),MAPNYE(1,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(MAPNYE(2,i),MAPNYE(1,i),MAPNYE(2,j),MAPNYE(1,j)) = w(i)*w(j)*m66(i,j)
  end do; end do

end function math_66toSym3333


!--------------------------------------------------------------------------------------------------
!> @brief Convert 6 Voigt stress vector into symmetric 3x3 tensor.
!--------------------------------------------------------------------------------------------------
pure function math_Voigt6to33_stress(sigma_tilde) result(sigma)

  real(pReal), dimension(3,3) :: sigma
  real(pReal), dimension(6), intent(in) :: sigma_tilde


  sigma = reshape([sigma_tilde(1), sigma_tilde(6), sigma_tilde(5), &
                   sigma_tilde(6), sigma_tilde(2), sigma_tilde(4), &
                   sigma_tilde(5), sigma_tilde(4), sigma_tilde(3)],[3,3])

end function math_Voigt6to33_stress


!--------------------------------------------------------------------------------------------------
!> @brief Convert 6 Voigt strain vector into symmetric 3x3 tensor.
!--------------------------------------------------------------------------------------------------
pure function math_Voigt6to33_strain(epsilon_tilde) result(epsilon)

  real(pReal), dimension(3,3) :: epsilon
  real(pReal), dimension(6), intent(in) :: epsilon_tilde


  epsilon = reshape([          epsilon_tilde(1), 0.5_pReal*epsilon_tilde(6), 0.5_pReal*epsilon_tilde(5), &
                     0.5_pReal*epsilon_tilde(6),           epsilon_tilde(2), 0.5_pReal*epsilon_tilde(4), &
                     0.5_pReal*epsilon_tilde(5), 0.5_pReal*epsilon_tilde(4),           epsilon_tilde(3)],[3,3])

end function math_Voigt6to33_strain


!--------------------------------------------------------------------------------------------------
!> @brief Convert 3x3 stress tensor into 6 Voigt vector.
!--------------------------------------------------------------------------------------------------
pure function math_33toVoigt6_stress(sigma) result(sigma_tilde)

  real(pReal), dimension(6) :: sigma_tilde
  real(pReal), dimension(3,3), intent(in) :: sigma


  sigma_tilde = [sigma(1,1), sigma(2,2), sigma(3,3), &
                 sigma(3,2), sigma(3,1), sigma(1,2)]

end function math_33toVoigt6_stress


!--------------------------------------------------------------------------------------------------
!> @brief Convert 3x3 strain tensor into 6 Voigt vector.
!--------------------------------------------------------------------------------------------------
pure function math_33toVoigt6_strain(epsilon) result(epsilon_tilde)

  real(pReal), dimension(6) :: epsilon_tilde
  real(pReal), dimension(3,3), intent(in) :: epsilon


  epsilon_tilde = [          epsilon(1,1),           epsilon(2,2),           epsilon(3,3), &
                   2.0_pReal*epsilon(3,2), 2.0_pReal*epsilon(3,1), 2.0_pReal*epsilon(1,2)]

end function math_33toVoigt6_strain



!--------------------------------------------------------------------------------------------------
!> @brief Convert 6x6 Voigt stiffness matrix into symmetric 3x3x3x3 tensor.
!--------------------------------------------------------------------------------------------------
pure function math_Voigt66to3333_stiffness(C_tilde) result(C)

  real(pReal), dimension(3,3,3,3) :: C
  real(pReal), dimension(6,6), intent(in) :: C_tilde

  integer :: i,j


  do i=1,6; do j=1,6
    C(MAPVOIGT(1,i),MAPVOIGT(2,i),MAPVOIGT(1,j),MAPVOIGT(2,j)) = C_tilde(i,j)
    C(MAPVOIGT(2,i),MAPVOIGT(1,i),MAPVOIGT(1,j),MAPVOIGT(2,j)) = C_tilde(i,j)
    C(MAPVOIGT(1,i),MAPVOIGT(2,i),MAPVOIGT(2,j),MAPVOIGT(1,j)) = C_tilde(i,j)
    C(MAPVOIGT(2,i),MAPVOIGT(1,i),MAPVOIGT(2,j),MAPVOIGT(1,j)) = C_tilde(i,j)
  end do; end do

end function math_Voigt66to3333_stiffness


!--------------------------------------------------------------------------------------------------
!> @brief Convert 3x3x3x3 stiffness tensor into 6x6 Voigt matrix.
!--------------------------------------------------------------------------------------------------
pure function math_3333toVoigt66_stiffness(C) result(C_tilde)

  real(pReal), dimension(6,6) :: C_tilde
  real(pReal), dimension(3,3,3,3), intent(in) :: C

  integer :: i,j


#ifndef __INTEL_COMPILER
  do concurrent(i=1:6, j=1:6)
    C_tilde(i,j) = C(MAPVOIGT(1,i),MAPVOIGT(2,i),MAPVOIGT(1,j),MAPVOIGT(2,j))
  end do
#else
  forall(i=1:6, j=1:6) C_tilde(i,j) = C(MAPVOIGT(1,i),MAPVOIGT(2,i),MAPVOIGT(1,j),MAPVOIGT(2,j))
#endif

end function math_3333toVoigt66_stiffness


!--------------------------------------------------------------------------------------------------
!> @brief Draw a sample from a normal distribution.
!> @details https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
!>          https://masuday.github.io/fortran_tutorial/random.html
!--------------------------------------------------------------------------------------------------
impure elemental subroutine math_normal(x,mu,sigma)

  real(pReal), intent(out) :: x
  real(pReal), intent(in), optional :: mu, sigma

  real(pReal) :: sigma_, mu_
  real(pReal), dimension(2) :: rnd


  if (present(mu)) then
    mu_ = mu
  else
    mu_ = 0.0_pReal
  end if

  if (present(sigma)) then
    sigma_ = sigma
  else
    sigma_ = 1.0_pReal
  end if

  call random_number(rnd)
  x = mu_ + sigma_ * sqrt(-2.0_pReal*log(1.0_pReal-rnd(1)))*cos(TAU*(1.0_pReal - rnd(2)))

end subroutine math_normal


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric matrix
!--------------------------------------------------------------------------------------------------
pure subroutine math_eigh(w,v,error,m)

  real(pReal), dimension(:,:),                  intent(in)  :: m                                    !< quadratic matrix to compute eigenvectors and values of
  real(pReal), dimension(size(m,1)),            intent(out) :: w                                    !< eigenvalues
  real(pReal), dimension(size(m,1),size(m,1)),  intent(out) :: v                                    !< eigenvectors
  logical,                                      intent(out) :: error

  integer :: ierr
  real(pReal), dimension(size(m,1)**2) :: work


  v = m                                                                                             ! copy matrix to input (doubles as output) array
  call dsyev('V','U',size(m,1),v,size(m,1),w,work,size(work,1),ierr)
  error = (ierr /= 0)

end subroutine math_eigh


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric 3x3 matrix using an analytical expression
!> and the general LAPACK powered version for arbritrary sized matrices as fallback
!> @author Joachim Kopp, Max-Planck-Institut für Kernphysik, Heidelberg (Copyright (C) 2006)
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details See http://arxiv.org/abs/physics/0610206 (DSYEVH3)
!--------------------------------------------------------------------------------------------------
pure subroutine math_eigh33(w,v,m)

  real(pReal), dimension(3,3),intent(in)  :: m                                                      !< 3x3 matrix to compute eigenvectors and values of
  real(pReal), dimension(3),  intent(out) :: w                                                      !< eigenvalues
  real(pReal), dimension(3,3),intent(out) :: v                                                      !< eigenvectors

  real(pReal) :: T, U, norm, threshold
  logical :: error


  w = math_eigvalsh33(m)

  v(1:3,2) = [ m(1,2) * m(2,3) - m(1,3) * m(2,2), &
               m(1,3) * m(1,2) - m(2,3) * m(1,1), &
               m(1,2)**2]

  T = maxval(abs(w))
  U = max(T, T**2)
  threshold = sqrt(5.68e-14_pReal * U**2)

#ifndef __INTEL_LLVM_COMPILER
  v(1:3,1) = [m(1,3)*w(1) + v(1,2), &
              m(2,3)*w(1) + v(2,2), &
#else
  v(1:3,1) = [IEEE_FMA(m(1,3),w(1),v(1,2)), &
              IEEE_FMA(m(2,3),w(1),v(2,2)), &
#endif
              (m(1,1) - w(1)) * (m(2,2) - w(1)) - v(3,2)]
  norm = norm2(v(1:3, 1))
  fallback1: if (norm < threshold) then
    call math_eigh(w,v,error,m)
  else fallback1
    v(1:3,1) = v(1:3, 1) / norm
#ifndef __INTEL_LLVM_COMPILER
    v(1:3,2) = [m(1,3)*w(2) + v(1,2), &
                m(2,3)*w(2) + v(2,2), &
#else
    v(1:3,2) = [IEEE_FMA(m(1,3),w(2),v(1,2)), &
                IEEE_FMA(m(2,3),w(2),v(2,2)), &
#endif
                (m(1,1) - w(2)) * (m(2,2) - w(2)) - v(3,2)]
    norm = norm2(v(1:3, 2))
    fallback2: if (norm < threshold) then
      call math_eigh(w,v,error,m)
    else fallback2
      v(1:3,2) = v(1:3, 2) / norm
      v(1:3,3) = math_cross(v(1:3,1),v(1:3,2))
     end if fallback2
  end if fallback1

end subroutine math_eigh33


!--------------------------------------------------------------------------------------------------
!> @brief Calculate rotational part of a deformation gradient
!> @details https://www.jstor.org/stable/43637254
!!          https://www.jstor.org/stable/43637372
!!          https://doi.org/10.1023/A:1007407802076
!--------------------------------------------------------------------------------------------------
pure function math_rotationalPart(F) result(R)

  real(pReal), dimension(3,3), intent(in) :: &
    F                                                                                               ! deformation gradient
  real(pReal), dimension(3,3) :: &
    C, &                                                                                            ! right Cauchy-Green tensor
    R                                                                                               ! rotational part
  real(pReal), dimension(3) :: &
    lambda, &                                                                                       ! principal stretches
    I_C, &                                                                                          ! invariants of C
    I_U                                                                                             ! invariants of U
  real(pReal), dimension(2) :: &
    I_F                                                                                             ! first two invariants of F
  real(pReal) :: x,Phi


  C = matmul(transpose(F),F)
  I_C = math_invariantsSym33(C)
  I_F = [math_trace33(F), 0.5_pReal*(math_trace33(F)**2 - math_trace33(matmul(F,F)))]

  x = math_clip(I_C(1)**2 -3.0_pReal*I_C(2),0.0_pReal)**(3.0_pReal/2.0_pReal)
  if (dNeq0(x)) then
    Phi = acos(math_clip((I_C(1)**3 -4.5_pReal*I_C(1)*I_C(2) +13.5_pReal*I_C(3))/x,-1.0_pReal,1.0_pReal))
    lambda = I_C(1) +(2.0_pReal * sqrt(math_clip(I_C(1)**2-3.0_pReal*I_C(2),0.0_pReal))) &
                    *cos((Phi-TAU*[1.0_pReal,2.0_pReal,3.0_pReal])/3.0_pReal)
    lambda = sqrt(math_clip(lambda,0.0_pReal)/3.0_pReal)
  else
    lambda = sqrt(I_C(1)/3.0_pReal)
  end if

  I_U = [sum(lambda), lambda(1)*lambda(2)+lambda(2)*lambda(3)+lambda(3)*lambda(1), product(lambda)]

  R = I_U(1)*I_F(2) * math_I3 &
    +(I_U(1)**2-I_U(2)) * F &
    - I_U(1)*I_F(1) * transpose(F) &
    + I_U(1) * transpose(matmul(F,F)) &
    - matmul(F,C)
  R = R*math_det33(R)**(-1.0_pReal/3.0_pReal)

end function math_rotationalPart


!--------------------------------------------------------------------------------------------------
!> @brief Eigenvalues of symmetric matrix
! will return NaN on error
!--------------------------------------------------------------------------------------------------
pure function math_eigvalsh(m)

  real(pReal), dimension(:,:),                  intent(in)  :: m                                    !< symmetric matrix to compute eigenvalues of
  real(pReal), dimension(size(m,1))                         :: math_eigvalsh

  real(pReal), dimension(size(m,1),size(m,1))               :: m_
  integer :: ierr
  real(pReal), dimension(size(m,1)**2) :: work


  m_= m                                                                                             ! copy matrix to input (will be destroyed)
  call dsyev('N','U',size(m,1),m_,size(m,1),math_eigvalsh,work,size(work,1),ierr)
  if (ierr /= 0) math_eigvalsh = IEEE_value(1.0_pReal,IEEE_quiet_NaN)

end function math_eigvalsh


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues of symmetric 3x3 matrix using an analytical expression
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details similar to http://arxiv.org/abs/physics/0610206 (DSYEVC3)
!> but apparently more stable solution and has general LAPACK powered version for arbritrary sized
!> matrices as fallback
!--------------------------------------------------------------------------------------------------
pure function math_eigvalsh33(m)

  real(pReal), intent(in), dimension(3,3) :: m                                                      !< 3x3 symmetric matrix to compute eigenvalues of
  real(pReal), dimension(3) :: math_eigvalsh33,I
  real(pReal) :: P, Q, rho, phi
  real(pReal), parameter :: TOL=1.e-14_pReal


  I = math_invariantsSym33(m)                                                                       ! invariants are coefficients in characteristic polynomial apart for the sign of c0 and c2 in http://arxiv.org/abs/physics/0610206

  P = I(2)-I(1)**2/3.0_pReal                                                                        ! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)
  Q = product(I(1:2))/3.0_pReal &
    - 2.0_pReal/27.0_pReal*I(1)**3 &
    - I(3)                                                                                          ! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)

  if (all(abs([P,Q]) < TOL)) then
    math_eigvalsh33 = math_eigvalsh(m)
  else
    rho=sqrt(-3.0_pReal*P**3)/9.0_pReal
    phi=acos(math_clip(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
    math_eigvalsh33 = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                                                            [cos( phi              /3.0_pReal), &
                                                             cos((phi+TAU)/3.0_pReal), &
                                                             cos((phi+2.0_pReal*TAU)/3.0_pReal) &
                                                            ] &
                    + I(1)/3.0_pReal
  end if

end function math_eigvalsh33


!--------------------------------------------------------------------------------------------------
!> @brief invariants of symmetrix 3x3 matrix
!--------------------------------------------------------------------------------------------------
pure function math_invariantsSym33(m)

  real(pReal), dimension(3,3), intent(in) :: m
  real(pReal), dimension(3) :: math_invariantsSym33


  math_invariantsSym33(1) = math_trace33(m)
  math_invariantsSym33(2) = m(1,1)*m(2,2) + m(1,1)*m(3,3) + m(2,2)*m(3,3) &
                          -(m(1,2)**2     + m(1,3)**2     + m(2,3)**2)
  math_invariantsSym33(3) = math_detSym33(m)

end function math_invariantsSym33


!--------------------------------------------------------------------------------------------------
!> @brief factorial
!--------------------------------------------------------------------------------------------------
integer pure function math_factorial(n)

  integer, intent(in) :: n


  math_factorial = product(math_range(n))

end function math_factorial


!--------------------------------------------------------------------------------------------------
!> @brief binomial coefficient
!--------------------------------------------------------------------------------------------------
integer pure function math_binomial(n,k)

  integer, intent(in) :: n, k
  integer :: i, k_, n_


  k_ = min(k,n-k)
  n_ = n
  math_binomial = merge(1,0,k_>-1)                                                                  ! handling special cases k < 0 or k > n
  do i = 1, k_
    math_binomial = (math_binomial * n_)/i
    n_ = n_ -1
  end do

end function math_binomial


!--------------------------------------------------------------------------------------------------
!> @brief multinomial coefficient
!--------------------------------------------------------------------------------------------------
integer pure function math_multinomial(k)

  integer, intent(in), dimension(:) :: k
  integer :: i


  math_multinomial = product([(math_binomial(sum(k(1:i)),k(i)),i=1,size(k))])

end function math_multinomial


!--------------------------------------------------------------------------------------------------
!> @brief volume of tetrahedron given by four vertices
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_volTetrahedron(v1,v2,v3,v4)

  real(pReal), dimension (3), intent(in) :: v1,v2,v3,v4
  real(pReal), dimension (3,3) :: m


  m(1:3,1) = v1-v2
  m(1:3,2) = v1-v3
  m(1:3,3) = v1-v4

  math_volTetrahedron = abs(math_det33(m))/6.0_pReal

end function math_volTetrahedron


!--------------------------------------------------------------------------------------------------
!> @brief area of triangle given by three vertices
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_areaTriangle(v1,v2,v3)

  real(pReal), dimension (3), intent(in) :: v1,v2,v3


  math_areaTriangle = 0.5_pReal * norm2(math_cross(v1-v2,v1-v3))

end function math_areaTriangle


!--------------------------------------------------------------------------------------------------
!> @brief Limit a scalar value to a certain range (either one or two sided).
!--------------------------------------------------------------------------------------------------
real(pReal) pure elemental function math_clip(a, left, right)

  real(pReal), intent(in) :: a
  real(pReal), intent(in), optional :: left, right


  math_clip = a
  if (present(left))  math_clip = max(left,math_clip)
  if (present(right)) math_clip = min(right,math_clip)
  if (present(left) .and. present(right)) then
    if (left>right) error stop 'left > right'
  end if

end function math_clip


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some math functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  integer, dimension(2,4) :: &
    sort_in_   = reshape([+1,+5,  +5,+6,  -1,-1,  +3,-2],[2,4])
  integer, dimension(2,4), parameter :: &
    sort_out_  = reshape([-1,-1,  +1,+5,  +5,+6,  +3,-2],[2,4])

  integer,     dimension(5) :: range_out_ = [1,2,3,4,5]
  integer,     dimension(3) :: ijk

  real(pReal)                 :: det
  real(pReal), dimension(3)   :: v3_1,v3_2,v3_3,v3_4
  real(pReal), dimension(6)   :: v6
  real(pReal), dimension(9)   :: v9
  real(pReal), dimension(3,3) :: t33,t33_2
  real(pReal), dimension(6,6) :: t66
  real(pReal), dimension(9,9) :: t99,t99_2
  real(pReal), dimension(:,:), &
               allocatable    :: txx,txx_2
  real(pReal)                 :: r
  integer                     :: d
  logical                     :: e


  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,3.0_pReal,3.0_pReal,3.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1,2,3,0])) > tol_math_check)) &
    error stop 'math_expand [1,2,3] by [1,2,3,0] => [1,2,2,3,3,3]'

  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1,2])) > tol_math_check)) &
    error stop 'math_expand [1,2,3] by [1,2] => [1,2,2]'

  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,1.0_pReal,1.0_pReal,1.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal],[1,2,3])) > tol_math_check)) &
    error stop 'math_expand [1,2] by [1,2,3] => [1,2,2,1,1,1]'

  call math_sort(sort_in_,1,3,2)
  if (any(sort_in_ /= sort_out_)) &
    error stop 'math_sort'

  if (any(math_range(5) /= range_out_)) &
    error stop 'math_range'

  if (any(dNeq(math_exp33(math_I3,0),math_I3))) &
    error stop 'math_exp33(math_I3,1)'
  if (any(dNeq(math_exp33(math_I3,128),exp(1.0_pReal)*math_I3))) &
    error stop 'math_exp33(math_I3,128)'

  call random_number(v9)
  if (any(dNeq(math_33to9(math_9to33(v9)),v9))) &
    error stop 'math_33to9/math_9to33'

  call random_number(t99)
  if (any(dNeq(math_3333to99(math_99to3333(t99)),t99))) &
    error stop 'math_3333to99/math_99to3333'

  call random_number(v6)
  if (any(dNeq(math_sym33to6(math_6toSym33(v6)),v6))) &
    error stop 'math_sym33to6/math_6toSym33'

  call random_number(t66)
  if (any(dNeq(math_sym3333to66(math_66toSym3333(t66)),t66,1.0e-15_pReal))) &
    error stop 'math_sym3333to66/math_66toSym3333'

  if (any(dNeq(math_3333toVoigt66_stiffness(math_Voigt66to3333_stiffness(t66)),t66,1.0e-15_pReal))) &
    error stop 'math_3333toVoigt66/math_Voigt66to3333'

  call random_number(v6)
  if (any(dNeq0(math_6toSym33(v6) - math_symmetric33(math_6toSym33(v6))))) &
    error stop 'math_symmetric33'

  call random_number(v3_1)
  call random_number(v3_2)
  call random_number(v3_3)
  call random_number(v3_4)

  if (dNeq(abs(dot_product(math_cross(v3_1-v3_4,v3_2-v3_4),v3_3-v3_4))/6.0_pReal, &
          math_volTetrahedron(v3_1,v3_2,v3_3,v3_4),tol=1.0e-12_pReal)) &
    error stop 'math_volTetrahedron'

  call random_number(t33)
  if (dNeq(math_det33(math_symmetric33(t33)),math_detSym33(math_symmetric33(t33)),tol=1.0e-12_pReal)) &
    error stop 'math_det33/math_detSym33'

  if (any(dNeq(t33+transpose(t33),math_mul3333xx33(math_identity4th(),t33+transpose(t33))))) &
    error stop 'math_mul3333xx33/math_identity4th'

  if (any(dNeq0(math_eye(3),math_inv33(math_I3)))) &
    error stop 'math_inv33(math_I3)'

  do while(abs(math_det33(t33))<1.0e-9_pReal)
    call random_number(t33)
  end do
  if (any(dNeq0(matmul(t33,math_inv33(t33)) - math_eye(3),tol=1.0e-9_pReal))) &
    error stop 'math_inv33'

  call math_invert33(t33_2,det,e,t33)
  if (any(dNeq0(matmul(t33,t33_2) - math_eye(3),tol=1.0e-9_pReal)) .or. e) &
    error stop 'math_invert33: T:T^-1 != I'
  if (dNeq(det,math_det33(t33),tol=1.0e-12_pReal)) &
    error stop 'math_invert33 (determinant)'

  call math_invert(t33_2,e,t33)
  if (any(dNeq0(matmul(t33,t33_2) - math_eye(3),tol=1.0e-9_pReal)) .or. e) &
    error stop 'math_invert t33'

  do while(math_det33(t33)<1.0e-2_pReal)                                                            ! O(det(F)) = 1
    call random_number(t33)
  end do
  t33_2 = math_rotationalPart(transpose(t33))
  t33   = math_rotationalPart(t33)
  if (any(dNeq0(matmul(t33_2,t33) - math_I3,tol=1.0e-10_pReal))) &
    error stop 'math_rotationalPart (forward-backward)'
  if (dNeq(1.0_pReal,math_det33(math_rotationalPart(t33)),tol=1.0e-10_pReal)) &
    error stop 'math_rotationalPart (determinant)'

  call random_number(r)
  d = int(r*5.0_pReal) + 1
  txx = math_eye(d)
  allocate(txx_2(d,d))
  call math_invert(txx_2,e,txx)
  if (any(dNeq0(txx_2,txx) .or. e)) &
    error stop 'math_invert(txx)/math_eye'

  call math_invert(t99_2,e,t99) ! not sure how likely it is that we get a singular matrix
  if (any(dNeq0(matmul(t99_2,t99)-math_eye(9),tol=1.0e-9_pReal)) .or. e) &
    error stop 'math_invert(t99)'

  if (any(dNeq(math_clip([4.0_pReal,9.0_pReal],5.0_pReal,6.5_pReal),[5.0_pReal,6.5_pReal]))) &
    error stop 'math_clip'

  if (math_factorial(10) /= 3628800) &
    error stop 'math_factorial'

  if (math_binomial(49,6) /= 13983816) &
    error stop 'math_binomial'

  if (math_multinomial([1,2,3,4]) /= 12600) &
    error stop 'math_multinomial'

  ijk = cshift([1,2,3],int(r*1.0e2_pReal))
  if (dNeq(math_LeviCivita(ijk(1),ijk(2),ijk(3)),+1.0_pReal)) &
    error stop 'math_LeviCivita(even)'
  ijk = cshift([3,2,1],int(r*2.0e2_pReal))
  if (dNeq(math_LeviCivita(ijk(1),ijk(2),ijk(3)),-1.0_pReal)) &
    error stop 'math_LeviCivita(odd)'
  ijk = cshift([2,2,1],int(r*2.0e2_pReal))
  if (dNeq0(math_LeviCivita(ijk(1),ijk(2),ijk(3)))) &
    error stop 'math_LeviCivita'

  normal_distribution: block
    integer, parameter :: N = 1000000
    real(pReal), dimension(:), allocatable :: r
    real(pReal) :: mu, sigma

    allocate(r(N))
    call random_number(mu)
    call random_number(sigma)

    sigma = 1.0_pReal + sigma*5.0_pReal
    mu = (mu-0.5_pReal)*10_pReal

    call math_normal(r,mu,sigma)

    if (abs(mu -sum(r)/real(N,pReal))>5.0e-2_pReal) &
      error stop 'math_normal(mu)'

    mu = sum(r)/real(N,pReal)
    if (abs(sigma**2 -1.0_pReal/real(N-1,pReal) * sum((r-mu)**2))/sigma > 5.0e-2_pReal) &
      error stop 'math_normal(sigma)'
  end block normal_distribution

end subroutine selfTest

end module math
