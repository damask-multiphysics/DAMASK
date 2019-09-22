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
  use debug
  use numerics
  implicit none
  public
#if __INTEL_COMPILER >= 1900
  ! do not make use associated entities available to other modules
  private :: &
    prec
#endif

  real(pReal),    parameter :: PI = acos(-1.0_pReal)                                                !< ratio of a circle's circumference to its diameter
  real(pReal),    parameter :: INDEG = 180.0_pReal/PI                                               !< conversion from radian into degree
  real(pReal),    parameter :: INRAD = PI/180.0_pReal                                               !< conversion from degree into radian
  complex(pReal), parameter :: TWOPIIMG = cmplx(0.0_pReal,2.0_pReal*PI)                             !< Re(0.0), Im(2xPi)

  real(pReal), dimension(3,3), parameter :: &
    MATH_I3 = reshape([&
      1.0_pReal,0.0_pReal,0.0_pReal, &
      0.0_pReal,1.0_pReal,0.0_pReal, &
      0.0_pReal,0.0_pReal,1.0_pReal  &
      ],[3,3])                                                                                      !< 3x3 Identity

  real(pReal), dimension(6), parameter, private :: &
    nrmMandel = [&
      1.0_pReal,       1.0_pReal,       1.0_pReal, &
      sqrt(2.0_pReal), sqrt(2.0_pReal), sqrt(2.0_pReal) ]                                           !< weighting for Mandel notation (forward)

  real(pReal), dimension(6), parameter , private :: &
    invnrmMandel = [&
      1.0_pReal,                 1.0_pReal,                 1.0_pReal, &
      1.0_pReal/sqrt(2.0_pReal), 1.0_pReal/sqrt(2.0_pReal), 1.0_pReal/sqrt(2.0_pReal) ]             !< weighting for Mandel notation (backward)

  integer, dimension (2,6), parameter, private :: &
    mapNye = reshape([&
      1,1, &
      2,2, &
      3,3, &
      1,2, &
      2,3, &
      1,3  &
      ],[2,6])                                                                                      !< arrangement in Nye notation.

  integer, dimension (2,6), parameter, private :: &
    mapVoigt = reshape([&
      1,1, &
      2,2, &
      3,3, &
      2,3, &
      1,3, &
      1,2  &
      ],[2,6])                                                                                      !< arrangement in Voigt notation

  integer, dimension (2,9), parameter, private :: &
    mapPlain = reshape([&
      1,1, &
      1,2, &
      1,3, &
      2,1, &
      2,2, &
      2,3, &
      3,1, &
      3,2, &
      3,3  &
      ],[2,9])                                                                                      !< arrangement in Plain notation
     
!---------------------------------------------------------------------------------------------------
 private :: &
   unitTest

contains

!--------------------------------------------------------------------------------------------------
!> @brief initialization of random seed generator and internal checks
!--------------------------------------------------------------------------------------------------
subroutine math_init

  integer :: i
  real(pReal), dimension(4) :: randTest
  integer :: randSize
  integer, dimension(:), allocatable :: randInit
  
  write(6,'(/,a)')   ' <<<+-  math init  -+>>>'

  call random_seed(size=randSize)
  allocate(randInit(randSize))
  if (randomSeed > 0) then
    randInit = randomSeed
    call random_seed(put=randInit)
  else
    call random_seed()
    call random_seed(get = randInit)
    randInit(2:randSize) = randInit(1)
    call random_seed(put = randInit)
  endif

  do i = 1, 4
    call random_number(randTest(i))
  enddo

  write(6,'(a,I2)') ' size  of random seed:    ', randSize
  do i = 1,randSize
    write(6,'(a,I2,I14)') ' value of random seed:    ', i, randInit(i)
  enddo
  write(6,'(a,4(/,26x,f17.14),/)') ' start of random sequence: ', randTest

  call random_seed(put = randInit)
  call unitTest

end subroutine math_init


!--------------------------------------------------------------------------------------------------
!> @brief Quicksort algorithm for two-dimensional integer arrays
! Sorting is done with respect to array(sort,:) and keeps array(/=sort,:) linked to it.
! default: sort=1
!--------------------------------------------------------------------------------------------------
recursive subroutine math_sort(a, istart, iend, sortDim)

  integer, dimension(:,:), intent(inout) :: a
  integer, intent(in),optional :: istart,iend, sortDim
  integer :: ipivot,s,e,d
  
  if(present(istart)) then
    s = istart
  else
    s = lbound(a,2)
  endif
  
  if(present(iend)) then
    e = iend
  else
    e = ubound(a,2)
  endif
  
  if(present(sortDim)) then
    d = sortDim
  else
    d = 1
  endif
   
  if (s < e) then
    ipivot = qsort_partition(a,s, e, d)
    call math_sort(a, s, ipivot-1, d)
    call math_sort(a, ipivot+1, e, d)
  endif


  contains

  !-------------------------------------------------------------------------------------------------
  !> @brief Partitioning required for quicksort
  !-------------------------------------------------------------------------------------------------
  integer function qsort_partition(a, istart, iend, sort)
  
    integer, dimension(:,:), intent(inout) :: a
    integer,                 intent(in)    :: istart,iend,sort
    integer, dimension(size(a,1))          :: tmp
    integer :: i,j
   
    do
      ! find the first element on the right side less than or equal to the pivot point
      do j = iend, istart, -1
        if (a(sort,j) <= a(sort,istart)) exit
      enddo
      ! find the first element on the left side greater than the pivot point
      do i = istart, iend
        if (a(sort,i) > a(sort,istart)) exit
      enddo
      cross: if (i >= j) then ! if the indices cross, exchange left value with pivot and return with the partition index
        tmp         = a(:,istart)
        a(:,istart) = a(:,j)
        a(:,j)      = tmp
        qsort_partition = j
        return
      else cross ! if they do not cross, exchange values
        tmp    = a(:,i)
        a(:,i) = a(:,j)
        a(:,j) = tmp
      endif cross
    enddo
  
  end function qsort_partition

end subroutine math_sort


!--------------------------------------------------------------------------------------------------
!> @brief vector expansion
!> @details takes a set of numbers (a,b,c,...) and corresponding multiples (x,y,z,...)
!> to return a vector of x times a, y times b, z times c, ... 
!--------------------------------------------------------------------------------------------------
pure function math_expand(what,how)

  real(pReal),   dimension(:), intent(in) :: what
  integer, dimension(:), intent(in) :: how
  real(pReal), dimension(sum(how)) ::  math_expand
  integer :: i
 
  if (sum(how) == 0) return
 
  do i = 1, size(how)
    math_expand(sum(how(1:i-1))+1:sum(how(1:i))) = what(mod(i-1,size(what))+1)
  enddo

end function math_expand


!--------------------------------------------------------------------------------------------------
!> @brief range of integers starting at one
!--------------------------------------------------------------------------------------------------
pure function math_range(N)

  integer, intent(in) :: N                                                                          !< length of range
  integer :: i
  integer, dimension(N) :: math_range

  math_range = [(i,i=1,N)]

end function math_range


!--------------------------------------------------------------------------------------------------
!> @brief second rank identity tensor of specified dimension
!--------------------------------------------------------------------------------------------------
pure function math_identity2nd(dimen)

  integer, intent(in) :: dimen                                                                      !< tensor dimension
  integer :: i
  real(pReal), dimension(dimen,dimen) :: math_identity2nd

  math_identity2nd = 0.0_pReal
  do i=1, dimen
    math_identity2nd(i,i) = 1.0_pReal
  enddo

end function math_identity2nd


!--------------------------------------------------------------------------------------------------
!> @brief symmetric fourth rank identity tensor of specified dimension
!  from http://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)#Derivative_of_a_second-order_tensor_with_respect_to_itself
!--------------------------------------------------------------------------------------------------
pure function math_identity4th(dimen)

  integer, intent(in) :: dimen                                                                      !< tensor dimension
  integer :: i,j,k,l
  real(pReal), dimension(dimen,dimen,dimen,dimen) ::  math_identity4th
  real(pReal), dimension(dimen,dimen)             ::  identity2nd

  identity2nd = math_identity2nd(dimen)
  forall(i=1:dimen,j=1:dimen,k=1:dimen,l=1:dimen) &
    math_identity4th(i,j,k,l) = 0.5_pReal*(identity2nd(i,k)*identity2nd(j,l)+identity2nd(i,l)*identity2nd(j,k))

end function math_identity4th


!--------------------------------------------------------------------------------------------------
!> @brief permutation tensor e_ijk used for computing cross product of two tensors
! e_ijk =  1 if even permutation of ijk
! e_ijk = -1 if odd permutation of ijk
! e_ijk =  0 otherwise
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_civita(i,j,k)

  integer, intent(in) :: i,j,k

  math_civita = 0.0_pReal
  if (((i == 1).and.(j == 2).and.(k == 3)) .or. &
      ((i == 2).and.(j == 3).and.(k == 1)) .or. &
      ((i == 3).and.(j == 1).and.(k == 2))) math_civita = 1.0_pReal
  if (((i == 1).and.(j == 3).and.(k == 2)) .or. &
      ((i == 2).and.(j == 1).and.(k == 3)) .or. &
      ((i == 3).and.(j == 2).and.(k == 1))) math_civita = -1.0_pReal

end function math_civita


!--------------------------------------------------------------------------------------------------
!> @brief kronecker delta function d_ij
! d_ij = 1 if i = j
! d_ij = 0 otherwise
! inspired by http://fortraninacworld.blogspot.de/2012/12/ternary-operator.html
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
!> @brief outer product A \otimes B of arbitrary sized vectors A and B
!--------------------------------------------------------------------------------------------------
pure function math_outer(A,B)

  real(pReal), dimension(:), intent(in) ::  A,B
  real(pReal), dimension(size(A,1),size(B,1)) ::  math_outer
  integer :: i,j

  forall(i=1:size(A,1),j=1:size(B,1)) math_outer(i,j) = A(i)*B(j)

end function math_outer


!--------------------------------------------------------------------------------------------------
!> @brief outer product A \otimes B of arbitrary sized vectors A and B
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_inner(A,B)

  real(pReal), dimension(:),         intent(in) :: A
  real(pReal), dimension(size(A,1)), intent(in) :: B

  math_inner = sum(A*B)

end function math_inner


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 33xx33 = 1 (double contraction --> ij * ij)
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_mul33xx33(A,B)

  real(pReal), dimension(3,3), intent(in) ::  A,B
  integer :: i,j
  real(pReal), dimension(3,3) ::              C

  forall(i=1:3,j=1:3) C(i,j) = A(i,j) * B(i,j)
  math_mul33xx33 = sum(C)

end function math_mul33xx33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3333x33 = 33 (double contraction --> ijkl *kl = ij)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx33(A,B)

  real(pReal), dimension(3,3) :: math_mul3333xx33
  real(pReal), dimension(3,3,3,3), intent(in) ::  A
  real(pReal), dimension(3,3), intent(in) ::  B
  integer :: i,j 

  forall(i = 1:3,j = 1:3) math_mul3333xx33(i,j) = sum(A(i,j,1:3,1:3)*B(1:3,1:3))

end function math_mul3333xx33


!--------------------------------------------------------------------------------------------------
!> @brief matrix multiplication 3333x3333 = 3333 (ijkl *klmn = ijmn)
!--------------------------------------------------------------------------------------------------
pure function math_mul3333xx3333(A,B)

  integer :: i,j,k,l
  real(pReal), dimension(3,3,3,3), intent(in) ::  A
  real(pReal), dimension(3,3,3,3), intent(in) ::  B
  real(pReal), dimension(3,3,3,3) :: math_mul3333xx3333

  forall(i = 1:3,j = 1:3, k = 1:3, l= 1:3) &
    math_mul3333xx3333(i,j,k,l) = sum(A(i,j,1:3,1:3)*B(1:3,1:3,k,l))

end function math_mul3333xx3333


!--------------------------------------------------------------------------------------------------
!> @brief 3x3 matrix exponential up to series approximation order n (default 5)
!--------------------------------------------------------------------------------------------------
pure function math_exp33(A,n)

  integer :: i
  integer, intent(in), optional :: n
  real(pReal), dimension(3,3), intent(in) :: A
  real(pReal), dimension(3,3) :: B, math_exp33
  real(pReal) :: invFac
  integer     :: order

  B = math_I3                                                                                       ! init
  invFac = 1.0_pReal                                                                                ! 0!
  math_exp33 = B                                                                                    ! A^0 = eye2

  if (present(n)) then
    order = n
  else
    order = 5
  endif
  
  do i = 1, order
    invFac = invFac/real(i,pReal)                                                                   ! invfac = 1/i!
    B = matmul(B,A)
    math_exp33 = math_exp33 + invFac*B                                                              ! exp = SUM (A^i)/i!
  enddo

end function math_exp33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 33 matrix (function)
!> @details Direct Cramer inversion of matrix A. Returns all zeroes if not possible, i.e. 
! if determinant is close to zero
!--------------------------------------------------------------------------------------------------
pure function math_inv33(A)

  real(pReal),dimension(3,3),intent(in)  :: A
  real(pReal)                :: DetA
  real(pReal),dimension(3,3) :: math_inv33
  logical                    :: error

  call math_invert33(math_inv33,DetA,error,A)
  if(error) math_inv33 = 0.0_pReal

end function math_inv33


!--------------------------------------------------------------------------------------------------
!> @brief Cramer inversion of 33 matrix (subroutine)
!> @details Direct Cramer inversion of matrix A. Also returns determinant
!  Returns an error if not possible, i.e. if determinant is close to zero
!--------------------------------------------------------------------------------------------------
pure subroutine math_invert33(InvA, DetA, error, A)

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

  real(pReal),dimension(3,3,3,3)            :: math_invSym3333

  real(pReal),dimension(3,3,3,3),intent(in) :: A

  integer :: ierr
  integer,     dimension(6)        :: ipiv6
  real(pReal), dimension(6,6)      :: temp66
  real(pReal), dimension(6*(64+2)) :: work
  logical                          :: error
  external :: &
   dgetrf, &
   dgetri

  temp66 = math_sym3333to66(A)
  call dgetrf(6,6,temp66,6,ipiv6,ierr)
  error = (ierr /= 0)
  call dgetri(6,temp66,6,ipiv6,work,size(work,1),ierr)
  error = error .or. (ierr /= 0)
  if (error) then
    call IO_error(400, ext_msg = 'math_invSym3333')
  else
    math_invSym3333 = math_66toSym3333(temp66)
  endif

end function math_invSym3333


!--------------------------------------------------------------------------------------------------
!> @brief invert quadratic matrix of arbitrary dimension
!--------------------------------------------------------------------------------------------------
subroutine math_invert(InvA, error, A)

  real(pReal), dimension(:,:),                 intent(in)  :: A
  real(pReal), dimension(size(A,1),size(A,1)), intent(out) :: invA
  logical,                                     intent(out) :: error

  integer,     dimension(size(A,1))        :: ipiv
  real(pReal), dimension(size(A,1)*(64+2)) :: work
  integer                                  :: ierr
  external :: &
   dgetrf, &
   dgetri

  invA = A 
  call dgetrf(size(A,1),size(A,1),invA,size(A,1),ipiv,ierr)
  error = (ierr /= 0)
  call dgetri(size(A,1),InvA,size(A,1),ipiv,work,size(work,1),ierr)
  error = error .or. (ierr /= 0)
 
end subroutine math_invert


!--------------------------------------------------------------------------------------------------
!> @brief symmetrize a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_symmetric33(m)

  real(pReal), dimension(3,3) :: math_symmetric33
  real(pReal), dimension(3,3), intent(in) :: m

  math_symmetric33 = 0.5_pReal * (m + transpose(m))

end function math_symmetric33


!--------------------------------------------------------------------------------------------------
!> @brief symmetrize a 66 matrix
!--------------------------------------------------------------------------------------------------
pure function math_symmetric66(m)

  real(pReal), dimension(6,6) :: math_symmetric66
  real(pReal), dimension(6,6), intent(in) :: m

  math_symmetric66 = 0.5_pReal * (m + transpose(m))

end function math_symmetric66


!--------------------------------------------------------------------------------------------------
!> @brief skew part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_skew33(m)

  real(pReal), dimension(3,3) :: math_skew33
  real(pReal), dimension(3,3), intent(in) :: m

  math_skew33 = m - math_symmetric33(m)

end function math_skew33


!--------------------------------------------------------------------------------------------------
!> @brief hydrostatic part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_spherical33(m)

  real(pReal), dimension(3,3) :: math_spherical33
  real(pReal), dimension(3,3), intent(in) :: m

  math_spherical33 = math_I3 * math_trace33(m)/3.0_pReal

end function math_spherical33


!--------------------------------------------------------------------------------------------------
!> @brief deviatoric part of a 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_deviatoric33(m)

  real(pReal), dimension(3,3) :: math_deviatoric33
  real(pReal), dimension(3,3), intent(in) :: m

  math_deviatoric33 = m - math_spherical33(m)

end function math_deviatoric33


!--------------------------------------------------------------------------------------------------
!> @brief trace of a 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_trace33(m)

  real(pReal), dimension(3,3), intent(in) :: m

  math_trace33 = m(1,1) + m(2,2) + m(3,3)

end function math_trace33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_det33(m)

  real(pReal), dimension(3,3), intent(in) :: m
  
  math_det33 = m(1,1)* (m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
             - m(1,2)* (m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
             + m(1,3)* (m(2,1)*m(3,2)-m(2,2)*m(3,1))

end function math_det33


!--------------------------------------------------------------------------------------------------
!> @brief determinant of a symmetric 33 matrix
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_detSym33(m)

  real(pReal), dimension(3,3), intent(in) :: m
 
  math_detSym33 = -(m(1,1)*m(2,3)**2 + m(2,2)*m(1,3)**2 + m(3,3)*m(1,2)**2) &
                  + m(1,1)*m(2,2)*m(3,3) + 2.0_pReal * m(1,2)*m(1,3)*m(2,3)

end function  math_detSym33


!--------------------------------------------------------------------------------------------------
!> @brief convert 33 matrix into vector 9
!--------------------------------------------------------------------------------------------------
pure function math_33to9(m33)

  real(pReal), dimension(9)               :: math_33to9
  real(pReal), dimension(3,3), intent(in) :: m33
  
  integer :: i
  
  do i = 1, 9
    math_33to9(i) = m33(mapPlain(1,i),mapPlain(2,i))
  enddo

end function math_33to9


!--------------------------------------------------------------------------------------------------
!> @brief convert 9 vector into 33 matrix
!--------------------------------------------------------------------------------------------------
pure function math_9to33(v9)

  real(pReal), dimension(3,3)           :: math_9to33
  real(pReal), dimension(9), intent(in) :: v9
  
  integer :: i

  do i = 1, 9
    math_9to33(mapPlain(1,i),mapPlain(2,i)) = v9(i)
  enddo

end function math_9to33


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 33 matrix into 6 vector
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_sym33to6(m33,weighted)

  real(pReal), dimension(6)               :: math_sym33to6
  real(pReal), dimension(3,3), intent(in) :: m33                                                     !< symmetric matrix (no internal check)
  logical,     optional,       intent(in) :: weighted                                                !< weight according to Mandel (.true. by default)
  
  real(pReal), dimension(6) :: w
  integer :: i
  
  if(present(weighted)) then
    w = merge(nrmMandel,1.0_pReal,weighted)
  else
    w = nrmMandel
  endif

  do i = 1, 6
    math_sym33to6(i) = w(i)*m33(mapNye(1,i),mapNye(2,i))
  enddo 

end function math_sym33to6


!--------------------------------------------------------------------------------------------------
!> @brief convert 6 vector into symmetric 33 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_6toSym33(v6,weighted)

  real(pReal), dimension(3,3)           :: math_6toSym33
  real(pReal), dimension(6), intent(in) :: v6
  logical,     optional,     intent(in) :: weighted                                                  !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i
  
  if(present(weighted)) then
    w = merge(invnrmMandel,1.0_pReal,weighted)
  else
    w = invnrmMandel
  endif

  do i=1,6
    math_6toSym33(mapNye(1,i),mapNye(2,i)) = w(i)*v6(i)
    math_6toSym33(mapNye(2,i),mapNye(1,i)) = w(i)*v6(i)
  enddo

end function math_6toSym33


!--------------------------------------------------------------------------------------------------
!> @brief convert 3333 matrix into 99 matrix
!--------------------------------------------------------------------------------------------------
pure function math_3333to99(m3333)

  real(pReal), dimension(9,9)                 :: math_3333to99
  real(pReal), dimension(3,3,3,3), intent(in) :: m3333

  integer :: i,j

  do i=1,9; do j=1,9
    math_3333to99(i,j) = m3333(mapPlain(1,i),mapPlain(2,i),mapPlain(1,j),mapPlain(2,j))
  enddo; enddo

end function math_3333to99


!--------------------------------------------------------------------------------------------------
!> @brief convert 99 matrix into 3333 matrix
!--------------------------------------------------------------------------------------------------
pure function math_99to3333(m99)

  real(pReal), dimension(3,3,3,3)         :: math_99to3333
  real(pReal), dimension(9,9), intent(in) :: m99

  integer :: i,j

  do i=1,9; do j=1,9
    math_99to3333(mapPlain(1,i),mapPlain(2,i),mapPlain(1,j),mapPlain(2,j)) = m99(i,j)
  enddo; enddo

end function math_99to3333


!--------------------------------------------------------------------------------------------------
!> @brief convert symmetric 3333 matrix into 66 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_sym3333to66(m3333,weighted)

  real(pReal), dimension(6,6)                 :: math_sym3333to66
  real(pReal), dimension(3,3,3,3), intent(in) :: m3333                                               !< symmetric matrix (no internal check)
  logical,     optional,           intent(in) :: weighted                                            !< weight according to Mandel (.true. by default)

  real(pReal), dimension(6) :: w
  integer :: i,j
  
  if(present(weighted)) then
    w = merge(nrmMandel,1.0_pReal,weighted)
  else
    w = nrmMandel
  endif

  do i=1,6; do j=1,6
    math_sym3333to66(i,j) = w(i)*w(j)*m3333(mapNye(1,i),mapNye(2,i),mapNye(1,j),mapNye(2,j))
  enddo; enddo

end function math_sym3333to66


!--------------------------------------------------------------------------------------------------
!> @brief convert 66 matrix into symmetric 3333 matrix
!> @details Weighted conversion (default) rearranges according to Nye and weights shear
! components according to Mandel. Advisable for matrix operations.
! Unweighted conversion only changes order according to Nye
!--------------------------------------------------------------------------------------------------
pure function math_66toSym3333(m66,weighted)

  real(pReal), dimension(3,3,3,3)            :: math_66toSym3333
  real(pReal), dimension(6,6),    intent(in) :: m66
  logical,     optional,          intent(in) :: weighted                                             !< weight according to Mandel (.true. by default)
  
  real(pReal), dimension(6) :: w
  integer :: i,j
  
  if(present(weighted)) then
    w = merge(invnrmMandel,1.0_pReal,weighted)
  else
    w = invnrmMandel
  endif

  do i=1,6; do j=1,6
    math_66toSym3333(mapNye(1,i),mapNye(2,i),mapNye(1,j),mapNye(2,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(mapNye(2,i),mapNye(1,i),mapNye(1,j),mapNye(2,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(mapNye(1,i),mapNye(2,i),mapNye(2,j),mapNye(1,j)) = w(i)*w(j)*m66(i,j)
    math_66toSym3333(mapNye(2,i),mapNye(1,i),mapNye(2,j),mapNye(1,j)) = w(i)*w(j)*m66(i,j)
  enddo; enddo

end function math_66toSym3333


!--------------------------------------------------------------------------------------------------
!> @brief convert 66 Voigt matrix into symmetric 3333 matrix
!--------------------------------------------------------------------------------------------------
pure function math_Voigt66to3333(m66)

  real(pReal), dimension(3,3,3,3) :: math_Voigt66to3333
  real(pReal), dimension(6,6), intent(in) :: m66
  integer :: i,j

  do i=1,6; do j=1, 6
    math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(1,j),mapVoigt(2,j)) = m66(i,j)
    math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(1,j),mapVoigt(2,j)) = m66(i,j)
    math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(2,j),mapVoigt(1,j)) = m66(i,j)
    math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(2,j),mapVoigt(1,j)) = m66(i,j)
  enddo; enddo

end function math_Voigt66to3333


!--------------------------------------------------------------------------------------------------
!> @brief action of a quaternion on a vector (rotate vector v with Q)
!--------------------------------------------------------------------------------------------------
pure function math_qRot(Q,v)

  real(pReal), dimension(4), intent(in) :: Q
  real(pReal), dimension(3), intent(in) :: v
  real(pReal), dimension(3) :: math_qRot
  real(pReal), dimension(4,4) :: T
  integer :: i, j

  do i = 1,4
    do j = 1,i
      T(i,j) = Q(i) * Q(j)
    enddo
  enddo

  math_qRot = [-v(1)*(T(3,3)+T(4,4)) + v(2)*(T(3,2)-T(4,1)) + v(3)*(T(4,2)+T(3,1)), &
                v(1)*(T(3,2)+T(4,1)) - v(2)*(T(2,2)+T(4,4)) + v(3)*(T(4,3)-T(2,1)), &
                v(1)*(T(4,2)-T(3,1)) + v(2)*(T(4,3)+T(2,1)) - v(3)*(T(2,2)+T(3,3))]

  math_qRot = 2.0_pReal * math_qRot + v

end function math_qRot


!--------------------------------------------------------------------------------------------------
!> @brief rotation matrix from Bunge-Euler (3-1-3) angles (in radians)
!> @details deprecated
!--------------------------------------------------------------------------------------------------
pure function math_EulerToR(Euler)

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
!> @brief draw a random sample from Gauss variable
!--------------------------------------------------------------------------------------------------
real(pReal) function math_sampleGaussVar(meanvalue, stddev, width)

  real(pReal), intent(in) ::            meanvalue, &                                                ! meanvalue of gauss distribution
                                        stddev                                                      ! standard deviation of gauss distribution
  real(pReal), intent(in), optional ::  width                                                       ! width of considered values as multiples of standard deviation
  real(pReal), dimension(2) ::          rnd                                                         ! random numbers
  real(pReal) ::                        scatter, &                                                  ! normalized scatter around meanvalue
                                        myWidth

  if (abs(stddev) < tol_math_check) then
    math_sampleGaussVar = meanvalue
  else
    myWidth = merge(width,3.0_pReal,present(width))                                                 ! use +-3*sigma as default value for scatter if not given
   
    do
      call random_number(rnd)
      scatter = myWidth * (2.0_pReal * rnd(1) - 1.0_pReal)
      if (rnd(2) <= exp(-0.5_pReal * scatter ** 2.0_pReal)) exit                                    ! test if scattered value is drawn
    enddo

    math_sampleGaussVar = scatter * stddev
  endif

end function math_sampleGaussVar


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric matrix m
! ToDo: has wrong oder of arguments
!--------------------------------------------------------------------------------------------------
subroutine math_eigenValuesVectorsSym(m,values,vectors,error)

  real(pReal), dimension(:,:),                  intent(in)  :: m
  real(pReal), dimension(size(m,1)),            intent(out) :: values
  real(pReal), dimension(size(m,1),size(m,1)),  intent(out) :: vectors
  logical, intent(out) :: error
  integer :: ierr
  real(pReal), dimension((64+2)*size(m,1)) :: work                                                  ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f
  external :: &
    dsyev

  vectors = m                                                                                       ! copy matrix to input (doubles as output) array
  call dsyev('V','U',size(m,1),vectors,size(m,1),values,work,size(work,1),ierr)
  error = (ierr /= 0)

end subroutine math_eigenValuesVectorsSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues and eigenvectors of symmetric 33 matrix m using an analytical expression
!> and the general LAPACK powered version for arbritrary sized matrices as fallback
!> @author Joachim Kopp, Max-Planck-Institut für Kernphysik, Heidelberg (Copyright (C) 2006)
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details See http://arxiv.org/abs/physics/0610206 (DSYEVH3)
! ToDo: has wrong oder of arguments
!--------------------------------------------------------------------------------------------------
subroutine math_eigenValuesVectorsSym33(m,values,vectors)
 
  real(pReal), dimension(3,3),intent(in)  :: m
  real(pReal), dimension(3),  intent(out) :: values
  real(pReal), dimension(3,3),intent(out) :: vectors
  real(pReal) :: T, U, norm, threshold
  logical :: error
 
  values = math_eigenvaluesSym33(m)
 
  vectors(1:3,2) = [ m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2), &
                     m(1, 3) * m(1, 2) - m(2, 3) * m(1, 1), &
                     m(1, 2)**2]
 
  T = maxval(abs(values))
  U = max(T, T**2)
  threshold = sqrt(5.68e-14_pReal * U**2)
 
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
  vectors(1:3,3) = math_cross(vectors(1:3,1),vectors(1:3,2))

end subroutine math_eigenValuesVectorsSym33


!--------------------------------------------------------------------------------------------------
!> @brief eigenvector basis of symmetric matrix m
!--------------------------------------------------------------------------------------------------
function math_eigenvectorBasisSym(m)

  real(pReal), dimension(:,:),      intent(in)  :: m
  real(pReal), dimension(size(m,1))             :: values
  real(pReal), dimension(size(m,1),size(m,1))   :: vectors
  real(pReal), dimension(size(m,1),size(m,1))   :: math_eigenvectorBasisSym
  logical :: error
  integer :: i
 
  math_eigenvectorBasisSym = 0.0_pReal
  call math_eigenValuesVectorsSym(m,values,vectors,error)
  if(error) return
  
  do i=1, size(m,1)
    math_eigenvectorBasisSym = math_eigenvectorBasisSym &
                             + sqrt(values(i)) * math_outer(vectors(:,i),vectors(:,i))
  enddo

end function math_eigenvectorBasisSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvector basis of symmetric 33 matrix m
!--------------------------------------------------------------------------------------------------
pure function math_eigenvectorBasisSym33(m)

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
  ! this is not really correct, but at least the basis is correct
    EB(1,1,1)=1.0_pReal
    EB(2,2,2)=1.0_pReal
    EB(3,3,3)=1.0_pReal
  else threeSimilarEigenvalues
    rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
    phi=acos(math_clip(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
    values = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                              [cos(phi/3.0_pReal), &
                               cos((phi+2.0_pReal*PI)/3.0_pReal), &
                               cos((phi+4.0_pReal*PI)/3.0_pReal) &
                              ] + invariants(1)/3.0_pReal
    N(1:3,1:3,1) = m-values(1)*math_I3
    N(1:3,1:3,2) = m-values(2)*math_I3
    N(1:3,1:3,3) = m-values(3)*math_I3
    twoSimilarEigenvalues: if(abs(values(1)-values(2)) < TOL) then 
      EB(1:3,1:3,3)=matmul(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                                ((values(3)-values(1))*(values(3)-values(2)))
      EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,3)
    elseif(abs(values(2)-values(3)) < TOL) then twoSimilarEigenvalues
      EB(1:3,1:3,1)=matmul(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                                ((values(1)-values(2))*(values(1)-values(3)))
      EB(1:3,1:3,2)=math_I3-EB(1:3,1:3,1)
    elseif(abs(values(3)-values(1)) < TOL) then twoSimilarEigenvalues 
      EB(1:3,1:3,2)=matmul(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                                ((values(2)-values(1))*(values(2)-values(3)))
      EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,2)
    else twoSimilarEigenvalues
      EB(1:3,1:3,1)=matmul(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                                ((values(1)-values(2))*(values(1)-values(3)))
      EB(1:3,1:3,2)=matmul(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                                ((values(2)-values(1))*(values(2)-values(3)))
      EB(1:3,1:3,3)=matmul(N(1:3,1:3,1),N(1:3,1:3,2))/ &
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
pure function math_eigenvectorBasisSym33_log(m)

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
  ! this is not really correct, but at least the basis is correct
    EB(1,1,1)=1.0_pReal
    EB(2,2,2)=1.0_pReal
    EB(3,3,3)=1.0_pReal
  else threeSimilarEigenvalues
    rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
    phi=acos(math_clip(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
    values = 2.0_pReal*rho**(1.0_pReal/3.0_pReal)* &
                              [cos(phi/3.0_pReal), &
                               cos((phi+2.0_pReal*PI)/3.0_pReal), &
                               cos((phi+4.0_pReal*PI)/3.0_pReal) &
                              ] + invariants(1)/3.0_pReal
    N(1:3,1:3,1) = m-values(1)*math_I3
    N(1:3,1:3,2) = m-values(2)*math_I3
    N(1:3,1:3,3) = m-values(3)*math_I3
    twoSimilarEigenvalues: if(abs(values(1)-values(2)) < TOL) then 
      EB(1:3,1:3,3)=matmul(N(1:3,1:3,1),N(1:3,1:3,2))/ &
                                                ((values(3)-values(1))*(values(3)-values(2)))
      EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,3)
    elseif(abs(values(2)-values(3)) < TOL) then twoSimilarEigenvalues
      EB(1:3,1:3,1)=matmul(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                                ((values(1)-values(2))*(values(1)-values(3)))
      EB(1:3,1:3,2)=math_I3-EB(1:3,1:3,1)
    elseif(abs(values(3)-values(1)) < TOL) then twoSimilarEigenvalues 
      EB(1:3,1:3,2)=matmul(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                                ((values(2)-values(1))*(values(2)-values(3)))
      EB(1:3,1:3,1)=math_I3-EB(1:3,1:3,2)
    else twoSimilarEigenvalues
      EB(1:3,1:3,1)=matmul(N(1:3,1:3,2),N(1:3,1:3,3))/ &
                                                ((values(1)-values(2))*(values(1)-values(3)))
      EB(1:3,1:3,2)=matmul(N(1:3,1:3,1),N(1:3,1:3,3))/ &
                                                ((values(2)-values(1))*(values(2)-values(3)))
      EB(1:3,1:3,3)=matmul(N(1:3,1:3,1),N(1:3,1:3,2))/ &
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

  real(pReal), intent(in), dimension(3,3) :: m
  real(pReal), dimension(3,3) :: math_rotationalPart33
  real(pReal), dimension(3,3) :: U , Uinv

  U = math_eigenvectorBasisSym33(matmul(transpose(m),m))
  Uinv = math_inv33(U)

  inversionFailed: if (all(dEq0(Uinv))) then
    math_rotationalPart33 = math_I3
    call IO_warning(650)
  else inversionFailed
    math_rotationalPart33 = matmul(m,Uinv)
  endif inversionFailed

end function math_rotationalPart33


!--------------------------------------------------------------------------------------------------
!> @brief Eigenvalues of symmetric matrix m
! will return NaN on error
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym(m)
 
  real(pReal), dimension(:,:),                  intent(in)  :: m
  real(pReal), dimension(size(m,1))                         :: math_eigenvaluesSym
  real(pReal), dimension(size(m,1),size(m,1))               :: m_
  integer :: ierr
  real(pReal), dimension((64+2)*size(m,1)) :: work                                                  ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f
  external :: &
   dsyev

  m_= m                                                                                             ! copy matrix to input (will be destroyed)
  call dsyev('N','U',size(m,1),m_,size(m,1),math_eigenvaluesSym,work,size(work,1),ierr)
  if (ierr /= 0) math_eigenvaluesSym = IEEE_value(1.0_pReal,IEEE_quiet_NaN)

end function math_eigenvaluesSym


!--------------------------------------------------------------------------------------------------
!> @brief eigenvalues of symmetric 33 matrix m using an analytical expression
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @details similar to http://arxiv.org/abs/physics/0610206 (DSYEVC3)
!> but apparently more stable solution and has general LAPACK powered version for arbritrary sized
!> matrices as fallback
!--------------------------------------------------------------------------------------------------
function math_eigenvaluesSym33(m)

  real(pReal), intent(in), dimension(3,3) :: m
  real(pReal), dimension(3) :: math_eigenvaluesSym33,invariants
  real(pReal) :: P, Q, rho, phi
  real(pReal), parameter :: TOL=1.e-14_pReal

  invariants = math_invariantsSym33(m)                                                              ! invariants are coefficients in characteristic polynomial apart for the sign of c0 and c2 in http://arxiv.org/abs/physics/0610206

  P = invariants(2)-invariants(1)**2.0_pReal/3.0_pReal                                              ! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)
  Q = -2.0_pReal/27.0_pReal*invariants(1)**3.0_pReal+product(invariants(1:2))/3.0_pReal-invariants(3)! different from http://arxiv.org/abs/physics/0610206 (this formulation was in DAMASK)

  if(all(abs([P,Q]) < TOL)) then
    math_eigenvaluesSym33 = math_eigenvaluesSym(m)
  else
    rho=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
    phi=acos(math_clip(-Q/rho*0.5_pReal,-1.0_pReal,1.0_pReal))
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
  integer :: i, j
  
  j = min(k,n-k)
  math_binomial = product([(i, i=n, n-j+1, -1)])/math_factorial(j)

end function math_binomial


!--------------------------------------------------------------------------------------------------
!> @brief multinomial coefficient
!--------------------------------------------------------------------------------------------------
integer pure function math_multinomial(alpha)

  integer, intent(in), dimension(:) :: alpha
  integer :: i
 
  math_multinomial = 1
  do i = 1, size(alpha)
    math_multinomial = math_multinomial*math_binomial(sum(alpha(1:i)),alpha(i))
  enddo  

end function math_multinomial


!--------------------------------------------------------------------------------------------------
!> @brief volume of tetrahedron given by four vertices
!--------------------------------------------------------------------------------------------------
real(pReal) pure function math_volTetrahedron(v1,v2,v3,v4)

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

  real(pReal), dimension (3), intent(in) :: v1,v2,v3

  math_areaTriangle = 0.5_pReal * norm2(math_cross(v1-v2,v1-v3))

end function math_areaTriangle


!--------------------------------------------------------------------------------------------------
!> @brief rotate 33 tensor forward
!--------------------------------------------------------------------------------------------------
pure function math_rotate_forward33(tensor,R)

  real(pReal), dimension(3,3)             ::  math_rotate_forward33
  real(pReal), dimension(3,3), intent(in) :: tensor, R

  math_rotate_forward33 = matmul(R,matmul(tensor,transpose(R)))

end function math_rotate_forward33


!--------------------------------------------------------------------------------------------------
!> @brief rotate 33 tensor backward
!--------------------------------------------------------------------------------------------------
pure function math_rotate_backward33(tensor,R)

  real(pReal), dimension(3,3)             ::  math_rotate_backward33
  real(pReal), dimension(3,3), intent(in) :: tensor, R

  math_rotate_backward33 = matmul(transpose(R),matmul(tensor,R))

end function math_rotate_backward33


!--------------------------------------------------------------------------------------------------
!> @brief rotate 3333 tensor C'_ijkl=g_im*g_jn*g_ko*g_lp*C_mnop
!--------------------------------------------------------------------------------------------------
pure function math_rotate_forward3333(tensor,R)

  real(pReal), dimension(3,3,3,3)             ::  math_rotate_forward3333
  real(pReal), dimension(3,3),     intent(in) :: R
  real(pReal), dimension(3,3,3,3), intent(in) :: tensor
  integer :: i,j,k,l,m,n,o,p

  math_rotate_forward3333 = 0.0_pReal
  do i = 1,3;do j = 1,3;do k = 1,3;do l = 1,3
  do m = 1,3;do n = 1,3;do o = 1,3;do p = 1,3
    math_rotate_forward3333(i,j,k,l) = math_rotate_forward3333(i,j,k,l) &
                                     + R(i,m) * R(j,n) * R(k,o) * R(l,p) * tensor(m,n,o,p)
  enddo; enddo; enddo; enddo; enddo; enddo; enddo; enddo

end function math_rotate_forward3333


!--------------------------------------------------------------------------------------------------
!> @brief limits a scalar value to a certain range (either one or two sided)
! Will return NaN if left > right
!--------------------------------------------------------------------------------------------------
real(pReal) pure elemental function math_clip(a, left, right)

  real(pReal), intent(in) :: a
  real(pReal), intent(in), optional :: left, right

  math_clip = a
  if (present(left))  math_clip = max(left,math_clip)
  if (present(right)) math_clip = min(right,math_clip)
  if (present(left) .and. present(right)) &
    math_clip = merge (IEEE_value(1.0_pReal,IEEE_quiet_NaN),math_clip, left>right)

end function math_clip


!--------------------------------------------------------------------------------------------------
!> @brief check correctness of (some) math functions
!--------------------------------------------------------------------------------------------------
subroutine unitTest
 
  integer, dimension(2,4) :: &
    sort_in_   = reshape([+1,+5,  +5,+6,  -1,-1,  +3,-2],[2,4])
  integer, dimension(2,4), parameter :: &
    sort_out_  = reshape([-1,-1,  +1,+5,  +5,+6,  +3,-2],[2,4])
    
  real(pReal), dimension(5) :: range_out_ = [1.0_pReal,2.0_pReal,3.0_pReal,4.0_pReal,5.0_pReal]

  real(pReal)                 :: det
  real(pReal), dimension(6)   :: v6
  real(pReal), dimension(9)   :: v9
  real(pReal), dimension(3,3) :: t33,t33_2
  real(pReal), dimension(6,6) :: t66
  real(pReal), dimension(9,9) :: t99,t99_2
  logical                     :: e
  

  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,3.0_pReal,3.0_pReal,3.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1,2,3,0])) > tol_math_check)) &
    call IO_error(401,ext_msg='math_expand [1,2,3] by [1,2,3,0] => [1,2,2,3,3,3]')

  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal,3.0_pReal],[1,2])) > tol_math_check)) &
    call IO_error(401,ext_msg='math_expand [1,2,3] by [1,2] => [1,2,2]')

  if (any(abs([1.0_pReal,2.0_pReal,2.0_pReal,1.0_pReal,1.0_pReal,1.0_pReal] - &
              math_expand([1.0_pReal,2.0_pReal],[1,2,3])) > tol_math_check)) &
    call IO_error(401,ext_msg='math_expand [1,2] by [1,2,3] => [1,2,2,1,1,1]')


  call math_sort(sort_in_,1,3,2)
  if(any(sort_in_ /= sort_out_)) &
    call IO_error(401,ext_msg='math_sort')
  
  if(any(math_range(5) /= range_out_)) &
    call IO_error(401,ext_msg='math_range')
 
 
  call random_number(v9)
  if(any(dNeq(math_33to9(math_9to33(v9)),v9))) &
    call IO_error(401,ext_msg='math_33to9/math_9to33')
    
  call random_number(t99)
  if(any(dNeq(math_3333to99(math_99to3333(t99)),t99))) &
    call IO_error(401,ext_msg='math_3333to99/math_99to3333')
  
  call random_number(v6)
  if(any(dNeq(math_sym33to6(math_6toSym33(v6)),v6))) &
    call IO_error(401,ext_msg='math_sym33to6/math_6toSym33')
  
  call random_number(t66)
  if(any(dNeq(math_sym3333to66(math_66toSym3333(t66)),t66))) &
    call IO_error(401,ext_msg='math_sym3333to66/math_66toSym3333')
    
  call random_number(v6)
  if(any(dNeq0(math_6toSym33(v6) - math_symmetric33(math_6toSym33(v6))))) &
    call IO_error(401,ext_msg='math_symmetric33')
 
  call random_number(t33)
  if(dNeq(math_det33(math_symmetric33(t33)),math_detSym33(math_symmetric33(t33)),tol=1.0e-12_pReal)) &
    call IO_error(401,ext_msg='math_det33/math_detSym33')
  
  do while(abs(math_det33(t33))<1.0e-9_pReal)
    call random_number(t33)
  enddo
  if(any(dNeq0(matmul(t33,math_inv33(t33)) - math_identity2nd(3),tol=1.0e-9_pReal))) &
    call IO_error(401,ext_msg='math_inv33')
    
  call math_invert33(t33_2,det,e,t33)
  if(any(dNeq0(matmul(t33,t33_2) - math_identity2nd(3),tol=1.0e-9_pReal)) .or. e) &
    call IO_error(401,ext_msg='math_invert33: T:T^-1 != I')
  if(dNeq(det,math_det33(t33),tol=1.0e-12_pReal)) &
    call IO_error(401,ext_msg='math_invert33 (determinant)')
    
  call math_invert(t33_2,e,t33)
  if(any(dNeq0(matmul(t33,t33_2) - math_identity2nd(3),tol=1.0e-9_pReal)) .or. e) &
    call IO_error(401,ext_msg='math_invert t33')
    
  t33_2 = transpose(math_rotationalPart33(t33))
  if(any(dNeq0(math_rotationalPart33(matmul(t33_2,t33)) - MATH_I3,tol=5.0e-4_pReal))) &
    call IO_error(401,ext_msg='math_rotationalPart33')


  call math_invert(t99_2,e,t99) ! not sure how likely it is that we get a singular matrix 
  if(any(dNeq0(matmul(t99_2,t99)-math_identity2nd(9),tol=1.0e-9_pReal)) .or. e) &
    call IO_error(401,ext_msg='math_invert t99')

  if(any(math_clip([4.0_pReal,9.0_pReal],5.0_pReal,6.5_pReal)/=[5.0_pReal,6.5_pReal])) &
     call IO_error(401,ext_msg='math_clip')  

end subroutine unitTest

end module math
