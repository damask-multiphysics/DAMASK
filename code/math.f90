! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
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
!##############################################################
!* $Id$
!##############################################################
 MODULE math   
!##############################################################


 use, intrinsic :: iso_c_binding
 use prec, only: pReal,pInt
 use IO, only: IO_error
 implicit none

 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal
 real(pReal), parameter :: inDeg = 180.0_pReal/pi
 real(pReal), parameter :: inRad = pi/180.0_pReal

! *** 3x3 Identity ***
 real(pReal), dimension(3,3), parameter :: math_I3 = &
 reshape( (/ &
 1.0_pReal,0.0_pReal,0.0_pReal, &
 0.0_pReal,1.0_pReal,0.0_pReal, &
 0.0_pReal,0.0_pReal,1.0_pReal /),(/3,3/))

! *** Mandel notation ***
 integer(pInt), dimension (2,6), parameter :: mapMandel = &
 reshape((/&
  1_pInt,1_pInt, &
  2_pInt,2_pInt, &
  3_pInt,3_pInt, &
  1_pInt,2_pInt, &
  2_pInt,3_pInt, &
  1_pInt,3_pInt  &
 /),(/2,6/))

 real(pReal), dimension(6), parameter :: nrmMandel = &
 (/1.0_pReal,1.0_pReal,1.0_pReal, 1.414213562373095_pReal, 1.414213562373095_pReal, 1.414213562373095_pReal/)
 real(pReal), dimension(6), parameter :: invnrmMandel = &
 (/1.0_pReal,1.0_pReal,1.0_pReal,0.7071067811865476_pReal,0.7071067811865476_pReal,0.7071067811865476_pReal/)

! *** Voigt notation ***
 integer(pInt), dimension (2,6), parameter :: mapVoigt = &
 reshape((/&
  1_pInt,1_pInt, &
  2_pInt,2_pInt, &
  3_pInt,3_pInt, &
  2_pInt,3_pInt, &
  1_pInt,3_pInt, &
  1_pInt,2_pInt  &
 /),(/2,6/))

 real(pReal), dimension(6), parameter :: nrmVoigt = &
 (/1.0_pReal,1.0_pReal,1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal/)
 real(pReal), dimension(6), parameter :: invnrmVoigt = &
 (/1.0_pReal,1.0_pReal,1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal/)

! *** Plain notation ***
 integer(pInt), dimension (2,9), parameter :: mapPlain = &
 reshape((/&
  1_pInt,1_pInt, &
  1_pInt,2_pInt, &
  1_pInt,3_pInt, &
  2_pInt,1_pInt, &
  2_pInt,2_pInt, &
  2_pInt,3_pInt, &
  3_pInt,1_pInt, &
  3_pInt,2_pInt, &
  3_pInt,3_pInt  &
 /),(/2,9/))

! Symmetry operations as quaternions
! 24 for cubic, 12 for hexagonal = 36
integer(pInt), dimension(2), parameter :: math_NsymOperations = (/24_pInt,12_pInt/)
real(pReal), dimension(4,36), parameter :: math_symOperations = &
  reshape((/&
     1.0_pReal,                 0.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      ! cubic symmetry operations
     0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal,  0.7071067811865476_pReal, &       !     2-fold symmetry
     0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, -0.7071067811865476_pReal, &
     0.0_pReal,                -0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.7071067811865476_pReal, -0.7071067811865476_pReal,  0.0_pReal, &
     0.5_pReal,                 0.5_pReal,                 0.5_pReal,                 0.5_pReal, &                      !     3-fold symmetry
    -0.5_pReal,                 0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                -0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
    -0.5_pReal,                -0.5_pReal,                 0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.5_pReal,                -0.5_pReal,                 0.5_pReal, &
    -0.5_pReal,                 0.5_pReal,                -0.5_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.5_pReal,                 0.5_pReal,                -0.5_pReal, &
    -0.5_pReal,                 0.5_pReal,                 0.5_pReal,                -0.5_pReal, &
     0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal, &                      !     4-fold symmetry
     0.0_pReal,                 1.0_pReal,                 0.0_pReal,                 0.0_pReal, &
    -0.7071067811865476_pReal,  0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal, &
     0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal, &
     0.0_pReal,                 0.0_pReal,                 1.0_pReal,                 0.0_pReal, &
    -0.7071067811865476_pReal,  0.0_pReal,                 0.7071067811865476_pReal,  0.0_pReal, &
     0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.0_pReal,                 1.0_pReal, &
    -0.7071067811865476_pReal,  0.0_pReal,                 0.0_pReal,                 0.7071067811865476_pReal, &
     1.0_pReal,                 0.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      ! hexagonal symmetry operations
     0.0_pReal,                 1.0_pReal,                 0.0_pReal,                 0.0_pReal, &                      !     2-fold symmetry
     0.0_pReal,                 0.0_pReal,                 1.0_pReal,                 0.0_pReal, &
     0.0_pReal,                 0.5_pReal,                 0.866025403784439_pReal,   0.0_pReal, &
     0.0_pReal,                -0.5_pReal,                 0.866025403784439_pReal,   0.0_pReal, &
     0.0_pReal,                 0.866025403784439_pReal,   0.5_pReal,                 0.0_pReal, &
     0.0_pReal,                -0.866025403784439_pReal,   0.5_pReal,                 0.0_pReal, &
     0.866025403784439_pReal,   0.0_pReal,                 0.0_pReal,                 0.5_pReal, &                      !     6-fold symmetry
    -0.866025403784439_pReal,   0.0_pReal,                 0.0_pReal,                 0.5_pReal, &
     0.5_pReal,                 0.0_pReal,                 0.0_pReal,                 0.866025403784439_pReal, &
    -0.5_pReal,                 0.0_pReal,                 0.0_pReal,                 0.866025403784439_pReal, &
     0.0_pReal,                 0.0_pReal,                 0.0_pReal,                 1.0_pReal &
  /),(/4,36/))

  include 'fftw3.f03'

 CONTAINS

!**************************************************************************
! initialization of module
!**************************************************************************
 SUBROUTINE math_init ()

 use prec,     only: tol_math_check
 use numerics, only: fixedSeed
 use IO,       only: IO_error
 use debug,    only: debug_verbosity
 implicit none

 integer(pInt) :: i
 real(pReal), dimension(3,3) :: R,R2
 real(pReal), dimension(3) ::   Eulers
 real(pReal), dimension(4) ::   q,q2,axisangle,randTest
! the following variables are system dependend and shound NOT be pInt
 integer :: randSize                                ! gfortran requires a variable length to compile 
 integer, dimension(:), allocatable :: randInit     ! if recalculations of former randomness (with given seed) is necessary
                                                    ! comment the first random_seed call out, set randSize to 1, and use ifort
 character(len=64) :: error_msg
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  math init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 !$OMP END CRITICAL (write2out)
 
 call random_seed(size=randSize)
 allocate(randInit(randSize))
 if (fixedSeed > 0_pInt) then
   randInit(1:randSize) = int(fixedSeed)                      ! fixedSeed is of type pInt, randInit not
   call random_seed(put=randInit)
 else
   call random_seed()
 endif

 call random_seed(get=randInit)

 do i = 1_pInt, 4_pInt
   call random_number(randTest(i))
 enddo

 !$OMP CRITICAL (write2out)
 ! this critical block did cause trouble at IWM
 write(6,*) 'value of random seed:    ', randInit(1)
 write(6,*) 'size of random seed:     ', randSize
 write(6,'(a,4(/,26x,f16.14))') ' start of random sequence: ', randTest
 write(6,*) ''
 !$OMP END CRITICAL (write2out)
  
 call random_seed(put=randInit)
 call random_seed(get=randInit)

 call halton_seed_set(randInit(1))
 call halton_ndim_set(3_pInt)

 ! --- check rotation dictionary ---

 ! +++ q -> a -> q  +++
 q = math_qRnd();
 axisangle = math_QuaternionToAxisAngle(q);
 q2 = math_AxisAngleToQuaternion(axisangle(1:3),axisangle(4))
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(670_pInt,ext_msg=error_msg)
 endif 
 
 ! +++ q -> R -> q  +++
 R = math_QuaternionToR(q);
 q2 = math_RToQuaternion(R)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(671_pInt,ext_msg=error_msg)
 endif 
 
 ! +++ q -> euler -> q  +++
 Eulers = math_QuaternionToEuler(q);
 q2 = math_EulerToQuaternion(Eulers)
 if ( any(abs( q-q2) > tol_math_check) .and. &
      any(abs(-q-q2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation',min(maxval(abs( q-q2)),maxval(abs(-q-q2)))
   call IO_error(672_pInt,ext_msg=error_msg)
 endif 

 ! +++ R -> euler -> R  +++
 Eulers = math_RToEuler(R);
 R2 = math_EulerToR(Eulers)
 if ( any(abs( R-R2) > tol_math_check) ) then
   write (error_msg, '(a,e14.6)' ) 'maximum deviation',maxval(abs( R-R2))
   call IO_error(673_pInt,ext_msg=error_msg)
 endif 
 
 ENDSUBROUTINE math_init
 

!**************************************************************************
! Quicksort algorithm for two-dimensional integer arrays
!
! Sorting is done with respect to array(1,:)
! and keeps array(2:N,:) linked to it.
!**************************************************************************
 RECURSIVE SUBROUTINE qsort(a, istart, iend)

 implicit none
 integer(pInt), dimension(:,:) :: a
 integer(pInt) :: istart,iend,ipivot

 if (istart < iend) then
   ipivot = math_partition(a,istart, iend)
   call qsort(a, istart, ipivot-1_pInt)
   call qsort(a, ipivot+1_pInt, iend)
 endif
  
 ENDSUBROUTINE qsort


!**************************************************************************
! Partitioning required for quicksort
!**************************************************************************
 integer(pInt) function math_partition(a, istart, iend)

 implicit none
 integer(pInt), dimension(:,:) :: a
 integer(pInt) :: istart,iend,d,i,j,k,x,tmp

 d = size(a,1_pInt) ! number of linked data
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

 endfunction math_partition
 

!**************************************************************************
! range of integers starting at one
!**************************************************************************
 pure function math_range(N)  

 implicit none

 integer(pInt), intent(in) :: N
 integer(pInt) :: i
 integer(pInt), dimension(N) :: math_range

 forall (i=1_pInt:N) math_range(i) = i

 endfunction math_range


!**************************************************************************
! second rank identity tensor of specified dimension
!**************************************************************************
 pure function math_identity2nd(dimen)  

 implicit none

 integer(pInt), intent(in) :: dimen
 integer(pInt) :: i
 real(pReal), dimension(dimen,dimen) :: math_identity2nd

 math_identity2nd = 0.0_pReal 
 forall (i=1_pInt:dimen) math_identity2nd(i,i) = 1.0_pReal 

 endfunction math_identity2nd


!**************************************************************************
! permutation tensor e_ijk used for computing cross product of two tensors
! e_ijk =  1 if even permutation of ijk
! e_ijk = -1 if odd permutation of ijk
! e_ijk =  0 otherwise
!**************************************************************************
 pure function math_civita(i,j,k)

 implicit none

 integer(pInt), intent(in) :: i,j,k
 real(pReal) math_civita

 math_civita = 0.0_pReal
 if (((i == 1_pInt).and.(j == 2_pInt).and.(k == 3_pInt)) .or. &
     ((i == 2_pInt).and.(j == 3_pInt).and.(k == 1_pInt)) .or. &
     ((i == 3_pInt).and.(j == 1_pInt).and.(k == 2_pInt))) math_civita = 1.0_pReal
 if (((i == 1_pInt).and.(j == 3_pInt).and.(k == 2_pInt)) .or. &
     ((i == 2_pInt).and.(j == 1_pInt).and.(k == 3_pInt)) .or. &
     ((i == 3_pInt).and.(j == 2_pInt).and.(k == 1_pInt))) math_civita = -1.0_pReal

 endfunction math_civita


!**************************************************************************
! kronecker delta function d_ij
! d_ij = 1 if i = j
! d_ij = 0 otherwise
!**************************************************************************
 pure function math_delta(i,j)

 implicit none

 integer(pInt), intent (in) :: i,j
 real(pReal) :: math_delta

 math_delta = 0.0_pReal
 if (i == j) math_delta = 1.0_pReal

 endfunction math_delta


!**************************************************************************
! fourth rank identity tensor of specified dimension
!**************************************************************************
 pure function math_identity4th(dimen)  

 implicit none

 integer(pInt), intent(in) :: dimen
 integer(pInt) :: i,j,k,l
 real(pReal), dimension(dimen,dimen,dimen,dimen) ::  math_identity4th

 forall (i=1_pInt:dimen,j=1_pInt:dimen,k=1_pInt:dimen,l=1_pInt:dimen) math_identity4th(i,j,k,l) = &
        0.5_pReal*(math_I3(i,k)*math_I3(j,k)+math_I3(i,l)*math_I3(j,k)) 

 endfunction math_identity4th
 

!**************************************************************************
! vector product a x b
!**************************************************************************
 pure function math_vectorproduct(A,B)  

 implicit none

 real(pReal), dimension(3), intent(in) ::  A,B
 real(pReal), dimension(3) ::  math_vectorproduct

 math_vectorproduct(1) = A(2)*B(3)-A(3)*B(2)
 math_vectorproduct(2) = A(3)*B(1)-A(1)*B(3)
 math_vectorproduct(3) = A(1)*B(2)-A(2)*B(1)

 endfunction math_vectorproduct


!**************************************************************************
! tensor product a \otimes b
!**************************************************************************
 pure function math_tensorproduct(A,B)  

 implicit none

 real(pReal), dimension(3), intent(in) ::  A,B
 real(pReal), dimension(3,3) ::  math_tensorproduct
 integer(pInt) :: i,j
 
 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_tensorproduct(i,j) = A(i)*B(j)

 endfunction math_tensorproduct


!**************************************************************************
! matrix multiplication 3x3 = 1
!**************************************************************************
 pure function math_mul3x3(A,B)  

 implicit none

 integer(pInt) :: i
 real(pReal), dimension(3), intent(in) ::  A,B
 real(pReal), dimension(3) ::              C
 real(pReal) :: math_mul3x3

 forall (i=1_pInt:3_pInt) C(i) = A(i)*B(i)
 math_mul3x3 = sum(C)

 endfunction math_mul3x3


!**************************************************************************
! matrix multiplication 6x6 = 1
!**************************************************************************
 pure function math_mul6x6(A,B)  

 implicit none

 integer(pInt) :: i
 real(pReal), dimension(6), intent(in) ::  A,B
 real(pReal), dimension(6) ::              C
 real(pReal) :: math_mul6x6

 forall (i=1_pInt:6_pInt) C(i) = A(i)*B(i)
 math_mul6x6 = sum(C)

 endfunction math_mul6x6

 
!**************************************************************************
! matrix multiplication 33x33 = 1 (double contraction --> ij * ij)
!**************************************************************************
 pure function math_mul33xx33(A,B)  

 implicit none

 integer(pInt) :: i,j
 real(pReal), dimension(3,3), intent(in) ::  A,B
 real(pReal), dimension(3,3) ::              C
 real(pReal) :: math_mul33xx33

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) C(i,j) = A(i,j) * B(i,j)
 math_mul33xx33 = sum(C)

 endfunction math_mul33xx33

 
!**************************************************************************
! matrix multiplication 3333x33 = 33 (double contraction --> ijkl *kl = ij)
!**************************************************************************
 pure function math_mul3333xx33(A,B)  

 implicit none

 integer(pInt) :: i,j
 real(pReal), dimension(3,3,3,3), intent(in) ::  A
 real(pReal), dimension(3,3), intent(in) ::  B
 real(pReal), dimension(3,3) :: math_mul3333xx33

 do i = 1_pInt,3_pInt
   do j = 1_pInt,3_pInt
     math_mul3333xx33(i,j) = sum(A(i,j,1:3,1:3)*B(1:3,1:3))
 enddo; enddo

 endfunction math_mul3333xx33


!**************************************************************************
! matrix multiplication 33x33 = 33
!**************************************************************************
 pure function math_mul33x33(A,B)  

 implicit none

 integer(pInt) :: i,j
 real(pReal), dimension(3,3), intent(in) ::  A,B
 real(pReal), dimension(3,3) ::  math_mul33x33

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_mul33x33(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j)

 endfunction math_mul33x33


!**************************************************************************
! matrix multiplication 66x66 = 66
!**************************************************************************
 pure function math_mul66x66(A,B)  

 implicit none

 integer(pInt) :: i,j
 real(pReal), dimension(6,6), intent(in) ::  A,B
 real(pReal), dimension(6,6) ::  math_mul66x66

 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) math_mul66x66(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) + &
   A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)

 endfunction math_mul66x66

 
!**************************************************************************
! matrix multiplication 99x99 = 99
!**************************************************************************
 pure function math_mul99x99(A,B)  

 use prec, only: pReal, pInt
 implicit none

 integer(pInt)  i,j
 real(pReal), dimension(9,9), intent(in) ::  A,B

 real(pReal), dimension(9,9) ::  math_mul99x99


 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_mul99x99(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) + &
   A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j) + &
   A(i,7)*B(7,j) + A(i,8)*B(8,j) + A(i,9)*B(9,j)

 endfunction math_mul99x99

 
!**************************************************************************
! matrix multiplication 33x3 = 3
!**************************************************************************
 pure function math_mul33x3(A,B)  

 implicit none

 integer(pInt) :: i
 real(pReal), dimension(3,3), intent(in) ::  A
 real(pReal), dimension(3),   intent(in) ::  B
 real(pReal), dimension(3) ::  math_mul33x3

 forall (i=1_pInt:3_pInt) math_mul33x3(i) = sum(A(i,1:3)*B)

 endfunction math_mul33x3
 
 !**************************************************************************
! matrix multiplication complex(33) x real(3) = complex(3)
!**************************************************************************
 pure function math_mul33x3_complex(A,B)  

 implicit none

 integer(pInt) :: i
 complex(pReal), dimension(3,3), intent(in) ::  A
 real(pReal),    dimension(3),   intent(in) ::  B
 complex(pReal), dimension(3) ::  math_mul33x3_complex

 forall (i=1_pInt:3_pInt) math_mul33x3_complex(i) = sum(A(i,1:3)*B)

 endfunction math_mul33x3_complex

 
!**************************************************************************
! matrix multiplication 66x6 = 6
!**************************************************************************
 pure function math_mul66x6(A,B)  

 implicit none

 integer(pInt) :: i
 real(pReal), dimension(6,6), intent(in) ::  A
 real(pReal), dimension(6),   intent(in) ::  B
 real(pReal), dimension(6) ::  math_mul66x6

 forall (i=1_pInt:6_pInt) math_mul66x6(i) = &
   A(i,1)*B(1) + A(i,2)*B(2) + A(i,3)*B(3) + &
   A(i,4)*B(4) + A(i,5)*B(5) + A(i,6)*B(6)

 endfunction math_mul66x6

 
!**************************************************************************
! random quaternion
!**************************************************************************
 function math_qRnd()  

 implicit none

 real(pReal), dimension(4) :: math_qRnd
 real(pReal), dimension(3) :: rnd
 
 call halton(3,rnd)
 math_qRnd(1) = cos(2.0_pReal*pi*rnd(1))*sqrt(rnd(3))
 math_qRnd(2) = sin(2.0_pReal*pi*rnd(2))*sqrt(1.0_pReal-rnd(3))
 math_qRnd(3) = cos(2.0_pReal*pi*rnd(2))*sqrt(1.0_pReal-rnd(3))
 math_qRnd(4) = sin(2.0_pReal*pi*rnd(1))*sqrt(rnd(3))

 endfunction math_qRnd

 
!**************************************************************************
! quaternion multiplication q1xq2 = q12
!**************************************************************************
 pure function math_qMul(A,B)  

 implicit none

 real(pReal), dimension(4), intent(in) ::  A, B
 real(pReal), dimension(4) ::  math_qMul

 math_qMul(1) = A(1)*B(1) - A(2)*B(2) - A(3)*B(3) - A(4)*B(4)
 math_qMul(2) = A(1)*B(2) + A(2)*B(1) + A(3)*B(4) - A(4)*B(3)
 math_qMul(3) = A(1)*B(3) - A(2)*B(4) + A(3)*B(1) + A(4)*B(2)
 math_qMul(4) = A(1)*B(4) + A(2)*B(3) - A(3)*B(2) + A(4)*B(1)

 endfunction math_qMul

 
!**************************************************************************
! quaternion dotproduct
!**************************************************************************
 pure function math_qDot(A,B)  

 implicit none

 real(pReal), dimension(4), intent(in) :: A, B
 real(pReal) :: math_qDot

 math_qDot = A(1)*B(1) + A(2)*B(2) + A(3)*B(3) + A(4)*B(4)

 endfunction math_qDot

 
!**************************************************************************
! quaternion conjugation
!**************************************************************************
 pure function math_qConj(Q)  

 implicit none

 real(pReal), dimension(4), intent(in) ::  Q
 real(pReal), dimension(4) ::  math_qConj

 math_qConj(1) = Q(1)
 math_qConj(2:4) = -Q(2:4)

 endfunction math_qConj

 
!**************************************************************************
! quaternion norm
!**************************************************************************
 pure function math_qNorm(Q)  

 implicit none

 real(pReal), dimension(4), intent(in) ::  Q
 real(pReal) :: math_qNorm
 
 math_qNorm = sqrt(max(0.0_pReal, Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4)))

 endfunction math_qNorm


!**************************************************************************
! quaternion inversion
!**************************************************************************
 pure function math_qInv(Q)  

 implicit none

 real(pReal), dimension(4), intent(in) ::  Q
 real(pReal), dimension(4) ::  math_qInv
 real(pReal) :: squareNorm
 
 math_qInv = 0.0_pReal
 
 squareNorm = math_qDot(Q,Q)
 if (squareNorm > tiny(squareNorm)) &
   math_qInv = math_qConj(Q) / squareNorm
 
 endfunction math_qInv

 
!**************************************************************************
! action of a quaternion on a vector (rotate vector v with Q)
!**************************************************************************
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
 
 math_qRot(1) = -v(1)*(T(3,3)+T(4,4)) + v(2)*(T(3,2)-T(4,1)) + v(3)*(T(4,2)+T(3,1))
 math_qRot(2) =  v(1)*(T(3,2)+T(4,1)) - v(2)*(T(2,2)+T(4,4)) + v(3)*(T(4,3)-T(2,1))
 math_qRot(3) =  v(1)*(T(4,2)-T(3,1)) + v(2)*(T(4,3)+T(2,1)) - v(3)*(T(2,2)+T(3,3))
 
 math_qRot = 2.0_pReal * math_qRot + v

 endfunction math_qRot

 
!**************************************************************************
! transposition of a 33 matrix
!**************************************************************************
pure function math_transpose33(A)

 implicit none

 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal),dimension(3,3) :: math_transpose33
 integer(pInt) :: i,j
 
 forall(i=1_pInt:3_pInt, j=1_pInt:3_pInt) math_transpose33(i,j) = A(j,i)

 endfunction math_transpose33
 

!**************************************************************************
! Cramer inversion of 33 matrix (function)
!**************************************************************************
 pure function math_inv33(A)

!   direct Cramer inversion of matrix A.
!   returns all zeroes if not possible, i.e. if det close to zero

 implicit none

 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal) :: DetA
 real(pReal),dimension(3,3) :: math_inv33
 
 math_inv33 = 0.0_pReal

 DetA =   A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2))&
        - A(1,2) * (A(2,1) * A(3,3) - A(2,3) * A(3,1))&
        + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))

 if (abs(DetA) > tiny(abs(DetA))) then
   math_inv33(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2)) / DetA
   math_inv33(2,1) = (-A(2,1) * A(3,3) + A(2,3) * A(3,1)) / DetA
   math_inv33(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1)) / DetA

   math_inv33(1,2) = (-A(1,2) * A(3,3) + A(1,3) * A(3,2)) / DetA
   math_inv33(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1)) / DetA
   math_inv33(3,2) = (-A(1,1) * A(3,2) + A(1,2) * A(3,1)) / DetA

   math_inv33(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2)) / DetA
   math_inv33(2,3) = (-A(1,1) * A(2,3) + A(1,3) * A(2,1)) / DetA
   math_inv33(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1)) / DetA
 endif

 endfunction math_inv33


!**************************************************************************
! Cramer inversion of 33 matrix (subroutine)
!**************************************************************************
 PURE SUBROUTINE math_invert33(A, InvA, DetA, error)

!   Bestimmung der Determinanten und Inversen einer 33-Matrix
!   A      = Matrix A
!   InvA   = Inverse of A
!   DetA   = Determinant of A
!   error  = logical

 implicit none

 logical, intent(out) :: error
 real(pReal),dimension(3,3),intent(in)  :: A
 real(pReal),dimension(3,3),intent(out) :: InvA
 real(pReal), intent(out) :: DetA

 DetA =   A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2))&
        - A(1,2) * (A(2,1) * A(3,3) - A(2,3) * A(3,1))&
        + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))

 if (abs(DetA) <= tiny(abs(DetA))) then
   error = .true.
 else
   InvA(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2)) / DetA
   InvA(2,1) = (-A(2,1) * A(3,3) + A(2,3) * A(3,1)) / DetA
   InvA(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1)) / DetA

   InvA(1,2) = (-A(1,2) * A(3,3) + A(1,3) * A(3,2)) / DetA
   InvA(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1)) / DetA
   InvA(3,2) = (-A(1,1) * A(3,2) + A(1,2) * A(3,1)) / DetA

   InvA(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2)) / DetA
   InvA(2,3) = (-A(1,1) * A(2,3) + A(1,3) * A(2,1)) / DetA
   InvA(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1)) / DetA
   
   error = .false.
 endif

 ENDSUBROUTINE math_invert33


!**************************************************************************
! Gauss elimination to invert matrix of arbitrary dimension
!**************************************************************************
 PURE SUBROUTINE math_invert(dimen,A, InvA, AnzNegEW, error)

!   Invertieren einer dimen x dimen - Matrix
!   A        = Matrix A
!   InvA     = Inverse von A
!   AnzNegEW = Anzahl der negativen Eigenwerte von A
!   error    = logical
!              = false: Inversion wurde durchgefuehrt.
!              = true:  Die Inversion in SymGauss wurde wegen eines verschwindenen
!                       Pivotelement abgebrochen.

 implicit none

 integer(pInt), intent(in) :: dimen
 real(pReal), dimension(dimen,dimen), intent(in)  :: A
 real(pReal), dimension(dimen,dimen), intent(out) :: InvA
 integer(pInt), intent(out) :: AnzNegEW
 logical, intent(out) :: error
 real(pReal) :: LogAbsDetA
 real(pReal), dimension(dimen,dimen) :: B

 InvA = math_identity2nd(dimen)
 B = A
 CALL Gauss(dimen,B,InvA,LogAbsDetA,AnzNegEW,error)

 ENDSUBROUTINE math_invert


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 PURE SUBROUTINE Gauss (dimen,A,B,LogAbsDetA,NegHDK,error)

!   Loesung eines linearen Gleichungsssystem A * X = B mit Hilfe des
!   GAUSS-Algorithmus
!   Zur numerischen Stabilisierung wird eine Zeilen- und Spaltenpivotsuche
!   durchgefuehrt.
!
!   Eingabeparameter:
!   A(dimen,dimen) = Koeffizientenmatrix A
!   B(dimen,dimen) = rechte Seiten B
!
!   Ausgabeparameter:
!   B(dimen,dimen)    = Matrix der Unbekanntenvektoren X
!   LogAbsDetA    = 10-Logarithmus des Betrages der Determinanten von A
!   NegHDK        = Anzahl der negativen Hauptdiagonalkoeffizienten nach der
!                   Vorwaertszerlegung
!   error         = logical
!                   = false: Das Gleichungssystem wurde geloest.
!                   = true : Matrix A ist singulaer.
!
!   A und B werden veraendert!

 implicit none

 logical, intent(out) :: error
 integer(pInt), intent(in) :: dimen
 integer(pInt), intent(out) :: NegHDK
 real(pReal), intent(out) :: LogAbsDetA
 real(pReal), intent(inout), dimension(dimen,dimen) ::  A, B
 logical :: SortX
 integer(pInt) :: PivotZeile, PivotSpalte, StoreI, I, IP1, J, K, L
 integer(pInt), dimension(dimen) :: XNr
 real(pReal) :: AbsA, PivotWert, EpsAbs, Quote
 real(pReal), dimension(dimen) :: StoreA, StoreB

 error = .true.; NegHDK = 1_pInt; SortX = .false.

!   Unbekanntennumerierung

 DO  I = 1_pInt, dimen
    XNr(I) = I
 ENDDO

!   Genauigkeitsschranke und Bestimmung des groessten Pivotelementes

 PivotWert   = ABS(A(1,1))
 PivotZeile  = 1_pInt
 PivotSpalte = 1_pInt

 do  I = 1_pInt, dimen; do  J = 1_pInt, dimen
        AbsA = ABS(A(I,J))
        IF (AbsA .GT. PivotWert) THEN
            PivotWert   = AbsA
            PivotZeile  = I
            PivotSpalte = J
        ENDIF
 enddo; enddo

 IF (PivotWert .LT. 0.0000001_pReal) RETURN   ! Pivotelement = 0?

 EpsAbs = PivotWert * 0.1_pReal ** PRECISION(1.0_pReal)

!   V O R W A E R T S T R I A N G U L A T I O N

 DO  I = 1_pInt, dimen - 1_pInt
!     Zeilentausch?
    IF (PivotZeile .NE. I) THEN
        StoreA(I:dimen)       = A(I,I:dimen)
        A(I,I:dimen)          = A(PivotZeile,I:dimen)
        A(PivotZeile,I:dimen) = StoreA(I:dimen)
        StoreB(1:dimen)        = B(I,1:dimen)
        B(I,1:dimen)           = B(PivotZeile,1:dimen)
        B(PivotZeile,1:dimen)  = StoreB(1:dimen)
        SortX                = .TRUE.
    ENDIF
!     Spaltentausch?
    IF (PivotSpalte .NE. I) THEN
        StoreA(1:dimen)        = A(1:dimen,I)
        A(1:dimen,I)           = A(1:dimen,PivotSpalte)
        A(1:dimen,PivotSpalte) = StoreA(1:dimen)
        StoreI                = XNr(I)
        XNr(I)                = XNr(PivotSpalte)
        XNr(PivotSpalte)      = StoreI
        SortX                 = .TRUE.
    ENDIF
!     Triangulation
    DO  J = I + 1_pInt, dimen
        Quote = A(J,I) / A(I,I)
        DO  K = I + 1_pInt, dimen
            A(J,K) = A(J,K) - Quote * A(I,K)
        ENDDO
        DO  K = 1_pInt, dimen
            B(J,K) = B(J,K) - Quote * B(I,K)
        ENDDO
    ENDDO
!     Bestimmung des groessten Pivotelementes
    IP1         = I + 1_pInt
    PivotWert   = ABS(A(IP1,IP1))
    PivotZeile  = IP1
    PivotSpalte = IP1
    DO  J = IP1, dimen
        DO  K = IP1, dimen
            AbsA = ABS(A(J,K))
            IF (AbsA .GT. PivotWert) THEN
                PivotWert   = AbsA
                PivotZeile  = J
                PivotSpalte = K
            ENDIF
        ENDDO
    ENDDO

    IF (PivotWert .LT. EpsAbs) RETURN   ! Pivotelement = 0?

 ENDDO

!   R U E C K W A E R T S A U F L O E S U N G

 DO  I = dimen, 1_pInt, -1_pInt
    DO  L = 1_pInt, dimen
        DO  J = I + 1_pInt, dimen
            B(I,L) = B(I,L) - A(I,J) * B(J,L)
        ENDDO
        B(I,L) = B(I,L) / A(I,I)
    ENDDO
 ENDDO

!   Sortieren der Unbekanntenvektoren?

 IF (SortX) THEN
    DO  L = 1_pInt, dimen
        StoreA(1:dimen) = B(1:dimen,L)
        DO  I = 1_pInt, dimen
            J      = XNr(I)
            B(J,L) = StoreA(I)
        ENDDO
    ENDDO
 ENDIF

!   Determinante

 LogAbsDetA = 0.0_pReal
 NegHDK     = 0_pInt

 DO  I = 1_pInt, dimen
    IF (A(I,I) .LT. 0.0_pReal) NegHDK = NegHDK + 1_pInt
    AbsA       = ABS(A(I,I))
    LogAbsDetA = LogAbsDetA + LOG10(AbsA)
 ENDDO

 error = .false.

 ENDSUBROUTINE Gauss


!********************************************************************
! symmetrize a 33 matrix
!********************************************************************
 function math_symmetric33(m)

 implicit none

 real(pReal), dimension(3,3) :: math_symmetric33
 real(pReal), dimension(3,3), intent(in) :: m
 integer(pInt) :: i,j
 
 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_symmetric33(i,j) = 0.5_pReal * (m(i,j) + m(j,i))

 endfunction math_symmetric33
 

!********************************************************************
! symmetrize a 66 matrix
!********************************************************************
 pure function math_symmetric66(m)

 implicit none

 integer(pInt) :: i,j
 real(pReal), dimension(6,6), intent(in) :: m
 real(pReal), dimension(6,6) :: math_symmetric66
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) math_symmetric66(i,j) = 0.5_pReal * (m(i,j) + m(j,i))

 endfunction math_symmetric66

 
!********************************************************************
! skew part of a 33 matrix
!********************************************************************
 function math_skew33(m)

 implicit none

 real(pReal), dimension(3,3) :: math_skew33
 real(pReal), dimension(3,3), intent(in) :: m
 integer(pInt) :: i,j
 
 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_skew33(i,j) = m(i,j) - 0.5_pReal * (m(i,j) + m(j,i))

 endfunction math_skew33

 
!********************************************************************
! deviatoric part of a 33 matrix
!********************************************************************
function math_deviatoric33(m)

implicit none

real(pReal), dimension(3,3) :: math_deviatoric33
real(pReal), dimension(3,3), intent(in) :: m
integer(pInt) :: i
real(pReal) :: hydrostatic

hydrostatic = (m(1,1) + m(2,2) + m(3,3)) / 3.0_pReal
math_deviatoric33 = m
forall (i=1_pInt:3_pInt) math_deviatoric33(i,i) = m(i,i) - hydrostatic

endfunction math_deviatoric33


!********************************************************************
! equivalent scalar quantity of a full strain tensor
!********************************************************************
 pure function math_equivStrain33(m)

 implicit none

 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: math_equivStrain33,e11,e22,e33,s12,s23,s31

 e11 = (2.0_pReal*m(1,1)-m(2,2)-m(3,3))/3.0_pReal
 e22 = (2.0_pReal*m(2,2)-m(3,3)-m(1,1))/3.0_pReal
 e33 = (2.0_pReal*m(3,3)-m(1,1)-m(2,2))/3.0_pReal
 s12 = 2.0_pReal*m(1,2)
 s23 = 2.0_pReal*m(2,3)
 s31 = 2.0_pReal*m(3,1)

 math_equivStrain33 = 2.0_pReal*(1.50_pReal*(e11**2.0_pReal+e22**2.0_pReal+e33**2.0_pReal) + &
                                 0.75_pReal*(s12**2.0_pReal+s23**2.0_pReal+s31**2.0_pReal))**(0.5_pReal)/3.0_pReal

 endfunction math_equivStrain33

!********************************************************************
 subroutine math_equivStrain33_field(res,tensor,vm)
!********************************************************************
!calculate von Mises equivalent of tensor field
!
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in),  dimension(res(1),res(2),res(3),3,3) :: tensor
 ! output variables
 real(pReal), intent(out),  dimension(res(1),res(2),res(3)) :: vm
 ! other variables
 integer(pInt) :: i, j, k
 real(pReal), dimension(3,3) :: deviator, delta = 0.0_pReal
 real(pReal) :: J_2
 
 delta(1,1) = 1.0_pReal
 delta(2,2) = 1.0_pReal
 delta(3,3) = 1.0_pReal
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   deviator = tensor(i,j,k,1:3,1:3) - 1.0_pReal/3.0_pReal*tensor(i,j,k,1,1)*tensor(i,j,k,2,2)*tensor(i,j,k,3,3)*delta
   J_2 = deviator(1,1)*deviator(2,2)&
       + deviator(2,2)*deviator(3,3)&
       + deviator(1,1)*deviator(3,3)&
       - (deviator(1,2))**2.0_pReal&
       - (deviator(2,3))**2.0_pReal&
       - (deviator(1,3))**2.0_pReal
   vm(i,j,k) = sqrt(3.0_pReal*J_2)
 enddo; enddo; enddo

 end subroutine math_equivStrain33_field


!********************************************************************
! determinant of a 33 matrix
!********************************************************************
 pure function math_det33(m)

 implicit none

 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: math_det33

 math_det33 = m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
              -m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
              +m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))

 endfunction math_det33

 
!********************************************************************
! norm of a 33 matrix
!********************************************************************
 pure function math_norm33(m)

 implicit none

 real(pReal), dimension(3,3), intent(in) :: m
 real(pReal) :: math_norm33

 math_norm33 = sqrt(sum(m**2.0_pReal))

 endfunction

 
!********************************************************************
! euclidic norm of a 3 vector
!********************************************************************
 pure function math_norm3(v)

 implicit none

 real(pReal), dimension(3), intent(in) :: v
 real(pReal) :: math_norm3

 math_norm3 = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
 
 endfunction math_norm3

 
!********************************************************************
! convert 33 matrix into vector 9
!********************************************************************
 pure function math_Plain33to9(m33)

 implicit none

 real(pReal), dimension(3,3), intent(in) :: m33
 real(pReal), dimension(9) :: math_Plain33to9
 integer(pInt) :: i
 
 forall (i=1_pInt:9_pInt) math_Plain33to9(i) = m33(mapPlain(1,i),mapPlain(2,i))

 endfunction math_Plain33to9
 
 
!********************************************************************
! convert Plain 9 back to 33 matrix
!********************************************************************
 pure function math_Plain9to33(v9)

 implicit none

 real(pReal), dimension(9), intent(in) :: v9
 real(pReal), dimension(3,3) :: math_Plain9to33
 integer(pInt) :: i
 
 forall (i=1_pInt:9_pInt) math_Plain9to33(mapPlain(1,i),mapPlain(2,i)) = v9(i)

 endfunction math_Plain9to33
 

!********************************************************************
! convert symmetric 33 matrix into Mandel vector 6
!********************************************************************
 pure function math_Mandel33to6(m33)

 implicit none

 real(pReal), dimension(3,3), intent(in) :: m33
 real(pReal), dimension(6) :: math_Mandel33to6
 integer(pInt) :: i
 
 forall (i=1_pInt:6_pInt) math_Mandel33to6(i) = nrmMandel(i)*m33(mapMandel(1,i),mapMandel(2,i))

 endfunction math_Mandel33to6


!********************************************************************
! convert Mandel 6 back to symmetric 33 matrix
!********************************************************************
 pure function math_Mandel6to33(v6)

 implicit none

 real(pReal), dimension(6), intent(in) :: v6
 real(pReal), dimension(3,3) :: math_Mandel6to33
 integer(pInt) :: i
 
 forall (i=1_pInt:6_pInt)
  math_Mandel6to33(mapMandel(1,i),mapMandel(2,i)) = invnrmMandel(i)*v6(i)
  math_Mandel6to33(mapMandel(2,i),mapMandel(1,i)) = invnrmMandel(i)*v6(i)
 end forall

 endfunction math_Mandel6to33


!********************************************************************
! convert 3333 tensor into plain matrix 99
!********************************************************************
 pure function math_Plain3333to99(m3333)

 implicit none

 real(pReal), dimension(3,3,3,3), intent(in) :: m3333
 real(pReal), dimension(9,9) :: math_Plain3333to99
 integer(pInt) :: i,j
 
 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_Plain3333to99(i,j) = &
   m3333(mapPlain(1,i),mapPlain(2,i),mapPlain(1,j),mapPlain(2,j))

 endfunction math_Plain3333to99
 
!********************************************************************
! plain matrix 99 into 3333 tensor
!********************************************************************
 pure function math_Plain99to3333(m99)

 implicit none

 real(pReal), dimension(9,9), intent(in) :: m99
 real(pReal), dimension(3,3,3,3) :: math_Plain99to3333
 integer(pInt) :: i,j
 
 forall (i=1_pInt:9_pInt,j=1_pInt:9_pInt) math_Plain99to3333(mapPlain(1,i),mapPlain(2,i),&
     mapPlain(1,j),mapPlain(2,j)) = m99(i,j)

 endfunction math_Plain99to3333


!********************************************************************
! convert Mandel matrix 66 into Plain matrix 66
!********************************************************************
 pure function math_Mandel66toPlain66(m66)

 implicit none

 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(6,6) :: math_Mandel66toPlain66
 integer(pInt) :: i,j
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) &
   math_Mandel66toPlain66(i,j) = invnrmMandel(i) * invnrmMandel(j) * m66(i,j)
 return

 endfunction


!********************************************************************
! convert Plain matrix 66 into Mandel matrix 66
!********************************************************************
 pure function math_Plain66toMandel66(m66)

 implicit none

 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(6,6) :: math_Plain66toMandel66
 integer(pInt) i,j
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) &
   math_Plain66toMandel66(i,j) = nrmMandel(i) * nrmMandel(j) * m66(i,j)
 return

 endfunction


!********************************************************************
! convert symmetric 3333 tensor into Mandel matrix 66
!********************************************************************
 pure function math_Mandel3333to66(m3333)

 implicit none

 real(pReal), dimension(3,3,3,3), intent(in) :: m3333
 real(pReal), dimension(6,6) :: math_Mandel3333to66
 integer(pInt) :: i,j
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) math_Mandel3333to66(i,j) = &
   nrmMandel(i)*nrmMandel(j)*m3333(mapMandel(1,i),mapMandel(2,i),mapMandel(1,j),mapMandel(2,j))

 endfunction math_Mandel3333to66


!********************************************************************
! convert Mandel matrix 66 back to symmetric 3333 tensor
!********************************************************************
 pure function math_Mandel66to3333(m66)

 implicit none

 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(3,3,3,3) :: math_Mandel66to3333
 integer(pInt) :: i,j
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) 
   math_Mandel66to3333(mapMandel(1,i),mapMandel(2,i),mapMandel(1,j),mapMandel(2,j)) = invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(2,i),mapMandel(1,i),mapMandel(1,j),mapMandel(2,j)) = invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(1,i),mapMandel(2,i),mapMandel(2,j),mapMandel(1,j)) = invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
   math_Mandel66to3333(mapMandel(2,i),mapMandel(1,i),mapMandel(2,j),mapMandel(1,j)) = invnrmMandel(i)*invnrmMandel(j)*m66(i,j)
 end forall

 endfunction math_Mandel66to3333


!********************************************************************
! convert Voigt matrix 66 back to symmetric 3333 tensor
!********************************************************************
 pure function math_Voigt66to3333(m66)

 implicit none

 real(pReal), dimension(6,6), intent(in) :: m66
 real(pReal), dimension(3,3,3,3) :: math_Voigt66to3333
 integer(pInt) :: i,j
 
 forall (i=1_pInt:6_pInt,j=1_pInt:6_pInt) 
   math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(1,j),mapVoigt(2,j)) = invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(1,j),mapVoigt(2,j)) = invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(1,i),mapVoigt(2,i),mapVoigt(2,j),mapVoigt(1,j)) = invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
   math_Voigt66to3333(mapVoigt(2,i),mapVoigt(1,i),mapVoigt(2,j),mapVoigt(1,j)) = invnrmVoigt(i)*invnrmVoigt(j)*m66(i,j)
 end forall

 endfunction math_Voigt66to3333


!********************************************************************
! Euler angles (in radians) from rotation matrix
!********************************************************************
 pure function math_RtoEuler(R)

 implicit none

 real(pReal), dimension (3,3), intent(in) :: R
 real(pReal), dimension(3) :: math_RtoEuler
 real(pReal) :: sqhkl, squvw, sqhk, val

 sqhkl=sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3)+R(3,3)*R(3,3))
 squvw=sqrt(R(1,1)*R(1,1)+R(2,1)*R(2,1)+R(3,1)*R(3,1))
 sqhk=sqrt(R(1,3)*R(1,3)+R(2,3)*R(2,3))
! calculate PHI
 val=R(3,3)/sqhkl
 
 if(val >  1.0_pReal) val =  1.0_pReal
 if(val < -1.0_pReal) val = -1.0_pReal
     
 math_RtoEuler(2) = acos(val)

 if(math_RtoEuler(2) < 1.0e-8_pReal) then
! calculate phi2
     math_RtoEuler(3) = 0.0_pReal
! calculate phi1
     val=R(1,1)/squvw
     if(val >  1.0_pReal) val =  1.0_pReal
     if(val < -1.0_pReal) val = -1.0_pReal
     
     math_RtoEuler(1) = acos(val)
     if(R(2,1) > 0.0_pReal) math_RtoEuler(1) = 2.0_pReal*pi-math_RtoEuler(1)
 else
! calculate phi2
     val=R(2,3)/sqhk
     if(val >  1.0_pReal) val =  1.0_pReal
     if(val < -1.0_pReal) val = -1.0_pReal
     
     math_RtoEuler(3) = acos(val)
     if(R(1,3) < 0.0) math_RtoEuler(3) = 2.0_pReal*pi-math_RtoEuler(3)
! calculate phi1
     val=-R(3,2)/sin(math_RtoEuler(2))
     if(val >  1.0_pReal) val =  1.0_pReal
     if(val < -1.0_pReal) val = -1.0_pReal
     
     math_RtoEuler(1) = acos(val)
     if(R(3,1) < 0.0) math_RtoEuler(1) = 2.0_pReal*pi-math_RtoEuler(1)
 end if
 
 endfunction math_RtoEuler


!********************************************************************
! quaternion (w+ix+jy+kz) from orientation matrix
!********************************************************************
! math adopted from http://code.google.com/p/mtex/source/browse/trunk/geometry/geometry_tools/mat2quat.m
 pure function math_RtoQuaternion(R)

 implicit none

 real(pReal), dimension (3,3), intent(in) :: R
 real(pReal), dimension(4)   :: absQ, math_RtoQuaternion
 real(pReal) :: max_absQ
 integer(pInt), dimension(1) :: largest

 absQ(1) = 1.0_pReal+R(1,1)+R(2,2)+R(3,3)
 absQ(2) = 1.0_pReal+R(1,1)-R(2,2)-R(3,3)
 absQ(3) = 1.0_pReal-R(1,1)+R(2,2)-R(3,3)
 absQ(4) = 1.0_pReal-R(1,1)-R(2,2)+R(3,3)
 math_RtoQuaternion = 0.0_pReal

 largest = maxloc(absQ)

 max_absQ=0.5_pReal * sqrt(absQ(largest(1))) 

 select case(largest(1))
   case (1_pInt)
      !1----------------------------------
      math_RtoQuaternion(2) = R(2,3)-R(3,2)
      math_RtoQuaternion(3) = R(3,1)-R(1,3)
      math_RtoQuaternion(4) = R(1,2)-R(2,1)
   
   case (2_pInt)
      math_RtoQuaternion(1) = R(2,3)-R(3,2)
      !2----------------------------------
      math_RtoQuaternion(3) = R(1,2)+R(2,1)
      math_RtoQuaternion(4) = R(3,1)+R(1,3)
   
   case (3_pInt)
      math_RtoQuaternion(1) = R(3,1)-R(1,3)
      math_RtoQuaternion(2) = R(1,2)+R(2,1)
      !3----------------------------------
      math_RtoQuaternion(4) = R(2,3)+R(3,2)
   
   case (4_pInt)
      math_RtoQuaternion (1) = R(1,2)-R(2,1)
      math_RtoQuaternion (2) = R(3,1)+R(1,3)
      math_RtoQuaternion (3) = R(3,2)+R(2,3)
      !4----------------------------------
 end select

 math_RtoQuaternion = math_RtoQuaternion*0.25_pReal/max_absQ
 math_RtoQuaternion(largest(1)) = max_absQ
 
 endfunction math_RtoQuaternion


!****************************************************************
! rotation matrix from Euler angles (in radians)
!****************************************************************
 pure function math_EulerToR(Euler)

 implicit none

 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_EulerToR
 real(pReal) c1, c, c2, s1, s, s2

 C1 = cos(Euler(1))
 C = cos(Euler(2))
 C2 = cos(Euler(3))
 S1 = sin(Euler(1))
 S = sin(Euler(2))
 S2 = sin(Euler(3))

 math_EulerToR(1,1)=C1*C2-S1*S2*C
 math_EulerToR(1,2)=S1*C2+C1*S2*C
 math_EulerToR(1,3)=S2*S
 math_EulerToR(2,1)=-C1*S2-S1*C2*C
 math_EulerToR(2,2)=-S1*S2+C1*C2*C
 math_EulerToR(2,3)=C2*S
 math_EulerToR(3,1)=S1*S
 math_EulerToR(3,2)=-C1*S
 math_EulerToR(3,3)=C
 
 endfunction math_EulerToR
 

!********************************************************************
! quaternion (w+ix+jy+kz) from 3-1-3 Euler angles (in radians)
!********************************************************************
 pure function math_EulerToQuaternion(eulerangles)

 implicit none

 real(pReal), dimension(3), intent(in) :: eulerangles
 real(pReal), dimension(4) :: math_EulerToQuaternion
 real(pReal), dimension(3) :: halfangles
 real(pReal) :: c, s
 
 halfangles = 0.5_pReal * eulerangles
 
 c = cos(halfangles(2))
 s = sin(halfangles(2))
 
 math_EulerToQuaternion(1) = cos(halfangles(1)+halfangles(3)) * c
 math_EulerToQuaternion(2) = cos(halfangles(1)-halfangles(3)) * s
 math_EulerToQuaternion(3) = sin(halfangles(1)-halfangles(3)) * s
 math_EulerToQuaternion(4) = sin(halfangles(1)+halfangles(3)) * c
  
 endfunction math_EulerToQuaternion


!****************************************************************
! rotation matrix from axis and angle (in radians)  
!****************************************************************
 pure function math_AxisAngleToR(axis,omega)

 implicit none

 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: axisNrm
 real(pReal), dimension(3,3) :: math_AxisAngleToR
 real(pReal) :: norm,s,c,c1
 integer(pInt) :: i

 norm = sqrt(math_mul3x3(axis,axis))
 if (norm > 1.0e-8_pReal) then                             ! non-zero rotation
   forall (i=1_pInt:3_pInt) axisNrm(i) = axis(i)/norm      ! normalize axis to be sure

   s = sin(omega)
   c = cos(omega)
   c1 = 1.0_pReal - c
  
   ! formula for active rotation taken from http://mathworld.wolfram.com/RodriguesRotationFormula.html
   ! below is transposed form to get passive rotation
  
   math_AxisAngleToR(1,1) =  c + c1*axisNrm(1)**2.0_pReal
   math_AxisAngleToR(2,1) = -s*axisNrm(3) + c1*axisNrm(1)*axisNrm(2) 
   math_AxisAngleToR(3,1) =  s*axisNrm(2) + c1*axisNrm(1)*axisNrm(3)
  
   math_AxisAngleToR(1,2) =  s*axisNrm(3) + c1*axisNrm(2)*axisNrm(1)
   math_AxisAngleToR(2,2) =  c + c1*axisNrm(2)**2.0_pReal
   math_AxisAngleToR(3,2) = -s*axisNrm(1) + c1*axisNrm(2)*axisNrm(3)
  
   math_AxisAngleToR(1,3) = -s*axisNrm(2) + c1*axisNrm(3)*axisNrm(1)
   math_AxisAngleToR(2,3) =  s*axisNrm(1) + c1*axisNrm(3)*axisNrm(2)
   math_AxisAngleToR(3,3) =  c + c1*axisNrm(3)**2.0_pReal
 else
   math_AxisAngleToR = math_I3
 endif
 

 endfunction math_AxisAngleToR


!****************************************************************
! quaternion (w+ix+jy+kz) from axis and angle (in radians)  
!****************************************************************
 pure function math_AxisAngleToQuaternion(axis,omega)

 implicit none

 real(pReal), dimension(3), intent(in) :: axis
 real(pReal), intent(in) :: omega
 real(pReal), dimension(3) :: axisNrm
 real(pReal), dimension(4) :: math_AxisAngleToQuaternion
 real(pReal) :: s,c,norm
 integer(pInt) :: i

 norm = sqrt(math_mul3x3(axis,axis))
 if (norm > 1.0e-8_pReal) then                       ! non-zero rotation
   forall (i=1_pInt:3_pInt) axisNrm(i) = axis(i)/norm          ! normalize axis to be sure
   ! formula taken from http://en.wikipedia.org/wiki/Rotation_representation_%28mathematics%29#Rodrigues_parameters
   s = sin(omega/2.0_pReal)
   c = cos(omega/2.0_pReal)
   math_AxisAngleToQuaternion(1) =   c
   math_AxisAngleToQuaternion(2:4) = s * axisNrm(1:3)
 else
   math_AxisAngleToQuaternion = (/1.0_pReal,0.0_pReal,0.0_pReal,0.0_pReal/)   ! no rotation
 endif

 endfunction math_AxisAngleToQuaternion


!********************************************************************
! orientation matrix from quaternion (w+ix+jy+kz)
!********************************************************************
 pure function math_QuaternionToR(Q)

 implicit none

 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3,3) :: math_QuaternionToR, T,S
 integer(pInt) :: i, j
 
 forall (i = 1_pInt:3_pInt, j = 1_pInt:3_pInt) &
   T(i,j) = Q(i+1_pInt) * Q(j+1_pInt)
 S = reshape( (/0.0_pReal,     Q(4),    -Q(3), &
                    -Q(4),0.0_pReal,    +Q(2), &
                     Q(3),    -Q(2),0.0_pReal/),(/3,3/))  ! notation is transposed!

 math_QuaternionToR = (2.0_pReal * Q(1)*Q(1) - 1.0_pReal) * math_I3 + &
                      2.0_pReal * T - &
                      2.0_pReal * Q(1) * S

 endfunction math_QuaternionToR


!********************************************************************
! 3-1-3 Euler angles (in radians) from quaternion (w+ix+jy+kz)
!********************************************************************
 pure function math_QuaternionToEuler(Q)

 implicit none

 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3) :: math_QuaternionToEuler
 real(pReal) :: acos_arg

 math_QuaternionToEuler(2) = acos(1.0_pReal-2.0_pReal*(Q(2)*Q(2)+Q(3)*Q(3)))

 if (abs(math_QuaternionToEuler(2)) < 1.0e-3_pReal) then
   acos_arg=Q(1)
   if(acos_arg > 1.0_pReal)acos_arg = 1.0_pReal 
   if(acos_arg < -1.0_pReal)acos_arg = -1.0_pReal 
   math_QuaternionToEuler(1) = 2.0_pReal*acos(acos_arg)
   math_QuaternionToEuler(3) = 0.0_pReal
 else
   math_QuaternionToEuler(1) = atan2(Q(1)*Q(3)+Q(2)*Q(4), Q(1)*Q(2)-Q(3)*Q(4))
   if (math_QuaternionToEuler(1) < 0.0_pReal) &
     math_QuaternionToEuler(1) = math_QuaternionToEuler(1) + 2.0_pReal * pi

   math_QuaternionToEuler(3) = atan2(-Q(1)*Q(3)+Q(2)*Q(4), Q(1)*Q(2)+Q(3)*Q(4))
   if (math_QuaternionToEuler(3) < 0.0_pReal) &
     math_QuaternionToEuler(3) = math_QuaternionToEuler(3) + 2.0_pReal * pi
 endif

 if (math_QuaternionToEuler(2) < 0.0_pReal) &
   math_QuaternionToEuler(2) = math_QuaternionToEuler(2) + pi

 endfunction math_QuaternionToEuler


!********************************************************************
! axis-angle (x, y, z, ang in radians) from quaternion (w+ix+jy+kz)
!********************************************************************
 pure function math_QuaternionToAxisAngle(Q)

 implicit none

 real(pReal), dimension(4), intent(in) :: Q
 real(pReal) :: halfAngle, sinHalfAngle
 real(pReal), dimension(4) :: math_QuaternionToAxisAngle  

 halfAngle = acos(max(-1.0_pReal, min(1.0_pReal, Q(1))))            ! limit to [-1,1] --> 0 to 180 deg
 sinHalfAngle = sin(halfAngle)
 
 if (sinHalfAngle <= 1.0e-4_pReal) then                              ! very small rotation angle?
   math_QuaternionToAxisAngle = 0.0_pReal
 else
   math_QuaternionToAxisAngle(1:3) = Q(2:4)/sinHalfAngle
   math_QuaternionToAxisAngle(4) = halfAngle*2.0_pReal
 endif

 endfunction math_QuaternionToAxisAngle


!********************************************************************
! Rodrigues vector (x, y, z) from unit quaternion (w+ix+jy+kz)
!********************************************************************
 pure function math_QuaternionToRodrig(Q)

 use prec, only: DAMASK_NaN
 implicit none

 real(pReal), dimension(4), intent(in) :: Q
 real(pReal), dimension(3) :: math_QuaternionToRodrig

 if (Q(1) /= 0.0_pReal) then                                   ! unless rotation by 180 deg
   math_QuaternionToRodrig = Q(2:4)/Q(1)
 else
   math_QuaternionToRodrig = DAMASK_NaN                        ! NaN since Rodrig is unbound for 180 deg...
 endif

 endfunction math_QuaternionToRodrig


!**************************************************************************
! misorientation angle between two sets of Euler angles
!**************************************************************************
 pure function math_EulerMisorientation(EulerA,EulerB)

 implicit none

 real(pReal), dimension(3), intent(in) :: EulerA,EulerB
 real(pReal), dimension(3,3) :: r
 real(pReal) :: math_EulerMisorientation, tr

 r = math_mul33x33(math_EulerToR(EulerB),transpose(math_EulerToR(EulerA)))

 tr = (r(1,1)+r(2,2)+r(3,3)-1.0_pReal)*0.4999999_pReal
 math_EulerMisorientation = abs(0.5_pReal*pi-asin(tr))

 endfunction math_EulerMisorientation


!**************************************************************************
! figures whether unit quat falls into stereographic standard triangle
!**************************************************************************
pure function math_QuaternionInSST(Q, symmetryType)

  implicit none

  !*** input variables 
  real(pReal), dimension(4), intent(in) ::      Q                           ! orientation
  integer(pInt), intent(in) ::                  symmetryType                ! Type of crystal symmetry; 1:cubic, 2:hexagonal

  !*** output variables
  logical           ::                          math_QuaternionInSST
  
  !*** local variables
  real(pReal), dimension(3) ::                  Rodrig                      ! Rodrigues vector of Q
 
  Rodrig = math_QuaternionToRodrig(Q)
  select case (symmetryType)
    case (1_pInt)
      math_QuaternionInSST = Rodrig(1) > Rodrig(2) .and. &
                             Rodrig(2) > Rodrig(3) .and. &
                             Rodrig(3) > 0.0_pReal
    case (2_pInt)
      math_QuaternionInSST = Rodrig(1) > sqrt(3.0_pReal)*Rodrig(2) .and. &
                             Rodrig(2) > 0.0_pReal .and. &
                             Rodrig(3) > 0.0_pReal
    case default
      math_QuaternionInSST = .true.
  end select
  
endfunction math_QuaternionInSST


!**************************************************************************
! calculates the disorientation for 2 unit quaternions
!**************************************************************************
function math_QuaternionDisorientation(Q1, Q2, symmetryType)

  use IO,   only: IO_error
  implicit none
  
  !*** input variables 
  real(pReal), dimension(4), intent(in) ::      Q1, &                       ! 1st orientation
                                                Q2                          ! 2nd orientation
  integer(pInt), intent(in) ::                  symmetryType                ! Type of crystal symmetry; 1:cubic, 2:hexagonal
  
  !*** output variables
  real(pReal), dimension(4) ::                  math_QuaternionDisorientation         ! disorientation
  
  !*** local variables
  real(pReal), dimension(4) ::                  dQ,dQsymA,mis
  integer(pInt)    ::                           i,j,k,s
  
  dQ = math_qMul(math_qConj(Q1),Q2)
  math_QuaternionDisorientation = dQ
    
  select case (symmetryType)
    case (0_pInt)
      if (math_QuaternionDisorientation(1) < 0.0_pReal) &
        math_QuaternionDisorientation = -math_QuaternionDisorientation          ! keep omega within 0 to 180 deg
    
    case (1_pInt,2_pInt)
      s = sum(math_NsymOperations(1:symmetryType-1_pInt))
      do i = 1_pInt,2_pInt
        dQ = math_qConj(dQ)                                     ! switch order of "from -- to"
        do j = 1_pInt,math_NsymOperations(symmetryType)              ! run through first crystal's symmetries
          dQsymA = math_qMul(math_symOperations(1:4,s+j),dQ)      ! apply sym
          do k = 1_pInt,math_NsymOperations(symmetryType)            ! run through 2nd crystal's symmetries
            mis = math_qMul(dQsymA,math_symOperations(1:4,s+k))   ! apply sym
            if (mis(1) < 0.0_pReal) &                           ! want positive angle
              mis = -mis
            if (mis(1)-math_QuaternionDisorientation(1) > -1e-8_pReal .and. &
                math_QuaternionInSST(mis,symmetryType)) &
              math_QuaternionDisorientation = mis               ! found better one
      enddo; enddo; enddo
  
    case default
      call IO_error(550_pInt,symmetryType)                           ! complain about unknown symmetry
  end select
  
endfunction math_QuaternionDisorientation


!********************************************************************
!   draw a random sample from Euler space
!********************************************************************
 function math_sampleRandomOri()

 implicit none

 real(pReal), dimension(3) :: math_sampleRandomOri, rnd

 call halton(3_pInt,rnd)
 math_sampleRandomOri(1) = rnd(1)*2.0_pReal*pi
 math_sampleRandomOri(2) = acos(2.0_pReal*rnd(2)-1.0_pReal)
 math_sampleRandomOri(3) = rnd(3)*2.0_pReal*pi

 endfunction math_sampleRandomOri


!********************************************************************
!   draw a random sample from Gauss component
!   with noise (in radians) half-width 
!********************************************************************
 function math_sampleGaussOri(center,noise)

 implicit none

 real(pReal), dimension(3) :: math_sampleGaussOri, center, disturb
 real(pReal), dimension(3), parameter :: origin = (/0.0_pReal,0.0_pReal,0.0_pReal/)
 real(pReal), dimension(5) :: rnd
 real(pReal) :: noise,scatter,cosScatter
 integer(pInt) i

if (noise==0.0_pReal) then
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
   disturb(1) = scatter * rnd(1)                                                      ! phi1
   disturb(2) = sign(1.0_pReal,rnd(2))*acos(cosScatter+(1.0_pReal-cosScatter)*rnd(4)) ! Phi
   disturb(3) = scatter * rnd(2)                                                      ! phi2
   if (rnd(5) <= exp(-1.0_pReal*(math_EulerMisorientation(origin,disturb)/scatter)**2_pReal)) exit   
 enddo

 math_sampleGaussOri = math_RtoEuler(math_mul33x33(math_EulerToR(disturb),math_EulerToR(center)))
 
 endfunction math_sampleGaussOri
 

!********************************************************************
!   draw a random sample from Fiber component
!   with noise (in radians)
!********************************************************************
 function math_sampleFiberOri(alpha,beta,noise)

 implicit none

 real(pReal), dimension(3) :: math_sampleFiberOri, fiberInC,fiberInS,axis
 real(pReal), dimension(2) :: alpha,beta, rnd
 real(pReal), dimension(3,3) :: oRot,fRot,pRot
 real(pReal) :: noise, scatter, cos2Scatter, angle
 integer(pInt), dimension(2,3), parameter :: rotMap = reshape((/2_pInt,3_pInt,&
                                                                3_pInt,1_pInt,&
                                                                1_pInt,2_pInt/),(/2,3/))
 integer(pInt) :: i

! Helming uses different distribution with Bessel functions
! therefore the gauss scatter width has to be scaled differently
 scatter = 0.95_pReal * noise
 cos2Scatter = cos(2.0_pReal*scatter)

! fiber axis in crystal coordinate system
 fiberInC(1)=sin(alpha(1))*cos(alpha(2))
 fiberInC(2)=sin(alpha(1))*sin(alpha(2))
 fiberInC(3)=cos(alpha(1))
! fiber axis in sample coordinate system
 fiberInS(1)=sin(beta(1))*cos(beta(2))
 fiberInS(2)=sin(beta(1))*sin(beta(2))
 fiberInS(3)=cos(beta(1))

! ---# rotation matrix from sample to crystal system #---
 angle = -acos(dot_product(fiberInC,fiberInS))
 if(angle /= 0.0_pReal) then
!   rotation axis between sample and crystal system (cross product)
   forall(i=1:3) axis(i) = fiberInC(rotMap(1,i))*fiberInS(rotMap(2,i))-fiberInC(rotMap(2,i))*fiberInS(rotMap(1,i))
   oRot = math_AxisAngleToR(math_vectorproduct(fiberInC,fiberInS),angle)
 else
   oRot = math_I3
 end if

! ---# rotation matrix about fiber axis (random angle) #---
 call halton(1_pInt,rnd)
 fRot = math_AxisAngleToR(fiberInS,rnd(1)*2.0_pReal*pi)

! ---# rotation about random axis perpend to fiber #---
! random axis pependicular to fiber axis 
 call halton(2_pInt,axis)
 if (fiberInS(3) /= 0.0_pReal) then
     axis(3)=-(axis(1)*fiberInS(1)+axis(2)*fiberInS(2))/fiberInS(3)
 else if(fiberInS(2) /= 0.0_pReal) then
     axis(3)=axis(2)
     axis(2)=-(axis(1)*fiberInS(1)+axis(3)*fiberInS(3))/fiberInS(2)
 else if(fiberInS(1) /= 0.0_pReal) then
     axis(3)=axis(1)
     axis(1)=-(axis(2)*fiberInS(2)+axis(3)*fiberInS(3))/fiberInS(1)
 end if

! scattered rotation angle 
 do
   call halton(2_pInt,rnd)
     angle = acos(cos2Scatter+(1.0_pReal-cos2Scatter)*rnd(1))
     if (rnd(2) <= exp(-1.0_pReal*(angle/scatter)**2.0_pReal)) exit
 enddo
 call halton(1_pInt,rnd)
 if (rnd(1) <= 0.5) angle = -angle
 pRot = math_AxisAngleToR(axis,angle)

! ---# apply the three rotations #---
 math_sampleFiberOri = math_RtoEuler(math_mul33x33(pRot,math_mul33x33(fRot,oRot))) 

 endfunction math_sampleFiberOri


!********************************************************************
!   symmetric Euler angles for given symmetry string
!   'triclinic' or '', 'monoclinic', 'orthotropic'
!********************************************************************
 pure function math_symmetricEulers(sym,Euler)

 implicit none

 integer(pInt), intent(in) :: sym
 real(pReal), dimension(3), intent(in) :: Euler
 real(pReal), dimension(3,3) :: math_symmetricEulers
 integer(pInt) :: i,j
 
 math_symmetricEulers(1,1) = pi+Euler(1)
 math_symmetricEulers(2,1) = Euler(2)
 math_symmetricEulers(3,1) = Euler(3)

 math_symmetricEulers(1,2) = pi-Euler(1)
 math_symmetricEulers(2,2) = pi-Euler(2)
 math_symmetricEulers(3,2) = pi+Euler(3)

 math_symmetricEulers(1,3) = 2.0_pReal*pi-Euler(1)
 math_symmetricEulers(2,3) = pi-Euler(2)
 math_symmetricEulers(3,3) = pi+Euler(3)

 forall (i=1_pInt:3_pInt,j=1_pInt:3_pInt) math_symmetricEulers(j,i) = modulo(math_symmetricEulers(j,i),2.0_pReal*pi)

 select case (sym)
   case (4_pInt) ! all done

   case (2_pInt)  ! return only first
     math_symmetricEulers(1:3,2:3) = 0.0_pReal

   case default         ! return blank
     math_symmetricEulers = 0.0_pReal
 end select

 endfunction math_symmetricEulers


!********************************************************************
!   draw a random sample from Gauss variable
!********************************************************************
function math_sampleGaussVar(meanvalue, stddev, width)

implicit none

!*** input variables
real(pReal), intent(in) ::            meanvalue, &      ! meanvalue of gauss distribution
                                      stddev            ! standard deviation of gauss distribution
real(pReal), intent(in), optional ::  width             ! width of considered values as multiples of standard deviation

!*** output variables
real(pReal) ::                        math_sampleGaussVar

!*** local variables
real(pReal), dimension(2) ::          rnd               ! random numbers
real(pReal) ::                        scatter, &        ! normalized scatter around meanvalue
                                      myWidth

if (stddev == 0.0_pReal) then
    math_sampleGaussVar = meanvalue
    return
endif

if (present(width)) then
  myWidth = width
else
  myWidth = 3.0_pReal                                         ! use +-3*sigma as default value for scatter
endif

do
  call halton(2_pInt, rnd)
  scatter = myWidth * (2.0_pReal * rnd(1) - 1.0_pReal)
  if (rnd(2) <= exp(-0.5_pReal * scatter ** 2.0_pReal)) &     ! test if scattered value is drawn
    exit
enddo

math_sampleGaussVar = scatter * stddev

endfunction math_sampleGaussVar


!****************************************************************
subroutine math_spectralDecompositionSym33(M,values,vectors,error)
!****************************************************************
 implicit none

 real(pReal), dimension(3,3), intent(in) :: M
 real(pReal), dimension(3),   intent(out) :: values
 real(pReal), dimension(3,3), intent(out) :: vectors
 logical, intent(out) :: error

 integer(pInt) :: info
 real(pReal), dimension((64+2)*3) :: work                          ! block size of 64 taken from http://www.netlib.org/lapack/double/dsyev.f
 
 vectors = M                                                       ! copy matrix to input (doubles as output) array
 call DSYEV('V','U',3,vectors,3,values,work,(64+2)*3,info)
 error = (info == 0_pInt)
 
 return
end subroutine


!****************************************************************
 pure subroutine math_pDecomposition(FE,U,R,error)
!-----FE = R.U 
!****************************************************************
 implicit none

 real(pReal), intent(in), dimension(3,3) :: FE
 real(pReal), intent(out), dimension(3,3) :: R, U
 logical, intent(out) :: error
 real(pReal), dimension(3,3) :: CE, EB1, EB2, EB3, UI
 real(pReal) :: EW1, EW2, EW3, det

 error = .false.
 ce = math_mul33x33(math_transpose33(FE),FE)

 CALL math_spectral1(CE,EW1,EW2,EW3,EB1,EB2,EB3)
 U=sqrt(EW1)*EB1+sqrt(EW2)*EB2+sqrt(EW3)*EB3
 call math_invert33(U,UI,det,error)
 if (.not. error) R = math_mul33x33(FE,UI)

 ENDSUBROUTINE math_pDecomposition


!**********************************************************************
 pure subroutine math_spectral1(M,EW1,EW2,EW3,EB1,EB2,EB3)
!**** EIGENWERTE UND EIGENWERTBASIS DER SYMMETRISCHEN 3X3 MATRIX M

 implicit none

 real(pReal), dimension(3,3), intent(in) :: M
 real(pReal), dimension(3,3), intent(out) :: EB1, EB2, EB3
 real(pReal), intent(out) :: EW1,EW2,EW3
 real(pReal) HI1M, HI2M, HI3M, R, S, T, P, Q, RHO, PHI, Y1, Y2, Y3, D1, D2, D3
 real(pReal), parameter :: TOL=1.e-14_pReal
 real(pReal), dimension(3,3) :: M1, M2, M3
 real(pReal) C1,C2,C3,arg

 CALL math_hi(M,HI1M,HI2M,HI3M)
 R=-HI1M
 S= HI2M
 T=-HI3M
 P=S-R**2.0_pReal/3.0_pReal
 Q=2.0_pReal/27.0_pReal*R**3.0_pReal-R*S/3.0_pReal+T
 EB1=0.0_pReal
 EB2=0.0_pReal
 EB3=0.0_pReal
 IF((ABS(P).LT.TOL).AND.(ABS(Q).LT.TOL))THEN
!   DREI GLEICHE EIGENWERTE
   EW1=HI1M/3.0_pReal
   EW2=EW1
   EW3=EW1
!   this is not really correct, but this way U is calculated
!   correctly in PDECOMPOSITION (correct is EB?=I)
   EB1(1,1)=1.0_pReal
   EB2(2,2)=1.0_pReal
   EB3(3,3)=1.0_pReal
 ELSE
   RHO=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   arg=-Q/RHO/2.0_pReal
   if(arg.GT.1.0_pReal) arg=1.0_pReal
   if(arg.LT.-1.0_pReal) arg=-1.0_pReal
   PHI=acos(arg)
   Y1=2.0_pReal*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal)
   Y2=2.0_pReal*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal+2.0_pReal/3.0_pReal*PI)
   Y3=2.0_pReal*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal+4.0_pReal/3.0_pReal*PI)
   EW1=Y1-R/3.0_pReal
   EW2=Y2-R/3.0_pReal
   EW3=Y3-R/3.0_pReal
   C1=ABS(EW1-EW2)
   C2=ABS(EW2-EW3) 
   C3=ABS(EW3-EW1)

   IF(C1.LT.TOL) THEN
!  EW1 is equal to EW2
  D3=1.0_pReal/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  EB3=math_mul33x33(M1,M2)*D3

  EB1=math_I3-EB3
!  both EB2 and EW2 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW2=0.0_pReal
   ELSE IF(C2.LT.TOL) THEN
!  EW2 is equal to EW3
  D1=1.0_pReal/(EW1-EW2)/(EW1-EW3)
  M2=M-math_I3*EW2
  M3=M-math_I3*EW3
  EB1=math_mul33x33(M2,M3)*D1
  EB2=math_I3-EB1
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0_pReal
   ELSE IF(C3.LT.TOL) THEN
!  EW1 is equal to EW3
  D2=1.0_pReal/(EW2-EW1)/(EW2-EW3) 
  M1=M-math_I3*EW1
  M3=M-math_I3*EW3
  EB2=math_mul33x33(M1,M3)*D2
  EB1=math_I3-EB2
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0_pReal
   ELSE
!  all three eigenvectors are different
  D1=1.0_pReal/(EW1-EW2)/(EW1-EW3)
  D2=1.0_pReal/(EW2-EW1)/(EW2-EW3) 
  D3=1.0_pReal/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  M3=M-EW3*math_I3
  EB1=math_mul33x33(M2,M3)*D1
  EB2=math_mul33x33(M1,M3)*D2
  EB3=math_mul33x33(M1,M2)*D3

   END IF
 END IF

 ENDSUBROUTINE math_spectral1


!**********************************************************************
 function math_eigenvalues33(M)
!**** Eigenvalues of symmetric 3X3 matrix M

 implicit none

 real(pReal), intent(in), dimension(3,3) :: M
 real(pReal), dimension(3,3) :: EB1 = 0.0_pReal, EB2 = 0.0_pReal, EB3 = 0.0_pReal
 real(pReal), dimension(3) :: math_eigenvalues33
 real(pReal) :: HI1M, HI2M, HI3M, R, S, T, P, Q, RHO, PHI, Y1, Y2, Y3, arg
 real(pReal), parameter :: TOL=1.e-14_pReal

 CALL math_hi(M,HI1M,HI2M,HI3M)
 R=-HI1M
 S= HI2M
 T=-HI3M
 P=S-R**2.0_pReal/3.0_pReal
 Q=2.0_pReal/27.0_pReal*R**3.0_pReal-R*S/3.0_pReal+T
 
 if((abs(P) < TOL) .and. (abs(Q) < TOL)) THEN
! three equivalent eigenvalues
   math_eigenvalues33(1) = HI1M/3.0_pReal
   math_eigenvalues33(2)=math_eigenvalues33(1)
   math_eigenvalues33(3)=math_eigenvalues33(1)
!   this is not really correct, but this way U is calculated
!   correctly in PDECOMPOSITION (correct is EB?=I)
   EB1(1,1)=1.0_pReal
   EB2(2,2)=1.0_pReal
   EB3(3,3)=1.0_pReal
 else
   RHO=sqrt(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   arg=-Q/RHO/2.0_pReal
   if(arg.GT.1.0_pReal) arg=1.0_pReal
   if(arg.LT.-1.0_pReal) arg=-1.0_pReal
   PHI=acos(arg)
   Y1=2*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal)
   Y2=2*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal+2.0_pReal/3.0_pReal*PI)
   Y3=2*RHO**(1.0_pReal/3.0_pReal)*cos(PHI/3.0_pReal+4.0_pReal/3.0_pReal*PI)
   math_eigenvalues33(1) = Y1-R/3.0_pReal
   math_eigenvalues33(2) = Y2-R/3.0_pReal
   math_eigenvalues33(3) = Y3-R/3.0_pReal
 endif
 endfunction  math_eigenvalues33


!********************************************************************** 
!**** HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M

 PURE SUBROUTINE math_hi(M,HI1M,HI2M,HI3M)
 
 implicit none

 real(pReal), intent(in) :: M(3,3) 
 real(pReal), intent(out) :: HI1M, HI2M, HI3M 

 HI1M=M(1,1)+M(2,2)+M(3,3)
 HI2M=HI1M**2.0_pReal/2.0_pReal-  (M(1,1)**2.0_pReal+M(2,2)**2.0_pReal+M(3,3)**2.0_pReal)&
                                     /2.0_pReal-M(1,2)*M(2,1)-M(1,3)*M(3,1)-M(2,3)*M(3,2)
 HI3M=math_det33(M)
! QUESTION: is 3rd equiv det(M) ?? if yes, use function math_det !agreed on YES

 ENDSUBROUTINE math_hi


!*******************************************************************************
! GET_SEED returns a seed for the random number generator.
!
!  The seed depends on the current time, and ought to be (slightly)
!  different every millisecond. Once the seed is obtained, a random
!  number generator should be called a few times to further process
!  the seed.
!
!  Parameters:
!  Output, integer SEED, a pseudorandom seed value.
!
!  Modified: 27 June 2000
!  Author:  John Burkardt
!
!  Modified: 29 April 2005
!  Author: Franz Roters
!
 SUBROUTINE get_seed(seed)
 implicit none

 integer(pInt) :: seed
 real(pReal) ::  temp = 0.0_pReal
 character(len = 10) :: time
 character(len = 8) :: today
 integer(pInt) :: values(8)
 character(len = 5) :: zone

 call date_and_time (today, time, zone, values)

 temp = temp + real(values(2)- 1_pInt, pReal) / 11.0_pReal
 temp = temp + real(values(3)- 1_pInt, pReal) / 30.0_pReal
 temp = temp + real(values(5),         pReal) / 23.0_pReal
 temp = temp + real(values(6),         pReal) / 59.0_pReal
 temp = temp + real(values(7),         pReal) / 59.0_pReal
 temp = temp + real(values(8),         pReal) / 999.0_pReal
 temp = temp / 6.0_pReal

 if (temp <= 0.0_pReal) then
   temp = 1.0_pReal / 3.0_pReal
 else if (1.0_pReal <= temp) then
   temp = 2.0_pReal / 3.0_pReal
 end if

 seed = int(real(huge(1_pInt),pReal)*temp, pInt)
!
!  Never use a seed of 0 or maximum integer.
!
 if (seed == 0_pInt) then
   seed = 1_pInt
 end if

 if (seed == huge(1_pInt)) then
   seed = seed -1_pInt
 end if

 ENDSUBROUTINE get_seed


!*******************************************************************************
! HALTON computes the next element in the Halton sequence.
!
!  Parameters:
!  Input, integer NDIM, the dimension of the element.
!  Output, real R(NDIM), the next element of the current Halton sequence.
!
!  Modified: 09 March 2003
!  Author: John Burkardt
!
!  Modified: 29 April 2005
!  Author: Franz Roters
!
 subroutine halton(ndim, r)
 implicit none

 integer(pInt), intent(in) :: ndim
 real(pReal), intent(out), dimension(ndim) :: r
 integer(pInt), dimension(ndim) :: base
 integer(pInt) :: seed
 integer(pInt), dimension(1) :: value_halton

 call halton_memory ('GET', 'SEED', 1_pInt, value_halton)
 seed = value_halton(1)

 call halton_memory ('GET', 'BASE', ndim, base)

 call i_to_halton (seed, base, ndim, r)

 value_halton(1) = 1_pInt
 call halton_memory ('INC', 'SEED', 1_pInt, value_halton)

 ENDSUBROUTINE halton


!*******************************************************************************
! HALTON_MEMORY sets or returns quantities associated with the Halton sequence.
!
!  Parameters:
!  Input, character (len = *) action_halton, the desired action.
!  'GET' means get the value of a particular quantity.
!  'SET' means set the value of a particular quantity.
!  'INC' means increment the value of a particular quantity.
!   (Only the SEED can be incremented.)
!
!  Input, character (len = *) name_halton, the name of the quantity.
!  'BASE' means the Halton base or bases.
!  'NDIM' means the spatial dimension.
!  'SEED' means the current Halton seed.
!
!  Input/output, integer NDIM, the dimension of the quantity.
!  If action_halton is 'SET' and action_halton is 'BASE', then NDIM is input, and
!  is the number of entries in value_halton to be put into BASE.
!
!  Input/output, integer value_halton(NDIM), contains a value.
!  If action_halton is 'SET', then on input, value_halton contains values to be assigned
!  to the internal variable.
!  If action_halton is 'GET', then on output, value_halton contains the values of
!  the specified internal variable.
!  If action_halton is 'INC', then on input, value_halton contains the increment to
!  be added to the specified internal variable.
!
!  Modified: 09 March 2003
!  Author: John Burkardt
!
!  Modified: 29 April 2005
!  Author:  Franz Roters

 subroutine halton_memory (action_halton, name_halton, ndim, value_halton)
 implicit none

 character(len = *), intent(in) :: action_halton, name_halton
 integer(pInt), dimension(*), intent(inout) :: value_halton
 integer(pInt), allocatable, save, dimension(:) :: base
 logical, save :: first_call = .true.
 integer(pInt), intent(in) :: ndim
 integer(pInt):: i
 integer(pInt), save :: ndim_save = 0_pInt, seed = 1_pInt


 if (first_call) then
   ndim_save = 1_pInt
   allocate(base(ndim_save))
   base(1) = 2_pInt
   first_call = .false.
 endif
!
!  Set
!
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
!
!  Get
!
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
!
!  Increment
!
 elseif(action_halton(1:1) == 'I' .or. action_halton(1:1) == 'i') then
   if(name_halton(1:1) == 'S' .or. name_halton(1:1) == 's') then
     seed = seed + value_halton(1)
   end if
 endif

 ENDSUBROUTINE halton_memory


!*******************************************************************************
! HALTON_NDIM_SET sets the dimension for a Halton sequence.
!
!  Parameters:
!  Input, integer NDIM, the dimension of the Halton vectors.
!
!  Modified: 26 February 2001
!  Author: John Burkardt
!
!  Modified: 29 April 2005
!  Author: Franz Roters
!
 subroutine halton_ndim_set (ndim)
 implicit none

 integer(pInt), intent(in) :: ndim
 integer(pInt) :: value_halton(1)

 value_halton(1) = ndim
 call halton_memory ('SET', 'NDIM', 1_pInt, value_halton)

 ENDSUBROUTINE halton_ndim_set


!*******************************************************************************
! HALTON_SEED_SET sets the "seed" for the Halton sequence.
!
!  Calling HALTON repeatedly returns the elements of the
!  Halton sequence in order, starting with element number 1.
!  An internal counter, called SEED, keeps track of the next element
!  to return. Each time the routine is called, the SEED-th element
!  is computed, and then SEED is incremented by 1.
! 
!  To restart the Halton sequence, it is only necessary to reset
!  SEED to 1. It might also be desirable to reset SEED to some other value.
!  This routine allows the user to specify any value of SEED.
! 
!  The default value of SEED is 1, which restarts the Halton sequence.
!
!  Parameters:
!  Input, integer SEED, the seed for the Halton sequence.
!
!  Modified: 26 February 2001
!  Author: John Burkardt
!
!  Modified: 29 April 2005
!  Author: Franz Roters
!
 subroutine halton_seed_set (seed)
 implicit none

 integer(pInt), parameter :: ndim = 1_pInt
 integer(pInt), intent(in) :: seed
 integer(pInt) :: value_halton(ndim)

 value_halton(1) = seed
 call halton_memory ('SET', 'SEED', ndim, value_halton)

 ENDSUBROUTINE halton_seed_set


!*******************************************************************************
! I_TO_HALTON computes an element of a Halton sequence.
!
!  Reference:
!  J H Halton: On the efficiency of certain quasi-random sequences of points
!  in evaluating multi-dimensional integrals, Numerische Mathematik, Volume 2, pages 84-90, 1960.
! 
!  Parameters:
!  Input, integer SEED, the index of the desired element.
!  Only the absolute value of SEED is considered. SEED = 0 is allowed,
!  and returns R = 0.
! 
!  Input, integer BASE(NDIM), the Halton bases, which should be
!  distinct prime numbers.  This routine only checks that each base
!  is greater than 1.
! 
!  Input, integer NDIM, the dimension of the sequence.
! 
!  Output, real R(NDIM), the SEED-th element of the Halton sequence
!  for the given bases.
!
!  Modified: 26 February 2001
!  Author: John Burkardt
!
!  Modified: 29 April 2005
!  Author: Franz RotersA

 subroutine i_to_halton (seed, base, ndim, r)
 implicit none

 integer(pInt), intent(in) :: ndim
 integer(pInt), intent(in), dimension(ndim) :: base
 real(pReal), dimension(ndim) ::  base_inv
 integer(pInt), dimension(ndim) :: digit
 integer(pInt) :: i
 real(pReal), dimension(ndim), intent(out) ::r
 integer(pInt) :: seed
 integer(pInt), dimension(ndim) :: seed2

 seed2(1:ndim) = abs(seed)

 r(1:ndim) = 0.0_pReal

 if (any (base(1:ndim) <= 1_pInt)) then
   !$OMP CRITICAL (write2out)
   write (*, '(a)') ' '
   write (*, '(a)') 'I_TO_HALTON - Fatal error!'
   write (*, '(a)') ' An input base BASE is <= 1!'
   do i = 1, ndim
     write (*, '(i6,i6)') i, base(i)
   enddo
   call flush(6)
   !$OMP END CRITICAL (write2out)
   stop
 end if

 base_inv(1:ndim) = 1.0_pReal / real (base(1:ndim), pReal)

 do while ( any ( seed2(1:ndim) /= 0_pInt) )
   digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim))
   r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), pReal) * base_inv(1:ndim)
   base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), pReal)
   seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
 enddo

 ENDSUBROUTINE i_to_halton


!*******************************************************************************
! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Note:
!  PRIME_MAX is 1500, and the largest prime stored is 12553.
!  Reference:
!  Milton Abramowitz and Irene Stegun: Handbook of Mathematical Functions,
!  US Department of Commerce, 1964, pages 870-873.
! 
!  Daniel Zwillinger: CRC Standard Mathematical Tables and Formulae,
!  30th Edition, CRC Press, 1996, pages 95-98.
! 
!  Parameters:
!  Input, integer N, the index of the desired prime number.
!  N = -1 returns PRIME_MAX, the index of the largest prime available.
!  N = 0 is legal, returning PRIME = 1.
!  It should generally be true that 0 <= N <= PRIME_MAX.
! 
!  Output, integer PRIME, the N-th prime.  If N is out of range, PRIME
!  is returned as 0.
! 
!  Modified:  21 June 2002
!  Author: John Burkardt
! 
!  Modified: 29 April 2005
!  Author: Franz Roters
!
 function prime(n)
 implicit none

 integer(pInt), parameter :: prime_max = 1500
 integer(pInt), save :: icall = 0_pInt
 integer(pInt), intent(in) :: n
 integer(pInt), save, dimension(prime_max) :: npvec
 integer(pInt)  prime

 if (icall == 0_pInt) then
   icall = 1_pInt
  
   npvec(1:100) = (/&
          2_pInt,    3_pInt,    5_pInt,    7_pInt,   11_pInt,   13_pInt,   17_pInt,   19_pInt,   23_pInt,   29_pInt, &
         31_pInt,   37_pInt,   41_pInt,   43_pInt,   47_pInt,   53_pInt,   59_pInt,   61_pInt,   67_pInt,   71_pInt, &
         73_pInt,   79_pInt,   83_pInt,   89_pInt,   97_pInt,  101_pInt,  103_pInt,  107_pInt,  109_pInt,  113_pInt, &
        127_pInt,  131_pInt,  137_pInt,  139_pInt,  149_pInt,  151_pInt,  157_pInt,  163_pInt,  167_pInt,  173_pInt, &
        179_pInt,  181_pInt,  191_pInt,  193_pInt,  197_pInt,  199_pInt,  211_pInt,  223_pInt,  227_pInt,  229_pInt, &
        233_pInt,  239_pInt,  241_pInt,  251_pInt,  257_pInt,  263_pInt,  269_pInt,  271_pInt,  277_pInt,  281_pInt, &
        283_pInt,  293_pInt,  307_pInt,  311_pInt,  313_pInt,  317_pInt,  331_pInt,  337_pInt,  347_pInt,  349_pInt, &
        353_pInt,  359_pInt,  367_pInt,  373_pInt,  379_pInt,  383_pInt,  389_pInt,  397_pInt,  401_pInt,  409_pInt, &
        419_pInt,  421_pInt,  431_pInt,  433_pInt,  439_pInt,  443_pInt,  449_pInt,  457_pInt,  461_pInt,  463_pInt, &
        467_pInt,  479_pInt,  487_pInt,  491_pInt,  499_pInt,  503_pInt,  509_pInt,  521_pInt,  523_pInt,  541_pInt/)
  
   npvec(101:200) = (/ &
         547_pInt,  557_pInt,  563_pInt,  569_pInt,  571_pInt,  577_pInt,  587_pInt,  593_pInt,  599_pInt,  601_pInt, &
         607_pInt,  613_pInt,  617_pInt,  619_pInt,  631_pInt,  641_pInt,  643_pInt,  647_pInt,  653_pInt,  659_pInt, &
         661_pInt,  673_pInt,  677_pInt,  683_pInt,  691_pInt,  701_pInt,  709_pInt,  719_pInt,  727_pInt,  733_pInt, &
         739_pInt,  743_pInt,  751_pInt,  757_pInt,  761_pInt,  769_pInt,  773_pInt,  787_pInt,  797_pInt,  809_pInt, &
         811_pInt,  821_pInt,  823_pInt,  827_pInt,  829_pInt,  839_pInt,  853_pInt,  857_pInt,  859_pInt,  863_pInt, &
         877_pInt,  881_pInt,  883_pInt,  887_pInt,  907_pInt,  911_pInt,  919_pInt,  929_pInt,  937_pInt,  941_pInt, &
         947_pInt,  953_pInt,  967_pInt,  971_pInt,  977_pInt,  983_pInt,  991_pInt,  997_pInt, 1009_pInt, 1013_pInt, &
        1019_pInt, 1021_pInt, 1031_pInt, 1033_pInt, 1039_pInt, 1049_pInt, 1051_pInt, 1061_pInt, 1063_pInt, 1069_pInt, &
        1087_pInt, 1091_pInt, 1093_pInt, 1097_pInt, 1103_pInt, 1109_pInt, 1117_pInt, 1123_pInt, 1129_pInt, 1151_pInt, &
        1153_pInt, 1163_pInt, 1171_pInt, 1181_pInt, 1187_pInt, 1193_pInt, 1201_pInt, 1213_pInt, 1217_pInt, 1223_pInt/)
  
   npvec(201:300) = (/ &
        1229_pInt, 1231_pInt, 1237_pInt, 1249_pInt, 1259_pInt, 1277_pInt, 1279_pInt, 1283_pInt, 1289_pInt, 1291_pInt, &
        1297_pInt, 1301_pInt, 1303_pInt, 1307_pInt, 1319_pInt, 1321_pInt, 1327_pInt, 1361_pInt, 1367_pInt, 1373_pInt, &
        1381_pInt, 1399_pInt, 1409_pInt, 1423_pInt, 1427_pInt, 1429_pInt, 1433_pInt, 1439_pInt, 1447_pInt, 1451_pInt, &
        1453_pInt, 1459_pInt, 1471_pInt, 1481_pInt, 1483_pInt, 1487_pInt, 1489_pInt, 1493_pInt, 1499_pInt, 1511_pInt, &
        1523_pInt, 1531_pInt, 1543_pInt, 1549_pInt, 1553_pInt, 1559_pInt, 1567_pInt, 1571_pInt, 1579_pInt, 1583_pInt, &
        1597_pInt, 1601_pInt, 1607_pInt, 1609_pInt, 1613_pInt, 1619_pInt, 1621_pInt, 1627_pInt, 1637_pInt, 1657_pInt, &
        1663_pInt, 1667_pInt, 1669_pInt, 1693_pInt, 1697_pInt, 1699_pInt, 1709_pInt, 1721_pInt, 1723_pInt, 1733_pInt, &
        1741_pInt, 1747_pInt, 1753_pInt, 1759_pInt, 1777_pInt, 1783_pInt, 1787_pInt, 1789_pInt, 1801_pInt, 1811_pInt, &
        1823_pInt, 1831_pInt, 1847_pInt, 1861_pInt, 1867_pInt, 1871_pInt, 1873_pInt, 1877_pInt, 1879_pInt, 1889_pInt, &
        1901_pInt, 1907_pInt, 1913_pInt, 1931_pInt, 1933_pInt, 1949_pInt, 1951_pInt, 1973_pInt, 1979_pInt, 1987_pInt/)
  
   npvec(301:400) = (/ &
        1993_pInt, 1997_pInt, 1999_pInt, 2003_pInt, 2011_pInt, 2017_pInt, 2027_pInt, 2029_pInt, 2039_pInt, 2053_pInt, &
        2063_pInt, 2069_pInt, 2081_pInt, 2083_pInt, 2087_pInt, 2089_pInt, 2099_pInt, 2111_pInt, 2113_pInt, 2129_pInt, &
        2131_pInt, 2137_pInt, 2141_pInt, 2143_pInt, 2153_pInt, 2161_pInt, 2179_pInt, 2203_pInt, 2207_pInt, 2213_pInt, &
        2221_pInt, 2237_pInt, 2239_pInt, 2243_pInt, 2251_pInt, 2267_pInt, 2269_pInt, 2273_pInt, 2281_pInt, 2287_pInt, &
        2293_pInt, 2297_pInt, 2309_pInt, 2311_pInt, 2333_pInt, 2339_pInt, 2341_pInt, 2347_pInt, 2351_pInt, 2357_pInt, &
        2371_pInt, 2377_pInt, 2381_pInt, 2383_pInt, 2389_pInt, 2393_pInt, 2399_pInt, 2411_pInt, 2417_pInt, 2423_pInt, &
        2437_pInt, 2441_pInt, 2447_pInt, 2459_pInt, 2467_pInt, 2473_pInt, 2477_pInt, 2503_pInt, 2521_pInt, 2531_pInt, &
        2539_pInt, 2543_pInt, 2549_pInt, 2551_pInt, 2557_pInt, 2579_pInt, 2591_pInt, 2593_pInt, 2609_pInt, 2617_pInt, &
        2621_pInt, 2633_pInt, 2647_pInt, 2657_pInt, 2659_pInt, 2663_pInt, 2671_pInt, 2677_pInt, 2683_pInt, 2687_pInt, &
        2689_pInt, 2693_pInt, 2699_pInt, 2707_pInt, 2711_pInt, 2713_pInt, 2719_pInt, 2729_pInt, 2731_pInt, 2741_pInt/)
  
   npvec(401:500) = (/ &
        2749_pInt, 2753_pInt, 2767_pInt, 2777_pInt, 2789_pInt, 2791_pInt, 2797_pInt, 2801_pInt, 2803_pInt, 2819_pInt, &
        2833_pInt, 2837_pInt, 2843_pInt, 2851_pInt, 2857_pInt, 2861_pInt, 2879_pInt, 2887_pInt, 2897_pInt, 2903_pInt, &
        2909_pInt, 2917_pInt, 2927_pInt, 2939_pInt, 2953_pInt, 2957_pInt, 2963_pInt, 2969_pInt, 2971_pInt, 2999_pInt, &
        3001_pInt, 3011_pInt, 3019_pInt, 3023_pInt, 3037_pInt, 3041_pInt, 3049_pInt, 3061_pInt, 3067_pInt, 3079_pInt, &
        3083_pInt, 3089_pInt, 3109_pInt, 3119_pInt, 3121_pInt, 3137_pInt, 3163_pInt, 3167_pInt, 3169_pInt, 3181_pInt, &
        3187_pInt, 3191_pInt, 3203_pInt, 3209_pInt, 3217_pInt, 3221_pInt, 3229_pInt, 3251_pInt, 3253_pInt, 3257_pInt, &
        3259_pInt, 3271_pInt, 3299_pInt, 3301_pInt, 3307_pInt, 3313_pInt, 3319_pInt, 3323_pInt, 3329_pInt, 3331_pInt, &
        3343_pInt, 3347_pInt, 3359_pInt, 3361_pInt, 3371_pInt, 3373_pInt, 3389_pInt, 3391_pInt, 3407_pInt, 3413_pInt, &
        3433_pInt, 3449_pInt, 3457_pInt, 3461_pInt, 3463_pInt, 3467_pInt, 3469_pInt, 3491_pInt, 3499_pInt, 3511_pInt, &
        3517_pInt, 3527_pInt, 3529_pInt, 3533_pInt, 3539_pInt, 3541_pInt, 3547_pInt, 3557_pInt, 3559_pInt, 3571_pInt/)
  
   npvec(501:600) = (/ &
        3581_pInt, 3583_pInt, 3593_pInt, 3607_pInt, 3613_pInt, 3617_pInt, 3623_pInt, 3631_pInt, 3637_pInt, 3643_pInt, &
        3659_pInt, 3671_pInt, 3673_pInt, 3677_pInt, 3691_pInt, 3697_pInt, 3701_pInt, 3709_pInt, 3719_pInt, 3727_pInt, &
        3733_pInt, 3739_pInt, 3761_pInt, 3767_pInt, 3769_pInt, 3779_pInt, 3793_pInt, 3797_pInt, 3803_pInt, 3821_pInt, &
        3823_pInt, 3833_pInt, 3847_pInt, 3851_pInt, 3853_pInt, 3863_pInt, 3877_pInt, 3881_pInt, 3889_pInt, 3907_pInt, &
        3911_pInt, 3917_pInt, 3919_pInt, 3923_pInt, 3929_pInt, 3931_pInt, 3943_pInt, 3947_pInt, 3967_pInt, 3989_pInt, &
        4001_pInt, 4003_pInt, 4007_pInt, 4013_pInt, 4019_pInt, 4021_pInt, 4027_pInt, 4049_pInt, 4051_pInt, 4057_pInt, &
        4073_pInt, 4079_pInt, 4091_pInt, 4093_pInt, 4099_pInt, 4111_pInt, 4127_pInt, 4129_pInt, 4133_pInt, 4139_pInt, &
        4153_pInt, 4157_pInt, 4159_pInt, 4177_pInt, 4201_pInt, 4211_pInt, 4217_pInt, 4219_pInt, 4229_pInt, 4231_pInt, &
        4241_pInt, 4243_pInt, 4253_pInt, 4259_pInt, 4261_pInt, 4271_pInt, 4273_pInt, 4283_pInt, 4289_pInt, 4297_pInt, &
        4327_pInt, 4337_pInt, 4339_pInt, 4349_pInt, 4357_pInt, 4363_pInt, 4373_pInt, 4391_pInt, 4397_pInt, 4409_pInt/)
  
   npvec(601:700) = (/ &
        4421_pInt, 4423_pInt, 4441_pInt, 4447_pInt, 4451_pInt, 4457_pInt, 4463_pInt, 4481_pInt, 4483_pInt, 4493_pInt, &
        4507_pInt, 4513_pInt, 4517_pInt, 4519_pInt, 4523_pInt, 4547_pInt, 4549_pInt, 4561_pInt, 4567_pInt, 4583_pInt, &
        4591_pInt, 4597_pInt, 4603_pInt, 4621_pInt, 4637_pInt, 4639_pInt, 4643_pInt, 4649_pInt, 4651_pInt, 4657_pInt, &
        4663_pInt, 4673_pInt, 4679_pInt, 4691_pInt, 4703_pInt, 4721_pInt, 4723_pInt, 4729_pInt, 4733_pInt, 4751_pInt, &
        4759_pInt, 4783_pInt, 4787_pInt, 4789_pInt, 4793_pInt, 4799_pInt, 4801_pInt, 4813_pInt, 4817_pInt, 4831_pInt, &
        4861_pInt, 4871_pInt, 4877_pInt, 4889_pInt, 4903_pInt, 4909_pInt, 4919_pInt, 4931_pInt, 4933_pInt, 4937_pInt, &
        4943_pInt, 4951_pInt, 4957_pInt, 4967_pInt, 4969_pInt, 4973_pInt, 4987_pInt, 4993_pInt, 4999_pInt, 5003_pInt, &
        5009_pInt, 5011_pInt, 5021_pInt, 5023_pInt, 5039_pInt, 5051_pInt, 5059_pInt, 5077_pInt, 5081_pInt, 5087_pInt, &
        5099_pInt, 5101_pInt, 5107_pInt, 5113_pInt, 5119_pInt, 5147_pInt, 5153_pInt, 5167_pInt, 5171_pInt, 5179_pInt, &
        5189_pInt, 5197_pInt, 5209_pInt, 5227_pInt, 5231_pInt, 5233_pInt, 5237_pInt, 5261_pInt, 5273_pInt, 5279_pInt/)
  
   npvec(701:800) = (/ &
        5281_pInt, 5297_pInt, 5303_pInt, 5309_pInt, 5323_pInt, 5333_pInt, 5347_pInt, 5351_pInt, 5381_pInt, 5387_pInt, &
        5393_pInt, 5399_pInt, 5407_pInt, 5413_pInt, 5417_pInt, 5419_pInt, 5431_pInt, 5437_pInt, 5441_pInt, 5443_pInt, &
        5449_pInt, 5471_pInt, 5477_pInt, 5479_pInt, 5483_pInt, 5501_pInt, 5503_pInt, 5507_pInt, 5519_pInt, 5521_pInt, &
        5527_pInt, 5531_pInt, 5557_pInt, 5563_pInt, 5569_pInt, 5573_pInt, 5581_pInt, 5591_pInt, 5623_pInt, 5639_pInt, &
        5641_pInt, 5647_pInt, 5651_pInt, 5653_pInt, 5657_pInt, 5659_pInt, 5669_pInt, 5683_pInt, 5689_pInt, 5693_pInt, &
        5701_pInt, 5711_pInt, 5717_pInt, 5737_pInt, 5741_pInt, 5743_pInt, 5749_pInt, 5779_pInt, 5783_pInt, 5791_pInt, &
        5801_pInt, 5807_pInt, 5813_pInt, 5821_pInt, 5827_pInt, 5839_pInt, 5843_pInt, 5849_pInt, 5851_pInt, 5857_pInt, &
        5861_pInt, 5867_pInt, 5869_pInt, 5879_pInt, 5881_pInt, 5897_pInt, 5903_pInt, 5923_pInt, 5927_pInt, 5939_pInt, &
        5953_pInt, 5981_pInt, 5987_pInt, 6007_pInt, 6011_pInt, 6029_pInt, 6037_pInt, 6043_pInt, 6047_pInt, 6053_pInt, &
        6067_pInt, 6073_pInt, 6079_pInt, 6089_pInt, 6091_pInt, 6101_pInt, 6113_pInt, 6121_pInt, 6131_pInt, 6133_pInt/)
  
   npvec(801:900) = (/ &
        6143_pInt, 6151_pInt, 6163_pInt, 6173_pInt, 6197_pInt, 6199_pInt, 6203_pInt, 6211_pInt, 6217_pInt, 6221_pInt, &
        6229_pInt, 6247_pInt, 6257_pInt, 6263_pInt, 6269_pInt, 6271_pInt, 6277_pInt, 6287_pInt, 6299_pInt, 6301_pInt, &
        6311_pInt, 6317_pInt, 6323_pInt, 6329_pInt, 6337_pInt, 6343_pInt, 6353_pInt, 6359_pInt, 6361_pInt, 6367_pInt, &
        6373_pInt, 6379_pInt, 6389_pInt, 6397_pInt, 6421_pInt, 6427_pInt, 6449_pInt, 6451_pInt, 6469_pInt, 6473_pInt, &
        6481_pInt, 6491_pInt, 6521_pInt, 6529_pInt, 6547_pInt, 6551_pInt, 6553_pInt, 6563_pInt, 6569_pInt, 6571_pInt, &
        6577_pInt, 6581_pInt, 6599_pInt, 6607_pInt, 6619_pInt, 6637_pInt, 6653_pInt, 6659_pInt, 6661_pInt, 6673_pInt, &
        6679_pInt, 6689_pInt, 6691_pInt, 6701_pInt, 6703_pInt, 6709_pInt, 6719_pInt, 6733_pInt, 6737_pInt, 6761_pInt, &
        6763_pInt, 6779_pInt, 6781_pInt, 6791_pInt, 6793_pInt, 6803_pInt, 6823_pInt, 6827_pInt, 6829_pInt, 6833_pInt, &
        6841_pInt, 6857_pInt, 6863_pInt, 6869_pInt, 6871_pInt, 6883_pInt, 6899_pInt, 6907_pInt, 6911_pInt, 6917_pInt, &
        6947_pInt, 6949_pInt, 6959_pInt, 6961_pInt, 6967_pInt, 6971_pInt, 6977_pInt, 6983_pInt, 6991_pInt, 6997_pInt/)
  
   npvec(901:1000) = (/ &
        7001_pInt, 7013_pInt, 7019_pInt, 7027_pInt, 7039_pInt, 7043_pInt, 7057_pInt, 7069_pInt, 7079_pInt, 7103_pInt, &
        7109_pInt, 7121_pInt, 7127_pInt, 7129_pInt, 7151_pInt, 7159_pInt, 7177_pInt, 7187_pInt, 7193_pInt, 7207_pInt, &
        7211_pInt, 7213_pInt, 7219_pInt, 7229_pInt, 7237_pInt, 7243_pInt, 7247_pInt, 7253_pInt, 7283_pInt, 7297_pInt, &
        7307_pInt, 7309_pInt, 7321_pInt, 7331_pInt, 7333_pInt, 7349_pInt, 7351_pInt, 7369_pInt, 7393_pInt, 7411_pInt, &
        7417_pInt, 7433_pInt, 7451_pInt, 7457_pInt, 7459_pInt, 7477_pInt, 7481_pInt, 7487_pInt, 7489_pInt, 7499_pInt, &
        7507_pInt, 7517_pInt, 7523_pInt, 7529_pInt, 7537_pInt, 7541_pInt, 7547_pInt, 7549_pInt, 7559_pInt, 7561_pInt, &
        7573_pInt, 7577_pInt, 7583_pInt, 7589_pInt, 7591_pInt, 7603_pInt, 7607_pInt, 7621_pInt, 7639_pInt, 7643_pInt, &
        7649_pInt, 7669_pInt, 7673_pInt, 7681_pInt, 7687_pInt, 7691_pInt, 7699_pInt, 7703_pInt, 7717_pInt, 7723_pInt, &
        7727_pInt, 7741_pInt, 7753_pInt, 7757_pInt, 7759_pInt, 7789_pInt, 7793_pInt, 7817_pInt, 7823_pInt, 7829_pInt, &
        7841_pInt, 7853_pInt, 7867_pInt, 7873_pInt, 7877_pInt, 7879_pInt, 7883_pInt, 7901_pInt, 7907_pInt, 7919_pInt/)
  
   npvec(1001:1100) = (/ &
        7927_pInt, 7933_pInt, 7937_pInt, 7949_pInt, 7951_pInt, 7963_pInt, 7993_pInt, 8009_pInt, 8011_pInt, 8017_pInt, &
        8039_pInt, 8053_pInt, 8059_pInt, 8069_pInt, 8081_pInt, 8087_pInt, 8089_pInt, 8093_pInt, 8101_pInt, 8111_pInt, &
        8117_pInt, 8123_pInt, 8147_pInt, 8161_pInt, 8167_pInt, 8171_pInt, 8179_pInt, 8191_pInt, 8209_pInt, 8219_pInt, &
        8221_pInt, 8231_pInt, 8233_pInt, 8237_pInt, 8243_pInt, 8263_pInt, 8269_pInt, 8273_pInt, 8287_pInt, 8291_pInt, &
        8293_pInt, 8297_pInt, 8311_pInt, 8317_pInt, 8329_pInt, 8353_pInt, 8363_pInt, 8369_pInt, 8377_pInt, 8387_pInt, &
        8389_pInt, 8419_pInt, 8423_pInt, 8429_pInt, 8431_pInt, 8443_pInt, 8447_pInt, 8461_pInt, 8467_pInt, 8501_pInt, &
        8513_pInt, 8521_pInt, 8527_pInt, 8537_pInt, 8539_pInt, 8543_pInt, 8563_pInt, 8573_pInt, 8581_pInt, 8597_pInt, &
        8599_pInt, 8609_pInt, 8623_pInt, 8627_pInt, 8629_pInt, 8641_pInt, 8647_pInt, 8663_pInt, 8669_pInt, 8677_pInt, &
        8681_pInt, 8689_pInt, 8693_pInt, 8699_pInt, 8707_pInt, 8713_pInt, 8719_pInt, 8731_pInt, 8737_pInt, 8741_pInt, &
        8747_pInt, 8753_pInt, 8761_pInt, 8779_pInt, 8783_pInt, 8803_pInt, 8807_pInt, 8819_pInt, 8821_pInt, 8831_pInt/)
  
   npvec(1101:1200) = (/ &
        8837_pInt, 8839_pInt, 8849_pInt, 8861_pInt, 8863_pInt, 8867_pInt, 8887_pInt, 8893_pInt, 8923_pInt, 8929_pInt, &
        8933_pInt, 8941_pInt, 8951_pInt, 8963_pInt, 8969_pInt, 8971_pInt, 8999_pInt, 9001_pInt, 9007_pInt, 9011_pInt, &
        9013_pInt, 9029_pInt, 9041_pInt, 9043_pInt, 9049_pInt, 9059_pInt, 9067_pInt, 9091_pInt, 9103_pInt, 9109_pInt, &
        9127_pInt, 9133_pInt, 9137_pInt, 9151_pInt, 9157_pInt, 9161_pInt, 9173_pInt, 9181_pInt, 9187_pInt, 9199_pInt, &
        9203_pInt, 9209_pInt, 9221_pInt, 9227_pInt, 9239_pInt, 9241_pInt, 9257_pInt, 9277_pInt, 9281_pInt, 9283_pInt, &
        9293_pInt, 9311_pInt, 9319_pInt, 9323_pInt, 9337_pInt, 9341_pInt, 9343_pInt, 9349_pInt, 9371_pInt, 9377_pInt, &
        9391_pInt, 9397_pInt, 9403_pInt, 9413_pInt, 9419_pInt, 9421_pInt, 9431_pInt, 9433_pInt, 9437_pInt, 9439_pInt, &
        9461_pInt, 9463_pInt, 9467_pInt, 9473_pInt, 9479_pInt, 9491_pInt, 9497_pInt, 9511_pInt, 9521_pInt, 9533_pInt, &
        9539_pInt, 9547_pInt, 9551_pInt, 9587_pInt, 9601_pInt, 9613_pInt, 9619_pInt, 9623_pInt, 9629_pInt, 9631_pInt, &
        9643_pInt, 9649_pInt, 9661_pInt, 9677_pInt, 9679_pInt, 9689_pInt, 9697_pInt, 9719_pInt, 9721_pInt, 9733_pInt/)
  
   npvec(1201:1300) = (/ &
         9739_pInt, 9743_pInt, 9749_pInt, 9767_pInt, 9769_pInt, 9781_pInt, 9787_pInt, 9791_pInt, 9803_pInt, 9811_pInt, &
         9817_pInt, 9829_pInt, 9833_pInt, 9839_pInt, 9851_pInt, 9857_pInt, 9859_pInt, 9871_pInt, 9883_pInt, 9887_pInt, &
         9901_pInt, 9907_pInt, 9923_pInt, 9929_pInt, 9931_pInt, 9941_pInt, 9949_pInt, 9967_pInt, 9973_pInt,10007_pInt, &
        10009_pInt,10037_pInt,10039_pInt,10061_pInt,10067_pInt,10069_pInt,10079_pInt,10091_pInt,10093_pInt,10099_pInt, &
        10103_pInt,10111_pInt,10133_pInt,10139_pInt,10141_pInt,10151_pInt,10159_pInt,10163_pInt,10169_pInt,10177_pInt, &
        10181_pInt,10193_pInt,10211_pInt,10223_pInt,10243_pInt,10247_pInt,10253_pInt,10259_pInt,10267_pInt,10271_pInt, &
        10273_pInt,10289_pInt,10301_pInt,10303_pInt,10313_pInt,10321_pInt,10331_pInt,10333_pInt,10337_pInt,10343_pInt, &
        10357_pInt,10369_pInt,10391_pInt,10399_pInt,10427_pInt,10429_pInt,10433_pInt,10453_pInt,10457_pInt,10459_pInt, &
        10463_pInt,10477_pInt,10487_pInt,10499_pInt,10501_pInt,10513_pInt,10529_pInt,10531_pInt,10559_pInt,10567_pInt, &
        10589_pInt,10597_pInt,10601_pInt,10607_pInt,10613_pInt,10627_pInt,10631_pInt,10639_pInt,10651_pInt,10657_pInt/)
  
   npvec(1301:1400) = (/ &
        10663_pInt,10667_pInt,10687_pInt,10691_pInt,10709_pInt,10711_pInt,10723_pInt,10729_pInt,10733_pInt,10739_pInt, &
        10753_pInt,10771_pInt,10781_pInt,10789_pInt,10799_pInt,10831_pInt,10837_pInt,10847_pInt,10853_pInt,10859_pInt, &
        10861_pInt,10867_pInt,10883_pInt,10889_pInt,10891_pInt,10903_pInt,10909_pInt,19037_pInt,10939_pInt,10949_pInt, &
        10957_pInt,10973_pInt,10979_pInt,10987_pInt,10993_pInt,11003_pInt,11027_pInt,11047_pInt,11057_pInt,11059_pInt, &
        11069_pInt,11071_pInt,11083_pInt,11087_pInt,11093_pInt,11113_pInt,11117_pInt,11119_pInt,11131_pInt,11149_pInt, &
        11159_pInt,11161_pInt,11171_pInt,11173_pInt,11177_pInt,11197_pInt,11213_pInt,11239_pInt,11243_pInt,11251_pInt, &
        11257_pInt,11261_pInt,11273_pInt,11279_pInt,11287_pInt,11299_pInt,11311_pInt,11317_pInt,11321_pInt,11329_pInt, &
        11351_pInt,11353_pInt,11369_pInt,11383_pInt,11393_pInt,11399_pInt,11411_pInt,11423_pInt,11437_pInt,11443_pInt, &
        11447_pInt,11467_pInt,11471_pInt,11483_pInt,11489_pInt,11491_pInt,11497_pInt,11503_pInt,11519_pInt,11527_pInt, &
        11549_pInt,11551_pInt,11579_pInt,11587_pInt,11593_pInt,11597_pInt,11617_pInt,11621_pInt,11633_pInt,11657_pInt/)
  
   npvec(1401:1500) = (/ &
        11677_pInt,11681_pInt,11689_pInt,11699_pInt,11701_pInt,11717_pInt,11719_pInt,11731_pInt,11743_pInt,11777_pInt, &
        11779_pInt,11783_pInt,11789_pInt,11801_pInt,11807_pInt,11813_pInt,11821_pInt,11827_pInt,11831_pInt,11833_pInt, &
        11839_pInt,11863_pInt,11867_pInt,11887_pInt,11897_pInt,11903_pInt,11909_pInt,11923_pInt,11927_pInt,11933_pInt, &
        11939_pInt,11941_pInt,11953_pInt,11959_pInt,11969_pInt,11971_pInt,11981_pInt,11987_pInt,12007_pInt,12011_pInt, &
        12037_pInt,12041_pInt,12043_pInt,12049_pInt,12071_pInt,12073_pInt,12097_pInt,12101_pInt,12107_pInt,12109_pInt, &
        12113_pInt,12119_pInt,12143_pInt,12149_pInt,12157_pInt,12161_pInt,12163_pInt,12197_pInt,12203_pInt,12211_pInt, &
        12227_pInt,12239_pInt,12241_pInt,12251_pInt,12253_pInt,12263_pInt,12269_pInt,12277_pInt,12281_pInt,12289_pInt, &
        12301_pInt,12323_pInt,12329_pInt,12343_pInt,12347_pInt,12373_pInt,12377_pInt,12379_pInt,12391_pInt,12401_pInt, &
        12409_pInt,12413_pInt,12421_pInt,12433_pInt,12437_pInt,12451_pInt,12457_pInt,12473_pInt,12479_pInt,12487_pInt, &
        12491_pInt,12497_pInt,12503_pInt,12511_pInt,12517_pInt,12527_pInt,12539_pInt,12541_pInt,12547_pInt,12553_pInt/)

 endif

 if(n == -1_pInt) then
   prime = prime_max
 else if (n == 0_pInt) then
   prime = 1_pInt
 else if (n <= prime_max) then
   prime = npvec(n)
 else                                 ! why not use io_error here?
   prime = 0_pInt
!$OMP CRITICAL (write2out)
   write (6, '(a)') ' '
   write (6, '(a)') 'PRIME - Fatal error!'
   write (6, '(a,i6)') '  Illegal prime index N = ', n
   write (6, '(a,i6)') '  N must be between 0 and PRIME_MAX = ', prime_max
   call flush(6)
!$OMP END CRITICAL (write2out)
   stop
 end if
 endfunction prime


!**************************************************************************
! volume of tetrahedron given by four vertices
!**************************************************************************
 pure function math_volTetrahedron(v1,v2,v3,v4)  

 implicit none

 real(pReal) math_volTetrahedron
 real(pReal), dimension (3), intent(in) :: v1,v2,v3,v4
 real(pReal), dimension (3,3) :: m

 m(1:3,1) = v1-v2
 m(1:3,2) = v2-v3
 m(1:3,3) = v3-v4

 math_volTetrahedron = math_det33(m)/6.0_pReal 

 endfunction math_volTetrahedron


!**************************************************************************
! rotate 33 tensor forward
!**************************************************************************
 pure function math_rotate_forward33(tensor,rot_tensor)

 implicit none

 real(pReal), dimension(3,3) ::  math_rotate_forward33
 real(pReal), dimension(3,3), intent(in) :: tensor, rot_tensor
 
 math_rotate_forward33 = math_mul33x33(rot_tensor,&
                         math_mul33x33(tensor,math_transpose33(rot_tensor)))
 
 endfunction math_rotate_forward33


!**************************************************************************
! rotate 33 tensor backward
!**************************************************************************
 pure function math_rotate_backward33(tensor,rot_tensor)

 implicit none

 real(pReal), dimension(3,3) ::  math_rotate_backward33
 real(pReal), dimension(3,3), intent(in) :: tensor, rot_tensor
 
 math_rotate_backward33 = math_mul33x33(math_transpose33(rot_tensor),&
                           math_mul33x33(tensor,rot_tensor))
 
 endfunction math_rotate_backward33


!**************************************************************************
! rotate 3333 tensor
! C'_ijkl=g_im*g_jn*g_ko*g_lp*C_mnop
!**************************************************************************
 pure function math_rotate_forward3333(tensor,rot_tensor)

 implicit none

 real(pReal), dimension(3,3,3,3) ::  math_rotate_forward3333
 real(pReal), dimension(3,3), intent(in) :: rot_tensor
 real(pReal), dimension(3,3,3,3), intent(in) :: tensor
 integer(pInt) :: i,j,k,l,m,n,o,p
 
 math_rotate_forward3333= 0.0_pReal

 do i = 1_pInt,3_pInt; do j = 1_pInt,3_pInt; do k = 1_pInt,3_pInt; do l = 1_pInt,3_pInt
   do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt; do o = 1_pInt,3_pInt; do p = 1_pInt,3_pInt
     math_rotate_forward3333(i,j,k,l) = tensor(i,j,k,l)+rot_tensor(m,i)*rot_tensor(n,j)*&
                                                           rot_tensor(o,k)*rot_tensor(p,l)*tensor(m,n,o,p)
 enddo; enddo; enddo; enddo; enddo; enddo; enddo; enddo
 
 endfunction math_rotate_forward3333


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Functions below are taken from the old postprocessingMath.f90
!  mostly they are used in combination with f2py to build fortran 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! put the next two funtions into mesh?
 function mesh_location(idx,resolution)
 ! small helper functions for indexing
 ! CAREFULL, index and location runs from 0 to N-1 (python style)

   integer(pInt), intent(in) :: idx
   integer(pInt), intent(in), dimension(3) :: resolution
   integer(pInt), dimension(3) :: mesh_location
   mesh_location = (/modulo(idx/ resolution(3) / resolution(2),resolution(1)), &
                     modulo(idx/ resolution(3),                resolution(2)), &
                     modulo(idx,                               resolution(3))/)

 end function mesh_location
 

 function mesh_index(location,resolution)
 ! small helper functions for indexing
 ! CAREFULL, index and location runs from 0 to N-1 (python style)
   integer(pInt), intent(in), dimension(3) :: resolution, location
   integer(pInt) :: mesh_index
   
   mesh_index = modulo(location(3), resolution(3))     +&
               (modulo(location(2), resolution(2)))*resolution(3) +&
               (modulo(location(1), resolution(1)))*resolution(3)*resolution(2)

 end function mesh_index


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 subroutine volume_compare(res,geomdim,defgrad,nodes,volume_mismatch)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to calculate the mismatch between volume of reconstructed (compatible
! cube and determinant of defgrad at the FP

 use debug, only: debug_verbosity
 implicit none
 
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal),   intent(in), dimension(3) :: geomdim
 real(pReal),   intent(in), dimension(res(1),       res(2),       res(3),       3,3) :: defgrad
 real(pReal),   intent(in), dimension(res(1)+1_pInt,res(2)+1_pInt,res(3)+1_pInt,3)   :: nodes
 ! output variables
 real(pReal),  intent(out), dimension(res(1),       res(2),       res(3))            :: volume_mismatch
 ! other variables
 real(pReal),   dimension(8,3) ::  coords
 integer(pInt) i,j,k
 real(pReal) vol_initial

 if (debug_verbosity > 0_pInt) then
   print*, 'Calculating volume mismatch'
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 vol_initial = geomdim(1)*geomdim(2)*geomdim(3)/(real(res(1)*res(2)*res(3), pReal))
 do k = 1_pInt,res(3)
   do j = 1_pInt,res(2)
     do i = 1_pInt,res(1) 
       coords(1,1:3) = nodes(i,       j,       k       ,1:3)
       coords(2,1:3) = nodes(i+1_pInt,j,       k       ,1:3)
       coords(3,1:3) = nodes(i+1_pInt,j+1_pInt,k       ,1:3)
       coords(4,1:3) = nodes(i,       j+1_pInt,k       ,1:3)
       coords(5,1:3) = nodes(i,       j,       k+1_pInt,1:3)
       coords(6,1:3) = nodes(i+1_pInt,j,       k+1_pInt,1:3)
       coords(7,1:3) = nodes(i+1_pInt,j+1_pInt,k+1_pInt,1:3)
       coords(8,1:3) = nodes(i,       j+1_pInt,k+1_pInt,1:3)
       volume_mismatch(i,j,k) = abs(math_volTetrahedron(coords(7,1:3),coords(1,1:3),coords(8,1:3),coords(4,1:3))) &
                              + abs(math_volTetrahedron(coords(7,1:3),coords(1,1:3),coords(8,1:3),coords(5,1:3))) &
                              + abs(math_volTetrahedron(coords(7,1:3),coords(1,1:3),coords(3,1:3),coords(4,1:3))) &
                              + abs(math_volTetrahedron(coords(7,1:3),coords(1,1:3),coords(3,1:3),coords(2,1:3))) &
                              + abs(math_volTetrahedron(coords(7,1:3),coords(5,1:3),coords(2,1:3),coords(6,1:3))) &
                              + abs(math_volTetrahedron(coords(7,1:3),coords(5,1:3),coords(2,1:3),coords(1,1:3)))
       volume_mismatch(i,j,k) = volume_mismatch(i,j,k)/math_det33(defgrad(i,j,k,1:3,1:3))
 enddo; enddo; enddo
 volume_mismatch = volume_mismatch/vol_initial

end subroutine volume_compare


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine shape_compare(res,geomdim,defgrad,nodes,centroids,shape_mismatch)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to calculate the mismatch between the vectors from the central point to
! the corners of reconstructed (combatible) volume element and the vectors calculated by deforming 
! the initial volume element with the  current deformation gradient 

 use debug, only: debug_verbosity
 implicit none

 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3)   :: geomdim
 real(pReal), intent(in),  dimension(res(1),       res(2),       res(3),       3,3) :: defgrad
 real(pReal), intent(in),  dimension(res(1)+1_pInt,res(2)+1_pInt,res(3)+1_pInt,3)   :: nodes
 real(pReal), intent(in),  dimension(res(1),       res(2),       res(3),       3)   :: centroids
 ! output variables
 real(pReal), intent(out), dimension(res(1),       res(2),       res(3))            :: shape_mismatch
 ! other variables
 real(pReal), dimension(8,3) :: coords_initial
 integer(pInt) i,j,k

 if (debug_verbosity > 0_pInt) then
   print*, 'Calculating shape mismatch'
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 coords_initial(1,1:3) = (/-geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           -geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           -geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(2,1:3) = (/+geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           -geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           -geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(3,1:3) = (/+geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           +geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           -geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(4,1:3) = (/-geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           +geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           -geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(5,1:3) = (/-geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           -geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           +geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(6,1:3) = (/+geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           -geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           +geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(7,1:3) = (/+geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           +geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           +geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 coords_initial(8,1:3) = (/-geomdim(1)/2.0_pReal/real(res(1),pReal),&
                           +geomdim(2)/2.0_pReal/real(res(2),pReal),&
                           +geomdim(3)/2.0_pReal/real(res(3),pReal)/)
 do i=1_pInt,8_pInt
 enddo
 do k = 1_pInt,res(3)
   do j = 1_pInt,res(2)
     do i = 1_pInt,res(1)
       shape_mismatch(i,j,k) = &
           sqrt(sum((nodes(i,       j,       k,       1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(1,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i+1_pInt,j,       k,       1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(2,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i+1_pInt,j+1_pInt,k,       1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(3,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i,       j+1_pInt,k,       1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(4,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i,       j,       k+1_pInt,1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(5,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i+1_pInt,j,       k+1_pInt,1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(6,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i+1_pInt,j+1_pInt,k+1_pInt,1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(7,1:3)))**2.0_pReal))&
         + sqrt(sum((nodes(i,       j+1_pInt,k+1_pInt,1:3) - centroids(i,j,k,1:3)&
                    - matmul(defgrad(i,j,k,1:3,1:3), coords_initial(8,1:3)))**2.0_pReal))
 enddo; enddo; enddo

 end subroutine shape_compare


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine mesh_regular_grid(res,geomdim,defgrad_av,centroids,nodes)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to build mesh of (distoreted) cubes for given coordinates (= center of the cubes)
!
 use debug, only: debug_verbosity
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3)   :: geomdim
 real(pReal), intent(in), dimension(3,3) :: defgrad_av
 real(pReal), intent(in), dimension(res(1),       res(2),       res(3),       3) :: centroids
 ! output variables
 real(pReal),intent(out), dimension(res(1)+1_pInt,res(2)+1_pInt,res(3)+1_pInt,3) :: nodes
 ! variables with dimension depending on input
 real(pReal),             dimension(res(1)+2_pInt,res(2)+2_pInt,res(3)+2_pInt,3) :: wrappedCentroids
 ! other variables
 integer(pInt) :: i,j,k,n
 integer(pInt), dimension(3), parameter :: diag = 1_pInt
 integer(pInt), dimension(3)            :: shift = 0_pInt, lookup = 0_pInt, me = 0_pInt
 integer(pInt), dimension(3,8) :: neighbor = reshape((/ &
                                     0_pInt, 0_pInt, 0_pInt, &
                                     1_pInt, 0_pInt, 0_pInt, &
                                     1_pInt, 1_pInt, 0_pInt, &
                                     0_pInt, 1_pInt, 0_pInt, &
                                     0_pInt, 0_pInt, 1_pInt, &
                                     1_pInt, 0_pInt, 1_pInt, &
                                     1_pInt, 1_pInt, 1_pInt, &
                                     0_pInt, 1_pInt, 1_pInt  &
                                    /), &
                                    (/3,8/))

 if (debug_verbosity > 0_pInt) then
   print*, 'Meshing cubes around centroids' 
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 nodes = 0.0_pReal
 wrappedCentroids = 0.0_pReal
 wrappedCentroids(2_pInt:res(1)+1_pInt,2_pInt:res(2)+1_pInt,2_pInt:res(3)+1_pInt,1:3) = centroids

 do k = 0_pInt,res(3)+1_pInt
   do j = 0_pInt,res(2)+1_pInt
     do i = 0_pInt,res(1)+1_pInt
       if (k==0_pInt .or. k==res(3)+1_pInt .or. &                               ! z skin
           j==0_pInt .or. j==res(2)+1_pInt .or. &                               ! y skin
           i==0_pInt .or. i==res(1)+1_pInt      ) then                          ! x skin
         me = (/i,j,k/)                                              ! me on skin
         shift = sign(abs(res+diag-2_pInt*me)/(res+diag),res+diag-2_pInt*me)
         lookup = me-diag+shift*res
   wrappedCentroids(i+1_pInt,j+1_pInt,k+1_pInt,1:3) = &
                                           centroids(lookup(1)+1_pInt,lookup(2)+1_pInt,lookup(3)+1_pInt,1:3) - &
                                           matmul(defgrad_av, shift*geomdim)
       endif
 enddo; enddo; enddo
 do k = 0_pInt,res(3)
   do j = 0_pInt,res(2)
     do i = 0_pInt,res(1)
       do n = 1_pInt,8_pInt
 nodes(i+1_pInt,j+1_pInt,k+1_pInt,1:3) = &
                                nodes(i+1_pInt,j+1_pInt,k+1_pInt,1:3) + wrappedCentroids(i+1_pInt+neighbor(1_pInt,n), &
                                                                                         j+1_pInt+neighbor(2,n), &
                                                                                         k+1_pInt+neighbor(3,n),1:3)
 enddo; enddo; enddo; enddo
 nodes = nodes/8.0_pReal

end subroutine mesh_regular_grid
 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine deformed_linear(res,geomdim,defgrad_av,defgrad,coord_avgCorner)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routine to calculate coordinates in current configuration for given defgrad
! using linear interpolation (blurres out high frequency defomation)
!
 use debug, only: debug_verbosity
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3)   :: geomdim
 real(pReal), intent(in), dimension(3,3) :: defgrad_av
 real(pReal), intent(in), dimension(     res(1),res(2),res(3),3,3) :: defgrad
 ! output variables
 real(pReal), intent(out), dimension(    res(1),res(2),res(3),3)   :: coord_avgCorner
 ! variables with dimension depending on input
 real(pReal),              dimension(8,6,res(1),res(2),res(3),3)   :: coord
 real(pReal),              dimension(  8,res(1),res(2),res(3),3)   :: coord_avgOrder
 ! other variables
 real(pReal), dimension(3) ::  myStep, fones = 1.0_pReal, parameter_coords, negative, positive
 integer(pInt), dimension(3) :: rear, init, ones = 1_pInt, oppo, me
 integer(pInt) i, j, k, s, o 
 integer(pInt), dimension(3,8) :: corner = reshape((/ &
                                              0_pInt, 0_pInt, 0_pInt,&
                                              1_pInt, 0_pInt, 0_pInt,&
                                              1_pInt, 1_pInt, 0_pInt,&
                                              0_pInt, 1_pInt, 0_pInt,&
                                              1_pInt, 1_pInt, 1_pInt,&
                                              0_pInt, 1_pInt, 1_pInt,&
                                              0_pInt, 0_pInt, 1_pInt,&
                                              1_pInt, 0_pInt, 1_pInt &
                                                  /), &
                                               (/3,8/))
 integer(pInt), dimension(3,8) :: step = reshape((/ &
                                            1_pInt, 1_pInt, 1_pInt,&
                                           -1_pInt, 1_pInt, 1_pInt,& 
                                           -1_pInt,-1_pInt, 1_pInt,&
                                            1_pInt,-1_pInt, 1_pInt,&
                                           -1_pInt,-1_pInt,-1_pInt,&
                                            1_pInt,-1_pInt,-1_pInt,&
                                            1_pInt, 1_pInt,-1_pInt,&
                                           -1_pInt, 1_pInt,-1_pInt & 
                                                /), &
                                             (/3,8/))
 integer(pInt), dimension(3,6) :: order = reshape((/ &
                                            1_pInt, 2_pInt, 3_pInt,&
                                            1_pInt, 3_pInt, 2_pInt,&
                                            2_pInt, 1_pInt, 3_pInt,&
                                            2_pInt, 3_pInt, 1_pInt,&
                                            3_pInt, 1_pInt, 2_pInt,&
                                            3_pInt, 2_pInt, 1_pInt &
                                                /), &
                                             (/3,6/))

 if (debug_verbosity > 0_pInt) then
   print*, 'Restore geometry using linear integration'
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 coord_avgOrder = 0.0_pReal
 
 do s = 0_pInt, 7_pInt                               ! corners (from 0 to 7)
   init = corner(:,s+1_pInt)*(res-ones) +ones
   oppo = corner(:,mod((s+4_pInt),8_pInt)+1_pInt)*(res-ones) +ones
   do o=1_pInt,6_pInt                                ! orders (from 1 to 6)
     do k = init(order(3,o)), oppo(order(3,o)), step(order(3,o),s+1_pInt)
       rear(order(2,o)) = init(order(2,o))
       do j = init(order(2,o)), oppo(order(2,o)), step(order(2,o),s+1_pInt)
         rear(order(1,o)) = init(order(1,o))
         do i = init(order(1,o)), oppo(order(1,o)), step(order(1,o),s+1_pInt)
           me(order(1,o)) = i
           me(order(2,o)) = j
           me(order(3,o)) = k
           if ( (me(1)==init(1)).and.(me(2)==init(2)).and. (me(3)==init(3)) ) then
             coord(s+1_pInt,o,me(1),me(2),me(3),1:3) = geomdim * (matmul(defgrad_av,corner(1:3,s+1)) + &
                                                     matmul(defgrad(me(1),me(2),me(3),1:3,1:3),0.5*step(1:3,s+1_pInt)/res))

           else
             myStep = (me-rear)*geomdim/res
             coord(s+1_pInt,o,me(1),me(2),me(3),1:3) = coord(s+1_pInt,o,rear(1),rear(2),rear(3),1:3) + &
                                             0.5*matmul(defgrad(me(1),me(2),me(3),1:3,1:3) + &
                                                        defgrad(rear(1),rear(2),rear(3),1:3,1:3),myStep)
           endif
           rear = me
   enddo; enddo; enddo; enddo
   do i = 1_pInt,6_pInt
     coord_avgOrder(s+1_pInt,1:res(1),1:res(2),1:res(3),1:3) = coord_avgOrder(s+1_pInt,  1:res(1),1:res(2),1:res(3),1:3)&
                                                                   + coord(s+1_pInt,i,1:res(1),1:res(2),1:res(3),1:3)/6.0
   enddo
 enddo

 do k = 0_pInt, res(3)-1_pInt
   do j = 0_pInt, res(2)-1_pInt
     do i = 0_pInt, res(1)-1_pInt
       parameter_coords = (2.0_pReal*(/real(i,pReal)+0.0_pReal,real(j,pReal)+0.0_pReal,real(k,pReal)+0.0_pReal/)&
                                      -real(res,pReal)+fones)/(real(res,pReal)-fones)
       positive = fones + parameter_coords
       negative = fones - parameter_coords
       coord_avgCorner(i+1_pInt,j+1_pInt,k+1_pInt,1:3)&
                       =(coord_avgOrder(1,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *negative(1)*negative(2)*negative(3)&
                       + coord_avgOrder(2,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *positive(1)*negative(2)*negative(3)&
                       + coord_avgOrder(3,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *positive(1)*positive(2)*negative(3)&
                       + coord_avgOrder(4,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *negative(1)*positive(2)*negative(3)&
                       + coord_avgOrder(5,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *positive(1)*positive(2)*positive(3)&
                       + coord_avgOrder(6,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *negative(1)*positive(2)*positive(3)&
                       + coord_avgOrder(7,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *negative(1)*negative(2)*positive(3)&
                       + coord_avgOrder(8,i+1_pInt,j+1_pInt,k+1_pInt,1:3) *positive(1)*negative(2)*positive(3))*0.125
 enddo; enddo; enddo

end subroutine deformed_linear


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine deformed_fft(res,geomdim,defgrad_av,scaling,defgrad,coords)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routine to calculate coordinates in current configuration for given defgrad
! using integration in Fourier space (more accurate than deformed(...))
!
 use numerics, only: fftw_timelimit, fftw_planner_flag  
 use debug, only: debug_verbosity
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3)   :: geomdim
 real(pReal), intent(in), dimension(3,3) :: defgrad_av
 real(pReal), intent(in)                 :: scaling
 real(pReal), intent(in),  dimension(res(1),              res(2),res(3),3,3) :: defgrad
 ! output variables
 real(pReal), intent(out), dimension(res(1),              res(2),res(3),3)   :: coords
! allocatable arrays for fftw c routines
 type(C_PTR) :: fftw_forth, fftw_back
 type(C_PTR) :: coords_fftw, defgrad_fftw
 real(pReal),    dimension(:,:,:,:,:), pointer :: defgrad_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: defgrad_complex
 real(pReal),    dimension(:,:,:,:),   pointer :: coords_real
 complex(pReal), dimension(:,:,:,:),   pointer :: coords_complex
 ! other variables
 integer(pInt) :: i, j, k, res1_red
 integer(pInt), dimension(3) :: k_s
 complex(pReal), parameter   :: integration_factor = cmplx(0.0_pReal,1.0_pReal)*pi*2.0_pReal
 real(pReal), dimension(3)   :: step, offset_coords
 
 if (debug_verbosity > 0_pInt) then
   print*, 'Restore geometry using FFT-based integration'
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 res1_red = res(1)/2_pInt + 1_pInt                                                                         ! size of complex array in first dimension (c2r, r2c)
 step = geomdim/real(res, pReal)

 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=102)
 call fftw_set_timelimit(fftw_timelimit)
 defgrad_fftw =         fftw_alloc_complex(int(res1_red     *res(2)*res(3)*9_pInt,C_SIZE_T)) !C_SIZE_T is of type integer(8)
 call c_f_pointer(defgrad_fftw, defgrad_real,   [res(1)+2_pInt,res(2),res(3),3,3])
 call c_f_pointer(defgrad_fftw, defgrad_complex,[res1_red     ,res(2),res(3),3,3])
 coords_fftw =          fftw_alloc_complex(int(res1_red     *res(2)*res(3)*3_pInt,C_SIZE_T))        !C_SIZE_T is of type integer(8)
 call c_f_pointer(coords_fftw, coords_real,     [res(1)+2_pInt,res(2),res(3),3])
 call c_f_pointer(coords_fftw, coords_complex,  [res1_red     ,res(2),res(3),3])

 fftw_forth = fftw_plan_many_dft_r2c(3,(/res(3),res(2) ,res(1)/),9_pInt,&                      ! dimensions , length in each dimension in reversed order
                          defgrad_real,(/res(3),res(2) ,res(1)+2_pInt/),&                               ! input data , physical length in each dimension in reversed order
                                     1,  res(3)*res(2)*(res(1)+2_pInt),&                                ! striding   , product of physical lenght in the 3 dimensions
                       defgrad_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 fftw_back  = fftw_plan_many_dft_c2r(3,(/res(3),res(2) ,res(1)/),3_pInt,&
                        coords_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,&
                           coords_real,(/res(3),res(2) ,res(1)+2_pInt/),&
                                     1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)


 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   defgrad_real(i,j,k,1:3,1:3) = defgrad(i,j,k,1:3,1:3)                                        ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo

 call fftw_execute_dft_r2c(fftw_forth, defgrad_real, defgrad_complex)

 coords_complex = 0.0
 do k = 1_pInt, res(3)
   k_s(3) = k-1_pInt
   if(k > res(3)/2_pInt+1_pInt) k_s(3) = k_s(3)-res(3)
   do j = 1_pInt, res(2)
     k_s(2) = j-1_pInt
     if(j > res(2)/2_pInt+1_pInt) k_s(2) = k_s(2)-res(2)
     do i = 1_pInt, res1_red
       k_s(1) = i-1_pInt
       if(i/=1_pInt) coords_complex(i,j,k,1:3) = coords_complex(i,j,k,1:3)&
                                + defgrad_complex(i,j,k,1:3,1)*geomdim(1)/(real(k_s(1),pReal)*integration_factor)
       if(j/=1_pInt) coords_complex(i,j,k,1:3) = coords_complex(i,j,k,1:3)&
                                + defgrad_complex(i,j,k,1:3,2)*geomdim(2)/(real(k_s(2),pReal)*integration_factor)
       if(k/=1_pInt) coords_complex(i,j,k,1:3) = coords_complex(i,j,k,1:3)&
                                + defgrad_complex(i,j,k,1:3,3)*geomdim(3)/(real(k_s(3),pReal)*integration_factor)
 enddo; enddo; enddo
 
 call fftw_execute_dft_c2r(fftw_back,coords_complex,coords_real)
 coords_real = coords_real/real(res(1)*res(2)*res(3))

 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   coords(i,j,k,1:3) = coords_real(i,j,k,1:3)                                        ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo

 offset_coords = matmul(defgrad(1,1,1,1:3,1:3),step/2.0_pReal) - scaling*coords(1,1,1,1:3)
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     coords(i,j,k,1:3) =  scaling*coords(i,j,k,1:3) + offset_coords + matmul(defgrad_av,&
                                                    (/step(1)*real(i-1_pInt,pReal),&
                                                      step(2)*real(j-1_pInt,pReal),&
                                                      step(3)*real(k-1_pInt,pReal)/))

 enddo; enddo; enddo

end subroutine deformed_fft


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine curl_fft(res,geomdim,vec_tens,field,curl)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates curl field using differentation in Fourier space
! use vec_tens to decide if tensor (3) or vector (1)

 use numerics, only: fftw_timelimit, fftw_planner_flag  
 use debug, only: debug_verbosity
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3) :: geomdim
 integer(pInt), intent(in) :: vec_tens
 real(pReal), intent(in),  dimension(res(1),              res(2),res(3),vec_tens,3) :: field
 ! output variables
 real(pReal), intent(out), dimension(res(1),              res(2),res(3),vec_tens,3) :: curl
 ! variables with dimension depending on input
 real(pReal),              dimension(res(1)/2_pInt+1_pInt,res(2),res(3),3) :: xi
 ! allocatable arrays for fftw c routines
 type(C_PTR) :: fftw_forth, fftw_back
 type(C_PTR) :: field_fftw, curl_fftw
 real(pReal),    dimension(:,:,:,:,:), pointer :: field_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: field_complex
 real(pReal),    dimension(:,:,:,:,:), pointer :: curl_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: curl_complex
 ! other variables
 integer(pInt) i, j, k, l, res1_red
 integer(pInt), dimension(3) :: k_s,cutting_freq
 real(pReal) :: wgt
 complex(pReal), parameter :: differentation_factor = cmplx(0.0_pReal,1.0_pReal)*2.0_pReal*pi       ! cmplx(0.0_pReal, 2.0_pReal*pi) gets huge rounding error (casting to single prec?)

 if (debug_verbosity > 0_pInt) then
   print*, 'Calculating curl of vector/tensor field'
   print '(a,e12.5,e12.5,e12.5)', ' Dimension: ', geomdim
   print '(a,i5,i5,i5)',          ' Resolution:', res
 endif

 wgt = 1.0_pReal/real(res(1)*res(2)*res(3),pReal)
 res1_red = res(1)/2_pInt + 1_pInt                                                                         ! size of complex array in first dimension (c2r, r2c)

 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=102)
 call fftw_set_timelimit(fftw_timelimit)
 field_fftw =         fftw_alloc_complex(int(res1_red     *res(2)*res(3)*vec_tens*3_pInt,C_SIZE_T)) !C_SIZE_T is of type integer(8)
 call c_f_pointer(field_fftw, field_real,   [res(1)+2_pInt,res(2),res(3),vec_tens,3])
 call c_f_pointer(field_fftw, field_complex,[res1_red     ,res(2),res(3),vec_tens,3])
 curl_fftw =          fftw_alloc_complex(int(res1_red     *res(2)*res(3)*vec_tens*3_pInt,C_SIZE_T))        !C_SIZE_T is of type integer(8)
 call c_f_pointer(curl_fftw, curl_real,     [res(1)+2_pInt,res(2),res(3),vec_tens,3])
 call c_f_pointer(curl_fftw, curl_complex,  [res1_red     ,res(2),res(3),vec_tens,3])

 fftw_forth = fftw_plan_many_dft_r2c(3,(/res(3),res(2) ,res(1)/),vec_tens*3_pInt,&                      ! dimensions , length in each dimension in reversed order
                            field_real,(/res(3),res(2) ,res(1)+2_pInt/),&                               ! input data , physical length in each dimension in reversed order
                                     1,  res(3)*res(2)*(res(1)+2_pInt),&                                ! striding   , product of physical lenght in the 3 dimensions
                         field_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 fftw_back  = fftw_plan_many_dft_c2r(3,(/res(3),res(2) ,res(1)/),vec_tens*3_pInt,&
                          curl_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,&
                             curl_real,(/res(3),res(2) ,res(1)+2_pInt/),&
                                     1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)


 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   field_real(i,j,k,1:vec_tens,1:3) = field(i,j,k,1:vec_tens,1:3)                                        ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo

 call fftw_execute_dft_r2c(fftw_forth, field_real, field_complex)

 do k = 1_pInt, res(3)                              ! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1, res1_red
           k_s(1) = i - 1_pInt
           xi(i,j,k,1:3) = real(k_s, pReal)/geomdim
 enddo; enddo; enddo
!remove the given highest frequencies
 cutting_freq = (/0_pInt,0_pInt,0_pInt/)                                                      ! for 0,0,0, just the highest freq. is removed
 xi( res(1)/2_pInt+1_pInt-cutting_freq(1):res1_red , 1:res(2) , 1:res(3) , 1 ) = 0.0_pReal
 xi(1:res1_red , res(2)/2_pInt+1_pInt-cutting_freq(2):res(2)/2_pInt+1_pInt+cutting_freq(2) , 1:res(3) , 2) = 0.0_pReal
 xi(1:res1_red , 1:res(2) , res(3)/2_pInt+1_pInt-cutting_freq(3):res(3)/2_pInt+1_pInt+cutting_freq(3) , 3) = 0.0_pReal

 do k = 1, res(3); do j = 1, res(2);do i = 1, res1_red
   do l = 1, vec_tens
     curl_complex(i,j,k,l,1) = ( field_complex(i,j,k,l,3)*xi(i,j,k,2)&
                                -field_complex(i,j,k,l,2)*xi(i,j,k,3) )*differentation_factor
     curl_complex(i,j,k,l,2) = (-field_complex(i,j,k,l,3)*xi(i,j,k,1)&
                                +field_complex(i,j,k,l,1)*xi(i,j,k,3) )*differentation_factor
     curl_complex(i,j,k,l,3) = ( field_complex(i,j,k,l,2)*xi(i,j,k,1)&
                                -field_complex(i,j,k,l,1)*xi(i,j,k,2) )*differentation_factor
   enddo
 enddo; enddo; enddo

 call fftw_execute_dft_c2r(fftw_back, curl_complex, curl_real)
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   curl(i,j,k,1:vec_tens,1:3) = curl_real(i,j,k,1:vec_tens,1:3)                                        ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo

 curl = curl * wgt
 
end subroutine curl_fft

 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine divergence_fft(res,geomdim,vec_tens,field,divergence)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates divergence field using integration in Fourier space
! use vec_tens to decide if tensor (3) or vector (1)

 use numerics, only: fftw_timelimit, fftw_planner_flag  
 use debug, only: debug_verbosity
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(3)   :: geomdim
 integer(pInt), intent(in) :: vec_tens
 real(pReal), intent(in),  dimension(res(1),    res(2),res(3),vec_tens,3) :: field
 ! output variables
 real(pReal), intent(out), dimension(res(1),    res(2),res(3),vec_tens) :: divergence
 ! variables with dimension depending on input
 real(pReal),    dimension(res(1)/2_pInt+1_pInt,res(2),res(3),3) :: xi
 ! allocatable arrays for fftw c routines
 type(C_PTR) :: fftw_forth, fftw_back
 type(C_PTR) :: field_fftw, divergence_fftw
 real(pReal),    dimension(:,:,:,:,:), pointer :: field_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: field_complex
 real(pReal),    dimension(:,:,:,:),   pointer :: divergence_real
 complex(pReal), dimension(:,:,:,:),   pointer :: divergence_complex
 ! other variables
 integer(pInt) :: i, j, k, l, res1_red
 real(pReal) :: wgt
 integer(pInt), dimension(3) :: k_s,cutting_freq
 complex(pReal), parameter :: differentation_factor = cmplx(0.0_pReal,1.0_pReal)*2.0_pReal*pi       ! cmplx(0.0_pReal, 2.0_pReal*pi) gets huge rounding error (casting to single prec?)

 if (debug_verbosity > 0_pInt) then 
   print '(a)', 'Calculating divergence of tensor/vector field using FFT'  
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)',          ' Resolution:', res
 endif

 res1_red = res(1)/2_pInt + 1_pInt                                                                         ! size of complex array in first dimension (c2r, r2c)
 wgt = 1.0_pReal/real(res(1)*res(2)*res(3),pReal)

 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=102)
 call fftw_set_timelimit(fftw_timelimit)
 field_fftw = fftw_alloc_complex(int(res1_red*res(2)*res(3)*vec_tens*3_pInt,C_SIZE_T))              !C_SIZE_T is of type integer(8)
 call c_f_pointer(field_fftw, field_real,             [res(1)+2_pInt,res(2),res(3),vec_tens,3])
 call c_f_pointer(field_fftw, field_complex,          [res1_red     ,res(2),res(3),vec_tens,3])
 divergence_fftw = fftw_alloc_complex(int(res1_red*res(2)*res(3)*vec_tens,C_SIZE_T))
 call c_f_pointer(divergence_fftw, divergence_real,   [res(1)+2_pInt,res(2),res(3),vec_tens])
 call c_f_pointer(divergence_fftw, divergence_complex,[res1_red     ,res(2),res(3),vec_tens])

 fftw_forth = fftw_plan_many_dft_r2c(3,(/res(3),res(2) ,res(1)/),vec_tens*3_pInt,&                  ! dimensions , length in each dimension in reversed order
                            field_real,(/res(3),res(2) ,res(1)+2_pInt/),&                           ! input data , physical length in each dimension in reversed order
                                     1,  res(3)*res(2)*(res(1)+2_pInt),&                            ! striding   , product of physical lenght in the 3 dimensions
                         field_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 fftw_back  = fftw_plan_many_dft_c2r(3,(/res(3),res(2) ,res(1)/),vec_tens,&
                    divergence_complex,(/res(3),res(2) ,res1_red/),&
                                     1,  res(3)*res(2)* res1_red,&
                       divergence_real,(/res(3),res(2) ,res(1)+2_pInt/),&
                                     1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 field_real = 0.0_pReal                                                                             ! padding 
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   field_real(i,j,k,1:vec_tens,1:3) = field(i,j,k,1:vec_tens,1:3)                                   ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo
 
 call fftw_execute_dft_r2c(fftw_forth, field_real, field_complex)
 do k = 1_pInt, res(3)                                                                              ! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1, res1_red
           k_s(1) = i - 1_pInt
           xi(i,j,k,1:3) = real(k_s, pReal)/geomdim
 enddo; enddo; enddo
 
!remove the given highest frequencies
 cutting_freq = 0_pInt                                                                              ! for 0,0,0, just the highest freq. is removed
 xi( res(1)/2_pInt+1_pInt-cutting_freq(1):res1_red , 1:res(2) , 1:res(3) , 1 ) = 0.0_pReal
 xi(1:res1_red , res(2)/2_pInt+1_pInt-cutting_freq(2):res(2)/2_pInt+1_pInt+cutting_freq(2) , 1:res(3) , 2) = 0.0_pReal
 xi(1:res1_red , 1:res(2) , res(3)/2_pInt+1_pInt-cutting_freq(3):res(3)/2_pInt+1_pInt+cutting_freq(3) , 3) = 0.0_pReal
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
   do l = 1_pInt, vec_tens
     divergence_complex(i,j,k,l) = sum(field_complex(i,j,k,l,1:3)*xi(i,j,k,1:3))&
                                   *differentation_factor
   enddo
 enddo; enddo; enddo
 call fftw_execute_dft_c2r(fftw_back, divergence_complex, divergence_real)

 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   divergence(i,j,k,1:vec_tens) = divergence_real(i,j,k,1:vec_tens)                                        ! ensure that data is aligned properly (fftw_alloc)
 enddo; enddo; enddo

 divergence = divergence * wgt
 
 end subroutine divergence_fft

          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 subroutine divergence_fdm(res,geomdim,vec_tens,order,field,divergence)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates divergence field using FDM with variable accuracy
! use vec_tes to decide if tensor (3) or vector (1)
 
 use debug, only: debug_verbosity
 implicit none
 integer(pInt), intent(in), dimension(3) :: res
 integer(pInt), intent(in)               :: vec_tens
 integer(pInt), intent(inout)            :: order
 real(pReal), intent(in), dimension(3)   :: geomdim
 real(pReal), intent(in), dimension(res(1),res(2),res(3),vec_tens,3) :: field
 ! output variables
 real(pReal), intent(out), dimension(res(1),res(2),res(3),vec_tens) :: divergence
 ! other variables
 integer(pInt), dimension(6,3) :: coordinates
 integer(pInt) i, j, k, m, l
 real(pReal), dimension(4,4), parameter :: FDcoefficient = reshape((/ & 
                    1.0_pReal/2.0_pReal,       0.0_pReal,           0.0_pReal,             0.0_pReal,& !from http://en.wikipedia.org/wiki/Finite_difference_coefficients
                    2.0_pReal/3.0_pReal,-1.0_pReal/12.0_pReal,      0.0_pReal,             0.0_pReal,&
                    3.0_pReal/4.0_pReal,-3.0_pReal/20.0_pReal,1.0_pReal/ 60.0_pReal,       0.0_pReal,&
                    4.0_pReal/5.0_pReal,-1.0_pReal/ 5.0_pReal,4.0_pReal/105.0_pReal,-1.0_pReal/280.0_pReal/),&
                               (/4,4/))
                               
 if (debug_verbosity > 0_pInt) then
   print*, 'Calculating divergence of tensor/vector field using FDM'
   print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
   print '(a,/,i5,i5,i5)', ' Resolution:', res
 endif

 divergence = 0.0_pReal
 order = order + 1_pInt
 do k = 0_pInt, res(3)-1_pInt; do j = 0_pInt, res(2)-1_pInt; do i = 0_pInt, res(1)-1_pInt
   do m = 1_pInt, order 
     coordinates(1,1:3) = mesh_location(mesh_index((/i+m,j,k/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     coordinates(2,1:3) = mesh_location(mesh_index((/i-m,j,k/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     coordinates(3,1:3) = mesh_location(mesh_index((/i,j+m,k/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     coordinates(4,1:3) = mesh_location(mesh_index((/i,j-m,k/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     coordinates(5,1:3) = mesh_location(mesh_index((/i,j,k+m/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     coordinates(6,1:3) = mesh_location(mesh_index((/i,j,k-m/),(/res(1),res(2),res(3)/)),(/res(1),res(2),res(3)/))&
                                                                                          + (/1_pInt,1_pInt,1_pInt/)
     do l = 1_pInt, vec_tens
       divergence(i+1_pInt,j+1_pInt,k+1_pInt,l) = divergence(i+1_pInt,j+1_pInt,k+1_pInt,l) + FDcoefficient(m,order) * &
                ((field(coordinates(1,1),coordinates(1,2),coordinates(1,3),l,1)- &
                  field(coordinates(2,1),coordinates(2,2),coordinates(2,3),l,1))*real(res(1),pReal)/geomdim(1) +&
                 (field(coordinates(3,1),coordinates(3,2),coordinates(3,3),l,2)- &
                  field(coordinates(4,1),coordinates(4,2),coordinates(4,3),l,2))*real(res(2),pReal)/geomdim(2) +&
                 (field(coordinates(5,1),coordinates(5,2),coordinates(5,3),l,3)- &
                  field(coordinates(6,1),coordinates(6,2),coordinates(6,3),l,3))*real(res(3),pReal)/geomdim(3))  
     enddo
   enddo
 enddo; enddo; enddo

 end subroutine divergence_fdm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 subroutine tensor_avg(res,tensor,avg)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate average of tensor field
!
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in), dimension(res(1),res(2),res(3),3,3) ::tensor
 ! output variables
 real(pReal), intent(out), dimension(3,3) :: avg
 ! other variables
 real(pReal) wgt
 integer(pInt) m,n
 
 wgt = 1.0_pReal/real(res(1)*res(2)*res(3), pReal)

 do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt
    avg(m,n) = sum(tensor(1:res(1),1:res(2),1:res(3),m,n)) * wgt   
 enddo; enddo

 end subroutine tensor_avg
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine logstrain_spat(res,defgrad,logstrain_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate logarithmic strain in spatial configuration for given defgrad field
!
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in),  dimension(res(1),res(2),res(3),3,3) :: defgrad
 ! output variables
 real(pReal), intent(out), dimension(res(1),res(2),res(3),3,3) :: logstrain_field
 ! other variables
 real(pReal), dimension(3,3) ::  temp33_Real, temp33_Real2
 real(pReal), dimension(3,3,3) :: eigenvectorbasis
 real(pReal), dimension(3) ::  eigenvalue
 integer(pInt) :: i, j, k
 logical :: errmatinv
 
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   call math_pDecomposition(defgrad(i,j,k,1:3,1:3),temp33_Real2,temp33_Real,errmatinv)  !store R in temp33_Real
   temp33_Real2 = math_inv33(temp33_Real)
   temp33_Real = math_mul33x33(defgrad(i,j,k,1:3,1:3),temp33_Real2)       ! v = F o inv(R), store in temp33_Real2
   call math_spectral1(temp33_Real,eigenvalue(1),              eigenvalue(2),              eigenvalue(3),&
                                   eigenvectorbasis(1,1:3,1:3),eigenvectorbasis(2,1:3,1:3),eigenvectorbasis(3,1:3,1:3))
   eigenvalue = log(sqrt(eigenvalue))
   logstrain_field(i,j,k,1:3,1:3) = eigenvalue(1)*eigenvectorbasis(1,1:3,1:3)+&
                                    eigenvalue(2)*eigenvectorbasis(2,1:3,1:3)+&
                                    eigenvalue(3)*eigenvectorbasis(3,1:3,1:3)
 enddo; enddo; enddo

 end subroutine logstrain_spat
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine logstrain_mat(res,defgrad,logstrain_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate logarithmic strain in material configuration for given defgrad field
!
 implicit none
  ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in),  dimension(res(1),res(2),res(3),3,3) :: defgrad
 ! output variables
 real(pReal), intent(out), dimension(res(1),res(2),res(3),3,3) :: logstrain_field
 ! other variables
 real(pReal), dimension(3,3) ::  temp33_Real, temp33_Real2
 real(pReal), dimension(3,3,3) :: eigenvectorbasis
 real(pReal), dimension(3) ::  eigenvalue
 integer(pInt) :: i, j, k
 logical :: errmatinv
 
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   call math_pDecomposition(defgrad(i,j,k,1:3,1:3),temp33_Real,temp33_Real2,errmatinv)  !store U in temp33_Real
   call math_spectral1(temp33_Real,eigenvalue(1),              eigenvalue(2),              eigenvalue(3),&
                                   eigenvectorbasis(1,1:3,1:3),eigenvectorbasis(2,1:3,1:3),eigenvectorbasis(3,1:3,1:3))
   eigenvalue = log(sqrt(eigenvalue))
   logstrain_field(i,j,k,1:3,1:3) = eigenvalue(1)*eigenvectorbasis(1,1:3,1:3)+&
                                    eigenvalue(2)*eigenvectorbasis(2,1:3,1:3)+&
                                    eigenvalue(3)*eigenvectorbasis(3,1:3,1:3)
 enddo; enddo; enddo

 end subroutine logstrain_mat
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine calculate_cauchy(res,defgrad,p_stress,c_stress)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate cauchy stress for given PK1 stress and defgrad field
!
 implicit none
 ! input variables
 integer(pInt), intent(in), dimension(3) :: res
 real(pReal), intent(in),  dimension(res(1),res(2),res(3),3,3) :: defgrad
 real(pReal), intent(in),  dimension(res(1),res(2),res(3),3,3) :: p_stress
 ! output variables
 real(pReal), intent(out),  dimension(res(1),res(2),res(3),3,3) :: c_stress
 ! other variables
 real(pReal) :: jacobi
 integer(pInt) :: i, j, k

 c_stress = 0.0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   jacobi = math_det33(defgrad(i,j,k,1:3,1:3))
   c_stress(i,j,k,1:3,1:3) = matmul(p_stress(i,j,k,1:3,1:3),transpose(defgrad(i,j,k,1:3,1:3)))/jacobi
 enddo; enddo; enddo

end subroutine calculate_cauchy

END MODULE math

!#############################################################################################################################
! BEGIN KDTREE2
!#############################################################################################################################
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
!#######################################################
! modifications: changed precision according to prec.f90
! k.komerla, m.diehl
!#######################################################

module kdtree2_priority_queue_module
  use prec
  !
  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload', 
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
  type kdtree2_result
      ! a pair of distances, indexes
      real(pReal)    :: dis        !=0.0
      integer(pInt) :: idx         !=-1   Initializers cause some bugs in compilers.
  end type kdtree2_result
  !
  ! A heap-based priority queue lets one efficiently implement the following
  ! operations, each in log(N) time, as opposed to linear time.
  !
  ! 1)  add a datum (push a datum onto the queue, increasing its length) 
  ! 2)  return the priority value of the maximum priority element 
  ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
  !     the size of the queue. 
  ! 4)  replace the datum with the maximum priority with a supplied datum
  !     (of either higher or lower priority), maintaining the size of the
  !     queue. 
  !
  !
  ! In the k-d tree case, the 'priority' is the square distance of a point in
  ! the data set to a reference point.   The goal is to keep the smallest M
  ! distances to a reference point.  The tree algorithm searches terminal
  ! nodes to decide whether to add points under consideration.
  !
  ! A priority queue is useful here because it lets one quickly return the
  ! largest distance currently existing in the list.  If a new candidate
  ! distance is smaller than this, then the new candidate ought to replace
  ! the old candidate.  In priority queue terms, this means removing the
  ! highest priority element, and inserting the new one.
  !
  ! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
  ! to Algorithms_, 1990, with further optimization by the author.
  !
  ! Originally informed by a C implementation by Sriranga Veeraraghavan.
  !
  ! This module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type pq
      !
      ! The priority queue consists of elements
      ! priority(1:heap_size), with associated payload(:).
      !
      ! There are heap_size active elements. 
      ! Assumes the allocation is always sufficient.  Will NOT increase it
      ! to match.
      integer(pInt) :: heap_size = 0
      type(kdtree2_result), pointer :: elems(:) 
  end type pq

  public :: kdtree2_result

  public :: pq
  public :: pq_create
  public :: pq_delete, pq_insert
  public :: pq_extract_max, pq_max, pq_replace_max, pq_maxpri
  private

contains


  function pq_create(results_in) result(res)
    !
    ! Create a priority queue from ALREADY allocated
    ! array pointers for storage.  NOTE! It will NOT
    ! add any alements to the heap, i.e. any existing
    ! data in the input arrays will NOT be used and may
    ! be overwritten.
    ! 
    ! usage:
    !    real(pReal), pointer :: x(:)
    !    integer(pInt), pointer :: k(:)
    !    allocate(x(1000),k(1000))
    !    pq => pq_create(x,k)
    !
    type(kdtree2_result), target:: results_in(:) 
    type(pq) :: res
    !
    !
    integer(pInt) :: nalloc

    nalloc = size(results_in,1)
    if (nalloc .lt. 1) then
       write (*,*) 'PQ_CREATE: error, input arrays must be allocated.'
    end if
    res%elems => results_in
    res%heap_size = 0
    return
  end function pq_create

  !
  ! operations for getting parents and left + right children
  ! of elements in a binary heap.
  !

!
! These are written inline for speed.
!    
!  integer(pInt) function parent(i)
!    integer(pInt), intent(in) :: i
!    parent = (i/2)
!    return
!  end function parent

!  integer(pInt) function left(i)
!    integer(pInt), intent(in) ::i
!    left = (2*i)
!    return
!  end function left

!  integer(pInt) function right(i)
!    integer(pInt), intent(in) :: i
!    right = (2*i)+1
!    return
!  end function right

!  logical function compare_priority(p1,p2)
!    real(pReal), intent(in) :: p1, p2
!
!    compare_priority = (p1 .gt. p2)
!    return
!  end function compare_priority

  subroutine heapify(a,i_in)
    !
    ! take a heap rooted at 'i' and force it to be in the
    ! heap canonical form.   This is performance critical 
    ! and has been tweaked a little to reflect this.
    !
    type(pq),pointer   :: a
    integer(pInt), intent(in) :: i_in
    !
    integer(pInt) :: i, l, r, largest

    real(pReal)    :: pri_i, pri_l, pri_r, pri_largest


    type(kdtree2_result) :: temp

    i = i_in

bigloop:  do
       l = 2*i ! left(i)
       r = l+1 ! right(i)
       ! 
       ! set 'largest' to the index of either i, l, r
       ! depending on whose priority is largest.
       !
       ! note that l or r can be larger than the heap size
       ! in which case they do not count.


       ! does left child have higher priority? 
       if (l .gt. a%heap_size) then
          ! we know that i is the largest as both l and r are invalid.
          exit 
       else
          pri_i = a%elems(i)%dis
          pri_l = a%elems(l)%dis 
          if (pri_l .gt. pri_i) then
             largest = l
             pri_largest = pri_l
          else
             largest = i
             pri_largest = pri_i
          endif

          !
          ! between i and l we have a winner
          ! now choose between that and r.
          !
          if (r .le. a%heap_size) then
             pri_r = a%elems(r)%dis
             if (pri_r .gt. pri_largest) then
                largest = r
             endif
          endif
       endif

       if (largest .ne. i) then
          ! swap data in nodes largest and i, then heapify

          temp = a%elems(i)
          a%elems(i) = a%elems(largest)
          a%elems(largest) = temp 
          ! 
          ! Canonical heapify() algorithm has tail-ecursive call: 
          !
          !        call heapify(a,largest)   
          ! we will simulate with cycle
          !
          i = largest
          cycle bigloop ! continue the loop 
       else
          return   ! break from the loop
       end if
    enddo bigloop
    return
  end subroutine heapify

  subroutine pq_max(a,e) 
    !
    ! return the priority and its payload of the maximum priority element
    ! on the queue, which should be the first one, if it is 
    ! in heapified form.
    !
    type(pq),pointer :: a
    type(kdtree2_result),intent(out)  :: e

    if (a%heap_size .gt. 0) then
       e = a%elems(1) 
    else
       write (*,*) 'PQ_MAX: ERROR, heap_size < 1'
       stop
    endif
    return
  end subroutine pq_max
  
  real(pReal) function pq_maxpri(a)
    type(pq), pointer :: a

    if (a%heap_size .gt. 0) then
       pq_maxpri = a%elems(1)%dis
    else
       write (*,*) 'PQ_MAX_PRI: ERROR, heapsize < 1'
       stop
    endif
    return
  end function pq_maxpri

  subroutine pq_extract_max(a,e)
    !
    ! return the priority and payload of maximum priority
    ! element, and remove it from the queue.
    ! (equivalent to 'pop()' on a stack)
    !
    type(pq),pointer :: a
    type(kdtree2_result), intent(out) :: e
    
    if (a%heap_size .ge. 1) then
       !
       ! return max as first element
       !
       e = a%elems(1) 
       
       !
       ! move last element to first
       !
       a%elems(1) = a%elems(a%heap_size) 
       a%heap_size = a%heap_size-1
       call heapify(a,1)
       return
    else
       write (*,*) 'PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ'
       stop
    end if
    
  end subroutine pq_extract_max


  real(pReal) function pq_insert(a,dis,idx) 
    !
    ! Insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    type(pq),pointer  :: a
    real(pReal), intent(in) :: dis
    integer(pInt), intent(in) :: idx
    !    type(kdtree2_result), intent(in) :: e
    !
    integer(pInt) :: i, isparent
    real(pReal)    :: parentdis
    !

    !    if (a%heap_size .ge. a%max_elems) then
    !       write (*,*) 'PQ_INSERT: error, attempt made to insert element on full PQ'
    !       stop
    !    else
    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while (i .gt. 1)
       isparent = int(i/2)
       parentdis = a%elems(isparent)%dis
       if (dis .gt. parentdis) then
          ! move what was in i's parent into i.
          a%elems(i)%dis = parentdis
          a%elems(i)%idx = a%elems(isparent)%idx
          i = isparent
       else
          exit
       endif
    end do

    ! insert the element at the determined position
    a%elems(i)%dis = dis
    a%elems(i)%idx = idx

    pq_insert = a%elems(1)%dis 
    return
    !    end if

  end function pq_insert

  subroutine pq_adjust_heap(a,i)
    type(pq),pointer  :: a
    integer(pInt), intent(in) :: i
    !
    ! nominally arguments (a,i), but specialize for a=1
    !
    ! This routine assumes that the trees with roots 2 and 3 are already heaps, i.e.
    ! the children of '1' are heaps.  When the procedure is completed, the
    ! tree rooted at 1 is a heap.
    real(pReal) :: prichild
    integer(pInt) :: parent, child, N

    type(kdtree2_result) :: e

    e = a%elems(i) 

    parent = i
    child = 2*i
    N = a%heap_size
    
    do while (child .le. N)
       if (child .lt. N) then
          if (a%elems(child)%dis .lt. a%elems(child+1)%dis) then
             child = child+1
          endif
       endif
       prichild = a%elems(child)%dis
       if (e%dis .ge. prichild) then
          exit 
       else
          ! move child into parent.
          a%elems(parent) = a%elems(child) 
          parent = child
          child = 2*parent
       end if
    end do
    a%elems(parent) = e
    return
  end subroutine pq_adjust_heap
    

  real(pReal) function pq_replace_max(a,dis,idx) 
    !
    ! Replace the extant maximum priority element
    ! in the PQ with (dis,idx).  Return
    ! the new maximum priority, which may be larger
    ! or smaller than the old one.
    !
    type(pq),pointer         :: a
    real(pReal), intent(in) :: dis
    integer(pInt), intent(in) :: idx
!    type(kdtree2_result), intent(in) :: e
    ! not tested as well!

    integer(pInt) :: parent, child, N
    real(pReal)    :: prichild, prichildp1

    type(kdtree2_result) :: etmp
    
    if (.true.) then
       N=a%heap_size
       if (N .ge. 1) then
          parent =1
          child=2

          loop: do while (child .le. N)
             prichild = a%elems(child)%dis

             !
             ! posibly child+1 has higher priority, and if
             ! so, get it, and increment child.
             !

             if (child .lt. N) then
                prichildp1 = a%elems(child+1)%dis
                if (prichild .lt. prichildp1) then
                   child = child+1
                   prichild = prichildp1
                endif
             endif

             if (dis .ge. prichild) then
                exit loop  
                ! we have a proper place for our new element, 
                ! bigger than either children's priority.
             else
                ! move child into parent.
                a%elems(parent) = a%elems(child) 
                parent = child
                child = 2*parent
             end if
          end do loop
          a%elems(parent)%dis = dis
          a%elems(parent)%idx = idx
          pq_replace_max = a%elems(1)%dis
       else
          a%elems(1)%dis = dis
          a%elems(1)%idx = idx
          pq_replace_max = dis
       endif
    else
       !
       ! slower version using elementary pop and push operations.
       !
       call pq_extract_max(a,etmp) 
       etmp%dis = dis
       etmp%idx = idx
       pq_replace_max = pq_insert(a,dis,idx)
    endif
    return
  end function pq_replace_max

  subroutine pq_delete(a,i)
    ! 
    ! delete item with index 'i'
    !
    type(pq),pointer :: a
    integer(pInt)           :: i

    if ((i .lt. 1) .or. (i .gt. a%heap_size)) then
       write (*,*) 'PQ_DELETE: error, attempt to remove out of bounds element.'
       stop
    endif

    ! swap the item to be deleted with the last element
    ! and shorten heap by one.
    a%elems(i) = a%elems(a%heap_size) 
    a%heap_size = a%heap_size - 1

    call heapify(a,i)

  end subroutine pq_delete

end module kdtree2_priority_queue_module


module kdtree2_module
  use prec
  use kdtree2_priority_queue_module
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported. 
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  public :: pReal
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy
  !---------------------------------------------------------------
  !-------------------SEARCH ROUTINES-----------------------------
  public :: kdtree2_n_nearest,kdtree2_n_nearest_around_point
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window. 
  !
  public :: kdtree2_r_nearest, kdtree2_r_nearest_around_point
  ! Return points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_sort_results
  ! Sort, in order of increasing distance, rseults from above.
  !
  public :: kdtree2_r_count, kdtree2_r_count_around_point 
  ! Count points within a fixed ball of arb vector/extant point 
  !
  public :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force
  ! brute force of kdtree2_[n|r]_nearest
  !----------------------------------------------------------------


  integer(pInt), parameter :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

  type interval
      real(pReal) :: lower,upper
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer(pInt) :: cut_dim
      ! the dimension to cut
      real(pReal) :: cut_val
      ! where to cut the dimension
      real(pReal) :: cut_val_left, cut_val_right  
      ! improved cutoffs knowing the spread in child boxes.
      integer(pInt) :: l, u
      type (tree_node), pointer :: left, right
      type(interval), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u] 


  end type tree_node

  type :: kdtree2
      ! Global information about the tree, one per tree
      integer(pInt) :: dimen=0, n=0
      ! dimensionality and total # of points
      real(pReal), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array 
      ! 
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable 
      !  vectorization on some processors and compilers. 

      integer(pInt), pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical       :: sort = .false.
      ! do we always sort output results?
      logical       :: rearrange = .false. 
      real(pReal), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer :: root => null()
      ! root pointer of the tree
  end type kdtree2


  type :: tree_search_record
      !
      ! One of these is created for each search.
      !
      private
      ! 
      ! Many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer(pInt)           :: dimen   
      integer(pInt)           :: nn, nfound
      real(pReal)      :: ballsize
      integer(pInt)           :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer(pInt)           :: nalloc  ! how much allocated for results(:)?
      logical           :: rearrange  ! are the data rearranged or original? 
      ! did the # of points found overflow the storage provided?
      logical           :: overflow
      real(pReal), pointer :: qv(:)  ! query vector
      type(kdtree2_result), pointer :: results(:) ! results
      type(pq) :: pq
      real(pReal), pointer :: data(:,:)  ! temp pointer to data
      integer(pInt), pointer      :: ind(:)     ! temp pointer to indexes
  end type tree_search_record

  private
  ! everything else is private.

  type(tree_search_record), save, target :: sr   ! A GLOBAL VARIABLE for search

contains

  function kdtree2_create(input_data,dim,sort,rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  If 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from SIZE(input_data,1).
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !                      
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness. 
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (kdtree2), pointer :: mr
    integer(pInt), intent(in), optional      :: dim
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    ! ..
    ! .. Array Arguments ..
    real(pReal), target :: input_data(:,:)
    !
    integer(pInt) :: i
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       write (*,*) 'KD_TREE_TRANS: likely user error.'
       write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
       write (*,*) 'KD_TREE_TRANS: and N=',mr%n
       write (*,*) 'KD_TREE_TRANS: note, that new format is data(1:D,1:N)'
       write (*,*) 'KD_TREE_TRANS: with usually N >> D.   If N =approx= D, then a k-d tree'
       write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
       stop
    end if

    call build_tree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (mr%rearrange) then
       allocate(mr%rearranged_data(mr%dimen,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:, &
           mr%ind(i))
       enddo
    else
       nullify(mr%rearranged_data)
    endif

  end function kdtree2_create

    subroutine build_tree(tp)
      type (kdtree2), pointer :: tp
      ! ..
      integer(pInt) :: j
      type(tree_node), pointer :: dummy => null()
      ! ..
      allocate (tp%ind(tp%n))
      forall (j=1:tp%n)
         tp%ind(j) = j
      end forall
      tp%root => build_tree_for_range(tp,1,tp%n, dummy)
    end subroutine build_tree

    recursive function build_tree_for_range(tp,l,u,parent) result (res)
      ! .. Function Return Cut_value ..
      type (tree_node), pointer :: res
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type (tree_node),pointer           :: parent
      ! ..
      ! .. Scalar Arguments ..
      integer(pInt), intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      integer(pInt) :: i, c, m, dimen
      logical :: recompute
      real(pReal)    :: average

!!$      If (.False.) Then 
!!$         If ((l .Lt. 1) .Or. (l .Gt. tp%n)) Then
!!$            Stop 'illegal L value in build_tree_for_range'
!!$         End If
!!$         If ((u .Lt. 1) .Or. (u .Gt. tp%n)) Then
!!$            Stop 'illegal u value in build_tree_for_range'
!!$         End If
!!$         If (u .Lt. l) Then
!!$            Stop 'U is less than L, thats illegal.'
!!$         End If
!!$      Endif
!!$      
      ! first compute min and max
      dimen = tp%dimen
      allocate (res)
      allocate(res%box(dimen))

      ! First, compute an APPROXIMATE bounding box of all points associated with this node.
      if ( u < l ) then
         ! no points in this box
         nullify(res)
         return
      end if

      if ((u-l)<=bucket_size) then
         !
         ! always compute true bounding box for terminal nodes.
         !
         do i=1,dimen
            call spread_in_coordinate(tp,i,l,u,res%box(i))
         end do
         res%cut_dim = 0
         res%cut_val = 0.0
         res%l = l
         res%u = u
         res%left =>null()
         res%right => null() 
      else
         ! 
         ! modify approximate bounding box.  This will be an
         ! overestimate of the true bounding box, as we are only recomputing 
         ! the bounding box for the dimension that the parent split on.
         !
         ! Going to a true bounding box computation would significantly
         ! increase the time necessary to build the tree, and usually
         ! has only a very small difference.  This box is not used
         ! for searching but only for deciding which coordinate to split on.
         !
         do i=1,dimen
            recompute=.true.
            if (associated(parent)) then
               if (i .ne. parent%cut_dim) then
                  recompute=.false.
               end if
            endif
            if (recompute) then
               call spread_in_coordinate(tp,i,l,u,res%box(i))
            else
               res%box(i) = parent%box(i)
            endif
         end do
         

         c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         !
         ! c is the identity of which coordinate has the greatest spread.
         !
         
         if (.false.) then
            ! select exact median to have fully balanced tree.
            m = (l+u)/2
            call select_on_coordinate(tp%the_data,tp%ind,c,m,l,u)
         else
            !
            ! select point halfway between min and max, as per A. Moore,
            ! who says this helps in some degenerate cases, or 
            ! actual arithmetic average. 
            !
            if (.true.) then
               ! actually compute average
               average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,pReal)
            else
               average = (res%box(c)%upper + res%box(c)%lower)/2.0
            endif
               
            res%cut_val = average
            m = select_on_coordinate_value(tp%the_data,tp%ind,c,average,l,u)
         endif
            
         ! moves indexes around
         res%cut_dim = c
         res%l = l
         res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

         res%left => build_tree_for_range(tp,l,m,res)
         res%right => build_tree_for_range(tp,m+1,u,res)

         if (associated(res%right) .eqv. .false.) then
            res%box = res%left%box
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = res%cut_val_left
         elseif (associated(res%left) .eqv. .false.) then
            res%box = res%right%box
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val = res%cut_val_right
         else
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = (res%cut_val_left + res%cut_val_right)/2


            ! now remake the true bounding box for self.  
            ! Since we are taking unions (in effect) of a tree structure,
            ! this is much faster than doing an exhaustive
            ! search over all points
            res%box%upper = max(res%left%box%upper,res%right%box%upper)
            res%box%lower = min(res%left%box%lower,res%right%box%lower) 
         endif
      end if
    end function build_tree_for_range

    integer(pInt) function select_on_coordinate_value(v,ind,c,alpha,li,ui) &
     result(res)
      ! Move elts of ind around between l and u, so that all points
      ! <= than alpha (in c cooordinate) are first, and then
      ! all points > alpha are second. 

      !
      ! Algorithm (matt kennel). 
      !
      ! Consider the list as having three parts: on the left,
      ! the points known to be <= alpha.  On the right, the points
      ! known to be > alpha, and in the middle, the currently unknown
      ! points.   The algorithm is to scan the unknown points, starting
      ! from the left, and swapping them so that they are added to
      ! the left stack or the right stack, as appropriate.
      ! 
      ! The algorithm finishes when the unknown stack is empty. 
      !
      ! .. Scalar Arguments ..
      integer(pInt), intent (In) :: c, li, ui
      real(pReal), intent(in) :: alpha
      ! ..
      real(pReal) :: v(1:,1:)
      integer(pInt) :: ind(1:)
      integer(pInt) :: tmp  
      ! ..
      integer(pInt) :: lb, rb
      !
      ! The points known to be <= alpha are in
      ! [l,lb-1]
      !
      ! The points known to be > alpha are in
      ! [rb+1,u].  
      !
      ! Therefore we add new points into lb or
      ! rb as appropriate.  When lb=rb
      ! we are done.  We return the location of the last point <= alpha.
      !
      ! 
      lb = li; rb = ui

      do while (lb < rb)
         if ( v(c,ind(lb)) <= alpha ) then
            ! it is good where it is.
            lb = lb+1
         else
            ! swap it with rb.
            tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
            rb = rb-1
         endif
      end do
      
      ! now lb .eq. ub 
      if (v(c,ind(lb)) <= alpha) then
         res = lb
      else
         res = lb-1
      endif
      
    end function select_on_coordinate_value

    subroutine select_on_coordinate(v,ind,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Scalar Arguments ..
      integer(pInt), intent (In) :: c, k, li, ui
      ! ..
      integer(pInt) :: i, l, m, s, t, u
      ! ..
      real(pReal) :: v(:,:)
      integer(pInt) :: ind(:)
      ! ..
      l = li
      u = ui
      do while (l<u)
         t = ind(l)
         m = l
         do i = l + 1, u
            if (v(c,ind(i))<v(c,t)) then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            end if
         end do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         if (m<=k) l = m + 1
         if (m>=k) u = m - 1
      end do
    end subroutine select_on_coordinate

   subroutine spread_in_coordinate(tp,c,l,u,interv) 
      ! the spread in coordinate 'c', between l and u. 
      !
      ! Return lower bound in 'smin', and upper in 'smax', 
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type(interval), intent(out) :: interv
      ! ..
      ! .. Scalar Arguments ..
      integer(pInt), intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(pReal) :: last, lmax, lmin, t, smin,smax
      integer(pInt) :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(pReal), pointer :: v(:,:)
      integer(pInt), pointer :: ind(:)
      ! ..
      v => tp%the_data(1:,1:)
      ind => tp%ind(1:)
      smin = v(c,ind(l))
      smax = smin

      ulocal = u

      do i = l + 2, ulocal, 2
         lmin = v(c,ind(i-1))
         lmax = v(c,ind(i))
         if (lmin>lmax) then
            t = lmin
            lmin = lmax
            lmax = t
         end if
         if (smin>lmin) smin = lmin
         if (smax<lmax) smax = lmax
      end do
      if (i==ulocal+1) then
         last = v(c,ind(ulocal))
         if (smin>last) smin = last
         if (smax<last) smax = last
      end if

      interv%lower = smin
      interv%upper = smax

    end subroutine spread_in_coordinate


  subroutine kdtree2_destroy(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    type (kdtree2), pointer :: tp
    ! ..
    call destroy_node(tp%root)

    deallocate (tp%ind)
    nullify (tp%ind)

    if (tp%rearrange) then
       deallocate(tp%rearranged_data)
       nullify(tp%rearranged_data)
    endif

    deallocate(tp)
    return

  contains
    recursive subroutine destroy_node(np)
      ! .. Structure Arguments ..
      type (tree_node), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      intrinsic ASSOCIATED
      ! ..
      if (associated(np%left)) then
         call destroy_node(np%left)
         nullify (np%left)
      end if
      if (associated(np%right)) then
         call destroy_node(np%right)
         nullify (np%right)
      end if
      if (associated(np%box)) deallocate(np%box)
      deallocate(np)
      return
      
    end subroutine destroy_node

  end subroutine kdtree2_destroy

  subroutine kdtree2_n_nearest(tp,qv,nn,results)
    ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (kdtree2), pointer      :: tp
    real(pReal), target, intent (In)    :: qv(:)
    integer(pInt), intent (In)         :: nn
    type(kdtree2_result), target :: results(:)


    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false. 

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    call validate_query_storage(nn) 
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
!    deallocate(sr%pqp)
    return
  end subroutine kdtree2_n_nearest

  subroutine kdtree2_n_nearest_around_point(tp,idxin,correltime,nn,results)
    ! Find the 'nn' vectors in the tree nearest to point 'idxin',
    ! with correlation window 'correltime', returing results in
    ! results(:), which must be pre-allocated upon entry.
    type (kdtree2), pointer        :: tp
    integer(pInt), intent (In)           :: idxin, correltime, nn
    type(kdtree2_result), target   :: results(:)

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = huge(1.0)       ! the largest real(pReal) number
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%nn = nn
    sr%nfound = 0

    sr%dimen = tp%dimen
    sr%nalloc = nn

    sr%results => results

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif

    call validate_query_storage(nn)
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
    deallocate (sr%qv)
    return
  end subroutine kdtree2_n_nearest_around_point

  subroutine kdtree2_r_nearest(tp,qv,r2,nfound,nalloc,results) 
    ! find the nearest neighbors to point 'idxin', within SQUARED
    ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
    ! size of memory allocated for results(1:nalloc).  Upon
    ! EXIT, nfound is the number actually found within the ball. 
    !
    !  Note that if nfound .gt. nalloc then more neighbors were found
    !  than there were storage to store.  The resulting list is NOT
    !  the smallest ball inside norm r^2 
    !
    ! Results are NOT sorted unless tree was created with sort option.
    type (kdtree2), pointer      :: tp
    real(pReal), target, intent (In)    :: qv(:)
    real(pReal), intent(in)             :: r2
    integer(pInt), intent(out)         :: nfound
    integer(pInt), intent (In)         :: nalloc
    type(kdtree2_result), target :: results(:)

    !
    sr%qv => qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    sr%results => results

    call validate_query_storage(nalloc)
    sr%nalloc = nalloc
    sr%overflow = .false. 
    sr%ind => tp%ind
    sr%rearrange= tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound, results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    return
  end subroutine kdtree2_r_nearest

  subroutine kdtree2_r_nearest_around_point(tp,idxin,correltime,r2,&
   nfound,nalloc,results)
    !
    ! Like kdtree2_r_nearest, but around a point 'idxin' already existing
    ! in the data set. 
    ! 
    ! Results are NOT sorted unless tree was created with sort option.
    !
    type (kdtree2), pointer      :: tp
    integer(pInt), intent (In)         :: idxin, correltime, nalloc
    real(pReal), intent(in)             :: r2
    integer(pInt), intent(out)         :: nfound
    type(kdtree2_result), target :: results(:)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%results => results

    sr%nalloc = nalloc
    sr%overflow = .false.

    call validate_query_storage(nalloc)

    !    sr%dsl = HUGE(sr%dsl)    ! set to huge positive values
    !    sr%il = -1               ! set to invalid indexes

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%rearrange = tp%rearrange
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call kdtree2_sort_results(nfound,results)
    endif

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from kdtree2_r_nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    endif

    deallocate (sr%qv)
    return
  end subroutine kdtree2_r_nearest_around_point

  function kdtree2_r_count(tp,qv,r2) result(nfound)
    ! Count the number of neighbors within square distance 'r2'. 
    type (kdtree2), pointer   :: tp
    real(pReal), target, intent (In) :: qv(:)
    real(pReal), intent(in)          :: r2
    integer(pInt)                   :: nfound
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    sr%qv => qv
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    
    nullify(sr%results) ! for some reason, FTN 95 chokes on '=> null()'

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.
    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count

  function kdtree2_r_count_around_point(tp,idxin,correltime,r2) &
   result(nfound)
    ! Count the number of neighbors within square distance 'r2' around
    ! point 'idxin' with decorrelation time 'correltime'.
    !
    type (kdtree2), pointer :: tp
    integer(pInt), intent (In)    :: correltime, idxin
    real(pReal), intent(in)        :: r2
    integer(pInt)                 :: nfound
    ! ..
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic HUGE
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin)
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    nullify(sr%results)

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif
    sr%dimen = tp%dimen

    !
    !sr%dsl = Huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call search(tp%root)

    nfound = sr%nfound

    return
  end function kdtree2_r_count_around_point


  subroutine validate_query_storage(n)
    !
    ! make sure we have enough storage for n
    !
    integer(pInt), intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'KD_TREE_TRANS:  you did not provide enough storage for results(1:n)'
       stop
       return
    endif

    return
  end subroutine validate_query_storage

  function square_distance(d, iv,qv) result (res)
    ! distance between iv[1:n] and qv[1:n] 
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(pReal) :: res
    ! ..
    ! ..
    ! .. Scalar Arguments ..
    integer(pInt) :: d
    ! ..
    ! .. Array Arguments ..
    real(pReal) :: iv(:),qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:d)-qv(1:d))**2 )
  end function square_distance
  
  recursive subroutine search(node)
    !
    ! This is the innermost core routine of the kd-tree search.  Along
    ! with "process_terminal_node", it is the performance bottleneck. 
    !
    ! This version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (Tree_node), pointer          :: node
    ! ..
    type(tree_node),pointer            :: ncloser, nfarther
    !
    integer(pInt)                            :: cut_dim, i
    ! ..
    real(pReal)                               :: qval, dis
    real(pReal)                               :: ballsize
    real(pReal), pointer           :: qv(:)
    type(interval), pointer :: box(:) 

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call process_terminal_node_fixedball(node)
       else
          call process_terminal_node(node)
       endif
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
       endif

       if (associated(ncloser)) call search(ncloser)

       ! we may need to search the second node. 
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false. 
             !
             box => node%box(1:)
             do i=1,sr%dimen
                if (i .ne. cut_dim) then
                   dis = dis + dis2_from_bnd(qv(i),box(i)%lower,box(i)%upper)
                   if (dis > ballsize) then
                      return
                   endif
                endif
             end do
             
             !
             ! if we are still here then we need to search mroe.
             !
             call search(nfarther)
          endif
       endif
    end if
  end subroutine search


  real(pReal) function dis2_from_bnd(x,amin,amax) result (res)
    real(pReal), intent(in) :: x, amin,amax

    if (x > amax) then
       res = (x-amax)**2;
       return
    else
       if (x < amin) then
          res = (amin-x)**2;
          return
       else
          res = 0.0
          return
       endif
    endif
    return
  end function dis2_from_bnd

  logical function box_in_search_range(node, sr) result(res)
    !
    ! Return the distance from 'qv' to the CLOSEST corner of node's
    ! bounding box
    ! for all coordinates outside the box.   Coordinates inside the box
    ! contribute nothing to the distance.
    !
    type (tree_node), pointer :: node
    type (tree_search_record), pointer :: sr

    integer(pInt) :: dimen, i
    real(pReal)    :: dis, ballsize
    real(pReal)    :: l, u

    dimen = sr%dimen
    ballsize = sr%ballsize
    dis = 0.0
    res = .true.
    do i=1,dimen
       l = node%box(i)%lower
       u = node%box(i)%upper
       dis = dis + (dis2_from_bnd(sr%qv(i),l,u))
       if (dis > ballsize) then
          res = .false.
          return
       endif
    end do
    res = .true.
    return
  end function box_in_search_range


  subroutine process_terminal_node(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (tree_node), pointer          :: node
    !
    real(pReal), pointer          :: qv(:)
    integer(pInt), pointer       :: ind(:)
    real(pReal), pointer          :: data(:,:)
    !
    integer(pInt)                :: dimen, i, indexofi, k, centeridx, correltime
    real(pReal)                   :: ballsize, sd, newpri
    logical                :: rearrange
    type(pq), pointer      :: pqp 
    !
    ! copy values from sr to local variables
    !
    !
    ! Notice, making local pointers with an EXPLICIT lower bound
    ! seems to generate faster code.
    ! why?  I don't know.
    qv => sr%qv(1:) 
    pqp => sr%pq
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)     
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window? 
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       endif


       ! 
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1 
          newpri = pq_insert(pqp,sd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the PQ. 

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list. 
          ! Hence we replace that with the current one.
          !
          ballsize = pq_replace_max(pqp,sd,indexofi)
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize 

  end subroutine process_terminal_node

  subroutine process_terminal_node_fixedball(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (tree_node), pointer          :: node
    !
    real(pReal), pointer          :: qv(:)
    integer(pInt), pointer       :: ind(:)
    real(pReal), pointer          :: data(:,:)
    !
    integer(pInt)                :: nfound
    integer(pInt)                :: dimen, i, indexofi, k
    integer(pInt)                :: centeridx, correltime, nn
    real(pReal)                   :: ballsize, sd
    logical                :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    dimen = sr%dimen
    ballsize = sr%ballsize 
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       ! 
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity. 
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor. 
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one. 

       ! which index to the point do we use? 

       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (data(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       endif

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       endif

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%idx = indexofi
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine process_terminal_node_fixedball

  subroutine kdtree2_n_nearest_brute_force(tp,qv,nn,results) 
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(pReal), intent (In)       :: qv(:)
    integer(pInt), intent (In)    :: nn
    type(kdtree2_result)    :: results(:) 

    integer(pInt) :: i, j, k
    real(pReal), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
       results(i)%dis =  huge(1.0)
       results(i)%idx = -1
    end do
    do i = 1, tp%n
       if (all_distances(i)<results(nn)%dis) then
          ! insert it somewhere on the list
          do j = 1, nn
             if (all_distances(i)<results(j)%dis) exit
          end do
          ! now we know 'j'
          do k = nn - 1, j, -1
             results(k+1) = results(k)
          end do
          results(j)%dis = all_distances(i)
          results(j)%idx = i
       end if
    end do
    deallocate (all_distances)
  end subroutine kdtree2_n_nearest_brute_force
  

  subroutine kdtree2_r_nearest_brute_force(tp,qv,r2,nfound,results) 
    ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type (kdtree2), pointer :: tp
    real(pReal), intent (In)       :: qv(:)
    real(pReal), intent (In)       :: r2
    integer(pInt), intent(out)    :: nfound
    type(kdtree2_result)    :: results(:) 

    integer(pInt) :: i, nalloc
    real(pReal), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = square_distance(tp%dimen,qv,tp%the_data(:,i))
    end do
    
    nfound = 0
    nalloc = size(results,1)

    do i = 1, tp%n
       if (all_distances(i)< r2) then
          ! insert it somewhere on the list
          if (nfound .lt. nalloc) then
             nfound = nfound+1
             results(nfound)%dis = all_distances(i)
             results(nfound)%idx = i
          endif
       end if
    enddo
    deallocate (all_distances)

    call kdtree2_sort_results(nfound,results)


  end subroutine kdtree2_r_nearest_brute_force

  subroutine kdtree2_sort_results(nfound,results)
    !  Use after search to sort results(1:nfound) in order of increasing 
    !  distance.
    integer(pInt), intent(in)          :: nfound
    type(kdtree2_result), target :: results(:) 
    !
    !

    !THIS IS BUGGY WITH INTEL FORTRAN
    !    If (nfound .Gt. 1) Call heapsort(results(1:nfound)%dis,results(1:nfound)%ind,nfound)
    !
    if (nfound .gt. 1) call heapsort_struct(results,nfound)

    return
  end subroutine kdtree2_sort_results

  subroutine heapsort(a,ind,n)
    !
    ! Sort a(1:n) in ascending order, permuting ind(1:n) similarly.
    ! 
    ! If ind(k) = k upon input, then it will give a sort index upon output.
    !
    integer(pInt),intent(in)          :: n
    real(pReal), intent(inout)         :: a(:) 
    integer(pInt), intent(inout)      :: ind(:)

    !
    !
    real(pReal)        :: value   ! temporary for a value from a()
    integer(pInt)     :: ivalue  ! temporary for a value from ind()

    integer(pInt)     :: i,j
    integer(pInt)     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft); ivalue=ind(ileft)
       else
          value=a(iright); ivalue=ind(iright)
          a(iright)=a(1); ind(iright)=ind(1)
          iright=iright-1
          if (iright == 1) then
             a(1)=value;ind(1)=ivalue
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j) < a(j+1)) j=j+1
          endif
          if(value < a(j)) then
             a(i)=a(j); ind(i)=ind(j)
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value; ind(i)=ivalue
    end do
  end subroutine heapsort

  subroutine heapsort_struct(a,n)
    !
    ! Sort a(1:n) in ascending order
    ! 
    !
    integer(pInt),intent(in)                 :: n
    type(kdtree2_result),intent(inout) :: a(:)

    !
    !
    type(kdtree2_result) :: value ! temporary value

    integer(pInt)     :: i,j
    integer(pInt)     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return                  

    do 
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright) 
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          endif
          if(value%dis < a(j)%dis) then
             a(i)=a(j); 
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value
    end do
  end subroutine heapsort_struct

end module kdtree2_module
!#############################################################################################################################
! END KDTREE2
!#############################################################################################################################
