!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Interpolation data used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_quadrature
  use prec

  implicit none(type,external)
  private

  integer, parameter :: &
    maxOrder = 5                                                                                    !< maximum integration order
  real(pREAL),   dimension(2,3),       parameter :: &
    triangle    = reshape([-1.0_pREAL, -1.0_pREAL, &
                            1.0_pREAL, -1.0_pREAL, &
                           -1.0_pREAL,  1.0_pREAL], shape=shape(triangle))
  real(pREAL),   dimension(3,4),       parameter :: &
    tetrahedron = reshape([-1.0_pREAL, -1.0_pREAL, -1.0_pREAL, &
                            1.0_pREAL, -1.0_pREAL, -1.0_pREAL, &
                           -1.0_pREAL,  1.0_pREAL, -1.0_pREAL, &
                           -1.0_pREAL, -1.0_pREAL,  1.0_pREAL], shape=shape(tetrahedron))

  type :: group_real                                                                                !< variable length datatype
    real(pREAL), dimension(:), allocatable :: p
  end type group_real

  integer,          dimension(2:3,maxOrder), public, protected :: &
    FEM_nQuadrature                                                                                 !< number of quadrature points for spatial dimension(2-3) and interpolation order (1-maxOrder)
  type(group_real), dimension(2:3,maxOrder), public, protected :: &
    FEM_quadrature_weights, &                                                                       !< quadrature weights for each quadrature rule
    FEM_quadrature_points                                                                           !< quadrature point coordinates (in simplical system) for each quadrature rule

  public :: &
    FEM_quadrature_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize FEM interpolation data.
!--------------------------------------------------------------------------------------------------
subroutine FEM_quadrature_init()

  print'(/,1x,a)', '<<<+-  FEM_quadrature init  -+>>>'; flush(6)

  print'(/,1x,a)', 'L. Zhang et al., Journal of Computational Mathematics 27(1):89-96, 2009'
  print'(  1x,a)', 'https://www.jstor.org/stable/43693493'

!--------------------------------------------------------------------------------------------------
! 2D linear
  FEM_nQuadrature(2,1) = 1

  allocate(FEM_quadrature_weights(2,1)%p(FEM_nQuadrature(2,1)))
  FEM_quadrature_weights(2,1)%p(1) = 1._pREAL

  FEM_quadrature_points (2,1)%p = permutationStar3([1._pREAL/3._pREAL])

!--------------------------------------------------------------------------------------------------
! 2D quadratic
  FEM_nQuadrature(2,2) = 3

  allocate(FEM_quadrature_weights(2,2)%p(FEM_nQuadrature(2,2)))
  FEM_quadrature_weights(2,2)%p(1:3) = 1._pREAL/3._pREAL

  FEM_quadrature_points (2,2)%p = permutationStar21([1._pREAL/6._pREAL])

!--------------------------------------------------------------------------------------------------
! 2D cubic
  FEM_nQuadrature(2,3) = 6

  allocate(FEM_quadrature_weights(2,3)%p(FEM_nQuadrature(2,3)))
  FEM_quadrature_weights(2,3)%p(1:3) = 2.8114980244097965e-1_pReal
  FEM_quadrature_weights(2,3)%p(4:6) = 5.218353089235369e-2_pReal

  FEM_quadrature_points (2,3)%p = [ &
    permutationStar21([1.628828503958919e-1_pReal]), &
    permutationStar21([4.7791988356756370e-1_pReal])]

!--------------------------------------------------------------------------------------------------
! 2D quartic
  FEM_nQuadrature(2,4) = 6

  allocate(FEM_quadrature_weights(2,4)%p(FEM_nQuadrature(2,4)))
  FEM_quadrature_weights(2,4)%p(1:3) = (620._pReal + sqrt(213125._pReal - 53320._pReal * sqrt(10._pReal))) / 3720._pReal
  FEM_quadrature_weights(2,4)%p(4:6) = (620._pReal - sqrt(213125._pReal - 53320._pReal * sqrt(10._pReal))) / 3720._pReal

  FEM_quadrature_points (2,4)%p = [ &
    permutationStar21([(8._pReal - sqrt(10._pReal) + sqrt(38._pReal - 44._pReal * sqrt(0.4_pReal))) / 18._pReal]), &
    permutationStar21([(8._pReal - sqrt(10._pReal) - sqrt(38._pReal - 44._pReal * sqrt(0.4_pReal))) / 18._pReal])]

!--------------------------------------------------------------------------------------------------
! 2D quintic
  FEM_nQuadrature(2,5) = 7

  allocate(FEM_quadrature_weights(2,5)%p(FEM_nQuadrature(2,5)))
  FEM_quadrature_weights(2,5)%p(1:3)   = (155._pReal - sqrt(15._pReal)) / 1200._pReal
  FEM_quadrature_weights(2,5)%p(4:6)   = (155._pReal + sqrt(15._pReal)) / 1200._pReal
  FEM_quadrature_weights(2,5)%p(7:7)   = 9._pReal/40._pReal

  FEM_quadrature_points (2,5)%p = [ &
    permutationStar21([(6._pReal - sqrt(15._pReal)) / 21._pReal]), &
    permutationStar21([(6._pReal + sqrt(15._pReal)) / 21._pReal]), &
    permutationStar3([1._pREAL/3._pREAL])]

!--------------------------------------------------------------------------------------------------
! 3D linear
  FEM_nQuadrature(3,1) = 1

  allocate(FEM_quadrature_weights(3,1)%p(FEM_nQuadrature(3,1)))
  FEM_quadrature_weights(3,1)%p(1) = 1.0_pREAL

  FEM_quadrature_points (3,1)%p = permutationStar4([0.25_pREAL])

!--------------------------------------------------------------------------------------------------
! 3D quadratic
  FEM_nQuadrature(3,2) = 4

  allocate(FEM_quadrature_weights(3,2)%p(FEM_nQuadrature(3,2)))
  FEM_quadrature_weights(3,2)%p(1:4) = 0.25_pREAL

  FEM_quadrature_points (3,2)%p = permutationStar31([(5._pReal - sqrt(5._pReal)) / 20._pReal])

!--------------------------------------------------------------------------------------------------
! 3D cubic
  FEM_nQuadrature(3,3) = 8

  allocate(FEM_quadrature_weights(3,3)%p(FEM_nQuadrature(3,3)))
  FEM_quadrature_weights(3,3)%p(1:4)  = .125_pReal &
                                      + sqrt((1715161837._pReal - 406006699._pReal * sqrt(17._pReal)) / 23101._pReal) / 3120._pReal
  FEM_quadrature_weights(3,3)%p(5:8)  = .125_pReal &
                                      - sqrt((1715161837._pReal - 406006699._pReal * sqrt(17._pReal)) / 23101._pReal) / 3120._pReal

  FEM_quadrature_points (3,3)%p = [ &
    permutationStar31([(55._pReal - 3._pReal * sqrt(17._pReal) + sqrt(1022._pReal - 134._pReal * sqrt(17._pReal))) / 196._pReal]), &
    permutationStar31([(55._pReal - 3._pReal * sqrt(17._pReal) - sqrt(1022._pReal - 134._pReal * sqrt(17._pReal))) / 196._pReal])]

!--------------------------------------------------------------------------------------------------
! 3D quartic
  FEM_nQuadrature(3,4) = 14

  allocate(FEM_quadrature_weights(3,4)%p(FEM_nQuadrature(3,4)))
  FEM_quadrature_weights(3,4)%p(1:4)  = 7.3493043116361949e-2_pReal
  FEM_quadrature_weights(3,4)%p(5:8)  = 1.1268792571801585e-1_pReal
  FEM_quadrature_weights(3,4)%p(9:14) = 4.2546020777081467e-2_pReal

  FEM_quadrature_points (3,4)%p = [ &
    permutationStar31([9.273525031089123e-2_pREAL]), &
    permutationStar31([3.108859192633006e-1_pREAL]), &
    permutationStar22([4.5503704125649650e-2_pReal])]

!--------------------------------------------------------------------------------------------------
! 3D quintic
  FEM_nQuadrature(3,5) = 14

  allocate(FEM_quadrature_weights(3,5)%p(FEM_nQuadrature(3,5)))
  FEM_quadrature_weights(3,5)%p(1:4)  = 1.1268792571801585e-1_pReal
  FEM_quadrature_weights(3,5)%p(5:8)  = 7.3493043116361950e-2_pReal
  FEM_quadrature_weights(3,5)%p(9:14) = 4.2546020777081466e-2_pReal

  FEM_quadrature_points (3,5)%p = [ &
    permutationStar31([3.108859192633006e-1_pReal]), &
    permutationStar31([9.273525031089123e-2_pReal]), &
    permutationStar22([4.5503704125649649e-2_pReal])]

  call selfTest()

end subroutine FEM_quadrature_init


!--------------------------------------------------------------------------------------------------
!> @brief Star 3 permutation.
!--------------------------------------------------------------------------------------------------
pure function permutationStar3(point) result(qPt)

  real(pREAL), dimension(2)             :: qPt
  real(pREAL), dimension(1), intent(in) :: point


  qPt = pack(matmul(triangle,reshape([ &
    point(1), point(1), point(1)],[3,1])),.true.)

end function permutationStar3


!--------------------------------------------------------------------------------------------------
!> @brief Star 21 permutation.
!--------------------------------------------------------------------------------------------------
pure function permutationStar21(point) result(qPt)

  real(pREAL), dimension(6)             :: qPt
  real(pREAL), dimension(1), intent(in) :: point


  qPt = pack(matmul(triangle,reshape([ &
    point(1), point(1), 1.0_pREAL - 2.0_pREAL*point(1), &
    point(1), 1.0_pREAL - 2.0_pREAL*point(1), point(1), &
    1.0_pREAL - 2.0_pREAL*point(1), point(1), point(1)],[3,3])),.true.)

end function permutationStar21


!--------------------------------------------------------------------------------------------------
!> @brief Star 4 permutation.
!--------------------------------------------------------------------------------------------------
pure function permutationStar4(point) result(qPt)

  real(pREAL), dimension(3)             :: qPt
  real(pREAL), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), point(1), point(1)],[4,1])),.true.)

end function permutationStar4


!--------------------------------------------------------------------------------------------------
!> @brief Star 31 permutation.
!--------------------------------------------------------------------------------------------------
pure function permutationStar31(point) result(qPt)

  real(pREAL), dimension(12)            :: qPt
  real(pREAL), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), point(1), 1.0_pREAL - 3.0_pREAL*point(1), &
    point(1), point(1), 1.0_pREAL - 3.0_pREAL*point(1), point(1), &
    point(1), 1.0_pREAL - 3.0_pREAL*point(1), point(1), point(1), &
    1.0_pREAL - 3.0_pREAL*point(1), point(1), point(1), point(1)],[4,4])),.true.)

end function permutationStar31


!--------------------------------------------------------------------------------------------------
!> @brief Star 22 permutation.
!--------------------------------------------------------------------------------------------------
function permutationStar22(point) result(qPt)

  real(pREAL), dimension(18)            :: qPt
  real(pREAL), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), 0.5_pREAL - point(1), 0.5_pREAL - point(1), &
    point(1), 0.5_pREAL - point(1), point(1), 0.5_pREAL - point(1), &
    0.5_pREAL - point(1), point(1), point(1), 0.5_pREAL - point(1), &
    0.5_pREAL - point(1), point(1), 0.5_pREAL - point(1), point(1), &
    0.5_pREAL - point(1), 0.5_pREAL - point(1), point(1), point(1), &
    point(1), 0.5_pREAL - point(1), 0.5_pREAL - point(1), point(1)],[4,6])),.true.)

end function permutationStar22


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of quadrature weights and points.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer :: o, d, n
  real(pREAL), dimension(2:3), parameter :: w = [3.0_pREAL,2.0_pREAL]


  do d = lbound(FEM_quadrature_weights,1), ubound(FEM_quadrature_weights,1)
    do o = lbound(FEM_quadrature_weights(d,:),1), ubound(FEM_quadrature_weights(d,:),1)
      if (dNeq(sum(FEM_quadrature_weights(d,o)%p),1.0_pREAL,5e-15_pREAL)) &
        error stop 'quadrature weights'
    end do
  end do

  do d = lbound(FEM_quadrature_points,1), ubound(FEM_quadrature_points,1)
    do o = lbound(FEM_quadrature_points(d,:),1), ubound(FEM_quadrature_points(d,:),1)
      n = size(FEM_quadrature_points(d,o)%p,1)/d
      if (any(dNeq(sum(reshape(FEM_quadrature_points(d,o)%p,[d,n]),2),-real(n,pREAL)/w(d),1.e-14_pREAL))) &
        error stop 'quadrature points'
    end do
  end do

end subroutine selfTest

end module FEM_quadrature
