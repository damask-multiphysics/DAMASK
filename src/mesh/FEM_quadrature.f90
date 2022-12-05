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
  real(pReal),   dimension(2,3),       parameter :: &
    triangle    = reshape([-1.0_pReal, -1.0_pReal, &
                            1.0_pReal, -1.0_pReal, &
                           -1.0_pReal,  1.0_pReal], shape=[2,3])
  real(pReal),   dimension(3,4),       parameter :: &
    tetrahedron = reshape([-1.0_pReal, -1.0_pReal, -1.0_pReal, &
                            1.0_pReal, -1.0_pReal, -1.0_pReal, &
                           -1.0_pReal,  1.0_pReal, -1.0_pReal, &
                           -1.0_pReal, -1.0_pReal,  1.0_pReal], shape=[3,4])

  type :: group_float                                                                               !< variable length datatype
    real(pReal), dimension(:), allocatable :: p
  end type group_float

  integer,             dimension(2:3,maxOrder), public, protected :: &
    FEM_nQuadrature                                                                                 !< number of quadrature points for spatial dimension(2-3) and interpolation order (1-maxOrder)
  type(group_float),   dimension(2:3,maxOrder), public, protected :: &
    FEM_quadrature_weights, &                                                                       !< quadrature weights for each quadrature rule
    FEM_quadrature_points                                                                           !< quadrature point coordinates (in simplical system) for each quadrature rule

  public :: &
    FEM_quadrature_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief initializes FEM interpolation data
!--------------------------------------------------------------------------------------------------
subroutine FEM_quadrature_init()

  print'(/,1x,a)', '<<<+-  FEM_quadrature init  -+>>>'; flush(6)

  print'(/,1x,a)', 'L. Zhang et al., Journal of Computational Mathematics 27(1):89-96, 2009'
  print'(  1x,a)', 'https://www.jstor.org/stable/43693493'

!--------------------------------------------------------------------------------------------------
! 2D linear
  FEM_nQuadrature(2,1) = 1

  allocate(FEM_quadrature_weights(2,1)%p(FEM_nQuadrature(2,1)))
  FEM_quadrature_weights(2,1)%p(1) = 1._pReal

  FEM_quadrature_points (2,1)%p = permutationStar3([1._pReal/3._pReal])

!--------------------------------------------------------------------------------------------------
! 2D quadratic
  FEM_nQuadrature(2,2) = 3

  allocate(FEM_quadrature_weights(2,2)%p(FEM_nQuadrature(2,2)))
  FEM_quadrature_weights(2,2)%p(1:3) = 1._pReal/3._pReal

  FEM_quadrature_points (2,2)%p = permutationStar21([1._pReal/6._pReal])

!--------------------------------------------------------------------------------------------------
! 2D cubic
  FEM_nQuadrature(2,3) = 6

  allocate(FEM_quadrature_weights(2,3)%p(FEM_nQuadrature(2,3)))
  FEM_quadrature_weights(2,3)%p(1:3) = 2.2338158967801147e-1_pReal
  FEM_quadrature_weights(2,3)%p(4:6) = 1.0995174365532187e-1_pReal

  FEM_quadrature_points (2,3)%p = [ &
    permutationStar21([4.4594849091596489e-1_pReal]), &
    permutationStar21([9.157621350977074e-2_pReal]) ]

!--------------------------------------------------------------------------------------------------
! 2D quartic
  FEM_nQuadrature(2,4) = 12

  allocate(FEM_quadrature_weights(2,4)%p(FEM_nQuadrature(2,4)))
  FEM_quadrature_weights(2,4)%p(1:3)  = 1.1678627572637937e-1_pReal
  FEM_quadrature_weights(2,4)%p(4:6)  = 5.0844906370206817e-2_pReal
  FEM_quadrature_weights(2,4)%p(7:12) = 8.285107561837358e-2_pReal

  FEM_quadrature_points (2,4)%p = [ &
    permutationStar21([2.4928674517091042e-1_pReal]), &
    permutationStar21([6.308901449150223e-2_pReal]), &
    permutationStar111([3.1035245103378440e-1_pReal, 5.3145049844816947e-2_pReal]) ]

!--------------------------------------------------------------------------------------------------
! 2D quintic
  FEM_nQuadrature(2,5) = 16

  allocate(FEM_quadrature_weights(2,5)%p(FEM_nQuadrature(2,5)))
  FEM_quadrature_weights(2,5)%p(1:1)   = 1.4431560767778717e-1_pReal
  FEM_quadrature_weights(2,5)%p(2:4)   = 9.509163426728463e-2_pReal
  FEM_quadrature_weights(2,5)%p(5:7)   = 1.0321737053471825e-1_pReal
  FEM_quadrature_weights(2,5)%p(8:10)  = 3.2458497623198080e-2_pReal
  FEM_quadrature_weights(2,5)%p(11:16) = 2.7230314174434994e-2_pReal

  FEM_quadrature_points (2,5)%p = [ &
    permutationStar3([1._pReal/3._pReal]), &
    permutationStar21([4.5929258829272316e-1_pReal]), &
    permutationStar21([1.705693077517602e-1_pReal]), &
    permutationStar21([5.0547228317030975e-2_pReal]), &
    permutationStar111([2.631128296346381e-1_pReal, 8.3947774099576053e-2_pReal]) ]

!--------------------------------------------------------------------------------------------------
! 3D linear
  FEM_nQuadrature(3,1) = 1

  allocate(FEM_quadrature_weights(3,1)%p(FEM_nQuadrature(3,1)))
  FEM_quadrature_weights(3,1)%p(1) = 1.0_pReal

  FEM_quadrature_points (3,1)%p = permutationStar4([0.25_pReal])

!--------------------------------------------------------------------------------------------------
! 3D quadratic
  FEM_nQuadrature(3,2) = 4

  allocate(FEM_quadrature_weights(3,2)%p(FEM_nQuadrature(3,2)))
  FEM_quadrature_weights(3,2)%p(1:4) = 0.25_pReal

  FEM_quadrature_points (3,2)%p = permutationStar31([1.3819660112501052e-1_pReal])

!--------------------------------------------------------------------------------------------------
! 3D cubic
  FEM_nQuadrature(3,3) = 14

  allocate(FEM_quadrature_weights(3,3)%p(FEM_nQuadrature(3,3)))
  FEM_quadrature_weights(3,3)%p(1:4)  = 7.3493043116361949e-2_pReal
  FEM_quadrature_weights(3,3)%p(5:8)  = 1.1268792571801585e-1_pReal
  FEM_quadrature_weights(3,3)%p(9:14) = 4.2546020777081467e-2_pReal

  FEM_quadrature_points (3,3)%p = [ &
    permutationStar31([9.273525031089123e-2_pReal]), &
    permutationStar31([3.108859192633006e-1_pReal]), &
    permutationStar22([4.5503704125649649e-2_pReal]) ]

!--------------------------------------------------------------------------------------------------
! 3D quartic (lower precision/unknown source)
  FEM_nQuadrature(3,4) = 35

  allocate(FEM_quadrature_weights(3,4)%p(FEM_nQuadrature(3,4)))
  FEM_quadrature_weights(3,4)%p(1:4)   = 0.0021900463965388_pReal
  FEM_quadrature_weights(3,4)%p(5:16)  = 0.0143395670177665_pReal
  FEM_quadrature_weights(3,4)%p(17:22) = 0.0250305395686746_pReal
  FEM_quadrature_weights(3,4)%p(23:34) = 0.0479839333057554_pReal
  FEM_quadrature_weights(3,4)%p(35)    = 0.0931745731195340_pReal

  FEM_quadrature_points (3,4)%p = [ &
    permutationStar31([0.0267367755543735_pReal]), &
    permutationStar211([0.0391022406356488_pReal, 0.7477598884818090_pReal]), &
    permutationStar22([0.4547545999844830_pReal]), &
    permutationStar211([0.2232010379623150_pReal, 0.0504792790607720_pReal]), &
    permutationStar4([0.25_pReal]) ]

!--------------------------------------------------------------------------------------------------
! 3D quintic (lower precision/unknown source)
  FEM_nQuadrature(3,5) = 56

  allocate(FEM_quadrature_weights(3,5)%p(FEM_nQuadrature(3,5)))
  FEM_quadrature_weights(3,5)%p(1:4)   = 0.0010373112336140_pReal
  FEM_quadrature_weights(3,5)%p(5:16)  = 0.0096016645399480_pReal
  FEM_quadrature_weights(3,5)%p(17:28) = 0.0164493976798232_pReal
  FEM_quadrature_weights(3,5)%p(29:40) = 0.0153747766513310_pReal
  FEM_quadrature_weights(3,5)%p(41:52) = 0.0293520118375230_pReal
  FEM_quadrature_weights(3,5)%p(53:56) = 0.0366291366405108_pReal

  FEM_quadrature_points (3,5)%p = [ &
    permutationStar31([0.0149520651530592_pReal]), &
    permutationStar211([0.0340960211962615_pReal, 0.1518319491659370_pReal]), &
    permutationStar211([0.0462051504150017_pReal, 0.3549340560639790_pReal]), &
    permutationStar211([0.2281904610687610_pReal, 0.0055147549744775_pReal]), &
    permutationStar211([0.3523052600879940_pReal, 0.0992057202494530_pReal]), &
    permutationStar31([0.1344783347929940_pReal]) ]

  call selfTest

end subroutine FEM_quadrature_init


!--------------------------------------------------------------------------------------------------
!> @brief star 3 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar3(point) result(qPt)

  real(pReal), dimension(2)             :: qPt
  real(pReal), dimension(1), intent(in) :: point


  qPt = pack(matmul(triangle,reshape([ &
    point(1), point(1), point(1)],[3,1])),.true.)

end function permutationStar3


!--------------------------------------------------------------------------------------------------
!> @brief star 21 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar21(point) result(qPt)

  real(pReal), dimension(6)             :: qPt
  real(pReal), dimension(1), intent(in) :: point


  qPt = pack(matmul(triangle,reshape([ &
    point(1), point(1), 1.0_pReal - 2.0_pReal*point(1), &
    point(1), 1.0_pReal - 2.0_pReal*point(1), point(1), &
    1.0_pReal - 2.0_pReal*point(1), point(1), point(1)],[3,3])),.true.)

end function permutationStar21


!--------------------------------------------------------------------------------------------------
!> @brief star 111 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar111(point) result(qPt)

  real(pReal), dimension(12)            :: qPt
  real(pReal), dimension(2), intent(in) :: point


  qPt = pack(matmul(triangle,reshape([ &
    point(1), point(2), 1.0_pReal - point(1) - point(2), &
    point(1), 1.0_pReal - point(1) - point(2), point(2), &
    point(2), point(1), 1.0_pReal - point(1) - point(2), &
    point(2), 1.0_pReal - point(1) - point(2), point(1), &
    1.0_pReal - point(1) - point(2), point(2), point(1), &
    1.0_pReal - point(1) - point(2), point(1), point(2)],[3,6])),.true.)

end function permutationStar111


!--------------------------------------------------------------------------------------------------
!> @brief star 4 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar4(point) result(qPt)

  real(pReal), dimension(3)             :: qPt
  real(pReal), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), point(1), point(1)],[4,1])),.true.)

end function permutationStar4


!--------------------------------------------------------------------------------------------------
!> @brief star 31 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar31(point) result(qPt)

  real(pReal), dimension(12)            :: qPt
  real(pReal), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), point(1), 1.0_pReal - 3.0_pReal*point(1), &
    point(1), point(1), 1.0_pReal - 3.0_pReal*point(1), point(1), &
    point(1), 1.0_pReal - 3.0_pReal*point(1), point(1), point(1), &
    1.0_pReal - 3.0_pReal*point(1), point(1), point(1), point(1)],[4,4])),.true.)

end function permutationStar31


!--------------------------------------------------------------------------------------------------
!> @brief star 22 permutation of input
!--------------------------------------------------------------------------------------------------
function permutationStar22(point) result(qPt)

  real(pReal), dimension(18)            :: qPt
  real(pReal), dimension(1), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), 0.5_pReal - point(1), 0.5_pReal - point(1), &
    point(1), 0.5_pReal - point(1), point(1), 0.5_pReal - point(1), &
    0.5_pReal - point(1), point(1), point(1), 0.5_pReal - point(1), &
    0.5_pReal - point(1), point(1), 0.5_pReal - point(1), point(1), &
    0.5_pReal - point(1), 0.5_pReal - point(1), point(1), point(1), &
    point(1), 0.5_pReal - point(1), 0.5_pReal - point(1), point(1)],[4,6])),.true.)

end function permutationStar22


!--------------------------------------------------------------------------------------------------
!> @brief star 211 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar211(point) result(qPt)

  real(pReal), dimension(36)            :: qPt
  real(pReal), dimension(2), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(1), point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2), &
    point(1), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(2), &
    point(1), point(2), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), &
    point(1), point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), &
    point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(2), &
    point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(2), point(1), &
    point(2), point(1), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), &
    point(2), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), &
    point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(1), &
    1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(1), point(2), &
    1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(2), point(1), &
    1.0_pReal - 2.0_pReal*point(1) - point(2), point(2), point(1), point(1)],[4,12])),.true.)

end function permutationStar211


!--------------------------------------------------------------------------------------------------
!> @brief star 1111 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar1111(point) result(qPt)

  real(pReal), dimension(72)            :: qPt
  real(pReal), dimension(3), intent(in) :: point


  qPt = pack(matmul(tetrahedron,reshape([ &
    point(1), point(2), point(3), 1.0_pReal - point(1) - point(2)- point(3), &
    point(1), point(2), 1.0_pReal - point(1) - point(2)- point(3), point(3), &
    point(1), point(3), point(2), 1.0_pReal - point(1) - point(2)- point(3), &
    point(1), point(3), 1.0_pReal - point(1) - point(2)- point(3), point(2), &
    point(1), 1.0_pReal - point(1) - point(2)- point(3), point(2), point(3), &
    point(1), 1.0_pReal - point(1) - point(2)- point(3), point(3), point(2), &
    point(2), point(1), point(3), 1.0_pReal - point(1) - point(2)- point(3), &
    point(2), point(1), 1.0_pReal - point(1) - point(2)- point(3), point(3), &
    point(2), point(3), point(1), 1.0_pReal - point(1) - point(2)- point(3), &
    point(2), point(3), 1.0_pReal - point(1) - point(2)- point(3), point(1), &
    point(2), 1.0_pReal - point(1) - point(2)- point(3), point(1), point(3), &
    point(2), 1.0_pReal - point(1) - point(2)- point(3), point(3), point(1), &
    point(3), point(1), point(2), 1.0_pReal - point(1) - point(2)- point(3), &
    point(3), point(1), 1.0_pReal - point(1) - point(2)- point(3), point(2), &
    point(3), point(2), point(1), 1.0_pReal - point(1) - point(2)- point(3), &
    point(3), point(2), 1.0_pReal - point(1) - point(2)- point(3), point(1), &
    point(3), 1.0_pReal - point(1) - point(2)- point(3), point(1), point(2), &
    point(3), 1.0_pReal - point(1) - point(2)- point(3), point(2), point(1), &
    1.0_pReal - point(1) - point(2)- point(3), point(1), point(2), point(3), &
    1.0_pReal - point(1) - point(2)- point(3), point(1), point(3), point(2), &
    1.0_pReal - point(1) - point(2)- point(3), point(2), point(1), point(3), &
    1.0_pReal - point(1) - point(2)- point(3), point(2), point(3), point(1), &
    1.0_pReal - point(1) - point(2)- point(3), point(3), point(1), point(2), &
    1.0_pReal - point(1) - point(2)- point(3), point(3), point(2), point(1)],[4,24])),.true.)

end function permutationStar1111


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of quadrature weights and points.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  integer :: o, d, n
  real(pReal), dimension(2:3), parameter :: w = [3.0_pReal,2.0_pReal]


  do d = lbound(FEM_quadrature_weights,1), ubound(FEM_quadrature_weights,1)
    do o = lbound(FEM_quadrature_weights(d,:),1), ubound(FEM_quadrature_weights(d,:),1)
      if (dNeq(sum(FEM_quadrature_weights(d,o)%p),1.0_pReal,5e-15_pReal)) &
        error stop 'quadrature weights'
    end do
  end do

  do d = lbound(FEM_quadrature_points,1), ubound(FEM_quadrature_points,1)
    do o = lbound(FEM_quadrature_points(d,:),1), ubound(FEM_quadrature_points(d,:),1)
      n = size(FEM_quadrature_points(d,o)%p,1)/d
      if (any(dNeq(sum(reshape(FEM_quadrature_points(d,o)%p,[d,n]),2),-real(n,pReal)/w(d),1.e-14_pReal))) &
        error stop 'quadrature points'
    end do
  end do

end subroutine selfTest

end module FEM_quadrature
