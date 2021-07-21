!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Interpolation data used by the FEM solver
!--------------------------------------------------------------------------------------------------
module FEM_quadrature
  use prec

  implicit none
  private

  integer, parameter :: &
    maxOrder = 5                                                                                    !< current max interpolation set at cubic (intended to be arbitrary)
  real(pReal),   dimension(2,3),       parameter :: &
    triangle    = reshape([-1.0_pReal, -1.0_pReal, &
                            1.0_pReal, -1.0_pReal, &
                           -1.0_pReal,  1.0_pReal], shape=[2,3])
  real(pReal),   dimension(3,4),       parameter :: &
    tetrahedron = reshape([-1.0_pReal, -1.0_pReal, -1.0_pReal, &
                            1.0_pReal, -1.0_pReal, -1.0_pReal, &
                           -1.0_pReal,  1.0_pReal, -1.0_pReal, &
                           -1.0_pReal, -1.0_pReal,  1.0_pReal], shape=[3,4])

  type :: group_float                                                                               !< variable length datatype used for storage of state
    real(pReal), dimension(:), pointer :: p
  end type group_float

  integer,             dimension(2:3,maxOrder), public, protected :: &
    FEM_nQuadrature                                                                                 !< number of quadrature points for a given spatial dimension(2-3) and interpolation order(1-maxOrder)
  type(group_float),   dimension(2:3,maxOrder), public, protected :: &
    FEM_quadrature_weights, &                                                                       !< quadrature weights for each quadrature rule
    FEM_quadrature_points                                                                           !< quadrature point coordinates (in simplical system) for each quadrature rule

  public :: &
    FEM_quadrature_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief initializes FEM interpolation data
!--------------------------------------------------------------------------------------------------
subroutine FEM_quadrature_init

  print'(/,a)', ' <<<+-  FEM_quadrature init  -+>>>'; flush(6)

  print*, 'L. Zhang et al., Journal of Computational Mathematics 27(1):89-96, 2009'
  print*, 'https://www.jstor.org/stable/43693493'

!--------------------------------------------------------------------------------------------------
! 2D linear
  FEM_nQuadrature(2,1) = 1

  allocate(FEM_quadrature_weights(2,1)%p(1))
  FEM_quadrature_weights(2,1)%p(1)   = 1.0_pReal

  allocate(FEM_quadrature_points (2,1)%p(2))
  FEM_quadrature_points (2,1)%p(1:2) = permutationStar3([1.0_pReal/3.0_pReal])

!--------------------------------------------------------------------------------------------------
! 2D quadratic
  FEM_nQuadrature(2,2) = 3

  allocate(FEM_quadrature_weights(2,2)%p(3))
  FEM_quadrature_weights(2,2)%p(1:3) = 1.0_pReal/3.0_pReal

  allocate(FEM_quadrature_points (2,2)%p(6))
  FEM_quadrature_points (2,2)%p(1:6) = permutationStar21([1.0_pReal/6.0_pReal])

!--------------------------------------------------------------------------------------------------
! 2D cubic
  FEM_nQuadrature(2,3) = 6

  allocate(FEM_quadrature_weights(2,3)%p(6))                                    
  FEM_quadrature_weights(2,3)%p(1:3) = 0.2233815896780115_pReal
  FEM_quadrature_weights(2,3)%p(4:6) = 0.1099517436553219_pReal

  allocate(FEM_quadrature_points (2,3)%p(12))
  FEM_quadrature_points (2,3)%p(1:6) = permutationStar21([0.4459484909159649_pReal])
  FEM_quadrature_points (2,3)%p(7:12)= permutationStar21([0.09157621350977074_pReal])

!--------------------------------------------------------------------------------------------------
! 2D quartic
  FEM_nQuadrature(2,4) = 12

  allocate(FEM_quadrature_weights(2,4)%p(12))                                     
  FEM_quadrature_weights(2,4)%p(1:3)  = 0.1167862757263794_pReal
  FEM_quadrature_weights(2,4)%p(4:6)  = 0.0508449063702068_pReal
  FEM_quadrature_weights(2,4)%p(7:12) = 0.08285107561837358_pReal

  allocate(FEM_quadrature_points (2,4)%p(24))                                                         
  FEM_quadrature_points (2,4)%p(1:6)  = permutationStar21([0.2492867451709104_pReal])
  FEM_quadrature_points (2,4)%p(7:12) = permutationStar21([0.06308901449150223_pReal])
  FEM_quadrature_points (2,4)%p(13:24)= permutationStar111([0.3103524510337844_pReal, 0.05314504984481695_pReal])

!--------------------------------------------------------------------------------------------------
! 2D quintic
  FEM_nQuadrature(2,5) = 16

  allocate(FEM_quadrature_weights(2,5)%p(16))
  FEM_quadrature_weights(2,5)%p(1  )  = 0.1443156076777871_pReal
  FEM_quadrature_weights(2,5)%p(2:4)  = 0.09509163426728462_pReal
  FEM_quadrature_weights(2,5)%p(5:7)  = 0.1032173705347183_pReal
  FEM_quadrature_weights(2,5)%p(8:10) = 0.03245849762319808_pReal
  FEM_quadrature_weights(2,5)%p(11:16)= 0.02723031417443499_pReal

  allocate(FEM_quadrature_points (2,5)%p(32))
                                                          
  FEM_quadrature_points (2,5)%p(1:2)  = permutationStar3([0.3333333333333333_pReal])
  FEM_quadrature_points (2,5)%p(3:8)  = permutationStar21([0.4592925882927231_pReal])
  FEM_quadrature_points (2,5)%p(9:14) = permutationStar21([0.1705693077517602_pReal])
  FEM_quadrature_points (2,5)%p(15:20)= permutationStar21([0.0505472283170310_pReal])
  FEM_quadrature_points (2,5)%p(21:32)= permutationStar111([0.2631128296346381_pReal, 0.008394777409957605_pReal])

!--------------------------------------------------------------------------------------------------
! 3D linear
  FEM_nQuadrature(3,1) = 1

  allocate(FEM_quadrature_weights(3,1)%p(1))
  FEM_quadrature_weights(3,1)%p(1)  = 1.0_pReal

  allocate(FEM_quadrature_points (3,1)%p(3))
  FEM_quadrature_points (3,1)%p(1:3)= permutationStar4([0.25_pReal])

!--------------------------------------------------------------------------------------------------
! 3D quadratic
  FEM_nQuadrature(3,2) = 4

  allocate(FEM_quadrature_weights(3,2)%p(4))
  FEM_quadrature_weights(3,2)%p(1:4) = 0.25_pReal

  allocate(FEM_quadrature_points (3,2)%p(12))
                                                          
  FEM_quadrature_points (3,2)%p(1:12)= permutationStar31([0.1381966011250105_pReal])

!--------------------------------------------------------------------------------------------------
! 3D cubic
  FEM_nQuadrature(3,3) = 14

  allocate(FEM_quadrature_weights(3,3)%p(14))
                                        
  FEM_quadrature_weights(3,3)%p(5:8)  = 0.1126879257180159_pReal
  FEM_quadrature_weights(3,3)%p(1:4)  = 0.0734930431163620_pReal
  FEM_quadrature_weights(3,3)%p(9:14) = 0.0425460207770815_pReal

  allocate(FEM_quadrature_points (3,3)%p(42))
                                                           
  FEM_quadrature_points (3,3)%p(1:12) = permutationStar31([0.09273525031089123_pReal])
  FEM_quadrature_points (3,3)%p(13:24)= permutationStar31([0.3108859192633006_pReal])
  FEM_quadrature_points (3,3)%p(25:42)= permutationStar22([0.04550370412564965_pReal])

!--------------------------------------------------------------------------------------------------
! 3D quartic
  FEM_nQuadrature(3,4) = 35

  allocate(FEM_quadrature_weights(3,4)%p(35))
                                          
  FEM_quadrature_weights(3,4)%p(1:4)    = 0.0021900463965388_pReal
  FEM_quadrature_weights(3,4)%p(5:16)   = 0.0143395670177665_pReal
  FEM_quadrature_weights(3,4)%p(17:22)  = 0.0250305395686746_pReal
  FEM_quadrature_weights(3,4)%p(23:34)  = 0.0479839333057554_pReal
  FEM_quadrature_weights(3,4)%p(35)     = 0.0931745731195340_pReal

  allocate(FEM_quadrature_points (3,4)%p(105))
                                                                                        
  FEM_quadrature_points (3,4)%p(1:12)   = permutationStar31([0.0267367755543735_pReal])
  FEM_quadrature_points (3,4)%p(13:48)  = permutationStar211([0.0391022406356488_pReal, 0.7477598884818090_pReal])
  FEM_quadrature_points (3,4)%p(49:66)  = permutationStar22([0.4547545999844830_pReal])
  FEM_quadrature_points (3,4)%p(67:102) = permutationStar211([0.2232010379623150_pReal, 0.0504792790607720_pReal])
  FEM_quadrature_points (3,4)%p(103:105)= permutationStar4([0.25_pReal])

!--------------------------------------------------------------------------------------------------
! 3D quintic
  FEM_nQuadrature(3,5) = 56

  allocate(FEM_quadrature_weights(3,5)%p(56))
                                          
  FEM_quadrature_weights(3,5)%p(1:4)    = 0.0010373112336140_pReal
  FEM_quadrature_weights(3,5)%p(5:16)   = 0.0096016645399480_pReal
  FEM_quadrature_weights(3,5)%p(17:28)  = 0.0164493976798232_pReal
  FEM_quadrature_weights(3,5)%p(29:40)  = 0.0153747766513310_pReal
  FEM_quadrature_weights(3,5)%p(41:52)  = 0.0293520118375230_pReal
  FEM_quadrature_weights(3,5)%p(53:56)  = 0.0366291366405108_pReal

  allocate(FEM_quadrature_points (3,5)%p(168))
                                                                                        
  FEM_quadrature_points (3,5)%p(1:12)   = permutationStar31([0.0149520651530592_pReal])
  FEM_quadrature_points (3,5)%p(13:48)  = permutationStar211([0.0340960211962615_pReal, 0.1518319491659370_pReal])
  FEM_quadrature_points (3,5)%p(49:84)  = permutationStar211([0.0462051504150017_pReal, 0.3549340560639790_pReal])
  FEM_quadrature_points (3,5)%p(85:120) = permutationStar211([0.2281904610687610_pReal, 0.0055147549744775_pReal])
  FEM_quadrature_points (3,5)%p(121:156)= permutationStar211([0.3523052600879940_pReal, 0.0992057202494530_pReal])
  FEM_quadrature_points (3,5)%p(157:168)= permutationStar31([0.1344783347929940_pReal])

end subroutine FEM_quadrature_init


!--------------------------------------------------------------------------------------------------
!> @brief star 3 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar3(point) result(qPt)

  real(pReal), dimension(2)             :: qPt
  real(pReal), dimension(1), intent(in) :: point

  real(pReal), dimension(3,1) :: temp

  temp(:,1) = [point(1), point(1), point(1)]

  qPt = reshape(matmul(triangle, temp),[2])

end function permutationStar3


!--------------------------------------------------------------------------------------------------
!> @brief star 21 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar21(point) result(qPt)

  real(pReal), dimension(6)             :: qPt
  real(pReal), dimension(1), intent(in) :: point

  real(pReal), dimension(3,3) :: temp

  temp(:,1) = [point(1), point(1), 1.0_pReal - 2.0_pReal*point(1)]
  temp(:,2) = [point(1), 1.0_pReal - 2.0_pReal*point(1), point(1)]
  temp(:,3) = [1.0_pReal - 2.0_pReal*point(1), point(1), point(1)]

 qPt = reshape(matmul(triangle, temp),[6])

end function permutationStar21


!--------------------------------------------------------------------------------------------------
!> @brief star 111 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar111(point) result(qPt)

  real(pReal), dimension(12)            :: qPt
  real(pReal), dimension(2), intent(in) :: point

  real(pReal), dimension(3,6) :: temp

  temp(:,1) = [point(1), point(2), 1.0_pReal - point(1) - point(2)]
  temp(:,2) = [point(1), 1.0_pReal - point(1) - point(2), point(2)]
  temp(:,3) = [point(2), point(1), 1.0_pReal - point(1) - point(2)]
  temp(:,4) = [point(2), 1.0_pReal - point(1) - point(2), point(1)]
  temp(:,5) = [1.0_pReal - point(1) - point(2), point(2), point(1)]
  temp(:,6) = [1.0_pReal - point(1) - point(2), point(1), point(2)]

  qPt = reshape(matmul(triangle, temp),[12])

end function permutationStar111


!--------------------------------------------------------------------------------------------------
!> @brief star 4 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar4(point) result(qPt)

  real(pReal), dimension(3)             :: qPt
  real(pReal), dimension(1), intent(in) :: point

  real(pReal), dimension(4,1) :: temp

  temp(:,1) = [point(1), point(1), point(1), point(1)]

  qPt = reshape(matmul(tetrahedron, temp),[3])

end function permutationStar4


!--------------------------------------------------------------------------------------------------
!> @brief star 31 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar31(point) result(qPt)

  real(pReal), dimension(12)            :: qPt
  real(pReal), dimension(1), intent(in) :: point

  real(pReal), dimension(4,4) :: temp

  temp(:,1) = [point(1), point(1), point(1), 1.0_pReal - 3.0_pReal*point(1)]
  temp(:,2) = [point(1), point(1), 1.0_pReal - 3.0_pReal*point(1), point(1)]
  temp(:,3) = [point(1), 1.0_pReal - 3.0_pReal*point(1), point(1), point(1)]
  temp(:,4) = [1.0_pReal - 3.0_pReal*point(1), point(1), point(1), point(1)]

  qPt = reshape(matmul(tetrahedron, temp),[12])

end function permutationStar31


!--------------------------------------------------------------------------------------------------
!> @brief star 22 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar22(point) result(qPt)

  real(pReal), dimension(18)            :: qPt
  real(pReal), dimension(1), intent(in) :: point

  real(pReal), dimension(4,6) :: temp

  temp(:,1) = [point(1), point(1), 0.5_pReal - point(1), 0.5_pReal - point(1)]
  temp(:,2) = [point(1), 0.5_pReal - point(1), point(1), 0.5_pReal - point(1)]
  temp(:,3) = [0.5_pReal - point(1), point(1), point(1), 0.5_pReal - point(1)]
  temp(:,4) = [0.5_pReal - point(1), point(1), 0.5_pReal - point(1), point(1)]
  temp(:,5) = [0.5_pReal - point(1), 0.5_pReal - point(1), point(1), point(1)]
  temp(:,6) = [point(1), 0.5_pReal - point(1), 0.5_pReal - point(1), point(1)]

  qPt = reshape(matmul(tetrahedron, temp),[18])

end function permutationStar22


!--------------------------------------------------------------------------------------------------
!> @brief star 211 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar211(point) result(qPt)

  real(pReal), dimension(36)            :: qPt
  real(pReal), dimension(2), intent(in) :: point

  real(pReal), dimension(4,12) :: temp

  temp(:,1 ) = [point(1), point(1), point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2)]
  temp(:,2 ) = [point(1), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(2)]
  temp(:,3 ) = [point(1), point(2), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2)]
  temp(:,4 ) = [point(1), point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1)]
  temp(:,5 ) = [point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(2)]
  temp(:,6 ) = [point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(2), point(1)]
  temp(:,7 ) = [point(2), point(1), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2)]
  temp(:,8 ) = [point(2), point(1), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1)]
  temp(:,9 ) = [point(2), 1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(1)]
  temp(:,10) = [1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(1), point(2)]
  temp(:,11) = [1.0_pReal - 2.0_pReal*point(1) - point(2), point(1), point(2), point(1)]
  temp(:,12) = [1.0_pReal - 2.0_pReal*point(1) - point(2), point(2), point(1), point(1)]

  qPt = reshape(matmul(tetrahedron, temp),[36])

end function permutationStar211


!--------------------------------------------------------------------------------------------------
!> @brief star 1111 permutation of input
!--------------------------------------------------------------------------------------------------
pure function permutationStar1111(point) result(qPt)

  real(pReal), dimension(72)            :: qPt
  real(pReal), dimension(3), intent(in) :: point

  real(pReal), dimension(4,24) :: temp

  temp(:,1 ) = [point(1), point(2), point(3), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,2 ) = [point(1), point(2), 1.0_pReal - point(1) - point(2)- point(3), point(3)]
  temp(:,3 ) = [point(1), point(3), point(2), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,4 ) = [point(1), point(3), 1.0_pReal - point(1) - point(2)- point(3), point(2)]
  temp(:,5 ) = [point(1), 1.0_pReal - point(1) - point(2)- point(3), point(2), point(3)]
  temp(:,6 ) = [point(1), 1.0_pReal - point(1) - point(2)- point(3), point(3), point(2)]
  temp(:,7 ) = [point(2), point(1), point(3), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,8 ) = [point(2), point(1), 1.0_pReal - point(1) - point(2)- point(3), point(3)]
  temp(:,9 ) = [point(2), point(3), point(1), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,10) = [point(2), point(3), 1.0_pReal - point(1) - point(2)- point(3), point(1)]
  temp(:,11) = [point(2), 1.0_pReal - point(1) - point(2)- point(3), point(1), point(3)]
  temp(:,12) = [point(2), 1.0_pReal - point(1) - point(2)- point(3), point(3), point(1)]
  temp(:,13) = [point(3), point(1), point(2), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,14) = [point(3), point(1), 1.0_pReal - point(1) - point(2)- point(3), point(2)]
  temp(:,15) = [point(3), point(2), point(1), 1.0_pReal - point(1) - point(2)- point(3)]
  temp(:,16) = [point(3), point(2), 1.0_pReal - point(1) - point(2)- point(3), point(1)]
  temp(:,17) = [point(3), 1.0_pReal - point(1) - point(2)- point(3), point(1), point(2)]
  temp(:,18) = [point(3), 1.0_pReal - point(1) - point(2)- point(3), point(2), point(1)]
  temp(:,19) = [1.0_pReal - point(1) - point(2)- point(3), point(1), point(2), point(3)]
  temp(:,20) = [1.0_pReal - point(1) - point(2)- point(3), point(1), point(3), point(2)]
  temp(:,21) = [1.0_pReal - point(1) - point(2)- point(3), point(2), point(1), point(3)]
  temp(:,22) = [1.0_pReal - point(1) - point(2)- point(3), point(2), point(3), point(1)]
  temp(:,23) = [1.0_pReal - point(1) - point(2)- point(3), point(3), point(1), point(2)]
  temp(:,24) = [1.0_pReal - point(1) - point(2)- point(3), point(3), point(2), point(1)]

  qPt = reshape(matmul(tetrahedron, temp),[72])

end function permutationStar1111

end module FEM_quadrature