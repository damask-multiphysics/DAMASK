submodule(phase:mechanical) eigen

  integer, dimension(:), allocatable :: &
    Nmodels

  integer(kind(EIGEN_UNDEFINED_ID)),  dimension(:,:), allocatable :: &
    model
  integer(kind(EIGEN_UNDEFINED_ID)),  dimension(:), allocatable :: &
    model_damage

  interface
    module function damage_anisobrittle_init() result(myKinematics)
      logical, dimension(:), allocatable :: myKinematics
    end function damage_anisobrittle_init

    module function thermalexpansion_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function thermalexpansion_init

    module subroutine thermalexpansion_LiAndItsTangent(Li, dLi_dTstar, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(out), dimension(3,3) :: &
        Li                                                                                          !< thermal velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLi_dTstar                                                                                  !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
    end subroutine thermalexpansion_LiAndItsTangent

  end interface


contains


module subroutine eigen_init(phases)

  type(tDict), pointer :: &
    phases

  integer :: &
    ph
  type(tDict), pointer :: &
    phase, &
    mechanics
  type(tList), pointer :: &
    kinematics

  print'(/,1x,a)', '<<<+-  phase:mechanical:eigen init  -+>>>'

!--------------------------------------------------------------------------------------------------
! explicit eigen mechanisms
  allocate(Nmodels(phases%length),source = 0)

  do ph = 1,phases%length
    phase => phases%get_dict(ph)
    mechanics => phase%get_dict('mechanical')
    kinematics => mechanics%get_list('eigen',defaultVal=emptyList)
    Nmodels(ph) = kinematics%length
  end do

  allocate(model(maxval(Nmodels),phases%length), source = EIGEN_undefined_ID)

  if (maxval(Nmodels) /= 0) then
    where(thermalexpansion_init(maxval(Nmodels))) model = EIGEN_thermal_expansion_ID
  end if

  allocate(model_damage(phases%length),  source = EIGEN_UNDEFINED_ID)

  where(damage_anisobrittle_init())  model_damage = EIGEN_cleavage_opening_ID


end subroutine eigen_init


!--------------------------------------------------------------------------------------------------
!> @brief checks if a kinematic mechanism is active or not
!--------------------------------------------------------------------------------------------------
function kinematics_active(kinematics_label,kinematics_length)  result(active_kinematics)

  character(len=*), intent(in)         :: kinematics_label                                          !< name of kinematic mechanism
  integer,          intent(in)         :: kinematics_length                                         !< max. number of kinematics in system
  logical, dimension(:,:), allocatable :: active_kinematics

  type(tDict), pointer :: &
    phases, &
    phase, &
    mechanics, &
    kinematic
  type(tList), pointer :: &
    kinematics
  integer :: ph,k


  phases => config_material%get_dict('phase')
  allocate(active_kinematics(kinematics_length,phases%length), source = .false. )
  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    mechanics => phase%get_dict('mechanical')
    kinematics => mechanics%get_list('eigen',defaultVal=emptyList)
    do k = 1, kinematics%length
      kinematic => kinematics%get_dict(k)
      active_kinematics(k,ph) = kinematic%get_asString('type') == kinematics_label
    end do
  end do

end function kinematics_active



!--------------------------------------------------------------------------------------------------
!> @brief checks if a kinematic mechanism is active or not
!--------------------------------------------------------------------------------------------------
function kinematics_active2(kinematics_label)  result(active_kinematics)

  character(len=*), intent(in)       :: kinematics_label                                            !< name of kinematic mechanism
  logical, dimension(:), allocatable :: active_kinematics

  type(tDict), pointer :: &
    phases, &
    phase, &
    kinematics_type
  type(tList), pointer :: &
    kinematics
  integer :: ph

  phases => config_material%get_dict('phase')
  allocate(active_kinematics(phases%length), source = .false.)
  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    kinematics => phase%get_list('damage',defaultVal=emptyList)
    if (kinematics%length < 1) return
    kinematics_type => kinematics%get_dict(1)
    if (.not. kinematics_type%contains('type')) continue
    active_kinematics(ph) = kinematics_type%get_asString('type',defaultVal='n/a') == kinematics_label
  end do


end function kinematics_active2


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
module subroutine phase_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                         S, Fi, ph,en)

  integer, intent(in) :: &
    ph,en
  real(pReal),   intent(in),  dimension(3,3) :: &
    S                                                                                               !< 2nd Piola-Kirchhoff stress
  real(pReal),   intent(in),  dimension(3,3) :: &
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),   intent(out), dimension(3,3) :: &
    Li                                                                                              !< intermediate velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLi_dS, &                                                                                       !< derivative of Li with respect to S
    dLi_dFi

  real(pReal), dimension(3,3) :: &
    my_Li, &                                                                                        !< intermediate velocity gradient
    FiInv, &
    temp_33
  real(pReal), dimension(3,3,3,3) :: &
    my_dLi_dS
  real(pReal) :: &
    detFi
  integer :: &
    k, i, j
  logical :: active

  active = .false.
  Li = 0.0_pReal
  dLi_dS  = 0.0_pReal
  dLi_dFi = 0.0_pReal


  plasticType: select case (phase_plasticity(ph))
    case (PLASTIC_isotropic_ID) plasticType
      call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, S ,ph,en)
      Li = Li + my_Li
      dLi_dS = dLi_dS + my_dLi_dS
      active = .true.
  end select plasticType

  KinematicsLoop: do k = 1, Nmodels(ph)
    kinematicsType: select case (model(k,ph))
      case (EIGEN_thermal_expansion_ID) kinematicsType
        call thermalexpansion_LiAndItsTangent(my_Li, my_dLi_dS, ph,en)
        Li = Li + my_Li
        dLi_dS = dLi_dS + my_dLi_dS
        active = .true.
    end select kinematicsType
  end do KinematicsLoop

  select case (model_damage(ph))
    case (EIGEN_cleavage_opening_ID)
      call damage_anisobrittle_LiAndItsTangent(my_Li, my_dLi_dS, S, ph, en)
      Li = Li + my_Li
      dLi_dS = dLi_dS + my_dLi_dS
      active = .true.
  end select

  if (.not. active) return

  FiInv = math_inv33(Fi)
  detFi = math_det33(Fi)
  Li = matmul(matmul(Fi,Li),FiInv)*detFi                                                            !< push forward to intermediate configuration
  temp_33 = matmul(FiInv,Li)

  do i = 1,3; do j = 1,3
    dLi_dS(1:3,1:3,i,j)  = matmul(matmul(Fi,dLi_dS(1:3,1:3,i,j)),FiInv)*detFi
    dLi_dFi(1:3,1:3,i,j) = dLi_dFi(1:3,1:3,i,j) + Li*FiInv(j,i)
    dLi_dFi(1:3,i,1:3,j) = dLi_dFi(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
  end do; end do

end subroutine phase_LiAndItsTangents


end submodule eigen
