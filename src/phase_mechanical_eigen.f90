submodule(phase:mechanics) eigendeformation

  interface
    module function kinematics_cleavage_opening_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_cleavage_opening_init

    module function kinematics_slipplane_opening_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_slipplane_opening_init

    module function kinematics_thermal_expansion_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_thermal_expansion_init

    module subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_cleavage_opening_LiAndItsTangent

    module subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_slipplane_opening_LiAndItsTangent

    module subroutine thermalexpansion_LiAndItsTangent(Li, dLi_dTstar, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(out), dimension(3,3) :: &
        Li                                                                                          !< thermal velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLi_dTstar                                                                                  !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
    end subroutine thermalexpansion_LiAndItsTangent

  end interface


contains


module subroutine eigendeformation_init(phases)

  class(tNode), pointer :: &
    phases

  integer :: &
    ph
  class(tNode), pointer :: &
    phase, &
    kinematics, &
    damage

  print'(/,a)', ' <<<+-  phase:mechanics:eigendeformation init  -+>>>'

!--------------------------------------------------------------------------------------------------
! initialize kinematic mechanisms
  allocate(phase_Nkinematics(phases%length),source = 0)
  do ph = 1,phases%length
    phase => phases%get(ph)
    kinematics => phase%get('kinematics',defaultVal=emptyList)
    phase_Nkinematics(ph) = kinematics%length
    kinematics => phase%get('damage',defaultVal=emptyList)
    if(kinematics%length >0) then
      damage => kinematics%get(1)
      if(damage%get_asString('type',defaultVal='n/a') == 'anisobrittle')  phase_Nkinematics(ph) =  phase_Nkinematics(ph) +1
      if(damage%get_asString('type',defaultVal='n/a') == 'isoductile'  )  phase_Nkinematics(ph) =  phase_Nkinematics(ph) +1
    endif
  enddo

  allocate(phase_kinematics(maxval(phase_Nkinematics),phases%length), source = KINEMATICS_undefined_ID)

  if(maxval(phase_Nkinematics) /= 0) then
    where(kinematics_cleavage_opening_init(maxval(phase_Nkinematics)))  phase_kinematics = KINEMATICS_cleavage_opening_ID
    where(kinematics_slipplane_opening_init(maxval(phase_Nkinematics))) phase_kinematics = KINEMATICS_slipplane_opening_ID
    where(kinematics_thermal_expansion_init(maxval(phase_Nkinematics))) phase_kinematics = KINEMATICS_thermal_expansion_ID
  endif

end subroutine eigendeformation_init


!--------------------------------------------------------------------------------------------------
!> @brief checks if a kinematic mechanism is active or not
!--------------------------------------------------------------------------------------------------
function kinematics_active(kinematics_label,kinematics_length)  result(active_kinematics)

  character(len=*), intent(in)         :: kinematics_label                                          !< name of kinematic mechanism
  integer,          intent(in)         :: kinematics_length                                         !< max. number of kinematics in system
  logical, dimension(:,:), allocatable :: active_kinematics

  class(tNode), pointer :: &
    phases, &
    phase, &
    kinematics, &
    kinematics_type
  integer :: p,k

  phases => config_material%get('phase')
  allocate(active_kinematics(kinematics_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    kinematics => phase%get('kinematics',defaultVal=emptyList)
    do k = 1, kinematics%length
      kinematics_type => kinematics%get(k)
      active_kinematics(k,p) = kinematics_type%get_asString('type') == kinematics_label
    enddo
  enddo


end function kinematics_active



!--------------------------------------------------------------------------------------------------
!> @brief checks if a kinematic mechanism is active or not
!--------------------------------------------------------------------------------------------------
function kinematics_active2(kinematics_label,kinematics_length)  result(active_kinematics)

  character(len=*), intent(in)         :: kinematics_label                                          !< name of kinematic mechanism
  integer,          intent(in)         :: kinematics_length                                         !< max. number of kinematics in system
  logical, dimension(:,:), allocatable :: active_kinematics

  class(tNode), pointer :: &
    phases, &
    phase, &
    kinematics, &
    kinematics_type
  integer :: p

  phases => config_material%get('phase')
  allocate(active_kinematics(kinematics_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    kinematics => phase%get('damage',defaultVal=emptyList)
    kinematics_type => kinematics%get(1)
    if (.not. kinematics_type%contains('type')) continue
    active_kinematics(1,p) = kinematics_type%get_asString('type',defaultVal='n/a') == kinematics_label
  enddo


end function kinematics_active2


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
module subroutine phase_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                         S, Fi, ph,me)

  integer, intent(in) :: &
    ph,me
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

  Li = 0.0_pReal
  dLi_dS  = 0.0_pReal
  dLi_dFi = 0.0_pReal


  plasticType: select case (phase_plasticity(ph))
    case (PLASTICITY_isotropic_ID) plasticType
      call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, S ,phase_plasticInstance(ph),me)
    case default plasticType
      my_Li = 0.0_pReal
      my_dLi_dS = 0.0_pReal
  end select plasticType

  Li = Li + my_Li
  dLi_dS = dLi_dS + my_dLi_dS

  KinematicsLoop: do k = 1, phase_Nkinematics(ph)
    kinematicsType: select case (phase_kinematics(k,ph))
      case (KINEMATICS_cleavage_opening_ID) kinematicsType
        call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ph, me)
      case (KINEMATICS_slipplane_opening_ID) kinematicsType
        call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ph, me)
      case (KINEMATICS_thermal_expansion_ID) kinematicsType
        call thermalexpansion_LiAndItsTangent(my_Li, my_dLi_dS, ph,me)
      case default kinematicsType
        my_Li = 0.0_pReal
        my_dLi_dS = 0.0_pReal
    end select kinematicsType
    Li = Li + my_Li
    dLi_dS = dLi_dS + my_dLi_dS
  enddo KinematicsLoop

  FiInv = math_inv33(Fi)
  detFi = math_det33(Fi)
  Li = matmul(matmul(Fi,Li),FiInv)*detFi                                                            !< push forward to intermediate configuration
  temp_33 = matmul(FiInv,Li)

  do i = 1,3; do j = 1,3
    dLi_dS(1:3,1:3,i,j)  = matmul(matmul(Fi,dLi_dS(1:3,1:3,i,j)),FiInv)*detFi
    dLi_dFi(1:3,1:3,i,j) = dLi_dFi(1:3,1:3,i,j) + Li*FiInv(j,i)
    dLi_dFi(1:3,i,1:3,j) = dLi_dFi(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
  enddo; enddo

end subroutine phase_LiAndItsTangents


end submodule eigendeformation
