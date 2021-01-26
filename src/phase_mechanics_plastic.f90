submodule(constitutive:constitutive_mech) plastic

  interface

    module subroutine isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        instance, &
        me
    end subroutine isotropic_LpAndItsTangent

    pure module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        instance, &
        me
    end subroutine phenopowerlaw_LpAndItsTangent

    pure module subroutine kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        instance, &
        me
    end subroutine kinehardening_LpAndItsTangent

    module subroutine dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,T,instance,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        instance, &
        me
    end subroutine dislotwin_LpAndItsTangent

    pure module subroutine dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp,T,instance,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        instance, &
        me
    end subroutine dislotungsten_LpAndItsTangent

    module subroutine nonlocal_LpAndItsTangent(Lp,dLp_dMp, &
                                                       Mp,Temperature,instance,me,ip,el)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                     intent(in) :: &
        Temperature
      integer,                         intent(in) :: &
        instance, &
        me, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine nonlocal_LpAndItsTangent

  end interface

contains

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: Discuss whether it makes sense if crystallite handles the configuration conversion, i.e.
! Mp in, dLp_dMp out
!--------------------------------------------------------------------------------------------------
module subroutine plastic_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                     S, Fi, co, ip, el)
  integer, intent(in) :: &
    co, &                                                                                           !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal),   intent(in),  dimension(3,3) :: &
    S, &                                                                                            !< 2nd Piola-Kirchhoff stress
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),   intent(out), dimension(3,3) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLp_dS, &
    dLp_dFi                                                                                         !< derivative me Lp with respect to Fi

  real(pReal), dimension(3,3,3,3) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to Mandel stress
  real(pReal), dimension(3,3) :: &
    Mp                                                                                              !< Mandel stress work conjugate with Lp
  integer :: &
    i, j, instance, me, ph


  Mp = matmul(matmul(transpose(Fi),Fi),S)
  me = material_phasememberAt(co,ip,el)
  ph = material_phaseAt(co,el)
  instance = phase_plasticityInstance(ph)

  plasticityType: select case (phase_plasticity(material_phaseAt(co,el)))

    case (PLASTICITY_NONE_ID) plasticityType
      Lp = 0.0_pReal
      dLp_dMp = 0.0_pReal

    case (PLASTICITY_ISOTROPIC_ID) plasticityType
      call isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
      call phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,me)

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call nonlocal_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),instance,me,ip,el)

    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),instance,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticityType
      call dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),instance,me)

  end select plasticityType

  do i=1,3; do j=1,3
    dLp_dFi(i,j,1:3,1:3) = matmul(matmul(Fi,S),transpose(dLp_dMp(i,j,1:3,1:3))) + &
                           matmul(matmul(Fi,dLp_dMp(i,j,1:3,1:3)),S)
    dLp_dS(i,j,1:3,1:3)  = matmul(matmul(transpose(Fi),Fi),dLp_dMp(i,j,1:3,1:3))                     ! ToDo: @PS: why not:   dLp_dMp:(FiT Fi)
  enddo; enddo

end subroutine plastic_LpAndItsTangents

end submodule plastic
