submodule(constitutive:mechanics) plastic

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


        module subroutine isotropic_dotState(Mp,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine isotropic_dotState

    module subroutine phenopowerlaw_dotState(Mp,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine phenopowerlaw_dotState

    module subroutine plastic_kinehardening_dotState(Mp,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine plastic_kinehardening_dotState

    module subroutine dislotwin_dotState(Mp,T,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine dislotwin_dotState

    module subroutine dislotungsten_dotState(Mp,T,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine dislotungsten_dotState

    module subroutine nonlocal_dotState(Mp,Temperature,timestep,instance,me,ip,el)
      real(pReal), dimension(3,3), intent(in) :: &
        Mp                                                                                          !< MandelStress
      real(pReal), intent(in) :: &
        Temperature, &                                                                              !< temperature
        timestep                                                                                    !< substepped crystallite time increment
      integer, intent(in) :: &
        instance, &
        me, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine nonlocal_dotState

        module subroutine dislotwin_dependentState(T,instance,me)
      integer,       intent(in) :: &
        instance, &
        me
      real(pReal),   intent(in) :: &
        T
    end subroutine dislotwin_dependentState

    module subroutine dislotungsten_dependentState(instance,me)
      integer,       intent(in) :: &
        instance, &
        me
    end subroutine dislotungsten_dependentState

    module subroutine nonlocal_dependentState(instance, me, ip, el)
      integer, intent(in) :: &
        instance, &
        me, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine nonlocal_dependentState

    module subroutine plastic_kinehardening_deltaState(Mp,instance,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        me
    end subroutine plastic_kinehardening_deltaState

    module subroutine plastic_nonlocal_deltaState(Mp,instance,me,ip,el)
      real(pReal), dimension(3,3), intent(in) :: &
        Mp
      integer, intent(in) :: &
        instance, &
        me, &
        ip, &
        el
    end subroutine plastic_nonlocal_deltaState

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


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
module function plastic_dotState(subdt,co,ip,el,ph,me) result(broken)

  integer, intent(in) :: &
    co, &                                                                                           !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    me
  real(pReal),  intent(in) :: &
    subdt                                                                                           !< timestep
  real(pReal),              dimension(3,3) :: &
    Mp
  integer :: &
    instance
  logical :: broken


  instance = phase_plasticityInstance(ph)

  Mp = matmul(matmul(transpose(constitutive_mech_Fi(ph)%data(1:3,1:3,me)),&
                     constitutive_mech_Fi(ph)%data(1:3,1:3,me)),constitutive_mech_S(ph)%data(1:3,1:3,me))

  plasticityType: select case (phase_plasticity(ph))

    case (PLASTICITY_ISOTROPIC_ID) plasticityType
      call isotropic_dotState(Mp,instance,me)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
      call phenopowerlaw_dotState(Mp,instance,me)

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call plastic_kinehardening_dotState(Mp,instance,me)

    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call dislotwin_dotState(Mp,thermal_T(ph,me),instance,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticityType
      call dislotungsten_dotState(Mp,thermal_T(ph,me),instance,me)

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call nonlocal_dotState(Mp,thermal_T(ph,me),subdt,instance,me,ip,el)
  end select plasticityType
  broken = any(IEEE_is_NaN(plasticState(ph)%dotState(:,me)))


end function plastic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different plasticity constitutive models
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dependentState(co, ip, el)

  integer, intent(in) :: &
    co, &                                                                                           !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

  integer :: &
    ph, &
    instance, me


  ph = material_phaseAt(co,el)
  me = material_phasememberAt(co,ip,el)
  instance = phase_plasticityInstance(ph)

  plasticityType: select case (phase_plasticity(material_phaseAt(co,el)))

    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call dislotwin_dependentState(thermal_T(ph,me),instance,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticityType
      call dislotungsten_dependentState(instance,me)

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call nonlocal_dependentState(instance,me,ip,el)

  end select plasticityType

end subroutine plastic_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
module function plastic_deltaState(co, ip, el, ph, of) result(broken)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    of
  logical :: &
    broken

  real(pReal),               dimension(3,3) :: &
    Mp
  integer :: &
    instance, &
    myOffset, &
    mySize


  Mp = matmul(matmul(transpose(constitutive_mech_Fi(ph)%data(1:3,1:3,of)),&
                     constitutive_mech_Fi(ph)%data(1:3,1:3,of)),constitutive_mech_S(ph)%data(1:3,1:3,of))
  instance = phase_plasticityInstance(ph)

  plasticityType: select case (phase_plasticity(ph))

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call plastic_kinehardening_deltaState(Mp,instance,of)
      broken = any(IEEE_is_NaN(plasticState(ph)%deltaState(:,of)))

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call plastic_nonlocal_deltaState(Mp,instance,of,ip,el)
      broken = any(IEEE_is_NaN(plasticState(ph)%deltaState(:,of)))

    case default
      broken = .false.

  end select plasticityType

  if(.not. broken) then
    select case(phase_plasticity(ph))
      case (PLASTICITY_NONLOCAL_ID,PLASTICITY_KINEHARDENING_ID)

        myOffset = plasticState(ph)%offsetDeltaState
        mySize   = plasticState(ph)%sizeDeltaState
        plasticState(ph)%state(myOffset + 1:myOffset + mySize,of) = &
        plasticState(ph)%state(myOffset + 1:myOffset + mySize,of) + plasticState(ph)%deltaState(1:mySize,of)
    end select
  endif

end function plastic_deltaState

end submodule plastic
