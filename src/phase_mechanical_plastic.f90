submodule(phase:mechanical) plastic

  interface

    module function plastic_none_init()          result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_none_init

    module function plastic_isotropic_init()     result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_isotropic_init

    module function plastic_phenopowerlaw_init() result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_phenopowerlaw_init

    module function plastic_kinehardening_init() result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_kinehardening_init

    module function plastic_dislotwin_init()     result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_dislotwin_init

    module function plastic_dislotungsten_init() result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_dislotungsten_init

    module function plastic_nonlocal_init()      result(myPlasticity)
      logical, dimension(:), allocatable :: &
        myPlasticity
    end function plastic_nonlocal_init

    module subroutine isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        me
    end subroutine isotropic_LpAndItsTangent

    pure module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        me
    end subroutine phenopowerlaw_LpAndItsTangent

    pure module subroutine kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        me
    end subroutine kinehardening_LpAndItsTangent

    module subroutine dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,T,ph,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        ph, &
        me
    end subroutine dislotwin_LpAndItsTangent

    pure module subroutine dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp,T,ph,me)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        ph, &
        me
    end subroutine dislotungsten_LpAndItsTangent

    module subroutine nonlocal_LpAndItsTangent(Lp,dLp_dMp, &
                                                       Mp,Temperature,ph,me,ip,el)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                     intent(in) :: &
        Temperature
      integer,                         intent(in) :: &
        ph, &
        me, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine nonlocal_LpAndItsTangent


    module subroutine isotropic_dotState(Mp,ph,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        me
    end subroutine isotropic_dotState

    module subroutine phenopowerlaw_dotState(Mp,ph,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        me
    end subroutine phenopowerlaw_dotState

    module subroutine plastic_kinehardening_dotState(Mp,ph,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        me
    end subroutine plastic_kinehardening_dotState

    module subroutine dislotwin_dotState(Mp,T,ph,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        ph, &
        me
    end subroutine dislotwin_dotState

    module subroutine dislotungsten_dotState(Mp,T,ph,me)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        ph, &
        me
    end subroutine dislotungsten_dotState

    module subroutine nonlocal_dotState(Mp,Temperature,timestep,ph,me,ip,el)
      real(pReal), dimension(3,3), intent(in) :: &
        Mp                                                                                          !< MandelStress
      real(pReal), intent(in) :: &
        Temperature, &                                                                              !< temperature
        timestep                                                                                    !< substepped crystallite time increment
      integer, intent(in) :: &
        ph, &
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

module subroutine plastic_init


  print'(/,a)', ' <<<+-  phase:mechanical:plastic init  -+>>>'

  where(plastic_none_init())              phase_plasticity = PLASTICITY_NONE_ID
  where(plastic_isotropic_init())         phase_plasticity = PLASTICITY_ISOTROPIC_ID
  where(plastic_phenopowerlaw_init())     phase_plasticity = PLASTICITY_PHENOPOWERLAW_ID
  where(plastic_kinehardening_init())     phase_plasticity = PLASTICITY_KINEHARDENING_ID
  where(plastic_dislotwin_init())         phase_plasticity = PLASTICITY_DISLOTWIN_ID
  where(plastic_dislotungsten_init())     phase_plasticity = PLASTICITY_DISLOTUNGSTEN_ID
  where(plastic_nonlocal_init())          phase_plasticity = PLASTICITY_NONLOCAL_ID

end subroutine plastic_init

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
    i, j, me, ph


  Mp = matmul(matmul(transpose(Fi),Fi),S)
  me = material_phasememberAt(co,ip,el)
  ph = material_phaseAt(co,el)

  plasticType: select case (phase_plasticity(material_phaseAt(co,el)))

    case (PLASTICITY_NONE_ID) plasticType
      Lp = 0.0_pReal
      dLp_dMp = 0.0_pReal

    case (PLASTICITY_ISOTROPIC_ID) plasticType
      call isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticType
      call phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)

    case (PLASTICITY_KINEHARDENING_ID) plasticType
      call kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,me)

    case (PLASTICITY_NONLOCAL_ID) plasticType
      call nonlocal_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),ph,me,ip,el)

    case (PLASTICITY_DISLOTWIN_ID) plasticType
      call dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),ph,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticType
      call dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp, thermal_T(ph,me),ph,me)

  end select plasticType

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
  logical :: broken


  Mp = matmul(matmul(transpose(phase_mechanical_Fi(ph)%data(1:3,1:3,me)),&
                     phase_mechanical_Fi(ph)%data(1:3,1:3,me)),phase_mechanical_S(ph)%data(1:3,1:3,me))

  plasticType: select case (phase_plasticity(ph))

    case (PLASTICITY_ISOTROPIC_ID) plasticType
      call isotropic_dotState(Mp,ph,me)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticType
      call phenopowerlaw_dotState(Mp,ph,me)

    case (PLASTICITY_KINEHARDENING_ID) plasticType
      call plastic_kinehardening_dotState(Mp,ph,me)

    case (PLASTICITY_DISLOTWIN_ID) plasticType
      call dislotwin_dotState(Mp,thermal_T(ph,me),ph,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticType
      call dislotungsten_dotState(Mp,thermal_T(ph,me),ph,me)

    case (PLASTICITY_NONLOCAL_ID) plasticType
      call nonlocal_dotState(Mp,thermal_T(ph,me),subdt,ph,me,ip,el)
  end select plasticType
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
  instance = phase_plasticInstance(ph)

  plasticType: select case (phase_plasticity(material_phaseAt(co,el)))

    case (PLASTICITY_DISLOTWIN_ID) plasticType
      call dislotwin_dependentState(thermal_T(ph,me),instance,me)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticType
      call dislotungsten_dependentState(instance,me)

    case (PLASTICITY_NONLOCAL_ID) plasticType
      call nonlocal_dependentState(instance,me,ip,el)

  end select plasticType

end subroutine plastic_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
module function plastic_deltaState(co, ip, el, ph, me) result(broken)

  integer, intent(in) :: &
    co, &                                                                                           !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    me
  logical :: &
    broken

  real(pReal),               dimension(3,3) :: &
    Mp
  integer :: &
    instance, &
    myOffset, &
    mySize


  Mp = matmul(matmul(transpose(phase_mechanical_Fi(ph)%data(1:3,1:3,me)),&
                     phase_mechanical_Fi(ph)%data(1:3,1:3,me)),phase_mechanical_S(ph)%data(1:3,1:3,me))
  instance = phase_plasticInstance(ph)

  plasticType: select case (phase_plasticity(ph))

    case (PLASTICITY_KINEHARDENING_ID) plasticType
      call plastic_kinehardening_deltaState(Mp,instance,me)
      broken = any(IEEE_is_NaN(plasticState(ph)%deltaState(:,me)))

    case (PLASTICITY_NONLOCAL_ID) plasticType
      call plastic_nonlocal_deltaState(Mp,instance,me,ip,el)
      broken = any(IEEE_is_NaN(plasticState(ph)%deltaState(:,me)))

    case default
      broken = .false.

  end select plasticType

  if(.not. broken) then
    select case(phase_plasticity(ph))
      case (PLASTICITY_NONLOCAL_ID,PLASTICITY_KINEHARDENING_ID)

        myOffset = plasticState(ph)%offsetDeltaState
        mySize   = plasticState(ph)%sizeDeltaState
        plasticState(ph)%state(myOffset + 1:myOffset + mySize,me) = &
        plasticState(ph)%state(myOffset + 1:myOffset + mySize,me) + plasticState(ph)%deltaState(1:mySize,me)
    end select
  endif

end function plastic_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief checks if a plastic module is active or not
!--------------------------------------------------------------------------------------------------
function plastic_active(plastic_label)  result(active_plastic)

  character(len=*), intent(in)       :: plastic_label                                               !< type of plasticity model
  logical, dimension(:), allocatable :: active_plastic

  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl
  integer :: ph

  phases => config_material%get('phase')
  allocate(active_plastic(phases%length), source = .false. )
  do ph = 1, phases%length
    phase => phases%get(ph)
    mech  => phase%get('mechanics')
    pl    => mech%get('plasticity')
    if(pl%get_asString('type') == plastic_label) active_plastic(ph) = .true.
  enddo

end function plastic_active

end submodule plastic
