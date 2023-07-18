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

    module subroutine isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine isotropic_LpAndItsTangent

    pure module subroutine phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine phenopowerlaw_LpAndItsTangent

    pure module subroutine kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine kinehardening_LpAndItsTangent

    module subroutine dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine dislotwin_LpAndItsTangent

    pure module subroutine dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine dislotungsten_LpAndItsTangent

    module subroutine nonlocal_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)
      real(pREAL), dimension(3,3),     intent(out) :: &
        Lp
      real(pREAL), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp
      real(pREAL), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        ph, &
        en
    end subroutine nonlocal_LpAndItsTangent


    module function isotropic_dotState(Mp,ph,en) result(dotState)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
        dotState
    end function isotropic_dotState

    module function phenopowerlaw_dotState(Mp,ph,en) result(dotState)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
        dotState
    end function phenopowerlaw_dotState

    module function plastic_kinehardening_dotState(Mp,ph,en) result(dotState)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
        dotState
    end function plastic_kinehardening_dotState

    module function dislotwin_dotState(Mp,ph,en) result(dotState)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
        dotState
    end function dislotwin_dotState

    module function dislotungsten_dotState(Mp,ph,en) result(dotState)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
      real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
        dotState
    end function dislotungsten_dotState

    module subroutine nonlocal_dotState(Mp,timestep,ph,en)
      real(pREAL), dimension(3,3), intent(in) :: &
        Mp                                                                                          !< MandelStress
      real(pREAL), intent(in) :: &
        timestep                                                                                    !< substepped crystallite time increment
      integer, intent(in) :: &
        ph, &
        en
    end subroutine nonlocal_dotState

    module subroutine dislotwin_dependentState(ph,en)
      integer,       intent(in) :: &
        ph, &
        en
    end subroutine dislotwin_dependentState

    module subroutine dislotungsten_dependentState(ph,en)
      integer,       intent(in) :: &
        ph, &
        en
    end subroutine dislotungsten_dependentState

    module subroutine nonlocal_dependentState(ph,en)
      integer, intent(in) :: &
        ph, &
        en
    end subroutine nonlocal_dependentState

    module subroutine plastic_kinehardening_deltaState(Mp,ph,en)
      real(pREAL), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        ph, &
        en
    end subroutine plastic_kinehardening_deltaState

    module subroutine plastic_nonlocal_deltaState(Mp,ph,en)
      real(pREAL), dimension(3,3), intent(in) :: &
        Mp
      integer, intent(in) :: &
        ph, &
        en
    end subroutine plastic_nonlocal_deltaState

  end interface

contains

module subroutine plastic_init


  print'(/,1x,a)', '<<<+-  phase:mechanical:plasticity init  -+>>>'

  where(plastic_none_init())              mechanical_plasticity_type = MECHANICAL_PLASTICITY_NONE
  where(plastic_isotropic_init())         mechanical_plasticity_type = MECHANICAL_PLASTICITY_ISOTROPIC
  where(plastic_phenopowerlaw_init())     mechanical_plasticity_type = MECHANICAL_PLASTICITY_PHENOPOWERLAW
  where(plastic_kinehardening_init())     mechanical_plasticity_type = MECHANICAL_PLASTICITY_KINEHARDENING
  where(plastic_dislotwin_init())         mechanical_plasticity_type = MECHANICAL_PLASTICITY_DISLOTWIN
  where(plastic_dislotungsten_init())     mechanical_plasticity_type = MECHANICAL_PLASTICITY_DISLOTUNGSTEN
  where(plastic_nonlocal_init())          mechanical_plasticity_type = MECHANICAL_PLASTICITY_NONLOCAL

  if (any(mechanical_plasticity_type == UNDEFINED)) call IO_error(201)

end subroutine plastic_init

!--------------------------------------------------------------------------------------------------
!> @brief constitutive equation for calculating the velocity gradient
! ToDo: Discuss whether it makes sense if crystallite handles the configuration conversion, i.e.
! Mp in, dLp_dMp out
!--------------------------------------------------------------------------------------------------
module subroutine plastic_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                           S, Fi, ph,en)
  integer, intent(in) :: &
    ph,en
  real(pREAL),   intent(in),  dimension(3,3) :: &
    S, &                                                                                            !< 2nd Piola-Kirchhoff stress
    Fi                                                                                              !< intermediate deformation gradient
  real(pREAL),   intent(out), dimension(3,3) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL),   intent(out), dimension(3,3,3,3) :: &
    dLp_dS, &
    dLp_dFi                                                                                         !< derivative en Lp with respect to Fi

  real(pREAL), dimension(3,3,3,3) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to Mandel stress
  real(pREAL), dimension(3,3) :: &
    Mp                                                                                              !< Mandel stress work conjugate with Lp
  integer :: &
    i, j


  if (mechanical_plasticity_type(ph) == MECHANICAL_PLASTICITY_NONE) then
    Lp      = 0.0_pREAL
    dLp_dFi = 0.0_pREAL
    dLp_dS  = 0.0_pREAL
  else

    Mp = matmul(matmul(transpose(Fi),Fi),S)

    plasticType: select case (mechanical_plasticity_type(ph))

      case (MECHANICAL_PLASTICITY_ISOTROPIC) plasticType
        call isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

      case (MECHANICAL_PLASTICITY_PHENOPOWERLAW) plasticType
        call phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

      case (MECHANICAL_PLASTICITY_KINEHARDENING) plasticType
        call kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

      case (MECHANICAL_PLASTICITY_NONLOCAL) plasticType
        call nonlocal_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

      case (MECHANICAL_PLASTICITY_DISLOTWIN) plasticType
        call dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

      case (MECHANICAL_PLASTICITY_DISLOTUNGSTEN) plasticType
        call dislotungsten_LpAndItsTangent(Lp,dLp_dMp,Mp,ph,en)

    end select plasticType

    do i=1,3; do j=1,3
      dLp_dFi(i,j,1:3,1:3) = matmul(matmul(Fi,S),transpose(dLp_dMp(i,j,1:3,1:3))) + &
                             matmul(matmul(Fi,dLp_dMp(i,j,1:3,1:3)),S)
      dLp_dS(i,j,1:3,1:3)  = matmul(matmul(transpose(Fi),Fi),dLp_dMp(i,j,1:3,1:3))                  ! ToDo: @PS: why not:   dLp_dMp:(FiT Fi)
    end do; end do

  end if

end subroutine plastic_LpAndItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
module function plastic_dotState(subdt,ph,en) result(dotState)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL),  intent(in) :: &
    subdt                                                                                           !< timestep
  real(pREAL),              dimension(3,3) :: &
    Mp
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState


  if (mechanical_plasticity_type(ph) /= MECHANICAL_PLASTICITY_NONE) then
    Mp = matmul(matmul(transpose(phase_mechanical_Fi(ph)%data(1:3,1:3,en)),&
                       phase_mechanical_Fi(ph)%data(1:3,1:3,en)),phase_mechanical_S(ph)%data(1:3,1:3,en))

    plasticType: select case (mechanical_plasticity_type(ph))

      case (MECHANICAL_PLASTICITY_ISOTROPIC) plasticType
        dotState = isotropic_dotState(Mp,ph,en)

      case (MECHANICAL_PLASTICITY_PHENOPOWERLAW) plasticType
        dotState = phenopowerlaw_dotState(Mp,ph,en)

      case (MECHANICAL_PLASTICITY_KINEHARDENING) plasticType
        dotState = plastic_kinehardening_dotState(Mp,ph,en)

      case (MECHANICAL_PLASTICITY_DISLOTWIN) plasticType
        dotState = dislotwin_dotState(Mp,ph,en)

      case (MECHANICAL_PLASTICITY_DISLOTUNGSTEN) plasticType
        dotState = dislotungsten_dotState(Mp,ph,en)

      case (MECHANICAL_PLASTICITY_NONLOCAL) plasticType
        call nonlocal_dotState(Mp,subdt,ph,en)
        dotState = plasticState(ph)%dotState(:,en)

    end select plasticType
  end if

end function plastic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different plasticity constitutive models
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dependentState(ph,en)

  integer, intent(in) :: &
    ph, &
    en


  plasticType: select case (mechanical_plasticity_type(ph))

    case (MECHANICAL_PLASTICITY_DISLOTWIN) plasticType
      call dislotwin_dependentState(ph,en)

    case (MECHANICAL_PLASTICITY_DISLOTUNGSTEN) plasticType
      call dislotungsten_dependentState(ph,en)

    case (MECHANICAL_PLASTICITY_NONLOCAL) plasticType
      call nonlocal_dependentState(ph,en)

  end select plasticType

end subroutine plastic_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models that have an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
module function plastic_deltaState(ph, en) result(broken)

  integer, intent(in) :: &
    ph, &
    en
  logical :: broken

  real(pREAL), dimension(3,3) :: &
    Mp
  integer :: &
    mySize


  broken = .false.

  select case (mechanical_plasticity_type(ph))
    case (MECHANICAL_PLASTICITY_NONLOCAL,MECHANICAL_PLASTICITY_KINEHARDENING)

      Mp = matmul(matmul(transpose(phase_mechanical_Fi(ph)%data(1:3,1:3,en)),&
                         phase_mechanical_Fi(ph)%data(1:3,1:3,en)),&
                  phase_mechanical_S(ph)%data(1:3,1:3,en))

      plasticType: select case (mechanical_plasticity_type(ph))

        case (MECHANICAL_PLASTICITY_KINEHARDENING) plasticType
          call plastic_kinehardening_deltaState(Mp,ph,en)

        case (MECHANICAL_PLASTICITY_NONLOCAL) plasticType
          call plastic_nonlocal_deltaState(Mp,ph,en)

      end select plasticType

      broken = any(IEEE_is_NaN(plasticState(ph)%deltaState(:,en)))
      if (.not. broken) then
        mySize   = plasticState(ph)%sizeDeltaState
        plasticState(ph)%deltaState2(1:mySize,en) = plasticState(ph)%deltaState2(1:mySize,en) &
                                                  + plasticState(ph)%deltaState(1:mySize,en)
      end if

  end select

end function plastic_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief checks if a plastic module is active or not
!--------------------------------------------------------------------------------------------------
function plastic_active(plastic_label) result(active_plastic)

  character(len=*), intent(in)       :: plastic_label                                               !< type of plasticity model
  logical, dimension(:), allocatable :: active_plastic

  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl
  integer :: ph

  phases => config_material%get_dict('phase')
  allocate(active_plastic(phases%length), source = .false. )
  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    mech  => phase%get_dict('mechanical')
    pl    => mech%get_dict('plastic',defaultVal = emptyDict)
    active_plastic(ph) = pl%get_asStr('type',defaultVal='none') == plastic_label
  end do

end function plastic_active

end submodule plastic
