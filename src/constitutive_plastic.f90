submodule(constitutive) constitutive_plastic

  implicit none
 
  interface

    module subroutine plastic_none_init
    end subroutine plastic_none_init

    module subroutine plastic_isotropic_init
    end subroutine plastic_isotropic_init

    module subroutine plastic_phenopowerlaw_init
    end subroutine plastic_phenopowerlaw_init

    module subroutine plastic_kinehardening_init
    end subroutine plastic_kinehardening_init

    module subroutine plastic_dislotwin_init
    end subroutine plastic_dislotwin_init

    module subroutine plastic_disloUCLA_init
    end subroutine plastic_disloUCLA_init

    module subroutine plastic_nonlocal_init
    end subroutine plastic_nonlocal_init


    module subroutine plastic_isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_isotropic_LpAndItsTangent

    pure module subroutine plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_phenopowerlaw_LpAndItsTangent

    pure module subroutine plastic_kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_kinehardening_LpAndItsTangent

    module subroutine plastic_dislotwin_LpAndItsTangent(Lp,dLp_dMp,Mp,T,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_dislotwin_LpAndItsTangent

    pure module subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp,Mp,T,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                     intent(in) :: &
        T
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_disloUCLA_LpAndItsTangent

    module subroutine plastic_nonlocal_LpAndItsTangent(Lp,dLp_dMp, &
                                                       Mp,Temperature,instance,of,ip,el)
      real(pReal), dimension(3,3),     intent(out) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out) :: &
        dLp_dMp                                                                                     !< derivative of Lp with respect to the Mandel stress

      real(pReal), dimension(3,3),     intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                     intent(in) :: &
        Temperature
      integer,                         intent(in) :: &
        instance, &
        of, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine plastic_nonlocal_LpAndItsTangent


    module subroutine plastic_dislotwin_dependentState(T,instance,of)
      integer,       intent(in) :: &
        instance, &
        of
      real(pReal),   intent(in) :: &
        T
    end subroutine plastic_dislotwin_dependentState

    module subroutine plastic_disloUCLA_dependentState(instance,of)
      integer,       intent(in) :: &
        instance, &
        of
    end subroutine plastic_disloUCLA_dependentState

    module subroutine plastic_nonlocal_dependentState(F, Fp, instance, of, ip, el)
      real(pReal), dimension(3,3), intent(in) :: &
        F, &
        Fp
      integer, intent(in) :: &
        instance, &
        of, &
        ip, &
        el
    end subroutine plastic_nonlocal_dependentState

  end interface


contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
module subroutine plastic_init

!--------------------------------------------------------------------------------------------------
! initialized plasticity
  if (any(phase_plasticity == PLASTICITY_NONE_ID))          call plastic_none_init
  if (any(phase_plasticity == PLASTICITY_ISOTROPIC_ID))     call plastic_isotropic_init
  if (any(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)) call plastic_phenopowerlaw_init
  if (any(phase_plasticity == PLASTICITY_KINEHARDENING_ID)) call plastic_kinehardening_init
  if (any(phase_plasticity == PLASTICITY_DISLOTWIN_ID))     call plastic_dislotwin_init
  if (any(phase_plasticity == PLASTICITY_DISLOUCLA_ID))     call plastic_disloucla_init
  if (any(phase_plasticity == PLASTICITY_NONLOCAL_ID)) then
    call plastic_nonlocal_init
  else
    call geometry_plastic_nonlocal_disable
  endif

end subroutine plastic_init



!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
module procedure constitutive_plastic_dependentState

  integer :: &
    ho, &                                                                                           !< homogenization
    tme, &                                                                                          !< thermal member position
    instance, of

  ho  = material_homogenizationAt(el)
  tme = thermalMapping(ho)%p(ip,el)
  of  = material_phasememberAt(ipc,ip,el)
  instance = phase_plasticityInstance(material_phaseAt(ipc,el))

  plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call plastic_dislotwin_dependentState(temperature(ho)%p(tme),instance,of)
    case (PLASTICITY_DISLOUCLA_ID) plasticityType
      call plastic_disloUCLA_dependentState(instance,of)
    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call plastic_nonlocal_dependentState (F,Fp,instance,of,ip,el)
  end select plasticityType

end procedure constitutive_plastic_dependentState

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: Discuss whether it makes sense if crystallite handles the configuration conversion, i.e.
! Mp in, dLp_dMp out
!--------------------------------------------------------------------------------------------------
module procedure constitutive_plastic_LpAndItsTangents

  real(pReal), dimension(3,3,3,3) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to Mandel stress
  real(pReal), dimension(3,3) :: &
    Mp                                                                                              !< Mandel stress work conjugate with Lp
  integer :: &
    ho, &                                                                                           !< homogenization
    tme                                                                                             !< thermal member position
  integer :: &
    i, j, instance, of

  ho = material_homogenizationAt(el)
  tme = thermalMapping(ho)%p(ip,el)

  Mp = matmul(matmul(transpose(Fi),Fi),S)
  of = material_phasememberAt(ipc,ip,el)
  instance = phase_plasticityInstance(material_phaseAt(ipc,el))

  plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))

    case (PLASTICITY_NONE_ID) plasticityType
      Lp = 0.0_pReal
      dLp_dMp = 0.0_pReal

    case (PLASTICITY_ISOTROPIC_ID) plasticityType
      call plastic_isotropic_LpAndItsTangent   (Lp,dLp_dMp,Mp,instance,of)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
      call plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call plastic_kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call plastic_nonlocal_LpAndItsTangent     (Lp,dLp_dMp,Mp, temperature(ho)%p(tme),instance,of,ip,el)

    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call plastic_dislotwin_LpAndItsTangent    (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

    case (PLASTICITY_DISLOUCLA_ID) plasticityType
      call plastic_disloucla_LpAndItsTangent    (Lp,dLp_dMp,Mp,temperature(ho)%p(tme),instance,of)

  end select plasticityType

  do i=1,3; do j=1,3
    dLp_dFi(i,j,1:3,1:3) = matmul(matmul(Fi,S),transpose(dLp_dMp(i,j,1:3,1:3))) + &
                           matmul(matmul(Fi,dLp_dMp(i,j,1:3,1:3)),S)
    dLp_dS(i,j,1:3,1:3)  = matmul(matmul(transpose(Fi),Fi),dLp_dMp(i,j,1:3,1:3))                     ! ToDo: @PS: why not:   dLp_dMp:(FiT Fi)
  enddo; enddo

end procedure constitutive_plastic_LpAndItsTangents

end submodule constitutive_plastic

