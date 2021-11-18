submodule(phase:mechanical) elastic

  type :: tParameters
    real(pReal) :: &
      C_11 = 0.0_pReal, &
      C_12 = 0.0_pReal, &
      C_13 = 0.0_pReal, &
      C_33 = 0.0_pReal, &
      C_44 = 0.0_pReal, &
      C_66 = 0.0_pReal
  end type tParameters

  type(tParameters), allocatable, dimension(:) :: param

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize elasticity.
!--------------------------------------------------------------------------------------------------
module subroutine elastic_init(phases)

  class(tNode), pointer :: &
    phases

  integer :: &
    ph
  class(tNode), pointer :: &
    phase, &
    mech, &
    elastic


  print'(/,1x,a)', '<<<+-  phase:mechanical:elastic init  -+>>>'
  print'(/,1x,a)', '<<<+-  phase:mechanical:elastic:Hooke init  -+>>>'

  print'(/,a,i0)', ' # phases: ',phases%length; flush(IO_STDOUT)

  allocate(param(phases%length))

  do ph = 1, phases%length
    phase   => phases%get(ph)
    mech    => phase%get('mechanical')
    elastic => mech%get('elastic')
    if (elastic%get_asString('type') /= 'Hooke') call IO_error(200,ext_msg=elastic%get_asString('type'))

    associate(prm => param(ph))

      prm%C_11 = elastic%get_asFloat('C_11')
      prm%C_12 = elastic%get_asFloat('C_12')
      prm%C_44 = elastic%get_asFloat('C_44')

      if (any(phase_lattice(ph) == ['hP','tI'])) then
        prm%C_13 = elastic%get_asFloat('C_13')
        prm%C_33 = elastic%get_asFloat('C_33')
      end if
      if (phase_lattice(ph) == 'tI') prm%C_66 = elastic%get_asFloat('C_66')

    end associate
  end do

end subroutine elastic_init


!--------------------------------------------------------------------------------------------------
!> @brief Return 6x6 elasticity tensor.
!--------------------------------------------------------------------------------------------------
module function elastic_C66(ph,en) result(C66)

  integer, intent(in) :: &
    ph, &
    en
  real(pReal), dimension(6,6) :: C66


  associate(prm => param(ph))
    C66 = 0.0_pReal
    C66(1,1) = prm%C_11
    C66(1,2) = prm%C_12
    C66(4,4) = prm%C_44

    if (any(phase_lattice(ph) == ['hP','tI'])) then
      C66(1,3) = prm%C_13
      C66(3,3) = prm%C_33
    end if

    if (phase_lattice(ph) == 'tI') C66(6,6) = prm%C_66

    C66 = lattice_symmetrize_C66(C66,phase_lattice(ph))

  end associate

end function elastic_C66


!--------------------------------------------------------------------------------------------------
!> @brief Return shear modulus.
!--------------------------------------------------------------------------------------------------
module function elastic_mu(ph,en) result(mu)

  integer, intent(in) :: &
    ph, &
    en
  real(pReal) :: &
    mu


  mu = lattice_equivalent_mu(elastic_C66(ph,en),'voigt')

end function elastic_mu


!--------------------------------------------------------------------------------------------------
!> @brief Return Poisson ratio.
!--------------------------------------------------------------------------------------------------
module function elastic_nu(ph,en) result(nu)

  integer, intent(in) :: &
    ph, &
    en
  real(pReal) :: &
    nu


  nu = lattice_equivalent_nu(elastic_C66(ph,en),'voigt')

end function elastic_nu



!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic and intermediate deformation gradients using Hooke's law
!--------------------------------------------------------------------------------------------------
module subroutine phase_hooke_SandItsTangents(S, dS_dFe, dS_dFi, &
                                              Fe, Fi, ph, en)

  integer, intent(in) :: &
    ph, &
    en
  real(pReal),   intent(in),  dimension(3,3) :: &
    Fe, &                                                                                           !< elastic deformation gradient
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),   intent(out), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dS_dFe, &                                                                                       !< derivative of 2nd P-K stress with respect to elastic deformation gradient
    dS_dFi                                                                                          !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

  real(pReal), dimension(3,3) :: E
  real(pReal), dimension(3,3,3,3) :: C
  integer :: &
    i, j


  C = math_66toSym3333(phase_damage_C66(phase_homogenizedC66(ph,en),ph,en))

  E = 0.5_pReal*(matmul(transpose(Fe),Fe)-math_I3)                                                  !< Green-Lagrange strain in unloaded configuration
  S = math_mul3333xx33(C,matmul(matmul(transpose(Fi),E),Fi))                                        !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

  do i =1, 3;do j=1,3
    dS_dFe(i,j,1:3,1:3) = matmul(Fe,matmul(matmul(Fi,C(i,j,1:3,1:3)),transpose(Fi)))                !< dS_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
    dS_dFi(i,j,1:3,1:3) = 2.0_pReal*matmul(matmul(E,Fi),C(i,j,1:3,1:3))                             !< dS_ij/dFi_kl = C_ijln * E_km * Fe_mn
  end do; end do

end subroutine phase_hooke_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!> ToDo: homogenizedC66 would be more consistent
!--------------------------------------------------------------------------------------------------
module function phase_homogenizedC66(ph,en) result(C)

  real(pReal), dimension(6,6) :: C
  integer,      intent(in)    :: ph, en


  plasticType: select case (phase_plasticity(ph))
    case (PLASTICITY_DISLOTWIN_ID) plasticType
     C = plastic_dislotwin_homogenizedC(ph,en)
    case default plasticType
     C = math_sym3333to66(math_Voigt66to3333(elastic_C66(ph,en)))
  end select plasticType

end function phase_homogenizedC66


end submodule elastic
