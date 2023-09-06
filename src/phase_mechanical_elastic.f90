submodule(phase:mechanical) elastic

  type :: tParameters
    type(tPolynomial) :: &
      C_11, &
      C_12, &
      C_13, &
      C_33, &
      C_44, &
      C_66
  end type tParameters

  type(tParameters), allocatable, dimension(:) :: param

contains

!--------------------------------------------------------------------------------------------------
!> @brief initialize elasticity
!--------------------------------------------------------------------------------------------------
module subroutine elastic_init(phases)

  type(tDict), pointer :: &
    phases

  integer :: &
    ph
  type(tDict), pointer :: &
    phase, &
    mech, &
    elastic
  character(len=:), allocatable :: refs


  print'(/,1x,a)', '<<<+-  phase:mechanical:elastic init  -+>>>'
  print'(/,1x,a)', '<<<+-  phase:mechanical:elastic:Hooke init  -+>>>'

  print'(/,1x,a,1x,i0)', '# phases:',phases%length; flush(IO_STDOUT)


  allocate(param(phases%length))

  do ph = 1, phases%length
    phase   => phases%get_dict(ph)
    mech    => phase%get_dict('mechanical')
    elastic => mech%get_dict('elastic')
    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(elastic,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs
    if (elastic%get_asStr('type') /= 'Hooke') call IO_error(200,ext_msg=elastic%get_asStr('type'))

    associate(prm => param(ph))

      prm%C_11 = polynomial(elastic,'C_11','T')
      prm%C_12 = polynomial(elastic,'C_12','T')
      prm%C_44 = polynomial(elastic,'C_44','T')

      if (any(phase_lattice(ph) == ['hP','tI'])) then
        prm%C_13 = polynomial(elastic,'C_13','T')
        prm%C_33 = polynomial(elastic,'C_33','T')
      end if

      if (phase_lattice(ph) == 'tI') &
        prm%C_66 = polynomial(elastic,'C_66','T')

    end associate
  end do

end subroutine elastic_init


!--------------------------------------------------------------------------------------------------
!> @brief return 6x6 elasticity tensor
!--------------------------------------------------------------------------------------------------
pure module function elastic_C66(ph,en) result(C66)

  integer, intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(6,6) :: C66
  real(pREAL) :: T


  associate(prm => param(ph))

    C66 = 0.0_pREAL
    T = thermal_T(ph,en)

    C66(1,1) = prm%C_11%at(T)
    C66(1,2) = prm%C_12%at(T)
    C66(4,4) = prm%C_44%at(T)

    if (any(phase_lattice(ph) == ['hP','tI'])) then
      C66(1,3) = prm%C_13%at(T)
      C66(3,3) = prm%C_33%at(T)
    end if

    if (phase_lattice(ph) == 'tI') C66(6,6) = prm%C_66%at(T)

    C66 = crystal_symmetrize_C66(C66,phase_lattice(ph))

  end associate

end function elastic_C66


!--------------------------------------------------------------------------------------------------
!> @brief return shear modulus
!--------------------------------------------------------------------------------------------------
pure module function elastic_mu(ph,en,isotropic_bound) result(mu)

  integer, intent(in) :: &
    ph, &
    en
  character(len=*), intent(in) :: isotropic_bound
  real(pREAL) :: &
    mu


  associate(prm => param(ph))

    mu = crystal_isotropic_mu(elastic_C66(ph,en),isotropic_bound,phase_lattice(ph))

  end associate

end function elastic_mu


!--------------------------------------------------------------------------------------------------
!> @brief return Poisson ratio
!--------------------------------------------------------------------------------------------------
pure module function elastic_nu(ph,en,isotropic_bound) result(nu)

  integer, intent(in) :: &
    ph, &
    en
  character(len=*), intent(in) :: isotropic_bound
  real(pREAL) :: &
    nu


  associate(prm => param(ph))

    nu = crystal_isotropic_nu(elastic_C66(ph,en),isotropic_bound,phase_lattice(ph))

  end associate

end function elastic_nu



!--------------------------------------------------------------------------------------------------
!> @brief return the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to
!> the elastic and intermediate deformation gradients using Hooke's law
! ToDo: Use Voigt matrix directly
!--------------------------------------------------------------------------------------------------
module subroutine phase_hooke_SandItsTangents(S, dS_dFe, dS_dFi, &
                                              Fe, Fi, ph, en)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL),   intent(in),  dimension(3,3) :: &
    Fe, &                                                                                           !< elastic deformation gradient
    Fi                                                                                              !< intermediate deformation gradient
  real(pREAL),   intent(out), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola-Kirchhoff stress tensor in lattice configuration
  real(pREAL),   intent(out), dimension(3,3,3,3) :: &
    dS_dFe, &                                                                                       !< derivative of 2nd P-K stress with respect to elastic deformation gradient
    dS_dFi                                                                                          !< derivative of 2nd P-K stress with respect to intermediate deformation gradient

  real(pREAL), dimension(3,3) :: E
  real(pREAL), dimension(6,6) :: C66
  real(pREAL), dimension(3,3,3,3) :: C
  integer :: &
    i, j


  C66 = phase_damage_C66(phase_homogenizedC66(ph,en),ph,en)
  C = math_Voigt66to3333_stiffness(C66)

  E = 0.5_pREAL*(matmul(transpose(Fe),Fe)-math_I3)                                                  !< Green-Lagrange strain in unloaded configuration
  S = math_Voigt6to33_stress(matmul(C66,math_33toVoigt6_strain(matmul(matmul(transpose(Fi),E),Fi))))!< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

  do i =1,3; do j=1,3
    dS_dFe(i,j,1:3,1:3) = matmul(Fe,matmul(matmul(Fi,C(i,j,1:3,1:3)),transpose(Fi)))                !< dS_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
    dS_dFi(i,j,1:3,1:3) = 2.0_pREAL*matmul(matmul(E,Fi),C(i,j,1:3,1:3))                             !< dS_ij/dFi_kl = C_ijln * E_km * Fe_mn
  end do; end do

end subroutine phase_hooke_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief Return the homogenized elasticity matrix.
!--------------------------------------------------------------------------------------------------
module function phase_homogenizedC66(ph,en) result(C)

  real(pREAL), dimension(6,6) :: C
  integer,      intent(in)    :: ph, en


  plasticType: select case (mechanical_plasticity_type(ph))
    case (MECHANICAL_PLASTICITY_DISLOTWIN) plasticType
     C = plastic_dislotwin_homogenizedC(ph,en)
    case default plasticType
     C = elastic_C66(ph,en)
  end select plasticType

end function phase_homogenizedC66


end submodule elastic
