submodule(phase:mechanical) elastic

  enum, bind(c); enumerator :: &
    ELASTICITY_UNDEFINED_ID, &
    ELASTICITY_HOOKE_ID, &
    STIFFNESS_DEGRADATION_UNDEFINED_ID, &
    STIFFNESS_DEGRADATION_DAMAGE_ID
  end enum

  integer, dimension(:), allocatable :: &
    phase_NstiffnessDegradations
 
  integer(kind(ELASTICITY_UNDEFINED_ID)), dimension(:),   allocatable :: &
    phase_elasticity                                                                                !< elasticity of each phase
  integer(kind(STIFFNESS_DEGRADATION_UNDEFINED_ID)),     dimension(:,:), allocatable :: &
    phase_stiffnessDegradation                                                                      !< active stiffness degradation mechanisms of each phase

contains


module subroutine elastic_init(phases)

  class(tNode), pointer :: &
    phases

  integer :: &
    ph, &
    stiffDegradationCtr
  class(tNode), pointer :: &
    phase, &
    mech, &
    elastic, &
    stiffDegradation
 
  print'(/,a)', ' <<<+-  phase:mechanical:elastic init  -+>>>'

  allocate(phase_elasticity(phases%length), source = ELASTICITY_undefined_ID)
  allocate(phase_NstiffnessDegradations(phases%length),source=0)
  
  do ph = 1, phases%length
    phase   => phases%get(ph)
    mech    => phase%get('mechanical')
    elastic => mech%get('elastic')
    if(IO_lc(elastic%get_asString('type')) == 'hooke') then ! accept small letter h for the moment
      phase_elasticity(ph) = ELASTICITY_HOOKE_ID
    else
      call IO_error(200,ext_msg=elastic%get_asString('type'))
    endif
    stiffDegradation => mech%get('stiffness_degradation',defaultVal=emptyList)                      ! check for stiffness degradation mechanisms
    phase_NstiffnessDegradations(ph) = stiffDegradation%length
  enddo

  allocate(phase_stiffnessDegradation(maxval(phase_NstiffnessDegradations),phases%length), &
                        source=STIFFNESS_DEGRADATION_undefined_ID)

  if(maxVal(phase_NstiffnessDegradations)/=0) then
    do ph = 1, phases%length
      phase => phases%get(ph)
      mech    => phase%get('mechanical')
      stiffDegradation => mech%get('stiffness_degradation',defaultVal=emptyList)
      do stiffDegradationCtr = 1, stiffDegradation%length
        if(stiffDegradation%get_asString(stiffDegradationCtr) == 'damage') &
            phase_stiffnessDegradation(stiffDegradationCtr,ph) = STIFFNESS_DEGRADATION_damage_ID
      enddo
    enddo
  endif

end subroutine elastic_init


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
    d, &                                                                                            !< counter in degradation loop
    i, j

  C = math_66toSym3333(phase_homogenizedC(ph,en))
  C = phase_damage_C(C,ph,en)

  E = 0.5_pReal*(matmul(transpose(Fe),Fe)-math_I3)                                                  !< Green-Lagrange strain in unloaded configuration
  S = math_mul3333xx33(C,matmul(matmul(transpose(Fi),E),Fi))                                        !< 2PK stress in lattice configuration in work conjugate with GL strain pulled back to lattice configuration

  do i =1, 3;do j=1,3
    dS_dFe(i,j,1:3,1:3) = matmul(Fe,matmul(matmul(Fi,C(i,j,1:3,1:3)),transpose(Fi)))                !< dS_ij/dFe_kl = C_ijmn * Fi_lm * Fi_on * Fe_ko
    dS_dFi(i,j,1:3,1:3) = 2.0_pReal*matmul(matmul(E,Fi),C(i,j,1:3,1:3))                             !< dS_ij/dFi_kl = C_ijln * E_km * Fe_mn
  enddo; enddo

end subroutine phase_hooke_SandItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenized elasticity matrix
!> ToDo: homogenizedC66 would be more consistent
!--------------------------------------------------------------------------------------------------
module function phase_homogenizedC(ph,en) result(C)

  real(pReal), dimension(6,6) :: C
  integer,      intent(in)    :: ph, en

  plasticType: select case (phase_plasticity(ph))
    case (PLASTICITY_DISLOTWIN_ID) plasticType
     C = plastic_dislotwin_homogenizedC(ph,en)
    case default plasticType
     C = lattice_C66(1:6,1:6,ph)
  end select plasticType

end function phase_homogenizedC


end submodule elastic
