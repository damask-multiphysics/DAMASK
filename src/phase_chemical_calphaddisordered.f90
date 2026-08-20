! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Aadhithyan Kannan, KU Leuven
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Disordered regular solution for free energy using SGTE representation (second generation)
!> @details Thermodynamic data from CALPHAD-based databases
!> @details no contribution from other physics. For example, magnetic and pressure dependance
!--------------------------------------------------------------------------------------------------
submodule(phase:chemical) calphaddisordered

  type :: tParameters                                                                       !< container type for internal constitutive parameters
    integer :: &
      N_components
    real(pREAL), dimension(:), allocatable :: &
      Mobility, &
      c_0
    real(pREAL), dimension(:,:,:), allocatable :: &
      G0_coeff                                                 ! for each component, for each temperature range, and third dimension: Temp range(2) + coefficients(8)
    real(pREAL), dimension(:, :, :, :), allocatable :: &
      Tranges_L0                                                ! for each component, for each interacting component, for each temperature range, and the range
    type(tpolynomial), dimension(:, :, :), allocatable :: L0_T  ! for each component, for each interacting component, for each temperature range
  end type tParameters

  type(tParameters), dimension(:),   allocatable :: param                                   !< containers of constitutive parameters (len Ninstances)

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details Read in material parameters, allocate arrays and does sanity checks
!--------------------------------------------------------------------------------------------------
module function calphaddisordered_init() result(myChemicalEnergy)

  logical, dimension(:), allocatable :: myChemicalEnergy

  type(tDict), pointer :: &
    phases, &
    phase, &
    chemical, components, component, &
    G_xs, &
    G0_range, Gxs_range

  type(tList), pointer :: &
    G_0, &
    Gxs_comp

  integer :: ph, Nmembers, com, xs_com, maxSize_G0lists, maxSize_Gxscomlists, i

  character(len=:), allocatable :: &
    refs, &
    extmsg

  myChemicalEnergy = chemical_active('calphaddisordered')
  if (count(myChemicalEnergy) == 0) return

  print'(/,a)', ' <<<+-  phase:chemical:calphaddisordered init  -+>>>'
  print'(a,i2)', ' # phases: ',count(myChemicalEnergy); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))

  extmsg = ''

  do ph = 1, size(phases)
    if ( .not. myChemicalEnergy(ph)) cycle

    phase => phases%get_dict(ph)
    chemical => phase%get_dict('chemical')
    components => chemical%get_dict('species',defaultVal=emptyDict)

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(chemical,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

    maxSize_G0lists = 0
    maxSize_Gxscomlists = 0

    ! parse to allocate
    do com = 1, size(components)
      component => components%get_dict(com)
      maxSize_G0lists = max(maxSize_G0lists, size(component%get_list('G_0', defaultVal=emptyList)))
      if (com > 1) then
        G_xs => component%get_dict('G_xs',defaultVal=emptyDict)
        do xs_com = 1, size(G_xs)
          maxSize_Gxscomlists = max(maxSize_Gxscomlists, size(G_xs%get_list(G_xs%key(xs_com), defaultVal=emptyList)))
        end do
      endif
    end do

    associate(prm  => param(ph))

      prm%N_components = size(components)
      allocate(prm%c_0(size(components)),        source=0.0_pREAL)
      allocate(prm%Mobility(size(components)),   source=0.0_pREAL)
      allocate(prm%G0_coeff(size(components),maxSize_G0lists,10), source=0.0_pREAL)
      allocate(prm%L0_T(size(components),size(components),maxSize_Gxscomlists))
      allocate(prm%Tranges_L0(size(components),size(components),maxSize_Gxscomlists,2), source=0.0_pREAL)

      ! parse and store
      do com = 1, size(components)
        component => components%get_dict(com)
        prm%Mobility(com)           = component%get_asReal('M')
        prm%c_0     (com)           = component%get_asReal('c_0')

        G_0 => component%get_list('G_0')
        do i = 1, size(G_0)
          G0_range  => G_0%get_dict(i)
          prm%G0_coeff (com, i, 1:2)  = G0_range%get_as1dReal('T_range', requiredSize=2)
          prm%G0_coeff (com, i, 3:)   = G0_range%get_as1dReal('SGTE_coeff', requiredSize=8)
        end do

        if (com > 1) then                  ! by default expect interactions for the second component onwards
          G_xs => component%get_dict('G_xs')
          do xs_com = 1, size(G_xs)
            Gxs_comp => G_xs%get_list(components%key(xs_com)) ! order corresponds to order in components dictionary
            do i = 1, size(Gxs_comp)
              Gxs_range => Gxs_comp%get_dict(i)
              prm%Tranges_L0 (com, xs_com, i, :) = Gxs_range%get_as1dReal('T_range', requiredSize=2)
              prm%L0_T(com, xs_com, i) = polynomial(Gxs_range, 'L_0', 'T')   ! only lower triangle of interaction matrix. symmetry applied later
            end do
          end do
        endif
      end do

      ! sanity checks
      if (any(prm%Mobility < 0.0_pREAL)) extmsg = trim(extmsg)//' M'
      if (any(prm%c_0      < 0.0_pREAL)) extmsg = trim(extmsg)//' c_0'
      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

      ! allocate fieldQuantities
      Nmembers = count(material_ID_phase == ph)
      ! call phase_allocateState(chemicalState(ph),Nmembers,1,1,0)
      allocate(current(ph)%C(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))
      allocate(current(ph)%dot_C(prm%N_components,Nmembers),source=0.0_pREAL)
      allocate(current(ph)%C0(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))

    end associate

  end do

end function calphaddisordered_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate chemical potential explicit.
!--------------------------------------------------------------------------------------------------
module function calphaddisordered_mu_explicit(ph,en,comp_prev) result(mu_explicit)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:), intent(in)    :: comp_prev
  real(pREAL), dimension(:), allocatable   :: mu_explicit

  real(pREAL), dimension(:), allocatable   :: G0
  real(pREAL), dimension(:, :), allocatable :: L0
  integer :: com_i, com_j

  associate(prm => param(ph))

    allocate(mu_explicit(prm%N_components - 1), source=0.0_pREAL)

    L0 = get_L0(ph, en)
    G0 = get_G0(ph, en)

    do com_i = 1, prm%N_components-1
      mu_explicit(com_i) = G0(com_i) - G0(prm%N_components) + &
                              L0(com_i, prm%N_components)*(comp_prev(prm%N_components) - comp_prev(com_i))
      do com_j = 1, prm%N_components-1
        if (com_j /= com_i) then
          mu_explicit(com_i) = mu_explicit(com_i) + (L0(com_i, com_j) - L0(prm%N_components, com_j))*(comp_prev(com_j))
        endif
      end do
    end do

 end associate

end function calphaddisordered_mu_explicit


!--------------------------------------------------------------------------------------------------
!> @brief Calculate G0: reference free energy for pure components evaluated at temperature T
!--------------------------------------------------------------------------------------------------
pure function get_G0(ph, en) result(G0)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:),allocatable :: G0

  real(pREAL) :: T
  integer :: com_i, i

  T = thermal_T(ph,en)

  associate(prm => param(ph))

    allocate(G0(prm%N_components), source=0.0_pREAL)

    do com_i = 1, prm%N_components
      do i = 1, size(prm%G0_coeff,2)
        if (prm%G0_coeff(com_i,i,1) <= T .and. T <= prm%G0_coeff(com_i,i,2)) then
          G0(com_i) = prm%G0_coeff(com_i, i, 3) + prm%G0_coeff(com_i, i, 4)*T + prm%G0_coeff(com_i, i, 5)*T*log(T) + &
                                                  prm%G0_coeff(com_i, i, 6)*T**2 + prm%G0_coeff(com_i, i, 7)*T**3 + &
                                                  prm%G0_coeff(com_i, i, 8)*T**(-1) + prm%G0_coeff(com_i, i, 9)*T**7 + &
                                                  prm%G0_coeff(com_i, i, 10)*T**(-9)
          exit
        endif
      end do
    end do

 end associate

end function get_G0


!--------------------------------------------------------------------------------------------------
!> @brief calculates L0: interaction coefficients(zero order) evaluated at temperature T
!--------------------------------------------------------------------------------------------------
pure function get_L0(ph, en) result(L0)

  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:, :),allocatable :: L0

  real(pREAL) :: T
  integer :: com_i, com_j, k


  T = thermal_T(ph,en)

  associate(prm => param(ph))
    allocate(L0(prm%N_components, prm%N_components), source=0.0_pREAL)

    do com_i = 1, prm%N_components-1
      do com_j = com_i+1, prm%N_components
        do k = 1, size(prm%Tranges_L0,3)
          if (prm%Tranges_L0(com_j, com_i, k, 1) <= T .and. T <= prm%Tranges_L0(com_j, com_i, k, 2)) then
            L0(com_j, com_i) = prm%L0_T(com_j, com_i, k)%at(T)
            L0(com_i, com_j) = L0(com_j, com_i)
            exit
          endif
        end do
      end do
    end do

  end associate

end function get_L0


!--------------------------------------------------------------------------------------------------
!> @brief calculate composition from diffusion potential.
!--------------------------------------------------------------------------------------------------
module function calphaddisordered_composition(mu_chemical,comp_prev,ph,en) result(comp)
  real(pREAL), dimension(:), intent(in) :: mu_chemical
  real(pREAL), dimension(:), intent(in) :: comp_prev
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:), allocatable :: comp

  real(pREAL) :: T
  real(pREAL), dimension(:), allocatable :: mu_explicit

  associate(prm => param(ph))

    allocate(comp(prm%N_components), source=0.0_pREAL)
    T = thermal_T(ph,en)

    mu_explicit = calphaddisordered_mu_explicit(ph,en,comp_prev)
    comp(1:prm%N_components-1) = exp((mu_chemical - mu_explicit)/R/T)/ &
          (1.0_pREAL+sum(exp((mu_chemical - mu_explicit)/R/T)))

    comp(prm%N_components) = 1.0_pREAL - sum(comp(1:prm%N_components-1))

  end associate

end function calphaddisordered_composition


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derivative of composition with respect to diffusion potential.
!--------------------------------------------------------------------------------------------------
module function calphaddisordered_compositionTangent(mu_chemical,comp_prev,ph,en) result(comp_tangent)

  real(pREAL), dimension(:), intent(in) :: mu_chemical
  real(pREAL), dimension(:), intent(in) :: comp_prev
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:,:),allocatable :: comp_tangent

  real(pREAL) :: T
  real(pREAL), dimension(:), allocatable :: comp

  integer :: com_i, com_j

  associate(prm => param(ph))

    allocate(comp_tangent(prm%N_components-1,prm%N_components-1), source= 0.0_pREAL)
    T = thermal_T(ph,en)

    comp = calphaddisordered_composition(mu_chemical,comp_prev,ph,en)
    do com_i = 1, prm%N_components-1
      comp_tangent(com_i,com_i) = comp(com_i)/R/T
      do com_j = 1, prm%N_components-1
        comp_tangent(com_i,com_j) = comp_tangent(com_i,com_j) - comp(com_i)*comp(com_j)/R/T
      end do
    end do

  end associate

end function calphaddisordered_compositionTangent


!--------------------------------------------------------------------------------------------------
!> @brief Return mobility. ToDo: from mobility database?
!--------------------------------------------------------------------------------------------------
module function calphaddisordered_mobility(ph,en) result(mobility)
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:,:),allocatable :: mobility

  integer :: com

  associate(prm => param(ph))

    allocate(mobility(prm%N_components-1,prm%N_components-1), source=0.0_pREAL)
    do com = 1, prm%N_components-1
      mobility(com,com) = prm%Mobility(com)
    end do

  end associate

end function calphaddisordered_mobility


end submodule calphaddisordered
