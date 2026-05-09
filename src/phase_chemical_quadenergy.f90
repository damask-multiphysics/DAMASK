! SPDX-License-Identifier: AGPL-3.0-or-later
submodule(phase:chemical) quadenergy

  real(pREAL), parameter :: &
    R = 8.314459848_pREAL                                                                           !< gas constant. Not sure where the exact value comes from, check https://en.wikipedia.org/wiki/Gas_constant
    ! R = N_A * K_B

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      N_components
    real(pREAL) :: &
      coeff_constant                                                                                !< constant energy
    real(pREAL), dimension(:), allocatable :: &
      Mobility, &
      c_0, &
      c_eq, &
      coeff_linear, &
      coeff_quadratic
  end type tParameters

  type(tParameters), dimension(:),   allocatable :: param                                           !< containers of constitutive parameters (len Ninstances)

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function quadenergy_init() result(myChemicalEnergy)

  logical, dimension(:), allocatable :: myChemicalEnergy

  type(tDict), pointer :: &
    phases, &
    phase, &
    chemical, components, component
  integer :: Nmembers,ph,com


  myChemicalEnergy = chemical_active('quadenergy')
  if (count(myChemicalEnergy) == 0) return
  print'(/,a)', ' <<<+-  phase:chemical:quadenergy init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(myChemicalEnergy); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')

  allocate(param(size(phases)))

  do ph = 1, size(phases)
    if ( .not. myChemicalEnergy(ph)) cycle

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)

    associate(prm  => param(ph))
      phase => phases%get_dict(ph)
      chemical => phase%get_dict('chemical')

      ! read params
      components => chemical%get_dict('components',defaultVal=emptyDict)

      prm%N_components = size(components)

      Nmembers = count(material_ID_phase == ph)
      allocate(prm%c_0(size(components)),                     source=0.0_pREAL)
      allocate(prm%Mobility(size(components)),                source=0.0_pREAL)
      allocate(prm%c_eq(size(components)),                    source=0.0_pREAL)
      allocate(prm%coeff_linear(size(components)),            source=0.0_pREAL)
      allocate(prm%coeff_quadratic(size(components)),         source=0.0_pREAL)

      do com = 1, size(components)
        component => components%get_dict(com)
        prm%Mobility(com)           = component%get_asReal('M')
        prm%c_0     (com)           = component%get_asReal('c_0')
        prm%c_eq    (com)           = component%get_asReal('c_eq', defaultVal=0.0_pREAL)
        prm%coeff_linear (com)      = component%get_asReal('G,c',  defaultVal=0.0_pREAL)
        prm%coeff_quadratic (com)   = component%get_asReal('G,c^2',defaultVal=0.0_pREAL)
      end do

      Nmembers = count(material_ID_phase == ph)
      allocate(current(ph)%C(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))
      allocate(current(ph)%dot_C(prm%N_components,Nmembers),source=0.0_pREAL)
      allocate(current(ph)%C0(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))

    end associate
  end do

end function quadenergy_init


module function quadenergy_composition(mu_chemical,ph,en) result(comp)
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:), intent(in) :: mu_chemical
  real(pREAL), dimension(:), allocatable :: comp

  integer :: com

  associate(prm => param(ph))
    allocate(comp(prm%N_components))

    do com = 1, prm%N_components - 1
      comp(com) = prm%c_eq(com) + 0.5_pREAL *(mu_chemical(com) - prm%coeff_linear(com))/prm%coeff_quadratic(com)
    end do
    comp(prm%N_components) = 1.0_pREAL - sum(comp(1:prm%N_components-1))

  end associate

end function quadenergy_composition


module function quadenergy_compositionTangent(mu_chemical,ph,en) result(comp_tangent)
  real(pREAL), dimension(:), intent(in) :: mu_chemical
  integer, intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(:,:),allocatable :: comp_tangent

  integer :: com

  associate(prm => param(ph))

    allocate(comp_tangent(prm%N_components-1,prm%N_components-1), source = 0.0_pREAL)

    do com = 1, prm%N_components-1
      comp_tangent(com,com) = 0.5_pREAL/prm%coeff_quadratic(com)
    end do

  end associate

end function quadenergy_compositionTangent


module function quadenergy_mobility(ph,en) result(mobility)
  integer, intent(in) :: ph, en
  real(pREAL), dimension(:,:),allocatable :: mobility

  integer :: com

  associate(prm => param(ph))

    allocate(mobility(prm%N_components-1,prm%N_components-1))
    mobility = 0.0_pREAL
    do com = 1, prm%N_components-1
      mobility(com,com) = prm%Mobility(com)
    end do

  end associate

end function quadenergy_mobility


end submodule quadenergy
