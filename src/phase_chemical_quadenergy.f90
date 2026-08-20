! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Sharan Roongta, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quadratic approximation for chemical free energy
!--------------------------------------------------------------------------------------------------
submodule(phase:chemical) quadenergy

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      N_components
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
  character(len=:), allocatable :: &
    refs, &
    extmsg


  myChemicalEnergy = chemical_active('quadenergy')
  if (count(myChemicalEnergy) == 0) return
  print'(/,a)', ' <<<+-  phase:chemical:quadenergy init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(myChemicalEnergy); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')

  allocate(param(size(phases)))

  extmsg = ''

  do ph = 1, size(phases)
    if ( .not. myChemicalEnergy(ph)) cycle

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
!    refs = config_listReferences(chemical,indent=3)
!    if (len(refs) > 0) print'(/,1x,a)', refs

    associate(prm  => param(ph))
      phase => phases%get_dict(ph)
      chemical => phase%get_dict('chemical')

      ! read params
      components => chemical%get_dict('species',defaultVal=emptyDict)

      prm%N_components = size(components)

      Nmembers = count(material_ID_phase == ph)
      allocate(prm%c_0(size(components)),                     source=0.0_pREAL)
      allocate(prm%Mobility(size(components) - 1),                source=0.0_pREAL)
      allocate(prm%c_eq(size(components) - 1),                    source=0.0_pREAL)
      allocate(prm%coeff_linear(size(components) - 1),            source=0.0_pREAL)
      allocate(prm%coeff_quadratic(size(components) - 1),         source=0.0_pREAL)

      do com = 1, size(components) - 1
        component => components%get_dict(com)
        prm%Mobility(com)           = component%get_asReal('M')
        prm%c_0     (com)           = component%get_asReal('c_0')
        prm%c_eq    (com)           = component%get_asReal('c_eq')
        prm%coeff_linear (com)      = component%get_asReal('G,c',  defaultVal=0.0_pREAL)
        prm%coeff_quadratic (com)   = component%get_asReal('G,c^2')
      end do
      com = size(components)
      component => components%get_dict(com)
      prm%c_0(com) = component%get_asReal('c_0')

      ! sanity checks
      if (any(prm%Mobility < 0.0_pREAL)) extmsg = trim(extmsg)//' M'
      if (any(prm%c_0      < 0.0_pREAL)) extmsg = trim(extmsg)//' c_0'
      if (any(prm%c_eq     < 0.0_pREAL)) extmsg = trim(extmsg)//' c_eq'
      if (any(prm%coeff_quadratic == 0.0_pREAL)) extmsg = trim(extmsg)//'G,c^2'
      if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

      ! allocate fieldQuantities
      Nmembers = count(material_ID_phase == ph)
      allocate(current(ph)%C(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))
      allocate(current(ph)%dot_C(prm%N_components,Nmembers),source=0.0_pREAL)
      allocate(current(ph)%C0(prm%N_components,Nmembers),source=spread(prm%c_0,2,Nmembers))

    end associate
  end do

end function quadenergy_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate composition from diffusion potential.
!--------------------------------------------------------------------------------------------------
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


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derivative of composition with respect to diffusion potential.
!--------------------------------------------------------------------------------------------------
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


!--------------------------------------------------------------------------------------------------
!> @brief Return mobility.
!--------------------------------------------------------------------------------------------------
module function quadenergy_mobility(ph,en) result(mobility)
  integer, intent(in) :: ph, en
  real(pREAL), dimension(:,:), allocatable :: mobility

  integer :: com

  associate(prm => param(ph))

    allocate(mobility(prm%N_components-1,prm%N_components-1), source=0.0_pREAL)
    do com = 1, prm%N_components-1
      mobility(com,com) = prm%Mobility(com)
    end do

  end associate

end function quadenergy_mobility


end submodule quadenergy
