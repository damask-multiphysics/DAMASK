!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) damage

  use lattice

  interface

    module subroutine pass_init
    end subroutine pass_init

  end interface

  type :: tDataContainer
    real(pReal), dimension(:), allocatable :: phi
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param

contains


!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine damage_init()

  class(tNode), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationDamage, &
    num_generic, &
    material_homogenization
  integer :: ho
  integer :: Ninstances,Nmaterialpoints,h


  print'(/,a)', ' <<<+-  homogenization:damage init  -+>>>'
  print'(/,a)', ' <<<+-  homogenization:damage:pass init  -+>>>'

  configHomogenizations => config_material%get('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%phi(count(material_homogenizationID==ho)), source=1.0_pReal)
    configHomogenization => configHomogenizations%get(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('damage')) then
        configHomogenizationDamage => configHomogenization%get('damage')
#if defined (__GFORTRAN__)
        prm%output = output_as1dString(configHomogenizationDamage)
#else
        prm%output = configHomogenizationDamage%get_as1dString('output',defaultVal=emptyStringArray)
#endif
      else
        prm%output = emptyStringArray
      endif
    end associate
  enddo

!------------------------------------------------------------------------------------
! read numerics parameter
  num_generic => config_numerics%get('generic',defaultVal= emptyDict)
  num_damage%charLength = num_generic%get_asFloat('charLength',defaultVal=1.0_pReal)

  Ninstances = count(damage_type == DAMAGE_nonlocal_ID)

  material_homogenization => config_material%get('homogenization')
  do h = 1, material_homogenization%length
    if (damage_type(h) /= DAMAGE_NONLOCAL_ID) cycle

    Nmaterialpoints = count(material_homogenizationAt == h)
    damageState_h(h)%sizeState = 1
    allocate(damageState_h(h)%state0   (1,Nmaterialpoints), source=1.0_pReal)
    allocate(damageState_h(h)%state    (1,Nmaterialpoints), source=1.0_pReal)

  enddo

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(ce)

  real(pReal) :: phi
  integer,     intent(in) :: ce


  if(damageState_h(material_homogenizationID(ce))%sizeState < 1) return
  phi = damagestate_h(material_homogenizationID(ce))%state(1,material_homogenizationEntry(ce))
  call phase_set_phi(phi,1,ce)

end subroutine damage_partition



!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized nonlocal damage mobility
!--------------------------------------------------------------------------------------------------
module function homogenization_mu_phi(ce) result(M)

  integer, intent(in) :: ce
  real(pReal) :: M

  M = lattice_M(material_phaseID(1,ce))

end function homogenization_mu_phi


!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces
!--------------------------------------------------------------------------------------------------
module function homogenization_f_phi(phi,ce) result(f_phi)

  integer, intent(in) :: ce
  real(pReal), intent(in) :: &
    phi
  real(pReal) :: f_phi

  f_phi = phase_f_phi(phi, 1, ce)

end function homogenization_f_phi


!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
module subroutine homogenization_set_phi(phi,ce)

  integer, intent(in) :: ce
  real(pReal),   intent(in) :: &
    phi
  integer :: &
    ho, &
    en

  ho = material_homogenizationID(ce)
  en = material_homogenizationEntry(ce)
  damagestate_h(ho)%state(1,en) = phi
  current(ho)%phi(en) = phi

end subroutine homogenization_set_phi


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine damage_results(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(ho))
  outputsLoop: do o = 1,size(prm%output)
    select case(prm%output(o))
      case ('phi')
        call results_writeDataset(group,damagestate_h(ho)%state(1,:),prm%output(o),&
                                  'damage indicator','-')
    end select
  enddo outputsLoop
  end associate

end subroutine damage_results

end submodule damage
