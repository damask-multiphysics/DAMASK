!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) homogenization_damage

  use lattice

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
    configHomogenizationDamage
  integer :: ho


  print'(/,a)', ' <<<+-  homogenization:damage init  -+>>>'
  print'(/,a)', ' <<<+-  homogenization:damage:isodamage init  -+>>>'

  configHomogenizations => config_material%get('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%phi(count(material_homogenizationAt2==ho)), source=1.0_pReal)
    configHomogenization => configHomogenizations%get(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('damage')) then
        configHomogenizationDamage => configHomogenization%get('damage')
#if defined (__GFORTRAN__)
        prm%output = output_asStrings(configHomogenizationDamage)
#else
        prm%output = configHomogenizationDamage%get_asStrings('output',defaultVal=emptyStringArray)
#endif
      else
        prm%output = emptyStringArray
      endif
    end associate
  enddo

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(ce)

  real(pReal) :: phi
  integer,     intent(in) :: ce

  integer :: co


  phi     = current(material_homogenizationAt2(ce))%phi(material_homogenizationMemberAt2(ce))
  do co = 1, homogenization_Nconstituents(material_homogenizationAt2(ce))
    call constitutive_damage_set_phi(phi,co,ce)
  enddo

end subroutine damage_partition



!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized nonlocal damage mobility
!--------------------------------------------------------------------------------------------------
module function damage_nonlocal_getMobility(ip,el) result(M)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  integer :: &
    co
  real(pReal) :: M

  M = 0.0_pReal

  do co = 1, homogenization_Nconstituents(material_homogenizationAt(el))
    M = M + lattice_M(material_phaseAt(co,el))
  enddo

  M = M/real(homogenization_Nconstituents(material_homogenizationAt(el)),pReal)

end function damage_nonlocal_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces
!--------------------------------------------------------------------------------------------------
module subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    phi
  real(pReal) :: &
    phiDot, dPhiDot_dPhi

  phiDot = 0.0_pReal
  dPhiDot_dPhi = 0.0_pReal

  call constitutive_damage_getRateAndItsTangents(phiDot, dPhiDot_dPhi, phi, ip, el)
  phiDot = phiDot/real(homogenization_Nconstituents(material_homogenizationAt(el)),pReal)
  dPhiDot_dPhi = dPhiDot_dPhi/real(homogenization_Nconstituents(material_homogenizationAt(el)),pReal)

end subroutine damage_nonlocal_getSourceAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
module subroutine damage_nonlocal_putNonLocalDamage(phi,ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    phi
  integer :: &
    homog, &
    offset

  homog  = material_homogenizationAt(el)
  offset = material_homogenizationMemberAt(ip,el)
  damage(homog)%p(offset) = phi

end subroutine damage_nonlocal_putNonLocalDamage


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine damage_nonlocal_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(damage_typeInstance(homog)))
  outputsLoop: do o = 1,size(prm%output)
    select case(prm%output(o))
      case ('phi')
        call results_writeDataset(group,damage(homog)%p,prm%output(o),&
                                  'damage indicator','-')
    end select
  enddo outputsLoop
  end associate

end subroutine damage_nonlocal_results

end submodule homogenization_damage
