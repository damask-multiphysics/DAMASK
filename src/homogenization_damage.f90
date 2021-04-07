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
    configHomogenizationDamage
  integer :: ho


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

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(ce)

  real(pReal) :: phi
  integer,     intent(in) :: ce

  integer :: co


  if(damageState_h(material_homogenizationID(ce))%sizeState < 1) return
  phi = damagestate_h(material_homogenizationID(ce))%state(1,material_homogenizationEntry(ce))
  do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
    call phase_damage_set_phi(phi,co,ce)
  enddo

end subroutine damage_partition



!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized nonlocal damage mobility
!--------------------------------------------------------------------------------------------------
module function damage_nonlocal_getMobility(ce) result(M)

  integer, intent(in) :: ce
  integer :: &
    co
  real(pReal) :: M

  M = 0.0_pReal

  do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
    M = M + lattice_M(material_phaseID(co,ce))
  enddo

  M = M/real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end function damage_nonlocal_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces
!--------------------------------------------------------------------------------------------------
module subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, phi, ce)

  integer, intent(in) :: ce
  real(pReal), intent(in) :: &
    phi
  real(pReal), intent(out) :: &
    phiDot

  real(pReal) :: &
    dPhiDot_dPhi

  call phase_damage_getRateAndItsTangents(phiDot, dPhiDot_dPhi, phi, ce)
  phiDot = phiDot/real(homogenization_Nconstituents(material_homogenizationID(ce)),pReal)

end subroutine damage_nonlocal_getSourceAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
module subroutine damage_nonlocal_putNonLocalDamage(phi,ce)

  integer, intent(in) :: ce
  real(pReal),   intent(in) :: &
    phi
  integer :: &
    ho, &
    en

  ho = material_homogenizationID(ce)
  en = material_homogenizationEntry(ce)
  damagestate_h(ho)%state(1,en) = phi

end subroutine damage_nonlocal_putNonLocalDamage


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
