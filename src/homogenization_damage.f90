!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) damage

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

  type(tDict), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationDamage
  integer :: ho,Nmembers


  print'(/,1x,a)', '<<<+-  homogenization:damage init  -+>>>'


  configHomogenizations => config_material%get_dict('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    Nmembers = count(material_homogenizationID == ho)
    allocate(current(ho)%phi(Nmembers), source=1.0_pReal)
    configHomogenization => configHomogenizations%get_dict(ho)
    associate(prm => param(ho))
      if (configHomogenization%contains('damage')) then
        configHomogenizationDamage => configHomogenization%get_dict('damage')
#if defined (__GFORTRAN__)
        prm%output = output_as1dString(configHomogenizationDamage)
#else
        prm%output = configHomogenizationDamage%get_as1dString('output',defaultVal=emptyStringArray)
#endif
        damageState_h(ho)%sizeState = 1
        allocate(damageState_h(ho)%state0(1,Nmembers), source=1.0_pReal)
        allocate(damageState_h(ho)%state (1,Nmembers), source=1.0_pReal)
      else
        prm%output = emptyStringArray
      end if
    end associate
  end do

  call pass_init()

end subroutine damage_init


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine damage_partition(ce)

  integer, intent(in) :: ce

  real(pReal) :: phi
  integer :: co


  if (damageState_h(material_homogenizationID(ce))%sizeState < 1) return
  phi = damagestate_h(material_homogenizationID(ce))%state(1,material_homogenizationEntry(ce))
  do co = 1, homogenization_Nconstituents(material_homogenizationID(ce))
    call phase_set_phi(phi,co,ce)
  end do


end subroutine damage_partition


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize damage viscosity.
!--------------------------------------------------------------------------------------------------
module function homogenization_mu_phi(ce) result(mu)

  integer, intent(in) :: ce
  real(pReal) :: mu


  mu = phase_mu_phi(1,ce)

end function homogenization_mu_phi


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize damage conductivity.
!--------------------------------------------------------------------------------------------------
module function homogenization_K_phi(ce) result(K)

  integer, intent(in) :: ce
  real(pReal), dimension(3,3) :: K


  K = phase_K_phi(1,ce)

end function homogenization_K_phi


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize damage driving force.
!--------------------------------------------------------------------------------------------------
module function homogenization_f_phi(phi,ce) result(f)

  integer, intent(in) :: ce
  real(pReal), intent(in) :: phi
  real(pReal) :: f


  f = phase_f_phi(phi, 1, ce)

end function homogenization_f_phi


!--------------------------------------------------------------------------------------------------
!> @brief Set damage field.
!--------------------------------------------------------------------------------------------------
module subroutine homogenization_set_phi(phi,ce)

  integer, intent(in) :: ce
  real(pReal), intent(in) :: phi

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
            call results_writeDataset(current(ho)%phi,group,prm%output(o),&
                                      'damage indicator','-')
        end select
      end do outputsLoop
  end associate

end subroutine damage_results

end submodule damage
