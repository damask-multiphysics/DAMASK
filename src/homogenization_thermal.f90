!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!--------------------------------------------------------------------------------------------------
submodule(homogenization) thermal

  interface

    module subroutine pass_init
    end subroutine pass_init

    module subroutine isotemperature_init
    end subroutine isotemperature_init

  end interface

  type :: tDataContainer
    real(pREAL), dimension(:), allocatable :: T, dot_T
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param


contains

!--------------------------------------------------------------------------------------------------
!> @brief Allocate variables and set parameters.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_init()

  type(tDict), pointer :: &
    configHomogenizations, &
    configHomogenization, &
    configHomogenizationThermal
  integer :: ho


  print'(/,1x,a)', '<<<+-  homogenization:thermal init  -+>>>'


  configHomogenizations => config_material%get_dict('homogenization')
  allocate(param(configHomogenizations%length))
  allocate(current(configHomogenizations%length))

  do ho = 1, configHomogenizations%length
    allocate(current(ho)%T(count(material_ID_homogenization==ho)), source=T_ROOM)
    allocate(current(ho)%dot_T(count(material_ID_homogenization==ho)), source=0.0_pREAL)
    configHomogenization => configHomogenizations%get_dict(ho)
    associate(prm => param(ho))

      if (configHomogenization%contains('thermal')) then
        configHomogenizationThermal => configHomogenization%get_dict('thermal')
#if defined (__GFORTRAN__)
        prm%output = output_as1dStr(configHomogenizationThermal)
#else
        prm%output = configHomogenizationThermal%get_as1dStr('output',defaultVal=emptyStrArray)
#endif
        select case (configHomogenizationThermal%get_asStr('type'))

          case ('pass')
            call pass_init()

          case ('isotemperature')
            call isotemperature_init()

        end select
      else
        prm%output = emptyStrArray
      end if

    end associate
  end do

end subroutine thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief Check if thermal homogemization description is present in the configuration file
!--------------------------------------------------------------------------------------------------
module function homogenization_thermal_active() result(active)

  logical :: active

  active = any(thermal_active(:))

end function homogenization_thermal_active


!--------------------------------------------------------------------------------------------------
!> @brief Partition temperature onto the individual constituents.
!--------------------------------------------------------------------------------------------------
module subroutine thermal_partition(ce)

  integer, intent(in) :: ce

  real(pREAL) :: T, dot_T
  integer :: co


  T     = current(material_ID_homogenization(ce))%T(material_entry_homogenization(ce))
  dot_T = current(material_ID_homogenization(ce))%dot_T(material_entry_homogenization(ce))
  do co = 1, homogenization_Nconstituents(material_ID_homogenization(ce))
    call phase_thermal_setField(T,dot_T,co,ce)
  end do

end subroutine thermal_partition


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize thermal viscosity.
!--------------------------------------------------------------------------------------------------
module function homogenization_mu_T(ce) result(mu)

  integer, intent(in) :: ce
  real(pREAL) :: mu

  integer :: co


  mu = phase_mu_T(1,ce)*material_v(1,ce)
  do co = 2, homogenization_Nconstituents(material_ID_homogenization(ce))
    mu = mu + phase_mu_T(co,ce)*material_v(co,ce)
  end do

end function homogenization_mu_T


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize thermal conductivity.
!--------------------------------------------------------------------------------------------------
module function homogenization_K_T(ce) result(K)

  integer, intent(in) :: ce
  real(pREAL), dimension(3,3) :: K

  integer :: co


  K = phase_K_T(1,ce)*material_v(1,ce)
  do co = 2, homogenization_Nconstituents(material_ID_homogenization(ce))
    K = K + phase_K_T(co,ce)*material_v(co,ce)
  end do

end function homogenization_K_T


!--------------------------------------------------------------------------------------------------
!> @brief Homogenize heat generation rate.
!--------------------------------------------------------------------------------------------------
module function homogenization_f_T(ce) result(f)

  integer, intent(in) :: ce
  real(pREAL) :: f

  integer :: co


  f = phase_f_T(material_ID_phase(1,ce),material_entry_phase(1,ce))*material_v(1,ce)
  do co = 2, homogenization_Nconstituents(material_ID_homogenization(ce))
    f = f + phase_f_T(material_ID_phase(co,ce),material_entry_phase(co,ce))*material_v(co,ce)
  end do

end function homogenization_f_T


!--------------------------------------------------------------------------------------------------
!> @brief Set thermal field and its rate (T and dot_T).
!--------------------------------------------------------------------------------------------------
module subroutine homogenization_thermal_setField(T,dot_T)

  real(pREAL), dimension(:), intent(in) :: T, dot_T

  integer :: ho, en, ce


  do ce=max(lbound(T,1),lbound(dot_T,1)), min(ubound(T,1),ubound(dot_T,1))
    ho = material_ID_homogenization(ce)
    en = material_entry_homogenization(ce)
    current(ho)%T(en) = T(ce)
    current(ho)%dot_T(en) = dot_T(ce)
    call thermal_partition(ce)
  end do

end subroutine homogenization_thermal_setField


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine thermal_result(ho,group)

  integer,          intent(in) :: ho
  character(len=*), intent(in) :: group

  integer :: o


  associate(prm => param(ho))
    outputsLoop: do o = 1,size(prm%output)
      select case(trim(prm%output(o)))
        case('T')
          call result_writeDataset(current(ho)%T,group,'T','temperature','K')
      end select
    end do outputsLoop
  end associate

end subroutine thermal_result

end submodule thermal
