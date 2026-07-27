!--------------------------------------------------------------------------------------------------
!> @brief Constitutive model for a thermal source due to Joule heating.
!> @details This submodule subscribes to the electrical fields and provides a heat source
!> to the main thermal manager based on the local current and electric field.
!--------------------------------------------------------------------------------------------------

submodule(phase:thermal) source_Joule

  ! This submodule has no internal parameters, but the type is kept for structural consistency.
  type :: tParameters
     ! No parameters needed for this simple model yet.
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param

contains

!--------------------------------------------------------------------------------------------------
!> @brief Module initialization.
!> @details Reads the YAML file to see if the 'Joule' source is active for any phase.
!--------------------------------------------------------------------------------------------------
module function source_Joule_init(maxNsources) result(isMySource)

  integer, intent(in) :: maxNsources

  logical, dimension(:,:), allocatable :: isMySource
  type(tDict), pointer :: phases
  integer :: ph, Nmembers, so, Nsources


  isMySource = thermal_active('Joule', maxNsources)
  if (count(isMySource) == 0) return

  print'(/,1x,a)', '<<<+-  phase:thermal:source_Joule init  -+>>>'

  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))

  do ph = 1, size(phases)                                                                           ! no state variables for this model (size zero)
    Nsources = count(isMySource(:,ph))
    if (Nsources == 0) cycle
    if (Nsources > 1) call IO_error(600, ext_msg='Joule')                                           ! max once per phase

    Nmembers = count(material_ID_phase == ph)
    do so = 1, maxNsources
      if (isMySource(so,ph)) then
        call phase_allocateState(thermalState(ph)%p(so), Nmembers, 0, 0, 0)
        exit
      end if
    end do
  end do

end function source_Joule_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate the heat generation rate due to Joule heating (q = J . E).
!--------------------------------------------------------------------------------------------------
module function source_Joule_f_T(ph,en) result(f_T)

  integer, intent(in) :: ph, en
  real(pREAL) :: f_T


  f_T = dot_product(electrical_E(ph,en), electrical_J(ph,en))                                       ! f_t = E . J

end function source_Joule_f_T

end submodule source_Joule
