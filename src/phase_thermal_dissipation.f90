!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:thermal) dissipation

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      kappa                                                                                         !< TAYLOR-QUINNEY factor
  end type tParameters

  type(tParameters), dimension(:),   allocatable :: param                                           !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function dissipation_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length
  logical, dimension(:,:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    thermal, &
    src
  class(tList), pointer :: &
    sources
  integer :: so,Nmembers,ph


  mySources = thermal_active('dissipation',source_length)
  if (count(mySources) == 0) return
  print'(/,1x,a)', '<<<+-  phase:thermal:dissipation init  -+>>>'
  print'(/,a,i2)', ' # phases: ',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))

  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    if (count(mySources(:,ph)) == 0) cycle !ToDo: error if > 1
    thermal => phase%get_dict('thermal')
    sources => thermal%get_list('source')
    do so = 1, sources%length
      if (mySources(so,ph)) then
        associate(prm  => param(ph))
          src => sources%get_dict(so)

          prm%kappa = src%get_asFloat('kappa')
          Nmembers = count(material_phaseID == ph)
          call phase_allocateState(thermalState(ph)%p(so),Nmembers,0,0,0)

        end associate
      end if
    end do
  end do


end function dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief Ninstancess dissipation rate
!--------------------------------------------------------------------------------------------------
module function dissipation_f_T(ph,en) result(f_T)

  integer, intent(in) :: ph, en
  real(pReal) :: &
    f_T


  associate(prm => param(ph))
    f_T = prm%kappa*sum(abs(mechanical_S(ph,en)*mechanical_L_p(ph,en)))
  end associate

end function dissipation_f_T

end submodule dissipation
