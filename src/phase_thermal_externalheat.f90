!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!> @brief material subroutine for variable heat source
!--------------------------------------------------------------------------------------------------
submodule(phase:thermal) externalheat


  integer,           dimension(:),   allocatable :: &
    source_thermal_externalheat_offset                                                              !< which source is my current thermal dissipation mechanism?

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    type(tTable) :: f
  end type tParameters

  type(tParameters), dimension(:), allocatable  :: param                                            !< containers of constitutive parameters (len Ninstances)


contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function externalheat_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length
  logical, dimension(:,:), allocatable :: mySources

  type(tDict), pointer :: &
    phases, &
    phase, &
    thermal, &
    src
  type(tList), pointer :: &
    sources
  integer :: so,Nmembers,ph


  mySources = thermal_active('externalheat',source_length)
  if (count(mySources) == 0) return
  print'(/,1x,a)', '<<<+-  phase:thermal:externalheat init  -+>>>'
  print'(/,a,i2)', ' # phases: ',count(mySources); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')
  allocate(param(phases%length))
  allocate(source_thermal_externalheat_offset  (phases%length), source=0)

  do ph = 1, phases%length
    phase => phases%get_dict(ph)
    if (count(mySources(:,ph)) == 0) cycle
    thermal => phase%get_dict('thermal')
    sources => thermal%get_list('source')
    do so = 1, sources%length
      if (mySources(so,ph)) then
        source_thermal_externalheat_offset(ph) = so
        associate(prm  => param(ph))
          src => sources%get_dict(so)

          prm%f = table(src,'t','f')

          Nmembers = count(material_phaseID == ph)
          call phase_allocateState(thermalState(ph)%p(so),Nmembers,1,1,0)
        end associate
      end if
    end do
  end do

end function externalheat_init


!--------------------------------------------------------------------------------------------------
!> @brief rate of change of state
!> @details state only contains current time to linearly interpolate given heat powers
!--------------------------------------------------------------------------------------------------
module subroutine externalheat_dotState(ph, en)

  integer, intent(in) :: &
    ph, &
    en

  integer :: &
    so

  so = source_thermal_externalheat_offset(ph)

  thermalState(ph)%p(so)%dotState(1,en) = 1.0_pReal                                                 ! state is current time

end subroutine externalheat_dotState


!--------------------------------------------------------------------------------------------------
!> @brief returns local heat generation rate
!--------------------------------------------------------------------------------------------------
module function externalheat_f_T(ph,en) result(f_T)

  integer, intent(in) :: &
    ph, &
    en
  real(pReal) :: &
    f_T

  integer :: &
    so


  so = source_thermal_externalheat_offset(ph)

  associate(prm => param(ph))
    f_T = prm%f%at(thermalState(ph)%p(so)%state(1,en))
  end associate

end function externalheat_f_T

end submodule externalheat
