!--------------------------------------------------------------------------------------------------
!> @author Tilak Raj Pant, Indian Institute of Science
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!--------------------------------------------------------------------------------------------------
submodule(phase) electrical


  type :: tDataContainer                                                                            !< state management for the electrical fields
    real(pREAL), dimension(:,:), allocatable :: E, J
  end type tDataContainer

  type(tDataContainer), dimension(:), allocatable :: current

  type :: tParameters
    type(tPolynomial) :: sigma_11, sigma_33
    character(len=pSTRLEN), allocatable, dimension(:) :: output
  end type tParameters

  type(tParameters), allocatable, dimension(:) :: param

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize electrical properties and state variables from the material file (E-based solver).
!--------------------------------------------------------------------------------------------------
module subroutine electrical_init(phases)

  type(tDict), pointer :: phases

  integer                       :: ph, Nmembers
  type(tDict),      pointer     :: phase_dict, electrical_dict
  character(len=:), allocatable :: electrical_type


  print'(/,1x,a)', '<<<+-  phase:electrical init (E-based) -+>>>'

  allocate(param(size(phases)))
  allocate(current(size(phases)))

  do ph = 1, size(phases)
    Nmembers = count(material_ID_phase == ph)

    allocate(current(ph)%E(3, Nmembers), source = 0.0_pREAL)
    allocate(current(ph)%J(3, Nmembers), source = 0.0_pREAL)

    phase_dict => phases%get_dict(ph)

    if (phase_dict%contains('electrical')) then                                                  ! initialize electrical parameters
      electrical_dict => phase_dict%get_dict('electrical')
      print'(/,1x,a,1x,i0,a)', 'phase', ph, ': '//phases%key(ph)//' (electrical)'

      electrical_type = electrical_dict%get_asStr('type')

      associate(prm => param(ph))
        select case(electrical_type)
          case('iso')
            prm%sigma_11 = polynomial(electrical_dict,'sigma','T')
            prm%sigma_33 = prm%sigma_11
          case('aniso')
            prm%sigma_11 = polynomial(electrical_dict,'sigma_11','T')
            if (any(phase_lattice(ph) == ['hP','tI'])) &
              prm%sigma_33 = polynomial(electrical_dict,'sigma_33','T')
          case default
            call IO_error(200, ext_msg='electrical: '//electrical_type)
        end select
        prm%output = electrical_dict%get_as1dStr('output', defaultVal=emptyStrArray)
      end associate
    end if
  end do

end subroutine electrical_init


!--------------------------------------------------------------------------------------------------
!> @brief Main constitutive update routine for electricals. (Placeholder)
!> @details This routine is called by homogenization at every time step.
!> For a quasi-static model with no internal state variables, it performs no update.
!--------------------------------------------------------------------------------------------------
module function phase_electrical_constitutive(Delta_t,ph,en) result(status)

  real(pREAL), intent(in)  :: Delta_t
  integer,     intent(in)  :: ph, en
  integer(kind(STATUS_OK)) :: status


  status = STATUS_OK

end function phase_electrical_constitutive


!--------------------------------------------------------------------------------------------------
!> @brief Return 3x3 electrical conductivity tensor based on phase properties.
!--------------------------------------------------------------------------------------------------
module function phase_sigma(co,ce) result(sigma)

  integer, intent(in) :: co, ce

  real(pREAL), dimension(3,3) :: sigma
  real(pREAL) :: T
  integer     :: ph, en


  ph = material_ID_phase(co,ce)
  en = material_entry_phase(co,ce)

  associate(prm => param(ph))
    sigma = 0.0_pREAL
    T = thermal_T(ph,en)                                                                            ! current temperature from the thermal submodule

    sigma(1,1) = prm%sigma_11%at(T)                                                                 ! build the tensor components
    if (any(phase_lattice(ph) == ['hP','tI'])) then
      sigma(3,3) = prm%sigma_33%at(T)
    else
      sigma(3,3) = sigma(1,1)
    end if

    sigma = crystal_symmetrize_33(sigma,phase_lattice(ph))
    sigma = crystallite_push33ToRef(co,ce,sigma)
  end associate

end function phase_sigma


!----------------------------------------------------------------------------------------------
!< @brief Set Electric field and compute J (current density)
!----------------------------------------------------------------------------------------------
module subroutine phase_electrical_setField(E, co, ce)

  real(pREAL), dimension(3), intent(in) :: E
  integer,                   intent(in) :: co, ce

  integer :: ph, en
  real(pREAL), dimension(3,3) :: sigma
  real(pREAL), dimension(3,3) :: O


  ph = material_ID_phase(co, ce)
  en = material_entry_phase(co, ce)

  current(ph)%E(:, en) = E
  sigma = phase_sigma(co, ce)
  current(ph)%J(:, en) = matmul(sigma, E)

end subroutine phase_electrical_setField


!--------------------------------------------------------------------------------------------------
!> @brief Public "getter" function to provide the electric field.
!--------------------------------------------------------------------------------------------------
module function electrical_E(ph,en) result(E)

  integer, intent(in)       :: ph,en
  real(pREAL), dimension(3) :: E


  E = current(ph)%E(:,en)

end function electrical_E


!--------------------------------------------------------------------------------------------------
!> @brief Public "getter" function to provide the current density.
!--------------------------------------------------------------------------------------------------
module function electrical_J(ph,en) result(J)

  integer, intent(in)       :: ph,en
  real(pREAL), dimension(3) :: J


  J = current(ph)%J(:,en)

end function electrical_J


!--------------------------------------------------------------------------------------------------
!> @brief Write restart data for the electrical model. (Placeholder)
!--------------------------------------------------------------------------------------------------
module subroutine electrical_restartWrite(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer,        intent(in) :: ph

  ! No internal state to write for this model yet.

end subroutine electrical_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read restart data for the electrical model. (Placeholder)
!--------------------------------------------------------------------------------------------------
module subroutine electrical_restartRead(groupHandle,ph)

  integer(HID_T), intent(in) :: groupHandle
  integer,        intent(in) :: ph

  ! No internal state to read for this model yet.

end subroutine electrical_restartRead


!--------------------------------------------------------------------------------------------------
!> @brief Forward the state for the electrical model. (Placeholder)
!--------------------------------------------------------------------------------------------------
module subroutine electrical_forward()

  ! No internal state to forward for this model yet.

end subroutine electrical_forward


!--------------------------------------------------------------------------------------------------
!> @brief Write electrical results to the HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine electrical_result(group,ph)

  character(len=*), intent(in) :: group
  integer,          intent(in) :: ph

  integer :: ou


  if (.not. allocated(param(ph)%output)) return

  call result_closeGroup(result_addGroup(group//'electrical'))

  do ou = 1, size(param(ph)%output)
    select case(trim(param(ph)%output(ou)))
      case('E')
        call result_writeDataset(current(ph)%E, group//'electrical', 'E', 'electric_field', 'V/m')
      case('J')
        call result_writeDataset(current(ph)%J, group//'electrical', 'J', 'current_density', 'A/m^2')
    end select
  end do

end subroutine electrical_result

end submodule electrical
