!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
module source_thermal_dissipation
  use prec
  use debug
  use discretization
  use material
  use config

  implicit none
  private

  integer,           dimension(:),   allocatable :: &
    source_thermal_dissipation_offset, &                                                            !< which source is my current thermal dissipation mechanism?
    source_thermal_dissipation_instance                                                             !< instance of thermal dissipation source mechanism

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    real(pReal) :: &
      kappa
  end type tParameters

  type(tParameters), dimension(:),   allocatable :: param                                           !< containers of constitutive parameters (len Ninstance)


  public :: &
    source_thermal_dissipation_init, &
    source_thermal_dissipation_getRateAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_dissipation_init

  integer :: Ninstance,sourceOffset,NofMyPhase,p

  write(6,'(/,a)') ' <<<+-  source_'//SOURCE_thermal_dissipation_label//' init  -+>>>'; flush(6)

  Ninstance = count(phase_source == SOURCE_THERMAL_DISSIPATION_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(source_thermal_dissipation_offset  (size(config_phase)), source=0)
  allocate(source_thermal_dissipation_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    source_thermal_dissipation_instance(p) = count(phase_source(:,1:p) == SOURCE_THERMAL_DISSIPATION_ID)
    do sourceOffset = 1, phase_Nsources(p)
      if (phase_source(sourceOffset,p) == SOURCE_THERMAL_DISSIPATION_ID) then
        source_thermal_dissipation_offset(p) = sourceOffset
        exit
      endif
    enddo

    if (all(phase_source(:,p) /= SOURCE_THERMAL_DISSIPATION_ID)) cycle
    associate(prm => param(source_thermal_dissipation_instance(p)), &
              config => config_phase(p))

    prm%kappa = config%getFloat('dissipation_coldworkcoeff')

    NofMyPhase = count(material_phaseAt==p) * discretization_nIP
    call material_allocateSourceState(p,sourceOffset,NofMyPhase,0,0,0)

    end associate
  enddo

end subroutine source_thermal_dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief Ninstances dissipation rate
!--------------------------------------------------------------------------------------------------
subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDot_dT, Tstar, Lp, phase)

  integer, intent(in) :: &
    phase
  real(pReal),  intent(in), dimension(3,3) :: &
    Tstar
  real(pReal),  intent(in), dimension(3,3) :: &
    Lp

  real(pReal),  intent(out) :: &
    TDot, &
    dTDot_dT

  associate(prm => param(source_thermal_dissipation_instance(phase)))
  TDot = prm%kappa*sum(abs(Tstar*Lp))
  dTDot_dT = 0.0_pReal
  end associate

end subroutine source_thermal_dissipation_getRateAndItsTangent

end module source_thermal_dissipation
