!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:thermal) dissipation

  integer,           dimension(:),   allocatable :: &
    source_thermal_dissipation_offset, &                                                            !< which source is my current thermal dissipation mechanism?
    source_thermal_dissipation_instance                                                             !< instance of thermal dissipation source mechanism

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

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, thermal, &
    src
  integer :: Ninstances,sourceOffset,Nconstituents,p

  print'(/,a)', ' <<<+-  thermal_dissipation init  -+>>>'

  mySources = thermal_active('dissipation',source_length)

  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstances))
  allocate(source_thermal_dissipation_offset  (phases%length), source=0)
  allocate(source_thermal_dissipation_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p)
    if(any(mySources(:,p))) source_thermal_dissipation_instance(p) = count(mySources(:,1:p))
    if(count(mySources(:,p)) == 0) cycle
    thermal => phase%get('thermal')
    sources => thermal%get('source')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_thermal_dissipation_offset(p) = sourceOffset
        associate(prm  => param(source_thermal_dissipation_instance(p)))
        src => sources%get(sourceOffset)

        prm%kappa = src%get_asFloat('kappa')
        Nconstituents = count(material_phaseAt==p) * discretization_nIPs
        call constitutive_allocateState(thermalState(p)%p(sourceOffset),Nconstituents,0,0,0)

        end associate
      endif
    enddo
  enddo


end function dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief Ninstancess dissipation rate
!--------------------------------------------------------------------------------------------------
module subroutine thermal_dissipation_getRate(TDot, Tstar, Lp, phase)

  integer, intent(in) :: &
    phase
  real(pReal),  intent(in), dimension(3,3) :: &
    Tstar
  real(pReal),  intent(in), dimension(3,3) :: &
    Lp

  real(pReal),  intent(out) :: &
    TDot

  associate(prm => param(source_thermal_dissipation_instance(phase)))
  TDot = prm%kappa*sum(abs(Tstar*Lp))
  end associate

end subroutine thermal_dissipation_getRate

end submodule dissipation
