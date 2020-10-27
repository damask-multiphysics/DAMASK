!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for thermal source due to plastic dissipation
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(constitutive:constitutive_thermal) source_thermal_dissipation

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
module function source_thermal_dissipation_init(source_length) result(mySources)

  integer, intent(in)                  :: source_length  
  logical, dimension(:,:), allocatable :: mySources

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src 
  integer :: Ninstances,sourceOffset,Nconstituents,p

  print'(/,a)', ' <<<+-  source_thermal_dissipation init  -+>>>'

  mySources = source_active('thermal_dissipation',source_length)
  Ninstances = count(mySources)
  print'(a,i2)', ' # instances: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstances))
  allocate(source_thermal_dissipation_offset  (phases%length), source=0)
  allocate(source_thermal_dissipation_instance(phases%length), source=0)

  do p = 1, phases%length
    phase => phases%get(p) 
    if(count(mySources(:,p)) == 0) cycle
    if(any(mySources(:,p))) source_thermal_dissipation_instance(p) = count(mySources(:,1:p))
    sources => phase%get('source')
    do sourceOffset = 1, sources%length
      if(mySources(sourceOffset,p)) then
        source_thermal_dissipation_offset(p) = sourceOffset
        associate(prm  => param(source_thermal_dissipation_instance(p)))

        src => sources%get(sourceOffset) 
        prm%kappa = src%get_asFloat('kappa')
        Nconstituents = count(material_phaseAt==p) * discretization_nIPs
        call constitutive_allocateState(sourceState(p)%p(sourceOffset),Nconstituents,0,0,0)

        end associate
      endif
    enddo
  enddo


end function source_thermal_dissipation_init


!--------------------------------------------------------------------------------------------------
!> @brief Ninstancess dissipation rate
!--------------------------------------------------------------------------------------------------
module subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDot_dT, Tstar, Lp, phase)

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

end submodule source_thermal_dissipation
