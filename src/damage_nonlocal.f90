!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for non-locally evolving damage field
!--------------------------------------------------------------------------------------------------
module damage_nonlocal
  use prec
  use material
  use config
  use numerics
  use crystallite
  use lattice
  use source_damage_isoBrittle
  use source_damage_isoDuctile
  use source_damage_anisoBrittle
  use source_damage_anisoDuctile
  use results

  implicit none
  private

  type :: tParameters
    character(len=pStringLen), allocatable, dimension(:) :: &
      output
  end type tParameters

  type(tparameters),             dimension(:), allocatable :: &
    param

  public :: &
    damage_nonlocal_init, &
    damage_nonlocal_getSourceAndItsTangent, &
    damage_nonlocal_getDiffusion, &
    damage_nonlocal_getMobility, &
    damage_nonlocal_putNonLocalDamage, &
    damage_nonlocal_results

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_init

  integer :: Ninstance,NofMyHomog,h

  write(6,'(/,a)') ' <<<+-  damage_'//DAMAGE_nonlocal_label//' init  -+>>>'; flush(6)

  Ninstance = count(damage_type == DAMAGE_nonlocal_ID)
  allocate(param(Ninstance))

  do h = 1, size(config_homogenization)
    if (damage_type(h) /= DAMAGE_NONLOCAL_ID) cycle
    associate(prm => param(damage_typeInstance(h)),config => config_homogenization(h))

    prm%output = config%getStrings('(output)',defaultVal=emptyStringArray)

    NofMyHomog = count(material_homogenizationAt == h)
    damageState(h)%sizeState = 1
    allocate(damageState(h)%state0   (1,NofMyHomog), source=damage_initialPhi(h))
    allocate(damageState(h)%subState0(1,NofMyHomog), source=damage_initialPhi(h))
    allocate(damageState(h)%state    (1,NofMyHomog), source=damage_initialPhi(h))

    nullify(damageMapping(h)%p)
    damageMapping(h)%p => material_homogenizationMemberAt
    deallocate(damage(h)%p)
    damage(h)%p => damageState(h)%state(1,:)

    end associate
  enddo

end subroutine damage_nonlocal_init


!--------------------------------------------------------------------------------------------------
!> @brief  calculates homogenized damage driving forces
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    phi
  integer :: &
    phase, &
    grain, &
    source, &
    constituent
  real(pReal) :: &
    phiDot, dPhiDot_dPhi, localphiDot, dLocalphiDot_dPhi

  phiDot = 0.0_pReal
  dPhiDot_dPhi = 0.0_pReal
  do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
    phase = material_phaseAt(grain,el)
    constituent = material_phasememberAt(grain,ip,el)
    do source = 1, phase_Nsources(phase)
      select case(phase_source(source,phase))
        case (SOURCE_damage_isoBrittle_ID)
         call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

        case (SOURCE_damage_isoDuctile_ID)
         call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

        case (SOURCE_damage_anisoBrittle_ID)
         call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

        case (SOURCE_damage_anisoDuctile_ID)
         call source_damage_anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

        case default
         localphiDot = 0.0_pReal
         dLocalphiDot_dPhi = 0.0_pReal

      end select
      phiDot = phiDot + localphiDot
      dPhiDot_dPhi = dPhiDot_dPhi + dLocalphiDot_dPhi
    enddo
  enddo

  phiDot = phiDot/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)
  dPhiDot_dPhi = dPhiDot_dPhi/real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)

end subroutine damage_nonlocal_getSourceAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns homogenized non local damage diffusion tensor in reference configuration
!--------------------------------------------------------------------------------------------------
function damage_nonlocal_getDiffusion(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal), dimension(3,3) :: &
    damage_nonlocal_getDiffusion
  integer :: &
    homog, &
    grain

  homog  = material_homogenizationAt(el)
  damage_nonlocal_getDiffusion = 0.0_pReal
  do grain = 1, homogenization_Ngrains(homog)
    damage_nonlocal_getDiffusion = damage_nonlocal_getDiffusion + &
      crystallite_push33ToRef(grain,ip,el,lattice_DamageDiffusion(1:3,1:3,material_phaseAt(grain,el)))
  enddo

  damage_nonlocal_getDiffusion = &
    charLength**2*damage_nonlocal_getDiffusion/real(homogenization_Ngrains(homog),pReal)

end function damage_nonlocal_getDiffusion


!--------------------------------------------------------------------------------------------------
!> @brief Returns homogenized nonlocal damage mobility
!--------------------------------------------------------------------------------------------------
real(pReal) function damage_nonlocal_getMobility(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  integer :: &
    ipc

  damage_nonlocal_getMobility = 0.0_pReal

  do ipc = 1, homogenization_Ngrains(material_homogenizationAt(el))
    damage_nonlocal_getMobility = damage_nonlocal_getMobility + lattice_DamageMobility(material_phaseAt(ipc,el))
  enddo

  damage_nonlocal_getMobility = damage_nonlocal_getMobility/&
                                real(homogenization_Ngrains(material_homogenizationAt(el)),pReal)

end function damage_nonlocal_getMobility


!--------------------------------------------------------------------------------------------------
!> @brief updated nonlocal damage field with solution from damage phase field PDE
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_putNonLocalDamage(phi,ip,el)

  integer, intent(in) :: &
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in) :: &
    phi
  integer :: &
    homog, &
    offset

  homog  = material_homogenizationAt(el)
  offset = damageMapping(homog)%p(ip,el)
  damage(homog)%p(offset) = phi

end subroutine damage_nonlocal_putNonLocalDamage


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_results(homog,group)

  integer,          intent(in) :: homog
  character(len=*), intent(in) :: group

  integer :: o

  associate(prm => param(damage_typeInstance(homog)))
  outputsLoop: do o = 1,size(prm%output)
    select case(prm%output(o))
      case ('damage')
        call results_writeDataset(group,damage(homog)%p,'phi',&
                                  'damage indicator','-')
    end select
  enddo outputsLoop
  end associate

end subroutine damage_nonlocal_results

end module damage_nonlocal
