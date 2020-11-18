!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Defines phase and homogenization
!--------------------------------------------------------------------------------------------------
module material
  use prec
  use math
  use config
  use results
  use IO
  use rotations
  use discretization

  implicit none
  private

  enum, bind(c); enumerator :: &
    ELASTICITY_UNDEFINED_ID, &
    ELASTICITY_HOOKE_ID, &
    PLASTICITY_UNDEFINED_ID, &
    PLASTICITY_NONE_ID, &
    PLASTICITY_ISOTROPIC_ID, &
    PLASTICITY_PHENOPOWERLAW_ID, &
    PLASTICITY_KINEHARDENING_ID, &
    PLASTICITY_DISLOTWIN_ID, &
    PLASTICITY_DISLOTUNGSTEN_ID, &
    PLASTICITY_NONLOCAL_ID, &
    SOURCE_UNDEFINED_ID ,&
    SOURCE_THERMAL_DISSIPATION_ID, &
    SOURCE_THERMAL_EXTERNALHEAT_ID, &
    SOURCE_DAMAGE_ISOBRITTLE_ID, &
    SOURCE_DAMAGE_ISODUCTILE_ID, &
    SOURCE_DAMAGE_ANISOBRITTLE_ID, &
    SOURCE_DAMAGE_ANISODUCTILE_ID, &
    KINEMATICS_UNDEFINED_ID ,&
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID, &
    STIFFNESS_DEGRADATION_UNDEFINED_ID, &
    STIFFNESS_DEGRADATION_DAMAGE_ID, &
    THERMAL_ISOTHERMAL_ID, &
    THERMAL_ADIABATIC_ID, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONE_ID, &
    DAMAGE_LOCAL_ID, &
    DAMAGE_NONLOCAL_ID, &
    HOMOGENIZATION_UNDEFINED_ID, &
    HOMOGENIZATION_NONE_ID, &
    HOMOGENIZATION_ISOSTRAIN_ID, &
    HOMOGENIZATION_RGC_ID
  end enum

  character(len=pStringLen), public, protected, allocatable, dimension(:) :: &
    material_name_phase, &                                                                          !< name of each phase
    material_name_homogenization                                                                    !< name of each homogenization

  integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
    thermal_type                                                                                    !< thermal transport model
  integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
    damage_type                                                                                     !< nonlocal damage model
  integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
    homogenization_type                                                                             !< type of each homogenization

  integer, public, protected :: &
    homogenization_maxNconstituents                                                                  !< max number of grains in any USED homogenization

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_Nconstituents, &                                                                  !< number of grains in each homogenization
    homogenization_typeInstance, &                                                                  !< instance of particular type of each homogenization
    thermal_typeInstance, &                                                                         !< instance of particular type of each thermal transport
    damage_typeInstance                                                                             !< instance of particular type of each nonlocal damage

  real(pReal), dimension(:), allocatable, public, protected :: &
    thermal_initialT, &                                                                             !< initial temperature per each homogenization
    damage_initialPhi                                                                               !< initial damage per each homogenization

  integer, dimension(:),     allocatable, public, protected :: &                                    ! (elem)
    material_homogenizationAt                                                                       !< homogenization ID of each element
  integer, dimension(:,:),   allocatable, public, target :: &                                       ! (ip,elem) ToDo: ugly target for mapping hack
    material_homogenizationMemberAt                                                                 !< position of the element within its homogenization instance
  integer, dimension(:,:),   allocatable, public, protected :: &                                    ! (constituent,elem)
    material_phaseAt                                                                                !< phase ID of each element
  integer, dimension(:,:,:), allocatable, public, protected :: &                                    ! (constituent,IP,elem)
    material_phaseMemberAt                                                                          !< position of the element within its phase instance

  type(tState),        allocatable, dimension(:), public :: &
    homogState, &
    thermalState, &
    damageState

  type(Rotation), dimension(:,:,:), allocatable, public, protected :: &
    material_orientation0                                                                           !< initial orientation of each grain,IP,element

! BEGIN DEPRECATED
  integer, dimension(:,:),   allocatable, private, target :: mappingHomogenizationConst             !< mapping from material points to offset in constant state/field
! END DEPRECATED

  type(tHomogMapping), allocatable, dimension(:), public :: &
    thermalMapping, &                                                                               !< mapping for thermal state/fields
    damageMapping                                                                                   !< mapping for damage state/fields

  type(group_float),  allocatable, dimension(:), public :: &
    temperature, &                                                                                  !< temperature field
    damage, &                                                                                       !< damage field
    temperatureRate                                                                                 !< temperature change rate field

  public :: &
    material_init, &
    ELASTICITY_UNDEFINED_ID, &
    ELASTICITY_HOOKE_ID, &
    PLASTICITY_UNDEFINED_ID, &
    PLASTICITY_NONE_ID, &
    PLASTICITY_ISOTROPIC_ID, &
    PLASTICITY_PHENOPOWERLAW_ID, &
    PLASTICITY_KINEHARDENING_ID, &
    PLASTICITY_DISLOTWIN_ID, &
    PLASTICITY_DISLOTUNGSTEN_ID, &
    PLASTICITY_NONLOCAL_ID, &
    SOURCE_UNDEFINED_ID ,&
    SOURCE_THERMAL_DISSIPATION_ID, &
    SOURCE_THERMAL_EXTERNALHEAT_ID, &
    SOURCE_DAMAGE_ISOBRITTLE_ID, &
    SOURCE_DAMAGE_ISODUCTILE_ID, &
    SOURCE_DAMAGE_ANISOBRITTLE_ID, &
    SOURCE_DAMAGE_ANISODUCTILE_ID, &
    KINEMATICS_UNDEFINED_ID ,&
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID, &
    STIFFNESS_DEGRADATION_UNDEFINED_ID, &
    STIFFNESS_DEGRADATION_DAMAGE_ID, &
    THERMAL_ISOTHERMAL_ID, &
    THERMAL_ADIABATIC_ID, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONE_ID, &
    DAMAGE_LOCAL_ID, &
    DAMAGE_NONLOCAL_ID, &
    HOMOGENIZATION_NONE_ID, &
    HOMOGENIZATION_ISOSTRAIN_ID, &
    HOMOGENIZATION_RGC_ID

contains

!--------------------------------------------------------------------------------------------------
!> @brief parses material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_init(restart)

  logical, intent(in) :: restart
  integer             :: myHomog

  print'(/,a)', ' <<<+-  material init  -+>>>'; flush(IO_STDOUT)


  call material_parseMaterial
  print*, 'Material parsed'

  call material_parseHomogenization
  print*, 'Homogenization parsed'


  allocate(homogState      (size(material_name_homogenization)))
  allocate(thermalState    (size(material_name_homogenization)))
  allocate(damageState     (size(material_name_homogenization)))

  allocate(thermalMapping  (size(material_name_homogenization)))
  allocate(damageMapping   (size(material_name_homogenization)))

  allocate(temperature     (size(material_name_homogenization)))
  allocate(damage          (size(material_name_homogenization)))

  allocate(temperatureRate (size(material_name_homogenization)))


  if (.not. restart) then
    call results_openJobFile
    call results_mapping_constituent(material_phaseAt,material_phaseMemberAt,material_name_phase)
    call results_mapping_homogenization(material_homogenizationAt,material_homogenizationMemberAt,material_name_homogenization)
    call results_closeJobFile
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN DEPRECATED
  allocate(mappingHomogenizationConst(  discretization_nIPs,discretization_Nelems),source=1)

! hack needed to initialize field values used during constitutive initialization
  do myHomog = 1, size(material_name_homogenization)
    thermalMapping     (myHomog)%p => mappingHomogenizationConst
    damageMapping      (myHomog)%p => mappingHomogenizationConst
    allocate(temperature     (myHomog)%p(1), source=thermal_initialT(myHomog))
    allocate(damage          (myHomog)%p(1), source=damage_initialPhi(myHomog))
    allocate(temperatureRate (myHomog)%p(1), source=0.0_pReal)
  enddo
! END DEPRECATED

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
! ToDo: This should be done in homogenization
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization

  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogMech, &
    homogThermal, &
    homogDamage

  integer :: h

  material_homogenization => config_material%get('homogenization')

  allocate(homogenization_type(size(material_name_homogenization)),           source=HOMOGENIZATION_undefined_ID)
  allocate(thermal_type(size(material_name_homogenization)),                  source=THERMAL_isothermal_ID)
  allocate(damage_type (size(material_name_homogenization)),                  source=DAMAGE_none_ID)
  allocate(homogenization_typeInstance(size(material_name_homogenization)),   source=0)
  allocate(thermal_typeInstance(size(material_name_homogenization)),          source=0)
  allocate(damage_typeInstance(size(material_name_homogenization)),           source=0)
  allocate(thermal_initialT(size(material_name_homogenization)),              source=300.0_pReal)
  allocate(damage_initialPhi(size(material_name_homogenization)),             source=1.0_pReal)

  do h=1, size(material_name_homogenization)
    homog => material_homogenization%get(h)
    homogMech => homog%get('mechanics')
    select case (homogMech%get_asString('type'))
      case('none')
        homogenization_type(h) = HOMOGENIZATION_NONE_ID
      case('isostrain')
        homogenization_type(h) = HOMOGENIZATION_ISOSTRAIN_ID
      case('RGC')
        homogenization_type(h) = HOMOGENIZATION_RGC_ID
      case default
        call IO_error(500,ext_msg=homogMech%get_asString('type'))
    end select

    homogenization_typeInstance(h) = count(homogenization_type==homogenization_type(h))

    if(homog%contains('thermal')) then
      homogThermal => homog%get('thermal')
        thermal_initialT(h) =  homogThermal%get_asFloat('T_0',defaultVal=300.0_pReal)

        select case (homogThermal%get_asString('type'))
          case('isothermal')
            thermal_type(h) = THERMAL_isothermal_ID
          case('adiabatic')
            thermal_type(h) = THERMAL_adiabatic_ID
          case('conduction')
            thermal_type(h) = THERMAL_conduction_ID
          case default
            call IO_error(500,ext_msg=homogThermal%get_asString('type'))
        end select
    endif

    if(homog%contains('damage')) then
      homogDamage => homog%get('damage')
        damage_initialPhi(h) =  homogDamage%get_asFloat('phi_0',defaultVal=1.0_pReal)
        select case (homogDamage%get_asString('type'))
          case('none')
            damage_type(h) = DAMAGE_none_ID
          case('local')
            damage_type(h) = DAMAGE_local_ID
          case('nonlocal')
            damage_type(h) = DAMAGE_nonlocal_ID
          case default
            call IO_error(500,ext_msg=homogDamage%get_asString('type'))
        end select
    endif
  enddo

  do h=1, size(material_name_homogenization)
    homogenization_typeInstance(h)  = count(homogenization_type(1:h) == homogenization_type(h))
    thermal_typeInstance(h)         = count(thermal_type       (1:h) == thermal_type       (h))
    damage_typeInstance(h)          = count(damage_type        (1:h) == damage_type        (h))
  enddo

end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the material part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMaterial

  class(tNode), pointer :: materials, &                                                             !> list of materials
                           material, &                                                              !> material definition
                           constituents, &                                                          !> list of constituents
                           constituent, &                                                           !> constituent definition
                           phases, &
                           homogenizations, &
                           homogenization

  integer, dimension(:), allocatable :: &
    counterPhase, &
    counterHomogenization

  real(pReal) :: &
    frac
  integer :: &
    e, i, c, &
    h

  materials       => config_material%get('material')
  phases          => config_material%get('phase')
  homogenizations => config_material%get('homogenization')

  call sanityCheck(materials, homogenizations)
  material_name_phase          = getKeys(phases)
  material_name_homogenization = getKeys(homogenizations)

  allocate(homogenization_Nconstituents(homogenizations%length))
  do h=1, homogenizations%length
    homogenization => homogenizations%get(h)
    homogenization_Nconstituents(h) = homogenization%get_asInt('N_constituents')
  enddo
  homogenization_maxNconstituents = maxval(homogenization_Nconstituents)

  allocate(counterPhase(phases%length),source=0)
  allocate(counterHomogenization(homogenizations%length),source=0)

  allocate(material_homogenizationAt(discretization_Nelems),source=0)
  allocate(material_homogenizationMemberAt(discretization_nIPs,discretization_Nelems),source=0)
  allocate(material_phaseAt(homogenization_maxNconstituents,discretization_Nelems),source=0)
  allocate(material_phaseMemberAt(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems),source=0)

  allocate(material_orientation0(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems))

  do e = 1, discretization_Nelems
    material     => materials%get(discretization_materialAt(e))
    constituents => material%get('constituents')

    material_homogenizationAt(e) = homogenizations%getIndex(material%get_asString('homogenization'))
    do i = 1, discretization_nIPs
      counterHomogenization(material_homogenizationAt(e)) = counterHomogenization(material_homogenizationAt(e)) + 1
      material_homogenizationMemberAt(i,e)                = counterHomogenization(material_homogenizationAt(e))
    enddo

    frac = 0.0_pReal
    do c = 1, constituents%length
      constituent => constituents%get(c)
      frac = frac + constituent%get_asFloat('fraction')

      material_phaseAt(c,e) = phases%getIndex(constituent%get_asString('phase'))
      do i = 1, discretization_nIPs
        counterPhase(material_phaseAt(c,e)) = counterPhase(material_phaseAt(c,e)) + 1
        material_phaseMemberAt(c,i,e)       = counterPhase(material_phaseAt(c,e))

        call material_orientation0(c,i,e)%fromQuaternion(constituent%get_asFloats('O',requiredSize=4)) ! should be done in crystallite
      enddo

    enddo
    if (dNeq(frac,1.0_pReal)) call IO_error(153,ext_msg='constituent')

  enddo

end subroutine material_parseMaterial


!--------------------------------------------------------------------------------------------------
!> @brief Check if material.yaml is consistent and contains sufficient # of materials
!--------------------------------------------------------------------------------------------------
subroutine sanityCheck(materials,homogenizations)

  class(tNode), intent(in) :: materials, &
                              homogenizations

  class(tNode), pointer :: material, &
                           homogenization, &
                           constituents
  integer :: m

  if(maxval(discretization_materialAt) > materials%length) &
    call IO_error(155,ext_msg='More materials requested than found in material.yaml')

  do m = 1, materials%length
    material => materials%get(m)
    constituents   => material%get('constituents')
    homogenization => homogenizations%get(material%get_asString('homogenization'))
    if(constituents%length /= homogenization%get_asInt('N_constituents')) call IO_error(148)
  enddo

end subroutine sanityCheck


!--------------------------------------------------------------------------------------------------
!> @brief Get all keys from a dictionary
!--------------------------------------------------------------------------------------------------
function getKeys(dict)

  class(tNode), intent(in) :: dict
  character(len=pStringLen), dimension(:), allocatable :: getKeys

  integer :: i

  allocate(getKeys(dict%length))
  do i=1, dict%length
    getKeys(i) = dict%getKey(i)
  enddo

end function getKeys

end module material
