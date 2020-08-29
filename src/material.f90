!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Parses material config file, either solverJobName.materialConfig or material.config
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

  character(len=pStringLen),    public, protected, allocatable, dimension(:) :: &
    material_name_phase, &                                                                          !< name of each phase
    material_name_homogenization                                                                    !< name of each homogenization
 
  integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
    thermal_type                                                                                    !< thermal transport model
  integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
    damage_type                                                                                     !< nonlocal damage model
  integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
    homogenization_type                                                                             !< type of each homogenization

  integer, public, protected :: &
    material_Nhomogenization                                                                        !< number of homogenizations

  integer, public, protected :: &
    homogenization_maxNgrains                                                                       !< max number of grains in any USED homogenization

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_Ngrains, &                                                                       !< number of grains in each homogenization
    homogenization_typeInstance, &                                                                  !< instance of particular type of each homogenization
    thermal_typeInstance, &                                                                         !< instance of particular type of each thermal transport
    damage_typeInstance                                                                             !< instance of particular type of each nonlocal damage

  real(pReal), dimension(:), allocatable, public, protected :: &
    thermal_initialT, &                                                                             !< initial temperature per each homogenization
    damage_initialPhi                                                                               !< initial damage per each homogenization

  integer, dimension(:),     allocatable, public, protected :: &                                    ! (elem)
    material_homogenizationAt                                                                       !< homogenization ID of each element (copy of discretization_homogenizationAt)
  integer, dimension(:,:),   allocatable, public, target :: &                                       ! (ip,elem) ToDo: ugly target for mapping hack
    material_homogenizationMemberAt                                                                 !< position of the element within its homogenization instance
  integer, dimension(:,:), allocatable, public, protected :: &                                      ! (constituent,elem)
    material_phaseAt                                                                                !< phase ID of each element
  integer, dimension(:,:,:), allocatable, public, protected :: &                                    ! (constituent,elem)
    material_phaseMemberAt                                                                          !< position of the element within its phase instance

  type(tState),        allocatable, dimension(:), public :: &
    homogState, &
    thermalState, &
    damageState

  type(Rotation), dimension(:,:,:), allocatable, public, protected :: &
    material_orientation0                                                                           !< initial orientation of each grain,IP,element

  integer, dimension(:), allocatable, private :: &
    microstructure_Nconstituents                                                                    !< number of constituents in each microstructure



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

  integer            :: ph, myHomog
  class(tNode), pointer :: &
    debug_material, &                                                                               ! pointer to material debug options
    phases, &
    material_homogenization
  character(len=pStringLen) :: sectionName
  
  write(6,'(/,a)') ' <<<+-  material init  -+>>>'; flush(6)

  phases => material_root%get('phase')
  allocate(material_name_phase(phases%length))
  do ph = 1, phases%length
    write(sectionName,'(i0,a)') ph,'_'
    material_name_phase(ph) = trim(adjustl(sectionName))//phases%getKey(ph)    !ToDO: No reason to do. Update damage tests 
  enddo
  
  material_homogenization => material_root%get('homogenization')
  allocate(material_name_homogenization(material_homogenization%length))
  do myHomog = 1, material_homogenization%length
    write(sectionName,'(i0,a)') myHomog,'_' 
    material_name_homogenization(myHomog) = trim(adjustl(sectionName))//material_homogenization%getKey(myHomog)
  enddo

  debug_material => debug_root%get('material',defaultVal=emptyList)

  call material_parseMicrostructure()
  if (debug_material%contains('basic')) write(6,'(a)') ' Microstructure parsed'; flush(6)

  call material_parseHomogenization()
  if (debug_material%contains('basic')) write(6,'(a)') ' Homogenization parsed'; flush(6)


  if(homogenization_maxNgrains > size(material_phaseAt,1)) call IO_error(148)

  allocate(homogState      (material_Nhomogenization))
  allocate(thermalState    (material_Nhomogenization))
  allocate(damageState     (material_Nhomogenization))

  allocate(thermalMapping  (material_Nhomogenization))
  allocate(damageMapping   (material_Nhomogenization))

  allocate(temperature     (material_Nhomogenization))
  allocate(damage          (material_Nhomogenization))

  allocate(temperatureRate (material_Nhomogenization))


  if (.not. restart) then
    call results_openJobFile
    call results_mapping_constituent(material_phaseAt,material_phaseMemberAt,material_name_phase)
    call results_mapping_materialpoint(material_homogenizationAt,material_homogenizationMemberAt,material_name_homogenization)
    call results_closeJobFile
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN DEPRECATED
  allocate(mappingHomogenizationConst(  discretization_nIP,discretization_nElem),source=1)

! hack needed to initialize field values used during constitutive initialization
  do myHomog = 1,material_Nhomogenization
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
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization

  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogMech, &
    homogThermal, &
    homogDamage

  integer :: h

  material_homogenization => material_root%get('homogenization')
  material_Nhomogenization = material_homogenization%length

  allocate(homogenization_type(material_Nhomogenization),           source=HOMOGENIZATION_undefined_ID)
  allocate(thermal_type(material_Nhomogenization),                  source=THERMAL_isothermal_ID)
  allocate(damage_type (material_Nhomogenization),                  source=DAMAGE_none_ID)
  allocate(homogenization_typeInstance(material_Nhomogenization),   source=0)
  allocate(thermal_typeInstance(material_Nhomogenization),          source=0)
  allocate(damage_typeInstance(material_Nhomogenization),           source=0)
  allocate(homogenization_Ngrains(material_Nhomogenization),        source=0)
  allocate(thermal_initialT(material_Nhomogenization),              source=300.0_pReal)
  allocate(damage_initialPhi(material_Nhomogenization),             source=1.0_pReal)

  do h=1, material_Nhomogenization
    homog => material_homogenization%get(h)
    homogMech => homog%get('mech')
    select case (homogMech%get_asString('type'))
      case('none')
        homogenization_type(h) = HOMOGENIZATION_NONE_ID
        homogenization_Ngrains(h) = 1
      case('isostrain')
        homogenization_type(h) = HOMOGENIZATION_ISOSTRAIN_ID
        homogenization_Ngrains(h) = homogMech%get_asInt('N_constituents')
      case('RGC')
        homogenization_type(h) = HOMOGENIZATION_RGC_ID
        homogenization_Ngrains(h) = homogMech%get_asInt('N_constituents')
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

  do h=1, material_Nhomogenization
    homogenization_typeInstance(h)  = count(homogenization_type(1:h) == homogenization_type(h))
    thermal_typeInstance(h)         = count(thermal_type       (1:h) == thermal_type       (h))
    damage_typeInstance(h)          = count(damage_type        (1:h) == damage_type        (h))
  enddo

  homogenization_maxNgrains = maxval(homogenization_Ngrains)


end subroutine material_parseHomogenization


!--------------------------------------------------------------------------------------------------
!> @brief parses the microstructure part in the material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_parseMicrostructure

  class(tNode), pointer :: microstructure, &                                                        !> pointer to microstructure list
                           constituentsInMicrostructure, &                                          !> pointer to a microstructure list item
                           constituents, &                                                          !> pointer to constituents list 
                           constituent, &                                                           !> pointer to each constituent 
                           phases, &
                           homogenization

  integer, dimension(:), allocatable :: &
   CounterPhase, &
   CounterHomogenization
 
 
  real(pReal), dimension(:,:), allocatable :: &
    microstructure_fraction                                                                         !< vol fraction of each constituent in microstrcuture 

  integer :: &
    e, &
    i, &
    m, &
    c, &
    microstructure_maxNconstituents

  real(pReal), dimension(4) :: phase_orientation

  homogenization => material_root%get('homogenization') 
  phases         => material_root%get('phase')
  microstructure => material_root%get('microstructure')
  allocate(microstructure_Nconstituents(microstructure%length), source = 0)
 
  if(any(discretization_microstructureAt > microstructure%length)) &
   call IO_error(155,ext_msg='More microstructures in geometry than sections in material.yaml')

  do m = 1, microstructure%length
    constituentsInMicrostructure => microstructure%get(m)
    constituents => constituentsInMicrostructure%get('constituents')
    microstructure_Nconstituents(m) = constituents%length
  enddo
 
  microstructure_maxNconstituents = maxval(microstructure_Nconstituents)
  allocate(microstructure_fraction(microstructure_maxNconstituents,microstructure%length), source =0.0_pReal)
  allocate(material_phaseAt(microstructure_maxNconstituents,discretization_nElem), source =0)
  allocate(material_orientation0(microstructure_maxNconstituents,discretization_nIP,discretization_nElem))
  allocate(material_homogenizationAt(discretization_nElem))
  allocate(material_homogenizationMemberAt(discretization_nIP,discretization_nElem),source=0)
  allocate(material_phaseMemberAt(microstructure_maxNconstituents,discretization_nIP,discretization_nElem),source=0)

  allocate(CounterPhase(phases%length),source=0)
  allocate(CounterHomogenization(homogenization%length),source=0)

  do m = 1, microstructure%length
    constituentsInMicrostructure => microstructure%get(m)
    constituents => constituentsInMicrostructure%get('constituents')
    do c = 1, constituents%length
      constituent => constituents%get(c)
      microstructure_fraction(c,m) = constituent%get_asFloat('fraction')
    enddo
   if (dNeq(sum(microstructure_fraction(:,m)),1.0_pReal)) call IO_error(153,ext_msg='constituent') 
  enddo

  do e = 1, discretization_nElem
    do i = 1, discretization_nIP
      constituentsInMicrostructure => microstructure%get(discretization_microstructureAt(e))
      constituents => constituentsInMicrostructure%get('constituents')
      do c = 1, constituents%length
        constituent => constituents%get(c)
        material_phaseAt(c,e) = phases%getIndex(constituent%get_asString('phase'))
        phase_orientation = constituent%get_asFloats('orientation')
        call material_orientation0(c,i,e)%fromQuaternion(phase_orientation)
      enddo
    enddo
  enddo
 
  do e = 1, discretization_nElem
    do i = 1, discretization_nIP
      constituentsInMicrostructure => microstructure%get(discretization_microstructureAt(e))  
      material_homogenizationAt(e) = homogenization%getIndex(constituentsInMicrostructure%get_asString('homogenization'))  
      CounterHomogenization(material_homogenizationAt(e)) = CounterHomogenization(material_homogenizationAt(e)) + 1
      material_homogenizationMemberAt(i,e) = CounterHomogenization(material_homogenizationAt(e))
    enddo
  enddo

  do e = 1, discretization_nElem
    do i = 1, discretization_nIP
      constituentsInMicrostructure => microstructure%get(discretization_microstructureAt(e))
      constituents => constituentsInMicrostructure%get('constituents')
      do c = 1, constituents%length
        CounterPhase(material_phaseAt(c,e)) = &
        CounterPhase(material_phaseAt(c,e)) + 1
        material_phaseMemberAt(c,i,e) = CounterPhase(material_phaseAt(c,e))
      enddo
    enddo
  enddo


end subroutine material_parseMicrostructure

 
end module material
