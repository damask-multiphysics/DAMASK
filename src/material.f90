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
    THERMAL_ISOTHERMAL_ID, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONE_ID, &
    DAMAGE_NONLOCAL_ID, &
    HOMOGENIZATION_UNDEFINED_ID, &
    HOMOGENIZATION_NONE_ID, &
    HOMOGENIZATION_ISOSTRAIN_ID, &
    HOMOGENIZATION_RGC_ID
  end enum

  character(len=:), public, protected, allocatable, dimension(:) :: &
    material_name_phase, &                                                                          !< name of each phase
    material_name_homogenization                                                                    !< name of each homogenization

  integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
    thermal_type                                                                                    !< thermal transport model
  integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
    damage_type                                                                                     !< nonlocal damage model
  integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
    homogenization_type                                                                             !< type of each homogenization

  integer, public, protected :: &
    homogenization_maxNconstituents                                                                 !< max number of grains in any USED homogenization

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_Nconstituents, &                                                                 !< number of grains in each homogenization
    homogenization_typeInstance, &                                                                  !< instance of particular type of each homogenization
    thermal_typeInstance, &                                                                         !< instance of particular type of each thermal transport
    damage_typeInstance                                                                             !< instance of particular type of each nonlocal damage

  real(pReal), dimension(:), allocatable, public, protected :: &
    thermal_initialT                                                                                !< initial temperature per each homogenization

  integer, dimension(:),     allocatable, public, protected :: &                                    ! (elem)
    material_homogenizationAt, &                                                                    !< homogenization ID of each element
    material_homogenizationAt2, &                                                                   !< per cell
    material_homogenizationMemberAt2                                                                !< cell
  integer, dimension(:,:),   allocatable, public, protected :: &                                    ! (ip,elem)
    material_homogenizationMemberAt                                                                 !< position of the element within its homogenization instance
  integer, dimension(:,:),   allocatable, public, protected :: &                                    ! (constituent,elem)
    material_phaseAt, &                                                                             !< phase ID of each element
    material_phaseAt2, &                                                                            !< per constituent,cell
    material_phaseMemberAt2                                                                         !< per constituent, cell
  integer, dimension(:,:,:), allocatable, public, protected :: &                                    ! (constituent,IP,elem)
    material_phaseMemberAt                                                                          !< position of the element within its phase instance

  type(tState),        allocatable, dimension(:), public :: &
    homogState, &
    damageState_h

  type(Rotation), dimension(:,:,:), allocatable, public, protected :: &
    material_orientation0                                                                           !< initial orientation of each grain,IP,element

  public :: &
    material_init, &
    THERMAL_ISOTHERMAL_ID, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONE_ID, &
    DAMAGE_NONLOCAL_ID, &
    HOMOGENIZATION_NONE_ID, &
    HOMOGENIZATION_ISOSTRAIN_ID, &
    HOMOGENIZATION_RGC_ID, &
    material_parseHomogenization

contains

!--------------------------------------------------------------------------------------------------
!> @brief parses material configuration file
!--------------------------------------------------------------------------------------------------
subroutine material_init(restart)

  logical, intent(in) :: restart

  print'(/,a)', ' <<<+-  material init  -+>>>'; flush(IO_STDOUT)


  call material_parseMaterial
  print*, 'Material parsed'

  allocate(homogState      (size(material_name_homogenization)))
  allocate(damageState_h   (size(material_name_homogenization)))

  if (.not. restart) then
    call results_openJobFile
    call results_mapping_phase(material_phaseAt,material_phaseMemberAt,material_name_phase)
    call results_mapping_homogenization(material_homogenizationAt,material_homogenizationMemberAt,material_name_homogenization)
    call results_closeJobFile
  endif

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

  do h=1, size(material_name_homogenization)
    homog => material_homogenization%get(h)
    homogMech => homog%get('mechanics')
    select case (homogMech%get_asString('type'))
      case('pass')
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
          case('conduction')
            thermal_type(h) = THERMAL_conduction_ID
          case default
            call IO_error(500,ext_msg=homogThermal%get_asString('type'))
        end select
    endif

    if(homog%contains('damage')) then
      homogDamage => homog%get('damage')
        select case (homogDamage%get_asString('type'))
          case('none')
            damage_type(h) = DAMAGE_none_ID
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
    el, ip, co, &
    h, ce

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


  allocate(material_homogenizationAt2(discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_homogenizationMemberAt2(discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_phaseAt2(homogenization_maxNconstituents,discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_phaseMemberAt2(homogenization_maxNconstituents,discretization_nIPs*discretization_Nelems),source=0)

  allocate(material_orientation0(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems))

  do el = 1, discretization_Nelems
    material     => materials%get(discretization_materialAt(el))
    constituents => material%get('constituents')

    material_homogenizationAt(el) = homogenizations%getIndex(material%get_asString('homogenization'))
    do ip = 1, discretization_nIPs
      ce = (el-1)*discretization_nIPs + ip
      counterHomogenization(material_homogenizationAt(el)) = counterHomogenization(material_homogenizationAt(el)) + 1
      material_homogenizationMemberAt(ip,el)                = counterHomogenization(material_homogenizationAt(el))
      material_homogenizationAt2(ce)       = material_homogenizationAt(el)
      material_homogenizationMemberAt2(ce) = material_homogenizationMemberAt(ip,el)
    enddo

    frac = 0.0_pReal
    do co = 1, constituents%length
      constituent => constituents%get(co)
      frac = frac + constituent%get_asFloat('v')

      material_phaseAt(co,el) = phases%getIndex(constituent%get_asString('phase'))
      do ip = 1, discretization_nIPs
        ce = (el-1)*discretization_nIPs + ip
        counterPhase(material_phaseAt(co,el)) = counterPhase(material_phaseAt(co,el)) + 1
        material_phaseMemberAt(co,ip,el)       = counterPhase(material_phaseAt(co,el))

        material_phaseAt2(co,ce)       = material_phaseAt(co,el)
        material_phaseMemberAt2(co,ce) = material_phaseMemberAt(co,ip,el)
        call material_orientation0(co,ip,el)%fromQuaternion(constituent%get_asFloats('O',requiredSize=4)) ! should be done in crystallite
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
  character(len=:),          dimension(:), allocatable :: getKeys
  character(len=pStringLen), dimension(:), allocatable :: temp

  integer :: i,l

  allocate(temp(dict%length))
  l = 0
  do i=1, dict%length
    temp(i) = dict%getKey(i)
    l = max(len_trim(temp(i)),l)
  enddo

  allocate(character(l)::getKeys(dict%length))
  do i=1, dict%length
    getKeys(i) = trim(temp(i))
  enddo

end function getKeys

end module material
