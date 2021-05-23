!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Defines phase and homogenization
!--------------------------------------------------------------------------------------------------
module material
  use prec
  use config
  use results
  use IO
  use rotations
  use discretization
  use YAML_types

  implicit none
  private

  type :: tRotationContainer
    type(Rotation), dimension(:),  allocatable :: data
  end type

  type(tRotationContainer), dimension(:), allocatable :: material_orientation0

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_Nconstituents                                                                    !< number of grains in each homogenization
  integer, public, protected :: &
    homogenization_maxNconstituents                                                                 !< max number of grains in any homogenization

  character(len=:), public, protected, allocatable, dimension(:) :: &
    material_name_phase, &                                                                          !< name of each phase
    material_name_homogenization                                                                    !< name of each homogenization

  integer, dimension(:),     allocatable, public, protected :: &                                    ! (elem)
    material_homogenizationID, &                                                                    !< per cell TODO: material_ID_homogenization
    material_homogenizationEntry                                                                    !< per cell TODO: material_entry_homogenization
  integer, dimension(:,:),   allocatable, public, protected :: &                                    ! (constituent,elem)
    material_phaseAt, &                                                                             !< phase ID of each element TODO: remove
    material_phaseID, &                                                                             !< per (constituent,cell) TODO: material_ID_phase
    material_phaseEntry                                                                             !< per (constituent,cell) TODO: material_entry_phase
  integer, dimension(:,:,:), allocatable, public, protected :: &                                    ! (constituent,IP,elem)
    material_phaseMemberAt                                                                          !TODO: remove
  public :: &
    tRotationContainer, &
    material_orientation0, &
    material_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief Parse material configuration file (material.yaml).
!--------------------------------------------------------------------------------------------------
subroutine material_init(restart)

  logical, intent(in) :: restart


  print'(/,a)', ' <<<+-  material init  -+>>>'; flush(IO_STDOUT)


  call parse
  print*, 'parsed material.yaml'


  if (.not. restart) then
    call results_openJobFile
    call results_mapping_phase(material_phaseID,material_phaseEntry,material_name_phase)
    call results_mapping_homogenization(material_homogenizationID,material_homogenizationEntry,material_name_homogenization)
    call results_closeJobFile
  endif

end subroutine material_init


!--------------------------------------------------------------------------------------------------
!> @brief Parse material.yaml to get the global structure
!--------------------------------------------------------------------------------------------------
subroutine parse()

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
    el, ip, co, ma, &
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

  allocate(material_phaseAt(homogenization_maxNconstituents,discretization_Nelems),source=0)
  allocate(material_phaseMemberAt(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems),source=0)


  allocate(material_homogenizationID(discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_homogenizationEntry(discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_phaseID(homogenization_maxNconstituents,discretization_nIPs*discretization_Nelems),source=0)
  allocate(material_phaseEntry(homogenization_maxNconstituents,discretization_nIPs*discretization_Nelems),source=0)

  do el = 1, discretization_Nelems
    material     => materials%get(discretization_materialAt(el))
    constituents => material%get('constituents')

    do ip = 1, discretization_nIPs
      ce = (el-1)*discretization_nIPs + ip
      material_homogenizationID(ce) = homogenizations%getIndex(material%get_asString('homogenization'))
      counterHomogenization(material_homogenizationID(ce)) = counterHomogenization(material_homogenizationID(ce)) + 1
      material_homogenizationEntry(ce) = counterHomogenization(material_homogenizationID(ce))
    enddo

    frac = 0.0_pReal
    do co = 1, constituents%length
      constituent => constituents%get(co)
      frac = frac + constituent%get_asFloat('v')

      material_phaseAt(co,el) = phases%getIndex(constituent%get_asString('phase'))
      do ip = 1, discretization_nIPs
        ce = (el-1)*discretization_nIPs + ip
        counterPhase(material_phaseAt(co,el)) = counterPhase(material_phaseAt(co,el)) + 1
        material_phaseMemberAt(co,ip,el)      = counterPhase(material_phaseAt(co,el))
        material_phaseEntry(co,ce) = counterPhase(material_phaseAt(co,el))
        material_phaseID(co,ce) = material_phaseAt(co,el)
      enddo

    enddo
    if (dNeq(frac,1.0_pReal)) call IO_error(153,ext_msg='constituent')

  enddo

  allocate(material_orientation0(materials%length))

  do ma = 1, materials%length
    material     => materials%get(ma)
    constituents => material%get('constituents')
    allocate(material_orientation0(ma)%data(constituents%length))
    do co = 1, constituents%length
      constituent => constituents%get(co)
      call material_orientation0(ma)%data(co)%fromQuaternion(constituent%get_as1dFloat('O',requiredSize=4))
    enddo
 enddo

end subroutine parse


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
