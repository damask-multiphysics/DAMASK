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
  use math
  use IO
  use rotations
  use discretization
  use YAML_types

  implicit none
  private

  type, public :: tRotationContainer
    type(tRotation), dimension(:), allocatable :: data
  end type tRotationContainer

  type, public :: tTensorContainer
    real(pReal), dimension(:,:,:), allocatable :: data
  end type tTensorContainer


  type(tRotationContainer), dimension(:), allocatable, public, protected :: material_O_0
  type(tTensorContainer),   dimension(:), allocatable, public, protected :: material_F_i_0

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_Nconstituents                                                                    !< number of grains in each homogenization
  integer, public, protected :: &
    homogenization_maxNconstituents                                                                 !< max number of grains in any homogenization

  character(len=:), public, protected, allocatable, dimension(:) :: &
    material_name_phase, &                                                                          !< name of each phase
    material_name_homogenization                                                                    !< name of each homogenization

  integer, dimension(:),   allocatable, public, protected :: &                                      ! (cell)
    material_homogenizationID, &                                                                    ! TODO: rename to material_ID_homogenization
    material_homogenizationEntry                                                                    ! TODO: rename to material_entry_homogenization
  integer, dimension(:,:), allocatable, public, protected :: &                                      ! (constituent,cell)
    material_phaseID, &                                                                             ! TODO: rename to material_ID_phase
    material_phaseEntry                                                                             ! TODO: rename to material_entry_phase

  real(pReal), dimension(:,:), allocatable, public, protected :: &
    material_v                                                                                      ! fraction

  public :: &
    material_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief Parse material configuration file (material.yaml).
!--------------------------------------------------------------------------------------------------
subroutine material_init(restart)

  logical, intent(in) :: restart


  print'(/,1x,a)', '<<<+-  material init  -+>>>'; flush(IO_STDOUT)


  call parse()
  print'(/,1x,a)', 'parsed material.yaml'


  if (.not. restart) then
    call results_openJobFile
    call results_mapping_phase(material_phaseID,material_phaseEntry,material_name_phase)
    call results_mapping_homogenization(material_homogenizationID,material_homogenizationEntry,material_name_homogenization)
    call results_closeJobFile
  end if

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

  real(pReal) :: v
  integer :: &
    el, ip, &
    ho, ph, &
    co, ce, &
    ma

  materials       => config_material%get('material')
  phases          => config_material%get('phase')
  homogenizations => config_material%get('homogenization')

  call sanityCheck(materials, homogenizations)

#if defined (__GFORTRAN__)
  material_name_phase          = getKeys(phases)
  material_name_homogenization = getKeys(homogenizations)
#else
  material_name_phase          = phases%Keys()
  material_name_homogenization = homogenizations%Keys()
#endif

  allocate(homogenization_Nconstituents(homogenizations%length))
  do ho=1, homogenizations%length
    homogenization => homogenizations%get(ho)
    homogenization_Nconstituents(ho) = homogenization%get_asInt('N_constituents')
  end do
  homogenization_maxNconstituents = maxval(homogenization_Nconstituents)

  allocate(counterPhase(phases%length),source=0)
  allocate(counterHomogenization(homogenizations%length),source=0)

  allocate(material_homogenizationID(discretization_Ncells),source=0)
  allocate(material_homogenizationEntry(discretization_Ncells),source=0)

  allocate(material_phaseID(homogenization_maxNconstituents,discretization_Ncells),source=0)
  allocate(material_phaseEntry(homogenization_maxNconstituents,discretization_Ncells),source=0)

  allocate(material_v(homogenization_maxNconstituents,discretization_Ncells),source=0.0_pReal)

  do el = 1, discretization_Nelems
    material => materials%get(discretization_materialAt(el))

    ho = homogenizations%getIndex(material%get_asString('homogenization'))
    do ip = 1, discretization_nIPs
      ce = (el-1)*discretization_nIPs + ip
      material_homogenizationID(ce) = ho
      counterHomogenization(ho) = counterHomogenization(ho) + 1
      material_homogenizationEntry(ce) = counterHomogenization(ho)
    end do

    constituents => material%get('constituents')
    do co = 1, constituents%length
      constituent => constituents%get(co)

      v = constituent%get_asFloat('v')

      ph = phases%getIndex(constituent%get_asString('phase'))
      do ip = 1, discretization_nIPs
        ce = (el-1)*discretization_nIPs + ip
        material_phaseID(co,ce) = ph
        counterPhase(ph) = counterPhase(ph) + 1
        material_phaseEntry(co,ce) = counterPhase(ph)
        material_v(co,ce) = v
      end do

    end do
    if (dNeq(sum(material_v(1:constituents%length,ce)),1.0_pReal,1.e-9_pReal)) &
      call IO_error(153,ext_msg='constituent')

  end do

  allocate(material_O_0(materials%length))
  allocate(material_F_i_0(materials%length))

  do ma = 1, materials%length
    material     => materials%get(ma)
    constituents => material%get('constituents')
    allocate(material_O_0(ma)%data(constituents%length))
    allocate(material_F_i_0(ma)%data(1:3,1:3,constituents%length))
    do co = 1, constituents%length
      constituent => constituents%get(co)
      call material_O_0(ma)%data(co)%fromQuaternion(constituent%get_as1dFloat('O',requiredSize=4))
      material_F_i_0(ma)%data(1:3,1:3,co) = constituent%get_as2dFloat('F_i',defaultVal=math_I3,requiredShape=[3,3])
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

  if (maxval(discretization_materialAt) > materials%length) &
    call IO_error(155,ext_msg='More materials requested than found in material.yaml')

  do m = 1, materials%length
    material => materials%get(m)
    constituents   => material%get('constituents')
    homogenization => homogenizations%get(material%get_asString('homogenization'))
    if (constituents%length /= homogenization%get_asInt('N_constituents')) call IO_error(148)
  end do

end subroutine sanityCheck

#if defined (__GFORTRAN__)
!--------------------------------------------------------------------------------------------------
!> @brief %keys() is broken on gfortran
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
  end do

  allocate(character(l)::getKeys(dict%length))
  do i=1, dict%length
    getKeys(i) = trim(temp(i))
  end do

end function getKeys
#endif

end module material
