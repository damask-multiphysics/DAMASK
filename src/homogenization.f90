! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization
!--------------------------------------------------------------------------------------------------
module homogenization
  use prec
  use math
  use constants
  use IO
  use config
  use material
  use phase
  use discretization
  use HDF5
  use HDF5_utilities
  use result
  use crystal

  implicit none(type,external)
  private

  type :: tState
    integer :: &
      sizeState        = 0                                                                          !< size of state
    ! http://stackoverflow.com/questions/3948210
    real(pREAL), pointer,     dimension(:,:), contiguous :: &                                       !< is basically an allocatable+target, but in a type needs to be pointer
      state0, &
      state
  end type

  enum, bind(c); enumerator :: &
    THERMAL_UNDEFINED_ID, &
    THERMAL_PASS_ID, &
    THERMAL_ISOTEMPERATURE_ID, &
    CHEMICAL_UNDEFINED_ID, &
    CHEMICAL_PASS_ID
  end enum
  integer(kind(THERMAL_UNDEFINED_ID)), dimension(:),   allocatable :: &
    thermal_type                                                                                    !< type of each homogenization

  type(tState),        allocatable, dimension(:), public :: &
    homogState, &
    damageState_h

  logical,             allocatable, dimension(:) :: &
    thermal_active, &
    chemical_active, &
    damage_active

!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
  real(pREAL),   dimension(:,:,:),     allocatable, public :: &
    homogenization_F                                                                                !< def grad of IP to be reached at end of FE increment
  real(pREAL),   dimension(:,:,:),     allocatable, public, protected :: &
    homogenization_P                                                                                !< first P--K stress of IP
  real(pREAL),   dimension(:,:,:,:,:), allocatable, public, protected :: &
    homogenization_dPdF                                                                             !< tangent of first P--K stress at IP

!--------------------------------------------------------------------------------------------------
  interface

    module subroutine mechanical_init()
    end subroutine mechanical_init

    module subroutine thermal_init()
    end subroutine thermal_init

    module subroutine chemical_init()
    end subroutine chemical_init

    module subroutine damage_init()
    end subroutine damage_init

    module subroutine mechanical_partition(subF,ce)
      real(pREAL), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ce
    end subroutine mechanical_partition

    module subroutine thermal_partition(ce)
      integer,     intent(in) :: ce
    end subroutine thermal_partition

    module subroutine damage_partition(ce)
      integer,     intent(in) :: ce
    end subroutine damage_partition

    module subroutine chemical_partition(Delta_t, ce)
      real(pREAL), intent(in) :: Delta_t
      integer,     intent(in) :: ce
    end subroutine chemical_partition

    module subroutine mechanical_homogenize(Delta_t,ce)
      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: ce                                                                                           !< cell
    end subroutine mechanical_homogenize

    module subroutine mechanical_result(group_base,ho)
      character(len=*), intent(in) :: group_base
      integer, intent(in)          :: ho
    end subroutine mechanical_result

    module subroutine damage_result(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine damage_result

    module subroutine thermal_result(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine thermal_result

    module subroutine chemical_result(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine chemical_result

    module function mechanical_updateState(subdt,subF,ce) result(doneAndHappy)
      real(pREAL), intent(in) :: &
        subdt                                                                                       !< current time step
      real(pREAL), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ce                                                                                          !< cell
      logical, dimension(2) :: doneAndHappy
    end function mechanical_updateState

    module function homogenization_thermal_active() result(active)
      logical :: active
    end function homogenization_thermal_active

    module function homogenization_mu_T(ce) result(mu)
      integer, intent(in) :: ce
      real(pREAL) :: mu
    end function homogenization_mu_T

    module function homogenization_K_T(ce) result(K)
      integer, intent(in) :: ce
      real(pREAL), dimension(3,3) :: K
    end function homogenization_K_T

    module function homogenization_f_T(ce) result(f)
      integer, intent(in) :: ce
      real(pREAL) :: f
    end function homogenization_f_T

    module subroutine homogenization_thermal_setField(T,dot_T)
      real(pREAL), dimension(:),  intent(in) :: T, dot_T
    end subroutine homogenization_thermal_setField

    module function homogenization_damage_active() result(active)
      logical :: active
    end function homogenization_damage_active

    module function homogenization_mu_phi(ce) result(mu)
      integer, intent(in) :: ce
      real(pREAL) :: mu
    end function homogenization_mu_phi

    module function homogenization_K_phi(ce) result(K)
      integer, intent(in) :: ce
      real(pREAL), dimension(3,3) :: K
    end function homogenization_K_phi

    module function homogenization_f_phi(phi,ce) result(f)
      integer, intent(in) :: ce
      real(pREAL), intent(in) :: phi
      real(pREAL) :: f
    end function homogenization_f_phi

    module subroutine homogenization_set_phi(phi)
      real(pREAL), dimension(:), intent(in) :: phi
    end subroutine homogenization_set_phi

    module function homogenization_chemical_active() result(active)
      logical :: active
    end function homogenization_chemical_active

    module function homogenization_composition(mu, Delta_t, ce) result(comp)
      real(pREAL), dimension(:), intent(in) :: mu
      real(pREAL), intent(in) :: Delta_t
      integer,     intent(in) :: ce
      real(pREAL), dimension(:), allocatable :: comp
    end function homogenization_composition

    module function homogenization_compositionTangent(mu, Delta_t, ce) result(comp_tangent)
      real(pREAL), dimension(:), intent(in) :: mu
      real(pREAL), intent(in) :: Delta_t
      integer,     intent(in) :: ce
      real(pREAL), dimension(:,:),allocatable :: comp_tangent
    end function homogenization_compositionTangent

    module function homogenization_mobility(ce) result(mobility)
      integer, intent(in) :: ce
      real(pREAL), dimension(:,:), allocatable :: mobility
    end function homogenization_mobility

    module subroutine homogenization_chemical_setField(mu, comp, Delta_t, ce)
      real(pREAL), dimension(:),  intent(in) :: mu
      real(pREAL), dimension(:),  intent(in) :: comp
      real(pREAL), intent(in) :: Delta_t
      integer, intent(in) :: ce
    end subroutine homogenization_chemical_setField

  end interface

  public ::  &
    homogenization_init, &
    homogenization_mechanical_response, &
    homogenization_thermal_response, &
    homogenization_thermal_active, &
    homogenization_chemical_active, &
    homogenization_mu_T, &
    homogenization_K_T, &
    homogenization_f_T, &
    homogenization_thermal_setfield, &
    homogenization_damage_active, &
    homogenization_mu_phi, &
    homogenization_K_phi, &
    homogenization_f_phi, &
    homogenization_set_phi, &
    homogenization_composition, &
    homogenization_compositionTangent, &
    homogenization_mobility, &
    homogenization_chemical_setField, &
    homogenization_forward, &
    homogenization_result, &
    homogenization_restartRead, &
    homogenization_restartWrite

contains


!--------------------------------------------------------------------------------------------------
!> @brief Module initialization.
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init()

  print'(/,1x,a)', '<<<+-  homogenization init  -+>>>'; flush(IO_STDOUT)


  allocate(homogState      (size(material_name_homogenization)))
  allocate(damageState_h   (size(material_name_homogenization)))
  call parseHomogenization()

  call mechanical_init()
  call thermal_init()
  call damage_init()
  call chemical_init()

end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
subroutine homogenization_mechanical_response(status,Delta_t,cell_start,cell_end)

  integer(kind(STATUS_OK)), intent(out) :: status
  real(pREAL), intent(in) :: Delta_t                                                                !< time increment
  integer, intent(in) :: &
    cell_start, cell_end
  integer :: &
    co, ce, ho, en
  logical :: &
    converged
  logical, dimension(2) :: &
    doneAndHappy


  status = STATUS_OK
  !$OMP PARALLEL DO PRIVATE(en,ho,co,converged,doneAndHappy)
  do ce = cell_start, cell_end

    en = material_entry_homogenization(ce)
    ho = material_ID_homogenization(ce)

    call phase_restore(ce,.false.) ! wrong name (is more a forward function)

    if (homogState(ho)%sizeState > 0)  homogState(ho)%state(:,en) = homogState(ho)%state0(:,en)
    if (damageState_h(ho)%sizeState > 0) damageState_h(ho)%state(:,en) = damageState_h(ho)%state0(:,en)
    call damage_partition(ce)

    doneAndHappy = [.false.,.true.]

    convergenceLooping: do while (status == STATUS_OK .and. .not. doneAndHappy(1))

      call mechanical_partition(homogenization_F(1:3,1:3,ce),ce)
      converged = all([(phase_mechanical_constitutive(Delta_t,co,ce) == STATUS_OK,co=1,homogenization_Nconstituents(ho))])
      if (converged) then
        doneAndHappy = mechanical_updateState(Delta_t,homogenization_F(1:3,1:3,ce),ce)
        converged = all(doneAndHappy)
      else
        doneAndHappy = [.true.,.false.]
      end if
    end do convergenceLooping
    if (.not. converged) then
      if (status == STATUS_OK) &
        call IO_warning(600,'mechanical response of cell', ce, 'on MPI rank', worldrank, emph=[2,4])
      status = STATUS_FAIL_PHASE_MECHANICAL
    end if
    converged = converged .and. all([(phase_damage_constitutive(Delta_t,co,ce)==STATUS_OK,co=1,homogenization_Nconstituents(ho))])

    if (.not. converged) then
      if (status == STATUS_OK) then
        call IO_warning(600,'damage response of cell', ce, 'on MPI rank', worldrank, emph=[2,4])
        status = STATUS_FAIL_PHASE_DAMAGE
      end if
    end if
  end do
  !$OMP END PARALLEL DO

  if (status /= STATUS_OK) return

  !$OMP PARALLEL DO PRIVATE(ho)
  do ce = cell_start, cell_end
      ho = material_ID_homogenization(ce)
      do co = 1, homogenization_Nconstituents(ho)
        call crystallite_orientations(co,ce)
      end do
      call mechanical_homogenize(Delta_t,ce)
  end do
  !$OMP END PARALLEL DO

end subroutine homogenization_mechanical_response


!--------------------------------------------------------------------------------------------------
!> @brief
!--------------------------------------------------------------------------------------------------
subroutine homogenization_thermal_response(status, &
                                           Delta_t,cell_start,cell_end)

  integer(kind(STATUS_OK)), intent(out) :: status
  real(pREAL), intent(in) :: Delta_t                                                                !< time increment
  integer, intent(in) :: &
    cell_start, cell_end

  integer :: &
    co, ce, ho


  status = STATUS_OK
  !$OMP PARALLEL DO PRIVATE(ho)
  do ce = cell_start, cell_end
    if (status /= STATUS_OK) continue
    ho = material_ID_homogenization(ce)
    do co = 1, homogenization_Nconstituents(ho)
      if (phase_thermal_constitutive(Delta_t,material_ID_phase(co,ce),material_entry_phase(co,ce)) /= STATUS_OK) then
        if (status == STATUS_OK) &
          call IO_warning(600,'thermal response of cell', ce, 'on MPI rank', worldrank, emph=[2,4])
        status = STATUS_FAIL_PHASE_THERMAL
      end if
    end do
  end do
  !$OMP END PARALLEL DO

end subroutine homogenization_thermal_response


!--------------------------------------------------------------------------------------------------
!> @brief writes homogenization results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_result()

  integer :: ho
  character(len=:), allocatable :: group_base,group


  call result_closeGroup(result_addGroup('current/homogenization/'))

  do ho=1,size(material_name_homogenization)
    group_base = 'current/homogenization/'//trim(material_name_homogenization(ho))
    call result_closeGroup(result_addGroup(group_base))

    call mechanical_result(group_base,ho)

    if (damage_active(ho)) then
      group = trim(group_base)//'/damage'
      call result_closeGroup(result_addGroup(group))
      call damage_result(ho,group)
    end if

    if (thermal_active(ho)) then
      group = trim(group_base)//'/thermal'
      call result_closeGroup(result_addGroup(group))
      call thermal_result(ho,group)
    end if

    if (chemical_active(ho)) then
      group = trim(group_base)//'/chemical'
      call result_closeGroup(result_addGroup(group))
      call chemical_result(ho,group)
    end if

 end do

end subroutine homogenization_result


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine homogenization_forward()

  integer :: ho


  do ho = 1, size(material_name_homogenization)
    homogState (ho)%state0 = homogState (ho)%state
    if (damageState_h(ho)%sizeState > 0) &
      damageState_h(ho)%state0 = damageState_h(ho)%state
  end do

end subroutine homogenization_forward


!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine homogenization_restartWrite(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ho


  groupHandle(1) = HDF5_addGroup(fileHandle,'homogenization')

  do ho = 1, size(material_name_homogenization)

    groupHandle(2) = HDF5_addGroup(groupHandle(1),material_name_homogenization(ho))

    call HDF5_write(homogState(ho)%state,groupHandle(2),'omega_mechanical') ! ToDo: should be done by mech

    if (damageState_h(ho)%sizeState > 0) &
      call HDF5_write(damageState_h(ho)%state,groupHandle(2),'omega_damage') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  end do

  call HDF5_closeGroup(groupHandle(1))

end subroutine homogenization_restartWrite


!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
subroutine homogenization_restartRead(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ho


  groupHandle(1) = HDF5_openGroup(fileHandle,'homogenization')

  do ho = 1, size(material_name_homogenization)

    groupHandle(2) = HDF5_openGroup(groupHandle(1),material_name_homogenization(ho))

    call HDF5_read(homogState(ho)%state0,groupHandle(2),'omega_mechanical') ! ToDo: should be done by mech

    if (damageState_h(ho)%sizeState > 0) &
      call HDF5_read(damageState_h(ho)%state0,groupHandle(2),'omega_damage') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  end do

  call HDF5_closeGroup(groupHandle(1))

end subroutine homogenization_restartRead


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine parseHomogenization

  type(tDict), pointer :: &
    material_homogenization, &
    homog, &
    homogThermal, &
    homogDamage, &
    homogChemical

  integer :: h

  material_homogenization => config_material%get_dict('homogenization')

  allocate(thermal_type(size(material_name_homogenization)),source=THERMAL_UNDEFINED_ID)
  allocate(thermal_active(size(material_name_homogenization)),source=.false.)
  allocate(damage_active(size(material_name_homogenization)),source=.false.)
  allocate(chemical_active(size(material_name_homogenization)),source=.false.)

  do h=1, size(material_name_homogenization)
    homog => material_homogenization%get_dict(h)

    if (homog%contains('thermal')) then
      homogThermal => homog%get_dict('thermal')
        select case (homogThermal%get_asStr('type'))
          case('pass')
            thermal_type(h) = THERMAL_PASS_ID
            thermal_active(h) = .true.
          case('isotemperature')
            thermal_type(h) = THERMAL_ISOTEMPERATURE_ID
            thermal_active(h) = .true.
          case default
            call IO_error(500,ext_msg=homogThermal%get_asStr('type'))
        end select
    end if

    if (homog%contains('damage')) then
      homogDamage => homog%get_dict('damage')
        select case (homogDamage%get_asStr('type'))
          case('pass')
            damage_active(h) = .true.
          case default
            call IO_error(500,ext_msg=homogDamage%get_asStr('type'))
        end select
    end if

    if (homog%contains('chemical')) then
      homogChemical => homog%get_dict('chemical')
        select case (homogChemical%get_asStr('type'))
          case('pass')
            chemical_active(h) = .true.
          case default
            call IO_error(500,ext_msg=homogChemical%get_asStr('type'))
        end select
    end if
  end do


end subroutine parseHomogenization


end module homogenization
