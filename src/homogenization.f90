!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization
!--------------------------------------------------------------------------------------------------
module homogenization
  use prec
  use IO
  use config
  use math
  use material
  use phase
  use discretization
  use HDF5_utilities
  use results
  use lattice

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

    type(tState),        allocatable, dimension(:), public :: &
    homogState, &
    damageState_h

  integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable :: &
    thermal_type                                                                                    !< thermal transport model
  integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable :: &
    damage_type                                                                                     !< nonlocal damage model

  type, private :: tNumerics_damage
    real(pReal) :: &
    charLength                                                                                      !< characteristic length scale for gradient problems
  end type tNumerics_damage

  type(tNumerics_damage), private :: &
    num_damage


  logical, public :: &
    terminallyIll = .false.                                                                         !< at least one material point is terminally ill

!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
  real(pReal),   dimension(:,:,:),     allocatable, public :: &
    homogenization_F0, &                                                                            !< def grad of IP at start of FE increment
    homogenization_F                                                                                !< def grad of IP to be reached at end of FE increment
  real(pReal),   dimension(:,:,:),     allocatable, public :: & !, protected :: &                   Issue with ifort
    homogenization_P                                                                                !< first P--K stress of IP
  real(pReal),   dimension(:,:,:,:,:), allocatable, public :: & !, protected ::  &
    homogenization_dPdF                                                                             !< tangent of first P--K stress at IP


!--------------------------------------------------------------------------------------------------
  type :: tNumerics
    integer :: &
      nMPstate                                                                                      !< materialpoint state loop limit
  end type tNumerics

  type(tNumerics) :: num

!--------------------------------------------------------------------------------------------------
  interface

    module subroutine mechanical_init(num_homog)
      class(tNode), pointer, intent(in) :: &
        num_homog                                                                                   !< pointer to mechanical homogenization numerics data
    end subroutine mechanical_init

    module subroutine thermal_init
    end subroutine thermal_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine mechanical_partition(subF,ce)
      real(pReal), intent(in), dimension(3,3) :: &
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

    module subroutine mechanical_homogenize(dt,ce)
     real(pReal), intent(in) :: dt
     integer, intent(in) :: &
       ce                                                                                           !< cell
    end subroutine mechanical_homogenize

    module subroutine mechanical_results(group_base,ho)
      character(len=*), intent(in) :: group_base
      integer, intent(in)          :: ho
    end subroutine mechanical_results

    module subroutine damage_results(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine damage_results

    module subroutine thermal_results(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine thermal_results

    module function mechanical_updateState(subdt,subF,ce) result(doneAndHappy)
      real(pReal), intent(in) :: &
        subdt                                                                                       !< current time step
      real(pReal), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ce                                                                                          !< cell
      logical, dimension(2) :: doneAndHappy
    end function mechanical_updateState

    module function homogenization_mu_T(ce) result(mu)
      integer, intent(in) :: ce
      real(pReal) :: mu
    end function homogenization_mu_T

    module function homogenization_K_T(ce) result(K)
      integer, intent(in) :: ce
      real(pReal), dimension(3,3) :: K
    end function homogenization_K_T

    module function homogenization_f_T(ce) result(f)
      integer, intent(in) :: ce
      real(pReal) :: f
    end function homogenization_f_T

    module subroutine homogenization_thermal_setField(T,dot_T, ce)
      integer, intent(in) :: ce
      real(pReal),   intent(in) :: T, dot_T
    end subroutine homogenization_thermal_setField

    module function homogenization_mu_phi(ce) result(mu)
      integer, intent(in) :: ce
      real(pReal) :: mu
    end function homogenization_mu_phi

    module function homogenization_K_phi(ce) result(K)
      integer, intent(in) :: ce
      real(pReal), dimension(3,3) :: K
    end function homogenization_K_phi

    module function homogenization_f_phi(phi,ce) result(f)
      integer, intent(in) :: ce
      real(pReal), intent(in) :: phi
      real(pReal) :: f
    end function homogenization_f_phi

    module subroutine homogenization_set_phi(phi,ce)
      integer, intent(in) :: ce
      real(pReal),   intent(in) :: &
        phi
    end subroutine homogenization_set_phi

  end interface

  public ::  &
    homogenization_init, &
    materialpoint_stressAndItsTangent, &
    homogenization_mu_T, &
    homogenization_K_T, &
    homogenization_f_T, &
    homogenization_thermal_setfield, &
    homogenization_mu_phi, &
    homogenization_K_phi, &
    homogenization_f_phi, &
    homogenization_set_phi, &
    homogenization_forward, &
    homogenization_results, &
    homogenization_restartRead, &
    homogenization_restartWrite, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONLOCAL_ID

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init()

  class (tNode) , pointer :: &
    num_homog, &
    num_homogGeneric

  print'(/,a)', ' <<<+-  homogenization init  -+>>>'; flush(IO_STDOUT)


  allocate(homogState      (size(material_name_homogenization)))
  allocate(damageState_h   (size(material_name_homogenization)))
  call material_parseHomogenization()

  num_homog        => config_numerics%get('homogenization',defaultVal=emptyDict)
  num_homogGeneric => num_homog%get('generic',defaultVal=emptyDict)

  num%nMPstate  = num_homogGeneric%get_asInt('nMPstate',defaultVal=10)
  if (num%nMPstate < 1) call IO_error(301,ext_msg='nMPstate')

  call mechanical_init(num_homog)
  call thermal_init()
  call damage_init()

end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(dt,FEsolving_execIP,FEsolving_execElem)

  real(pReal), intent(in) :: dt                                                                     !< time increment
  integer, dimension(2), intent(in) :: FEsolving_execElem, FEsolving_execIP
  integer :: &
    NiterationMPstate, &
    ip, &                                                                                            !< integration point number
    el, &                                                                                            !< element number
    co, ce, ho, en, ph
  logical :: &
    converged
  logical, dimension(2) :: &
    doneAndHappy

  !$OMP PARALLEL
  !$OMP DO PRIVATE(ce,en,ho,NiterationMPstate,converged,doneAndHappy)
  do el = FEsolving_execElem(1),FEsolving_execElem(2)

    do ip = FEsolving_execIP(1),FEsolving_execIP(2)
      ce = (el-1)*discretization_nIPs + ip
      en = material_homogenizationEntry(ce)
      ho = material_homogenizationID(ce)

      call phase_restore(ce,.false.) ! wrong name (is more a forward function)

      if(homogState(ho)%sizeState > 0)  homogState(ho)%state(:,en) = homogState(ho)%state0(:,en)
      if(damageState_h(ho)%sizeState > 0) damageState_h(ho)%state(:,en) = damageState_h(ho)%state0(:,en)
      call damage_partition(ce)

      doneAndHappy = [.false.,.true.]

      NiterationMPstate = 0
      convergenceLooping: do while (.not. (terminallyIll .or. doneAndHappy(1)) &
                                    .and. NiterationMPstate < num%nMPstate)
        NiterationMPstate = NiterationMPstate + 1

        call mechanical_partition(homogenization_F(1:3,1:3,ce),ce)
        converged = .true.
        do co = 1, homogenization_Nconstituents(ho)
          converged = converged .and. crystallite_stress(dt,co,ip,el)
        enddo

        if (converged) then
          doneAndHappy = mechanical_updateState(dt,homogenization_F(1:3,1:3,ce),ce)
          converged = all(doneAndHappy)
        else
          doneAndHappy = [.true.,.false.]
        endif

      enddo convergenceLooping
      if (.not. converged) then
        if (.not. terminallyIll) print*, ' Integration point ', ip,' at element ', el, ' terminally ill'
        terminallyIll = .true.
      endif
    enddo
  enddo
  !$OMP END DO

  if (.not. terminallyIll) then
    !$OMP DO PRIVATE(ho,ph,ce)
    do el = FEsolving_execElem(1),FEsolving_execElem(2)
      if (terminallyIll) continue
      do ip = FEsolving_execIP(1),FEsolving_execIP(2)
        ce = (el-1)*discretization_nIPs + ip
        ho = material_homogenizationID(ce)
        call thermal_partition(ce)
        do co = 1, homogenization_Nconstituents(ho)
          ph = material_phaseID(co,ce)
          if (.not. thermal_stress(dt,ph,material_phaseMemberAt(co,ip,el))) then
            if (.not. terminallyIll) &                                                              ! so first signals terminally ill...
              print*, ' Integration point ', ip,' at element ', el, ' terminally ill'
            terminallyIll = .true.                                                                  ! ...and kills all others
          endif
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO PRIVATE(ho,ce)
    elementLooping3: do el = FEsolving_execElem(1),FEsolving_execElem(2)
      IpLooping3: do ip = FEsolving_execIP(1),FEsolving_execIP(2)
        ce = (el-1)*discretization_nIPs + ip
        ho = material_homogenizationID(ce)
        do co = 1, homogenization_Nconstituents(ho)
          call crystallite_orientations(co,ip,el)
        enddo
        call mechanical_homogenize(dt,ce)
      enddo IpLooping3
    enddo elementLooping3
    !$OMP END DO
  else
    print'(/,a,/)', ' << HOMOG >> Material Point terminally ill'
  endif
  !$OMP END PARALLEL

end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief writes homogenization results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_results

  integer :: ho
  character(len=:), allocatable :: group_base,group


  call results_closeGroup(results_addGroup('current/homogenization/'))

  do ho=1,size(material_name_homogenization)
    group_base = 'current/homogenization/'//trim(material_name_homogenization(ho))
    call results_closeGroup(results_addGroup(group_base))

    call mechanical_results(group_base,ho)

    select case(damage_type(ho))
      case(DAMAGE_NONLOCAL_ID)
        group = trim(group_base)//'/damage'
        call results_closeGroup(results_addGroup(group))
        call damage_results(ho,group)
    end select

    select case(thermal_type(ho))
      case(THERMAL_CONDUCTION_ID)
        group = trim(group_base)//'/thermal'
        call results_closeGroup(results_addGroup(group))
        call thermal_results(ho,group)
    end select

 enddo

end subroutine homogenization_results


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine homogenization_forward

  integer :: ho


  do ho = 1, size(material_name_homogenization)
    homogState (ho)%state0 = homogState (ho)%state
    if(damageState_h(ho)%sizeState > 0) &
      damageState_h(ho)%state0 = damageState_h(ho)%state
  enddo

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

    call HDF5_write(homogState(ho)%state,groupHandle(2),'omega') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  enddo

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

    call HDF5_read(homogState(ho)%state0,groupHandle(2),'omega') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine homogenization_restartRead


!--------------------------------------------------------------------------------------------------
!> @brief parses the homogenization part from the material configuration
!--------------------------------------------------------------------------------------------------
subroutine material_parseHomogenization

  class(tNode), pointer :: &
    material_homogenization, &
    homog, &
    homogThermal, &
    homogDamage

  integer :: h

  material_homogenization => config_material%get('homogenization')

  allocate(thermal_type(size(material_name_homogenization)),        source=THERMAL_isothermal_ID)
  allocate(damage_type (size(material_name_homogenization)),        source=DAMAGE_none_ID)

  do h=1, size(material_name_homogenization)
    homog => material_homogenization%get(h)

    if (homog%contains('thermal')) then
      homogThermal => homog%get('thermal')
        select case (homogThermal%get_asString('type'))
          case('pass')
            thermal_type(h) = THERMAL_conduction_ID
          case default
            call IO_error(500,ext_msg=homogThermal%get_asString('type'))
        end select
    endif

    if (homog%contains('damage')) then
      homogDamage => homog%get('damage')
        select case (homogDamage%get_asString('type'))
          case('pass')
            damage_type(h) = DAMAGE_nonlocal_ID
          case default
            call IO_error(500,ext_msg=homogDamage%get_asString('type'))
        end select
    endif
  enddo

end subroutine material_parseHomogenization


end module homogenization
