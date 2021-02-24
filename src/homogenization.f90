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

  integer, dimension(:), allocatable, public, protected :: &
    homogenization_typeInstance, &                                                                  !< instance of particular type of each homogenization
    thermal_typeInstance, &                                                                         !< instance of particular type of each thermal transport
    damage_typeInstance                                                                             !< instance of particular type of each nonlocal damage

  real(pReal), dimension(:), allocatable, public, protected :: &
    thermal_initialT

  integer(kind(THERMAL_isothermal_ID)),       dimension(:),   allocatable, public, protected :: &
    thermal_type                                                                                    !< thermal transport model
  integer(kind(DAMAGE_none_ID)),              dimension(:),   allocatable, public, protected :: &
    damage_type                                                                                     !< nonlocal damage model
  integer(kind(HOMOGENIZATION_undefined_ID)), dimension(:),   allocatable, public, protected :: &
    homogenization_type                                                                             !< type of each homogenization

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

    module subroutine mechanical_partition(subF,ip,el)
      real(pReal), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ip, &                                                                                       !< integration point
        el                                                                                          !< element number
    end subroutine mechanical_partition

    module subroutine thermal_partition(ce)
      integer,     intent(in) :: ce
    end subroutine thermal_partition

    module subroutine damage_partition(ce)
      integer,     intent(in) :: ce
    end subroutine damage_partition

    module subroutine thermal_homogenize(ip,el)
      integer, intent(in) :: ip,el
    end subroutine thermal_homogenize

    module subroutine mechanical_homogenize(dt,ip,el)
     real(pReal), intent(in) :: dt
     integer, intent(in) :: &
       ip, &                                                                                        !< integration point
       el                                                                                           !< element number
    end subroutine mechanical_homogenize

    module subroutine mechanical_results(group_base,h)
      character(len=*), intent(in) :: group_base
      integer, intent(in)          :: h
    end subroutine mechanical_results

    module function mechanical_updateState(subdt,subF,ip,el) result(doneAndHappy)
      real(pReal), intent(in) :: &
        subdt                                                                                           !< current time step
      real(pReal), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ip, &                                                                                           !< integration point
        el                                                                                              !< element number
      logical, dimension(2) :: doneAndHappy
    end function mechanical_updateState


    module function thermal_conduction_getConductivity(ip,el) result(K)

      integer, intent(in) :: &
        ip, &                                                                                           !< integration point number
        el                                                                                              !< element number
      real(pReal), dimension(3,3) :: K

    end function thermal_conduction_getConductivity

    module function thermal_conduction_getSpecificHeat(ce) result(c_P)
      integer, intent(in) :: ce
      real(pReal) :: c_P
    end function thermal_conduction_getSpecificHeat

    module function thermal_conduction_getMassDensity(ce) result(rho)
      integer, intent(in) :: ce
      real(pReal) :: rho
    end function thermal_conduction_getMassDensity

    module subroutine homogenization_thermal_setField(T,dot_T, ce)
      integer, intent(in) :: ce
      real(pReal),   intent(in) :: T, dot_T
    end subroutine homogenization_thermal_setField

    module subroutine thermal_conduction_results(ho,group)
      integer,          intent(in) :: ho
      character(len=*), intent(in) :: group
    end subroutine thermal_conduction_results

    module function homogenization_thermal_T(ce) result(T)
      integer, intent(in) :: ce
      real(pReal) :: T
    end function homogenization_thermal_T

    module subroutine thermal_conduction_getSource(Tdot, ip,el)
      integer, intent(in) :: &
        ip, &                                                                                           !< integration point number
        el                                                                                              !< element number
      real(pReal), intent(out) :: Tdot
    end subroutine thermal_conduction_getSource

   module function damage_nonlocal_getMobility(ip,el) result(M)
    integer, intent(in) :: &
      ip, &                                                                                           !< integration point number
      el                                                                                              !< element number
    real(pReal) :: M
    end function damage_nonlocal_getMobility

    module subroutine damage_nonlocal_getSourceAndItsTangent(phiDot, dPhiDot_dPhi, phi, ip, el)

      integer, intent(in) :: &
        ip, &                                                                                           !< integration point number
        el                                                                                              !< element number
      real(pReal),   intent(in) :: &
        phi
      real(pReal) :: &
        phiDot, dPhiDot_dPhi
    end subroutine damage_nonlocal_getSourceAndItsTangent


    module subroutine damage_nonlocal_putNonLocalDamage(phi,ip,el)

      integer, intent(in) :: &
        ip, &                                                                                           !< integration point number
        el                                                                                              !< element number
      real(pReal),   intent(in) :: &
        phi

    end subroutine damage_nonlocal_putNonLocalDamage

    module subroutine damage_nonlocal_results(homog,group)

      integer,          intent(in) :: homog
      character(len=*), intent(in) :: group

    end subroutine damage_nonlocal_results
  end interface

  public ::  &
    homogenization_init, &
    materialpoint_stressAndItsTangent, &
    thermal_conduction_getSpecificHeat, &
    thermal_conduction_getConductivity, &
    thermal_conduction_getMassDensity, &
    thermal_conduction_getSource, &
    damage_nonlocal_getMobility, &
    damage_nonlocal_getSourceAndItsTangent, &
    damage_nonlocal_putNonLocalDamage, &
    homogenization_thermal_setfield, &
    homogenization_thermal_T, &
    homogenization_forward, &
    homogenization_results, &
    homogenization_restartRead, &
    homogenization_restartWrite, &
    THERMAL_CONDUCTION_ID, &
    DAMAGE_NONLOCAL_ID

  public :: &
    damage_nonlocal_init, &
    damage_nonlocal_getDiffusion

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
  call material_parseHomogenization


  num_homog        => config_numerics%get('homogenization',defaultVal=emptyDict)
  num_homogGeneric => num_homog%get('generic',defaultVal=emptyDict)

  num%nMPstate          = num_homogGeneric%get_asInt  ('nMPstate',     defaultVal=10)
  if (num%nMPstate < 1)                   call IO_error(301,ext_msg='nMPstate')


  call mechanical_init(num_homog)
  call thermal_init()
  call damage_init()

  call damage_nonlocal_init


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
    myNgrains, co, ce, ho, me, ph
  logical :: &
    converged
  logical, dimension(2) :: &
    doneAndHappy

  !$OMP PARALLEL
  !$OMP DO PRIVATE(ce,me,ho,myNgrains,NiterationMPstate,converged,doneAndHappy)
  do el = FEsolving_execElem(1),FEsolving_execElem(2)
    ho = material_homogenizationAt(el)
    myNgrains = homogenization_Nconstituents(ho)
    do ip = FEsolving_execIP(1),FEsolving_execIP(2)
      ce = (el-1)*discretization_nIPs + ip
      me = material_homogenizationMemberAt2(ce)

      call phase_restore(ce,.false.) ! wrong name (is more a forward function)

      if(homogState(ho)%sizeState > 0)  homogState(ho)%state(:,me) = homogState(ho)%state0(:,me)
      if(damageState_h(ho)%sizeState > 0) damageState_h(ho)%state(:,me) = damageState_h(ho)%state0(:,me)
      call damage_partition(ce)

      doneAndHappy = [.false.,.true.]

      NiterationMPstate = 0
      convergenceLooping: do while (.not. (terminallyIll .or. doneAndHappy(1)) &
                                    .and. NiterationMPstate < num%nMPstate)
        NiterationMPstate = NiterationMPstate + 1


        if (.not. doneAndHappy(1)) then
          call mechanical_partition(homogenization_F(1:3,1:3,ce),ip,el)
          converged = .true.
          do co = 1, myNgrains
            converged = converged .and. crystallite_stress(dt,co,ip,el)
          enddo

          if (.not. converged) then
            doneAndHappy = [.true.,.false.]
          else
            doneAndHappy = mechanical_updateState(dt,homogenization_F(1:3,1:3,ce),ip,el)
            converged = all(doneAndHappy)
          endif
        endif

      enddo convergenceLooping
      if (.not. converged) then
        if (.not. terminallyIll) print*, ' Integration point ', ip,' at element ', el, ' terminally ill'
        terminallyIll = .true.
      endif
    enddo
  enddo
  !$OMP END DO

  if (.not. terminallyIll ) then
    !$OMP DO PRIVATE(ho,ph,ce)
    do el = FEsolving_execElem(1),FEsolving_execElem(2)
      if (terminallyIll) continue
      ho = material_homogenizationAt(el)
      do ip = FEsolving_execIP(1),FEsolving_execIP(2)
        ce = (el-1)*discretization_nIPs + ip
        call thermal_partition(ce)
        do co = 1, homogenization_Nconstituents(ho)
          ph = material_phaseAt(co,el)
          if (.not. thermal_stress(dt,ph,material_phaseMemberAt(co,ip,el))) then
            if (.not. terminallyIll) &                                                           ! so first signals terminally ill...
              print*, ' Integration point ', ip,' at element ', el, ' terminally ill'
            terminallyIll = .true.                                                                  ! ...and kills all others
         endif
         call thermal_homogenize(ip,el)
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO PRIVATE(ho)
    elementLooping3: do el = FEsolving_execElem(1),FEsolving_execElem(2)
      ho = material_homogenizationAt(el)
      IpLooping3: do ip = FEsolving_execIP(1),FEsolving_execIP(2)
        do co = 1, homogenization_Nconstituents(ho)
          call crystallite_orientations(co,ip,el)
        enddo
        call mechanical_homogenize(dt,ip,el)
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

    group = trim(group_base)//'/damage'
    call results_closeGroup(results_addGroup(group))
    select case(damage_type(ho))
      case(DAMAGE_NONLOCAL_ID)
        call damage_nonlocal_results(ho,group)
    end select

    group = trim(group_base)//'/thermal'
    call results_closeGroup(results_addGroup(group))
    select case(thermal_type(ho))
      case(THERMAL_CONDUCTION_ID)
        call thermal_conduction_results(ho,group)
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

    call HDF5_write(groupHandle(2),homogState(ho)%state,'omega') ! ToDo: should be done by mech

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

    call HDF5_read(groupHandle(2),homogState(ho)%state0,'omega') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine homogenization_restartRead



!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine damage_nonlocal_init

  integer :: Ninstances,Nmaterialpoints,h
  class(tNode), pointer :: &
    num_generic, &
    material_homogenization

  print'(/,a)', ' <<<+-  damage_nonlocal init  -+>>>'; flush(6)

!------------------------------------------------------------------------------------
! read numerics parameter
  num_generic => config_numerics%get('generic',defaultVal= emptyDict)
  num_damage%charLength = num_generic%get_asFloat('charLength',defaultVal=1.0_pReal)

  Ninstances = count(damage_type == DAMAGE_nonlocal_ID)

  material_homogenization => config_material%get('homogenization')
  do h = 1, material_homogenization%length
    if (damage_type(h) /= DAMAGE_NONLOCAL_ID) cycle

    Nmaterialpoints = count(material_homogenizationAt == h)
    damageState_h(h)%sizeState = 1
    allocate(damageState_h(h)%state0   (1,Nmaterialpoints), source=1.0_pReal)
    allocate(damageState_h(h)%state    (1,Nmaterialpoints), source=1.0_pReal)

  enddo

end subroutine damage_nonlocal_init


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
  do grain = 1, homogenization_Nconstituents(homog)
    damage_nonlocal_getDiffusion = damage_nonlocal_getDiffusion + &
      crystallite_push33ToRef(grain,ip,el,lattice_D(1:3,1:3,material_phaseAt(grain,el)))
  enddo

  damage_nonlocal_getDiffusion = &
    num_damage%charLength**2*damage_nonlocal_getDiffusion/real(homogenization_Nconstituents(homog),pReal)

end function damage_nonlocal_getDiffusion


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


end module homogenization
