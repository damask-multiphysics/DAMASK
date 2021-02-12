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
  use damage_nonlocal
  use HDF5_utilities
  use results

  implicit none
  private

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
    homogenization_restartWrite

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init()

  class (tNode) , pointer :: &
    num_homog, &
    num_homogGeneric

  print'(/,a)', ' <<<+-  homogenization init  -+>>>'; flush(IO_STDOUT)

  num_homog        => config_numerics%get('homogenization',defaultVal=emptyDict)
  num_homogGeneric => num_homog%get('generic',defaultVal=emptyDict)

  num%nMPstate          = num_homogGeneric%get_asInt  ('nMPstate',     defaultVal=10)
  if (num%nMPstate < 1)                   call IO_error(301,ext_msg='nMPstate')


  call mechanical_init(num_homog)
  call thermal_init()
  call damage_init()

  if (any(damage_type == DAMAGE_nonlocal_ID))  call damage_nonlocal_init


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

      if(homogState(ho)%sizeState > 0)  homogState(ho)%State(:,me) = homogState(ho)%State0(:,me)
      if(damageState_h(ho)%sizeState > 0) damageState_h(ho)%State(:,me) = damageState_h(ho)%State0(:,me)

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

    call HDF5_read(groupHandle(2),homogState(ho)%state,'omega') ! ToDo: should be done by mech

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

    call HDF5_write(groupHandle(2),homogState(ho)%state,'omega') ! ToDo: should be done by mech

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine homogenization_restartRead


end module homogenization
