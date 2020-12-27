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
  use constitutive
  use discretization
  use thermal_isothermal
  use thermal_conduction
  use damage_none
  use damage_nonlocal
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
    real(pReal) :: &
      subStepMinHomog, &                                                                            !< minimum (relative) size of sub-step allowed during cutback in homogenization
      subStepSizeHomog, &                                                                           !< size of first substep when cutback in homogenization
      stepIncreaseHomog                                                                             !< increase of next substep size when previous substep converged in homogenization
  end type tNumerics

  type(tNumerics) :: num

!--------------------------------------------------------------------------------------------------
  interface

    module subroutine mech_init(num_homog)
      class(tNode), pointer, intent(in) :: &
        num_homog                                                                                   !< pointer to mechanical homogenization numerics data
    end subroutine mech_init

    module subroutine mech_partition(subF,ip,el)
      real(pReal), intent(in), dimension(3,3) :: &
        subF
      integer,     intent(in) :: &
        ip, &                                                                                       !< integration point
        el                                                                                          !< element number
    end subroutine mech_partition

    module subroutine mech_homogenize(ip,el)
     integer, intent(in) :: &
       ip, &                                                                                        !< integration point
       el                                                                                           !< element number
    end subroutine mech_homogenize

    module subroutine mech_results(group_base,h)

      character(len=*), intent(in) :: group_base
      integer, intent(in)          :: h

    end subroutine mech_results

! -------- ToDo ---------------------------------------------------------
    module function mech_RGC_updateState(P,F,F0,avgF,dt,dPdF,ip,el)
      logical, dimension(2) :: mech_RGC_updateState
      real(pReal), dimension(:,:,:),     intent(in)    :: &
        P,&                                                                                         !< partitioned stresses
        F,&                                                                                         !< partitioned deformation gradients
        F0                                                                                          !< partitioned initial deformation gradients
      real(pReal), dimension(:,:,:,:,:), intent(in) :: dPdF                                         !< partitioned stiffnesses
      real(pReal), dimension(3,3),       intent(in) :: avgF                                         !< average F
      real(pReal),                       intent(in) :: dt                                           !< time increment
      integer,                           intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
    end function mech_RGC_updateState

  end interface
! -----------------------------------------------------------------------

  public ::  &
    homogenization_init, &
    materialpoint_stressAndItsTangent, &
    homogenization_forward, &
    homogenization_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init

  class (tNode) , pointer :: &
    num_homog, &
    num_homogGeneric

  print'(/,a)', ' <<<+-  homogenization init  -+>>>'; flush(IO_STDOUT)

  num_homog        => config_numerics%get('homogenization',defaultVal=emptyDict)
  num_homogGeneric => num_homog%get('generic',defaultVal=emptyDict)

  num%nMPstate          = num_homogGeneric%get_asInt  ('nMPstate',     defaultVal=10)
  num%subStepMinHomog   = num_homogGeneric%get_asFloat('subStepMin',   defaultVal=1.0e-3_pReal)
  num%subStepSizeHomog  = num_homogGeneric%get_asFloat('subStepSize',  defaultVal=0.25_pReal)
  num%stepIncreaseHomog = num_homogGeneric%get_asFloat('stepIncrease', defaultVal=1.5_pReal)

  if (num%nMPstate < 1)                   call IO_error(301,ext_msg='nMPstate')
  if (num%subStepMinHomog <= 0.0_pReal)   call IO_error(301,ext_msg='subStepMinHomog')
  if (num%subStepSizeHomog <= 0.0_pReal)  call IO_error(301,ext_msg='subStepSizeHomog')
  if (num%stepIncreaseHomog <= 0.0_pReal) call IO_error(301,ext_msg='stepIncreaseHomog')


  call mech_init(num_homog)

  if (any(thermal_type == THERMAL_isothermal_ID)) call thermal_isothermal_init
  if (any(thermal_type == THERMAL_conduction_ID)) call thermal_conduction_init

  if (any(damage_type == DAMAGE_none_ID))      call damage_none_init
  if (any(damage_type == DAMAGE_nonlocal_ID))  call damage_nonlocal_init


end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(dt,FEsolving_execIP,FEsolving_execElem)

  real(pReal), intent(in) :: dt                                                                     !< time increment
  integer, dimension(2), intent(in) :: FEsolving_execElem, FEsolving_execIP
  integer :: &
    NiterationHomog, &
    NiterationMPstate, &
    ip, &                                                                                            !< integration point number
    el, &                                                                                            !< element number
    myNgrains, co, ce, ho
  real(pReal) :: &
    subFrac, &
    subStep
  logical :: &
    converged
  logical, dimension(2) :: &
    doneAndHappy


!$OMP PARALLEL DO PRIVATE(ce,ho,myNgrains,NiterationMPstate,NiterationHomog,subFrac,converged,subStep,doneAndHappy)
  do el = FEsolving_execElem(1),FEsolving_execElem(2)
    ho = material_homogenizationAt(el)
    myNgrains = homogenization_Nconstituents(ho)
    do ip = FEsolving_execIP(1),FEsolving_execIP(2)

!--------------------------------------------------------------------------------------------------
! initialize restoration points
      call constitutive_initializeRestorationPoints(ip,el)

      subFrac = 0.0_pReal
      converged = .false.                                                                           ! pretend failed step ...
      subStep = 1.0_pReal/num%subStepSizeHomog                                                      ! ... larger then the requested calculation

      if (homogState(ho)%sizeState > 0) &
          homogState(ho)%subState0(:,material_homogenizationMemberAt(ip,el)) = &
          homogState(ho)%State0(   :,material_homogenizationMemberAt(ip,el))

      if (damageState(ho)%sizeState > 0) &
          damageState(ho)%subState0(:,material_homogenizationMemberAt(ip,el)) = &
          damageState(ho)%State0(   :,material_homogenizationMemberAt(ip,el))


      NiterationHomog = 0
      cutBackLooping: do while (.not. terminallyIll .and. subStep  > num%subStepMinHomog)

        if (converged) then
          subFrac = subFrac + subStep
          subStep = min(1.0_pReal-subFrac,num%stepIncreaseHomog*subStep)             ! introduce flexibility for step increase/acceleration

          steppingNeeded: if (subStep > num%subStepMinHomog) then

            ! wind forward grain starting point
            call constitutive_windForward(ip,el)

            if(homogState(ho)%sizeState > 0) &
                homogState(ho)%subState0(:,material_homogenizationMemberAt(ip,el)) = &
                homogState(ho)%State    (:,material_homogenizationMemberAt(ip,el))
            if(damageState(ho)%sizeState > 0) &
                damageState(ho)%subState0(:,material_homogenizationMemberAt(ip,el)) = &
                damageState(ho)%State    (:,material_homogenizationMemberAt(ip,el))

          endif steppingNeeded
        elseif ( (myNgrains == 1 .and. subStep <= 1.0 ) .or. &                                   ! single grain already tried internal subStepping in crystallite
               num%subStepSizeHomog * subStep <=  num%subStepMinHomog ) then                   ! would require too small subStep
                                                                                                    ! cutback makes no sense
          if (.not. terminallyIll) &                                                           ! so first signals terminally ill...
            print*, ' Integration point ', ip,' at element ', el, ' terminally ill'
          terminallyIll = .true.                                                                  ! ...and kills all others
        else                                                                                      ! cutback makes sense
          subStep = num%subStepSizeHomog * subStep                                      ! crystallite had severe trouble, so do a significant cutback

          call crystallite_restore(ip,el,subStep < 1.0_pReal)
          call constitutive_restore(ip,el)

          if(homogState(ho)%sizeState > 0) &
              homogState(ho)%State(    :,material_homogenizationMemberAt(ip,el)) = &
              homogState(ho)%subState0(:,material_homogenizationMemberAt(ip,el))
          if(damageState(ho)%sizeState > 0) &
              damageState(ho)%State(    :,material_homogenizationMemberAt(ip,el)) = &
              damageState(ho)%subState0(:,material_homogenizationMemberAt(ip,el))
        endif

        if (subStep > num%subStepMinHomog) doneAndHappy = [.false.,.true.]

        NiterationMPstate = 0
        convergenceLooping: do while (.not. terminallyIll &
                       .and. .not. doneAndHappy(1) &
                       .and. NiterationMPstate < num%nMPstate)
          NiterationMPstate = NiterationMPstate + 1

!--------------------------------------------------------------------------------------------------
! deformation partitioning

          if (.not. doneAndHappy(1)) then
            ce = (el-1)*discretization_nIPs + ip
            call mech_partition(homogenization_F0(1:3,1:3,ce) &
                                      + (homogenization_F(1:3,1:3,ce)-homogenization_F0(1:3,1:3,ce))&
                                         *(subStep+subFrac), &
                                      ip,el)
            converged = .true.
            do co = 1, myNgrains
              converged = converged .and. crystallite_stress(dt*subStep,co,ip,el)
            enddo

            if (.not. converged) then
              doneAndHappy = [.true.,.false.]
            else
              ce = (el-1)*discretization_nIPs + ip
              doneAndHappy = updateState(dt*subStep, &
                                                  homogenization_F0(1:3,1:3,ce) &
                                                  + (homogenization_F(1:3,1:3,ce)-homogenization_F0(1:3,1:3,ce)) &
                                                     *(subStep+subFrac), &
                                                 ip,el)
              converged = all(doneAndHappy)
            endif
          endif

        enddo convergenceLooping
        NiterationHomog = NiterationHomog + 1

      enddo cutBackLooping
    enddo
  enddo
  !$OMP END PARALLEL DO

  if (.not. terminallyIll ) then
    !$OMP PARALLEL DO PRIVATE(ho,myNgrains)
    elementLooping3: do el = FEsolving_execElem(1),FEsolving_execElem(2)
      ho = material_homogenizationAt(el)
      myNgrains = homogenization_Nconstituents(ho)
      IpLooping3: do ip = FEsolving_execIP(1),FEsolving_execIP(2)
        do co = 1, myNgrains
          call crystallite_orientations(co,ip,el)
        enddo
        call mech_homogenize(ip,el)
      enddo IpLooping3
    enddo elementLooping3
    !$OMP END PARALLEL DO
  else
    print'(/,a,/)', ' << HOMOG >> Material Point terminally ill'
  endif

end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and
!> "happy" with result
!--------------------------------------------------------------------------------------------------
function updateState(subdt,subF,ip,el)

  real(pReal), intent(in) :: &
    subdt                                                                                           !< current time step
  real(pReal), intent(in), dimension(3,3) :: &
    subF
  integer,     intent(in) :: &
    ip, &                                                                                           !< integration point
    el                                                                                              !< element number
  integer :: c
  logical, dimension(2) :: updateState
  real(pReal) :: dPdFs(3,3,3,3,homogenization_Nconstituents(material_homogenizationAt(el)))

  updateState = .true.
  chosenHomogenization: select case(homogenization_type(material_homogenizationAt(el)))
    case (HOMOGENIZATION_RGC_ID) chosenHomogenization
      do c=1,homogenization_Nconstituents(material_homogenizationAt(el))
        dPdFs(:,:,:,:,c) = crystallite_stressTangent(c,ip,el)
      enddo
      updateState = &
        updateState .and. &
          mech_RGC_updateState(crystallite_P(1:3,1:3,1:homogenization_Nconstituents(material_homogenizationAt(el)),ip,el), &
                        crystallite_partitionedF(1:3,1:3,1:homogenization_Nconstituents(material_homogenizationAt(el)),ip,el), &
                        crystallite_partitionedF0(1:3,1:3,1:homogenization_Nconstituents(material_homogenizationAt(el)),ip,el),&
                               subF,&
                               subdt, &
                               dPdFs, &
                               ip, &
                               el)
  end select chosenHomogenization

end function updateState


!--------------------------------------------------------------------------------------------------
!> @brief writes homogenization results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_results
  use material, only: &
    material_homogenization_type => homogenization_type

  integer :: p
  character(len=:), allocatable :: group_base,group

  call results_closeGroup(results_addGroup('current/homogenization/'))

  do p=1,size(material_name_homogenization)
    group_base = 'current/homogenization/'//trim(material_name_homogenization(p))
    call results_closeGroup(results_addGroup(group_base))

    call mech_results(group_base,p)

    group = trim(group_base)//'/damage'
    call results_closeGroup(results_addGroup(group))
    select case(damage_type(p))
      case(DAMAGE_NONLOCAL_ID)
        call damage_nonlocal_results(p,group)
    end select

    group = trim(group_base)//'/thermal'
    call results_closeGroup(results_addGroup(group))
    select case(thermal_type(p))
      case(THERMAL_CONDUCTION_ID)
        call thermal_conduction_results(p,group)
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
    damageState(ho)%state0 = damageState(ho)%state
  enddo

end subroutine homogenization_forward

end module homogenization
