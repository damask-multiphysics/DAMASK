!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, damage & thermal internal microstructure state
!--------------------------------------------------------------------------------------------------
module phase
  use prec
  use math
  use rotations
  use IO
  use config
  use material
  use results
  use lattice
  use discretization
  use parallelization
  use HDF5_utilities

  implicit none
  private

  type(rotation),            dimension(:,:,:),        allocatable :: &
    crystallite_orientation                                                                         !< current orientation

  type :: tTensorContainer
    real(pReal), dimension(:,:,:), allocatable :: data
  end type

  type :: tNumerics
    integer :: &
      iJacoLpresiduum, &                                                                            !< frequency of Jacobian update of residuum in Lp
      nState, &                                                                                     !< state loop limit
      nStress                                                                                       !< stress loop limit
    real(pReal) :: &
      subStepMinCryst, &                                                                            !< minimum (relative) size of sub-step allowed during cutback
      subStepSizeCryst, &                                                                           !< size of first substep when cutback
      subStepSizeLp, &                                                                              !< size of first substep when cutback in Lp calculation
      subStepSizeLi, &                                                                              !< size of first substep when cutback in Li calculation
      stepIncreaseCryst, &                                                                          !< increase of next substep size when previous substep converged
      rtol_crystalliteState, &                                                                      !< relative tolerance in state loop
      rtol_crystalliteStress, &                                                                     !< relative tolerance in stress loop
      atol_crystalliteStress                                                                        !< absolute tolerance in stress loop
  end type tNumerics

  type(tNumerics) :: num                                                                            ! numerics parameters. Better name?

  type :: tDebugOptions
    logical :: &
      basic, &
      extensive, &
      selective
    integer :: &
      element, &
      ip, &
      grain
  end type tDebugOptions

  type(tDebugOptions) :: debugCrystallite

  integer, dimension(:), allocatable, public :: &                                                   !< ToDo: should be protected (bug in Intel compiler)
    thermal_Nsources, &
    phase_Nsources, &                                                                               !< number of source mechanisms active in each phase
    phase_Nkinematics, &                                                                            !< number of kinematic mechanisms active in each phase
    phase_NstiffnessDegradations, &                                                                 !< number of stiffness degradation mechanisms active in each phase
    phase_plasticInstance, &                                                                     !< instance of particular plasticity of each phase
    phase_elasticityInstance                                                                        !< instance of particular elasticity of each phase

  logical, dimension(:), allocatable, public :: &                                                   ! ToDo: should be protected (bug in Intel Compiler)
    phase_localPlasticity                                                                           !< flags phases with local constitutive law

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tSourceState),  allocatable, dimension(:), public :: &
    damageState, thermalState


  integer, public, protected :: &
    phase_plasticity_maxSizeDotState, &
    phase_source_maxSizeDotState

  interface

! == cleaned:begin =================================================================================
    module subroutine mechanical_init(phases)
      class(tNode), pointer :: phases
    end subroutine mechanical_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init(phases)
      class(tNode), pointer :: phases
    end subroutine thermal_init


    module subroutine mechanical_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine mechanical_results

    module subroutine damage_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine damage_results

    module subroutine mechanical_windForward(ph,me)
      integer, intent(in) :: ph, me
    end subroutine mechanical_windForward


    module subroutine mechanical_forward()
    end subroutine mechanical_forward

    module subroutine thermal_forward()
    end subroutine thermal_forward


    module subroutine mechanical_restore(ce,includeL)
      integer, intent(in) :: ce
      logical, intent(in) :: includeL
    end subroutine mechanical_restore


    module function phase_mechanical_dPdF(dt,co,ip,el) result(dPdF)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: &
        co, &                                                                                       !< counter in constituent loop
        ip, &                                                                                       !< counter in integration point loop
        el                                                                                          !< counter in element loop
      real(pReal), dimension(3,3,3,3) :: dPdF
    end function phase_mechanical_dPdF

    module subroutine mechanical_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartWrite

    module subroutine mechanical_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mechanical_restartRead


    module function mechanical_S(ph,me) result(S)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: S
    end function mechanical_S

    module function mechanical_L_p(ph,me) result(L_p)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: L_p
    end function mechanical_L_p

    module function phase_mechanical_getF(co,ip,el) result(F)
      integer, intent(in) :: co, ip, el
      real(pReal), dimension(3,3) :: F
    end function phase_mechanical_getF

    module function mechanical_F_e(ph,me) result(F_e)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: F_e
    end function mechanical_F_e

    module function phase_mechanical_getP(co,ip,el) result(P)
      integer, intent(in) :: co, ip, el
      real(pReal), dimension(3,3) :: P
    end function phase_mechanical_getP

    module function phase_damage_get_phi(co,ip,el) result(phi)
      integer, intent(in) :: co, ip, el
      real(pReal) :: phi
    end function phase_damage_get_phi

    module function thermal_T(ph,me) result(T)
      integer, intent(in) :: ph,me
      real(pReal) :: T
    end function thermal_T

    module function thermal_dot_T(ph,me) result(dot_T)
      integer, intent(in) :: ph,me
      real(pReal) :: dot_T
    end function thermal_dot_T


    module subroutine phase_mechanical_setF(F,co,ip,el)
      real(pReal), dimension(3,3), intent(in) :: F
      integer, intent(in) :: co, ip, el
    end subroutine phase_mechanical_setF

    module subroutine phase_thermal_setField(T,dot_T, co,ce)
      real(pReal), intent(in) :: T, dot_T
      integer, intent(in) :: ce, co
    end subroutine phase_thermal_setField

    module subroutine phase_damage_set_phi(phi,co,ce)
      real(pReal), intent(in) :: phi
      integer, intent(in) :: co, ce
    end subroutine phase_damage_set_phi

! == cleaned:end ===================================================================================

    module function thermal_stress(Delta_t,ph,me) result(converged_)

      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: ph, me
      logical :: converged_

    end function thermal_stress

    module function integrateDamageState(dt,co,ip,el) result(broken)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: &
        el, &                                                                                            !< element index in element loop
        ip, &                                                                                            !< integration point index in ip loop
        co                                                                                               !< grain index in grain loop
      logical :: broken
    end function integrateDamageState

    module function crystallite_stress(dt,co,ip,el) result(converged_)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: co, ip, el
      logical :: converged_
    end function crystallite_stress

    module function phase_homogenizedC(ph,me) result(C)
      integer, intent(in) :: ph, me
      real(pReal), dimension(6,6) :: C
    end function phase_homogenizedC


    module subroutine phase_damage_getRateAndItsTangents(phiDot, dPhiDot_dPhi, phi, ip, el)
      integer, intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal), intent(in) :: &
        phi                                                                                         !< damage parameter
      real(pReal), intent(inout) :: &
        phiDot, &
        dPhiDot_dPhi
    end subroutine phase_damage_getRateAndItsTangents

    module subroutine phase_thermal_getRate(TDot, ph,me)
      integer, intent(in) :: ph, me
      real(pReal), intent(out) :: &
        TDot
    end subroutine phase_thermal_getRate

    module subroutine plastic_nonlocal_updateCompatibility(orientation,instance,i,e)
      integer, intent(in) :: &
        instance, &
        i, &
        e
      type(rotation), dimension(1,discretization_nIPs,discretization_Nelems), intent(in) :: &
        orientation                                                                                 !< crystal orientation
    end subroutine plastic_nonlocal_updateCompatibility

    module subroutine plastic_dependentState(co,ip,el)
      integer, intent(in) :: &
        co, &                                                                                       !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine plastic_dependentState

  end interface



  type(tDebugOptions) :: debugConstitutive
#if __INTEL_COMPILER >= 1900
  public :: &
    prec, &
    math, &
    rotations, &
    IO, &
    config, &
    material, &
    results, &
    lattice, &
    discretization, &
    HDF5_utilities
#endif

  public :: &
    phase_init, &
    phase_homogenizedC, &
    phase_damage_getRateAndItsTangents, &
    phase_thermal_getRate, &
    phase_results, &
    phase_allocateState, &
    phase_forward, &
    phase_restore, &
    plastic_nonlocal_updateCompatibility, &
    converged, &
    crystallite_init, &
    crystallite_stress, &
    thermal_stress, &
    phase_mechanical_dPdF, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    phase_restartWrite, &
    phase_restartRead, &
    integrateDamageState, &
    phase_thermal_setField, &
    phase_damage_set_phi, &
    phase_damage_get_phi, &
    phase_mechanical_getP, &
    phase_mechanical_setF, &
    phase_mechanical_getF, &
    phase_windForward

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialze constitutive models for individual physics
!--------------------------------------------------------------------------------------------------
subroutine phase_init

  integer :: &
    ph, &                                                                                            !< counter in phase loop
    so                                                                                               !< counter in source loop
  class (tNode), pointer :: &
    debug_constitutive, &
    phases


  print'(/,a)', ' <<<+-  phase init  -+>>>'; flush(IO_STDOUT)

  debug_constitutive => config_debug%get('constitutive', defaultVal=emptyList)
  debugConstitutive%basic      =  debug_constitutive%contains('basic')
  debugConstitutive%extensive  =  debug_constitutive%contains('extensive')
  debugConstitutive%selective  =  debug_constitutive%contains('selective')
  debugConstitutive%element    =  config_debug%get_asInt('element',defaultVal = 1)
  debugConstitutive%ip         =  config_debug%get_asInt('integrationpoint',defaultVal = 1)
  debugConstitutive%grain      =  config_debug%get_asInt('grain',defaultVal = 1)


  phases => config_material%get('phase')

  call mechanical_init(phases)
  call damage_init
  call thermal_init(phases)


  phase_source_maxSizeDotState = 0
  PhaseLoop2:do ph = 1,phases%length
!--------------------------------------------------------------------------------------------------
! partition and initialize state
    plasticState(ph)%state             = plasticState(ph)%state0
    forall(so = 1:phase_Nsources(ph))
      damageState(ph)%p(so)%state = damageState(ph)%p(so)%state0
    end forall

    phase_source_maxSizeDotState   = max(phase_source_maxSizeDotState, &
                                                maxval(damageState(ph)%p%sizeDotState))
  enddo PhaseLoop2
  phase_plasticity_maxSizeDotState = maxval(plasticState%sizeDotState)

end subroutine phase_init


!--------------------------------------------------------------------------------------------------
!> @brief Allocate the components of the state structure for a given phase
!--------------------------------------------------------------------------------------------------
subroutine phase_allocateState(state, &
                                  Nconstituents,sizeState,sizeDotState,sizeDeltaState)

  class(tState), intent(out) :: &
    state
  integer, intent(in) :: &
    Nconstituents, &
    sizeState, &
    sizeDotState, &
    sizeDeltaState


  state%sizeState        = sizeState
  state%sizeDotState     = sizeDotState
  state%sizeDeltaState   = sizeDeltaState
  state%offsetDeltaState = sizeState-sizeDeltaState                                                 ! deltaState occupies latter part of state by definition

  allocate(state%atol             (sizeState),               source=0.0_pReal)
  allocate(state%state0           (sizeState,Nconstituents), source=0.0_pReal)
  allocate(state%state            (sizeState,Nconstituents), source=0.0_pReal)

  allocate(state%dotState      (sizeDotState,Nconstituents), source=0.0_pReal)

  allocate(state%deltaState  (sizeDeltaState,Nconstituents), source=0.0_pReal)


end subroutine phase_allocateState


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine phase_restore(ce,includeL)

  logical, intent(in) :: includeL
  integer, intent(in) :: ce

  integer :: &
    co, &                                                                                            !< constituent number
    so


  do co = 1,homogenization_Nconstituents(material_homogenizationAt2(ce))
    do so = 1, phase_Nsources(material_phaseAt2(co,ce))
      damageState(material_phaseAt2(co,ce))%p(so)%state( :,material_phasememberAt2(co,ce)) = &
      damageState(material_phaseAt2(co,ce))%p(so)%state0(:,material_phasememberAt2(co,ce))
    enddo
  enddo

  call mechanical_restore(ce,includeL)

end subroutine phase_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine phase_forward()

  integer :: ph, so


  call mechanical_forward()
  call thermal_forward()

  do ph = 1, size(damageState)
    do so = 1,phase_Nsources(ph)
      damageState(ph)%p(so)%state0 = damageState(ph)%p(so)%state
  enddo; enddo

end subroutine phase_forward


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine phase_results()

  integer :: ph
  character(len=:), allocatable :: group


  call results_closeGroup(results_addGroup('/current/phase/'))

  do ph = 1, size(material_name_phase)

    group = '/current/phase/'//trim(material_name_phase(ph))//'/'
    call results_closeGroup(results_addGroup(group))

    call mechanical_results(group,ph)
    call damage_results(group,ph)

  enddo

end subroutine phase_results


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init()

  integer :: &
    ph, &
    me, &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el, &                                                                                           !< counter in element loop
    so, &
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax                                                                                            !< maximum number of elements

  class(tNode), pointer :: &
    num_crystallite, &
    debug_crystallite, &                                                                            ! pointer to debug options for crystallite
    phases


  print'(/,a)', ' <<<+-  crystallite init  -+>>>'

  debug_crystallite => config_debug%get('crystallite', defaultVal=emptyList)
  debugCrystallite%extensive = debug_crystallite%contains('extensive')

  cMax = homogenization_maxNconstituents
  iMax = discretization_nIPs
  eMax = discretization_Nelems

  allocate(crystallite_orientation(cMax,iMax,eMax))

  num_crystallite => config_numerics%get('crystallite',defaultVal=emptyDict)

  num%subStepMinCryst        = num_crystallite%get_asFloat ('subStepMin',       defaultVal=1.0e-3_pReal)
  num%subStepSizeCryst       = num_crystallite%get_asFloat ('subStepSize',      defaultVal=0.25_pReal)
  num%stepIncreaseCryst      = num_crystallite%get_asFloat ('stepIncrease',     defaultVal=1.5_pReal)
  num%subStepSizeLp          = num_crystallite%get_asFloat ('subStepSizeLp',    defaultVal=0.5_pReal)
  num%subStepSizeLi          = num_crystallite%get_asFloat ('subStepSizeLi',    defaultVal=0.5_pReal)
  num%rtol_crystalliteState  = num_crystallite%get_asFloat ('rtol_State',       defaultVal=1.0e-6_pReal)
  num%rtol_crystalliteStress = num_crystallite%get_asFloat ('rtol_Stress',      defaultVal=1.0e-6_pReal)
  num%atol_crystalliteStress = num_crystallite%get_asFloat ('atol_Stress',      defaultVal=1.0e-8_pReal)
  num%iJacoLpresiduum        = num_crystallite%get_asInt   ('iJacoLpresiduum',  defaultVal=1)
  num%nState                 = num_crystallite%get_asInt   ('nState',           defaultVal=20)
  num%nStress                = num_crystallite%get_asInt   ('nStress',          defaultVal=40)

  if(num%subStepMinCryst   <= 0.0_pReal)      call IO_error(301,ext_msg='subStepMinCryst')
  if(num%subStepSizeCryst  <= 0.0_pReal)      call IO_error(301,ext_msg='subStepSizeCryst')
  if(num%stepIncreaseCryst <= 0.0_pReal)      call IO_error(301,ext_msg='stepIncreaseCryst')

  if(num%subStepSizeLp <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLp')
  if(num%subStepSizeLi <= 0.0_pReal)          call IO_error(301,ext_msg='subStepSizeLi')

  if(num%rtol_crystalliteState  <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteState')
  if(num%rtol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='rtol_crystalliteStress')
  if(num%atol_crystalliteStress <= 0.0_pReal) call IO_error(301,ext_msg='atol_crystalliteStress')

  if(num%iJacoLpresiduum < 1)                 call IO_error(301,ext_msg='iJacoLpresiduum')

  if(num%nState < 1)                          call IO_error(301,ext_msg='nState')
  if(num%nStress< 1)                          call IO_error(301,ext_msg='nStress')


  phases => config_material%get('phase')

  do ph = 1, phases%length
    do so = 1, phase_Nsources(ph)
      allocate(damageState(ph)%p(so)%subState0,source=damageState(ph)%p(so)%state0)                 ! ToDo: hack
    enddo
  enddo

  print'(a42,1x,i10)', '    # of elements:                       ', eMax
  print'(a42,1x,i10)', '    # of integration points/element:     ', iMax
  print'(a42,1x,i10)', 'max # of constituents/integration point: ', cMax
  flush(IO_STDOUT)


  !$OMP PARALLEL DO PRIVATE(ph,me)
  do el = 1, size(material_phaseMemberAt,3)
    do ip = 1, size(material_phaseMemberAt,2)
      do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
        ph = material_phaseAt(co,el)
        me = material_phaseMemberAt(co,ip,el)
        call crystallite_orientations(co,ip,el)
        call plastic_dependentState(co,ip,el)                                          ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief Wind homog inc forward.
!--------------------------------------------------------------------------------------------------
subroutine phase_windForward(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                            !< integration point number
    el                                                                                               !< element number

  integer :: &
    co, &                                                                                            !< constituent number
    so, ph, me


  do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
    ph = material_phaseAt(co,el)
    me = material_phaseMemberAt(co,ip,el)

    call mechanical_windForward(ph,me)

    do so = 1, phase_Nsources(material_phaseAt(co,el))
      damageState(ph)%p(so)%state0(:,me) = damageState(ph)%p(so)%state(:,me)
    enddo

  enddo

end subroutine phase_windForward


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations(co,ip,el)

  integer, intent(in) :: &
    co, &                                                                                            !< counter in integration point component loop
    ip, &                                                                                            !< counter in integration point loop
    el                                                                                               !< counter in element loop


  call crystallite_orientation(co,ip,el)%fromMatrix(transpose(math_rotationalPart(&
    mechanical_F_e(material_phaseAt(co,el),material_phaseMemberAt(co,ip,el)))))

  if (plasticState(material_phaseAt(1,el))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(crystallite_orientation, &
                                              phase_plasticInstance(material_phaseAt(1,el)),ip,el)


end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(co,ip,el, tensor33)

  real(pReal), dimension(3,3), intent(in) :: tensor33
  integer, intent(in):: &
    el, &
    ip, &
    co
  real(pReal), dimension(3,3) :: crystallite_push33ToRef

  real(pReal), dimension(3,3)             :: T


  T = matmul(material_orientation0(co,ip,el)%asMatrix(),transpose(math_inv33(phase_mechanical_getF(co,ip,el)))) ! ToDo: initial orientation correct?

  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief determines whether a point is converged
!--------------------------------------------------------------------------------------------------
logical pure function converged(residuum,state,atol)

  real(pReal), intent(in), dimension(:) ::&
    residuum, state, atol
  real(pReal) :: &
    rTol

  rTol = num%rTol_crystalliteState

  converged = all(abs(residuum) <= max(atol, rtol*abs(state)))

end function converged


!--------------------------------------------------------------------------------------------------
!> @brief Write current  restart information (Field and constitutive data) to file.
! ToDo: Merge data into one file for MPI
!--------------------------------------------------------------------------------------------------
subroutine phase_restartWrite(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_addGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_addGroup(groupHandle(1),material_name_phase(ph))

    call mechanical_restartWrite(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine phase_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read data for restart
! ToDo: Merge data into one file for MPI
!--------------------------------------------------------------------------------------------------
subroutine phase_restartRead(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_openGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_openGroup(groupHandle(1),material_name_phase(ph))

    call mechanical_restartRead(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine phase_restartRead


end module phase
