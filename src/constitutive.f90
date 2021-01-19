!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, damage & thermal internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive
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
  use results

  implicit none
  private

  enum, bind(c); enumerator :: &
    KINEMATICS_UNDEFINED_ID ,&
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID
  end enum

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


  integer(kind(KINEMATICS_UNDEFINED_ID)),     dimension(:,:), allocatable :: &
    phase_kinematics                                                                                !< active kinematic mechanisms of each phase

  integer, dimension(:), allocatable, public :: &                                                   !< ToDo: should be protected (bug in Intel compiler)
    thermal_Nsources, &
    phase_Nsources, &                                                                               !< number of source mechanisms active in each phase
    phase_Nkinematics, &                                                                            !< number of kinematic mechanisms active in each phase
    phase_NstiffnessDegradations, &                                                                 !< number of stiffness degradation mechanisms active in each phase
    phase_plasticityInstance, &                                                                     !< instance of particular plasticity of each phase
    phase_elasticityInstance                                                                        !< instance of particular elasticity of each phase

  logical, dimension(:), allocatable, public :: &                                                   ! ToDo: should be protected (bug in Intel Compiler)
    phase_localPlasticity                                                                           !< flags phases with local constitutive law

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tSourceState),  allocatable, dimension(:), public :: &
    damageState, thermalState


  integer, public, protected :: &
    constitutive_plasticity_maxSizeDotState, &
    constitutive_source_maxSizeDotState

  interface

! == cleaned:begin =================================================================================
    module subroutine mech_init(phases)
      class(tNode), pointer :: phases
    end subroutine mech_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init(phases)
      class(tNode), pointer :: phases
    end subroutine thermal_init


    module subroutine mech_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine mech_results

    module subroutine damage_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine damage_results


    module subroutine mech_initializeRestorationPoints(ph,me)
      integer, intent(in) :: ph, me
    end subroutine mech_initializeRestorationPoints

    module subroutine constitutive_thermal_initializeRestorationPoints(ph,me)
      integer, intent(in) :: ph, me
    end subroutine constitutive_thermal_initializeRestorationPoints


    module subroutine mech_windForward(ph,me)
      integer, intent(in) :: ph, me
    end subroutine mech_windForward


    module subroutine mech_forward()
    end subroutine mech_forward

    module subroutine thermal_forward()
    end subroutine thermal_forward


    module subroutine mech_restore(ip,el,includeL)
      integer, intent(in) :: ip, el
      logical, intent(in) :: includeL
    end subroutine mech_restore


    module function constitutive_mech_dPdF(dt,co,ip,el) result(dPdF)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: &
        co, &                                                                                       !< counter in constituent loop
        ip, &                                                                                       !< counter in integration point loop
        el                                                                                          !< counter in element loop
      real(pReal), dimension(3,3,3,3) :: dPdF
    end function constitutive_mech_dPdF

    module subroutine mech_restartWrite(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mech_restartWrite

    module subroutine mech_restartRead(groupHandle,ph)
      integer(HID_T), intent(in) :: groupHandle
      integer, intent(in) :: ph
    end subroutine mech_restartRead


    module function mech_S(ph,me) result(S)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: S
    end function mech_S

    module function mech_L_p(ph,me) result(L_p)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: L_p
    end function mech_L_p

    module function constitutive_mech_getF(co,ip,el) result(F)
      integer, intent(in) :: co, ip, el
      real(pReal), dimension(3,3) :: F
    end function constitutive_mech_getF

    module function mech_F_e(ph,me) result(F_e)
      integer, intent(in) :: ph,me
      real(pReal), dimension(3,3) :: F_e
    end function mech_F_e

    module function constitutive_mech_getP(co,ip,el) result(P)
      integer, intent(in) :: co, ip, el
      real(pReal), dimension(3,3) :: P
    end function constitutive_mech_getP

    module function thermal_T(ph,me) result(T)
      integer, intent(in) :: ph,me
      real(pReal) :: T
    end function thermal_T


    module subroutine constitutive_mech_setF(F,co,ip,el)
      real(pReal), dimension(3,3), intent(in) :: F
      integer, intent(in) :: co, ip, el
    end subroutine constitutive_mech_setF

    module subroutine constitutive_thermal_setT(T,co,ce)
      real(pReal), intent(in) :: T
      integer, intent(in) :: co, ce
    end subroutine constitutive_thermal_setT

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

    module function constitutive_homogenizedC(ph,me) result(C)
      integer, intent(in) :: ph, me
      real(pReal), dimension(6,6) :: C
    end function constitutive_homogenizedC

    module subroutine source_damage_anisoBrittle_dotState(S, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),  intent(in), dimension(3,3) :: &
        S
    end subroutine source_damage_anisoBrittle_dotState

    module subroutine source_damage_anisoDuctile_dotState(co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine source_damage_anisoDuctile_dotState

    module subroutine source_damage_isoDuctile_dotState(co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine source_damage_isoDuctile_dotState

    module subroutine source_thermal_externalheat_dotState(phase, of)
      integer, intent(in) :: &
        phase, &
        of
    end subroutine source_thermal_externalheat_dotState

    module subroutine constitutive_damage_getRateAndItsTangents(phiDot, dPhiDot_dPhi, phi, ip, el)
      integer, intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal), intent(in) :: &
        phi                                                                                         !< damage parameter
      real(pReal), intent(inout) :: &
        phiDot, &
        dPhiDot_dPhi
    end subroutine constitutive_damage_getRateAndItsTangents

    module subroutine constitutive_thermal_getRate(TDot, ip,el)
      integer, intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal), intent(out) :: &
        TDot
    end subroutine constitutive_thermal_getRate



    module subroutine plastic_nonlocal_updateCompatibility(orientation,instance,i,e)
      integer, intent(in) :: &
        instance, &
        i, &
        e
      type(rotation), dimension(1,discretization_nIPs,discretization_Nelems), intent(in) :: &
        orientation                                                                                 !< crystal orientation
    end subroutine plastic_nonlocal_updateCompatibility


    module subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dMi,Mi,instance,of)
      real(pReal), dimension(3,3),     intent(out) :: &
        Li                                                                                          !< inleastic velocity gradient
      real(pReal), dimension(3,3,3,3), intent(out)  :: &
        dLi_dMi                                                                                     !< derivative of Li with respect to Mandel stress
      real(pReal), dimension(3,3),     intent(in) :: &
        Mi                                                                                          !< Mandel stress
      integer,                         intent(in) :: &
        instance, &
        of
    end subroutine plastic_isotropic_LiAndItsTangent

    module subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar, S, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_cleavage_opening_LiAndItsTangent

    module subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_slipplane_opening_LiAndItsTangent

    module subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(out), dimension(3,3) :: &
        Li                                                                                          !< thermal velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLi_dTstar                                                                                  !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
    end subroutine kinematics_thermal_expansion_LiAndItsTangent

    module subroutine constitutive_plastic_dependentState(co,ip,el)
      integer, intent(in) :: &
        co, &                                                                                       !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine constitutive_plastic_dependentState

  end interface



  type(tDebugOptions) :: debugConstitutive

  public :: &
    constitutive_init, &
    constitutive_homogenizedC, &
    constitutive_damage_getRateAndItsTangents, &
    constitutive_thermal_getRate, &
    constitutive_results, &
    constitutive_allocateState, &
    constitutive_forward, &
    constitutive_restore, &
    plastic_nonlocal_updateCompatibility, &
    kinematics_active, &
    converged, &
    crystallite_init, &
    crystallite_stress, &
    thermal_stress, &
    constitutive_mech_dPdF, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    constitutive_restartWrite, &
    constitutive_restartRead, &
    integrateDamageState, &
    constitutive_thermal_setT, &
    constitutive_mech_getP, &
    constitutive_mech_setF, &
    constitutive_mech_getF, &
    constitutive_initializeRestorationPoints, &
    constitutive_thermal_initializeRestorationPoints, &
    constitutive_windForward, &
    KINEMATICS_UNDEFINED_ID ,&
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialze constitutive models for individual physics
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init

  integer :: &
    ph, &                                                                                            !< counter in phase loop
    so                                                                                               !< counter in source loop
  class (tNode), pointer :: &
    debug_constitutive, &
    phases


  print'(/,a)', ' <<<+-  constitutive init  -+>>>'; flush(IO_STDOUT)

  debug_constitutive => config_debug%get('constitutive', defaultVal=emptyList)
  debugConstitutive%basic      =  debug_constitutive%contains('basic')
  debugConstitutive%extensive  =  debug_constitutive%contains('extensive')
  debugConstitutive%selective  =  debug_constitutive%contains('selective')
  debugConstitutive%element    =  config_debug%get_asInt('element',defaultVal = 1)
  debugConstitutive%ip         =  config_debug%get_asInt('integrationpoint',defaultVal = 1)
  debugConstitutive%grain      =  config_debug%get_asInt('grain',defaultVal = 1)


  phases => config_material%get('phase')

  call mech_init(phases)
  call damage_init
  call thermal_init(phases)


  constitutive_source_maxSizeDotState = 0
  PhaseLoop2:do ph = 1,phases%length
!--------------------------------------------------------------------------------------------------
! partition and initialize state
    plasticState(ph)%partitionedState0 = plasticState(ph)%state0
    plasticState(ph)%state             = plasticState(ph)%partitionedState0
    forall(so = 1:phase_Nsources(ph))
      damageState(ph)%p(so)%partitionedState0 = damageState(ph)%p(so)%state0
      damageState(ph)%p(so)%state             = damageState(ph)%p(so)%partitionedState0
    end forall

    constitutive_source_maxSizeDotState   = max(constitutive_source_maxSizeDotState, &
                                                maxval(damageState(ph)%p%sizeDotState))
  enddo PhaseLoop2
  constitutive_plasticity_maxSizeDotState = maxval(plasticState%sizeDotState)

end subroutine constitutive_init



!--------------------------------------------------------------------------------------------------
!> @brief checks if a kinematic mechanism is active or not
!--------------------------------------------------------------------------------------------------
function kinematics_active(kinematics_label,kinematics_length)  result(active_kinematics)

  character(len=*), intent(in)         :: kinematics_label                                          !< name of kinematic mechanism
  integer,          intent(in)         :: kinematics_length                                         !< max. number of kinematics in system
  logical, dimension(:,:), allocatable :: active_kinematics

  class(tNode), pointer :: &
    phases, &
    phase, &
    kinematics, &
    kinematics_type
  integer :: p,k

  phases => config_material%get('phase')
  allocate(active_kinematics(kinematics_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    kinematics => phase%get('kinematics',defaultVal=emptyList)
    do k = 1, kinematics%length
      kinematics_type => kinematics%get(k)
      if(kinematics_type%get_asString('type') == kinematics_label) active_kinematics(k,p) = .true.
    enddo
  enddo


end function kinematics_active


!--------------------------------------------------------------------------------------------------
!> @brief Allocate the components of the state structure for a given phase
!--------------------------------------------------------------------------------------------------
subroutine constitutive_allocateState(state, &
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
  allocate(state%partitionedState0(sizeState,Nconstituents), source=0.0_pReal)
  allocate(state%state            (sizeState,Nconstituents), source=0.0_pReal)

  allocate(state%dotState      (sizeDotState,Nconstituents), source=0.0_pReal)

  allocate(state%deltaState  (sizeDeltaState,Nconstituents), source=0.0_pReal)


end subroutine constitutive_allocateState


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_restore(ip,el,includeL)

  logical, intent(in) :: includeL
  integer, intent(in) :: &
    ip, &                                                                                            !< integration point number
    el                                                                                               !< element number

  integer :: &
    co, &                                                                                            !< constituent number
    so


  do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
    do so = 1, phase_Nsources(material_phaseAt(co,el))
      damageState(material_phaseAt(co,el))%p(so)%state(          :,material_phasememberAt(co,ip,el)) = &
      damageState(material_phaseAt(co,el))%p(so)%partitionedState0(:,material_phasememberAt(co,ip,el))
    enddo
  enddo

  call mech_restore(ip,el,includeL)

end subroutine constitutive_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_forward()

  integer :: ph, so


  call mech_forward()
  call thermal_forward()

  do ph = 1, size(damageState)
    do so = 1,phase_Nsources(ph)
      damageState(ph)%p(so)%state0 = damageState(ph)%p(so)%state
  enddo; enddo

end subroutine constitutive_forward


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine constitutive_results()

  integer :: ph
  character(len=:), allocatable :: group


  call results_closeGroup(results_addGroup('/current/phase/'))

  do ph = 1, size(material_name_phase)

    group = '/current/phase/'//trim(material_name_phase(ph))//'/'
    call results_closeGroup(results_addGroup(group))

    call mech_results(group,ph)
    call damage_results(group,ph)

  enddo

end subroutine constitutive_results


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
    phases, &
    phase, &
    mech


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
    do so = 1, thermal_Nsources(ph)
      allocate(thermalState(ph)%p(so)%subState0,source=thermalState(ph)%p(so)%state0)               ! ToDo: hack
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
        call constitutive_plastic_dependentState(co,ip,el)                                          ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief Backup data for homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_initializeRestorationPoints(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                            !< integration point number
    el                                                                                               !< element number
  integer :: &
    co, &                                                                                            !< constituent number
    so,ph, me

  do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
    ph = material_phaseAt(co,el)
    me = material_phaseMemberAt(co,ip,el)

    call mech_initializeRestorationPoints(ph,me)

    do so = 1, size(damageState(ph)%p)
      damageState(ph)%p(so)%partitionedState0(:,me) = damageState(ph)%p(so)%state0(:,me)
    enddo

  enddo

end subroutine constitutive_initializeRestorationPoints


!--------------------------------------------------------------------------------------------------
!> @brief Wind homog inc forward.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_windForward(ip,el)

  integer, intent(in) :: &
    ip, &                                                                                            !< integration point number
    el                                                                                               !< element number

  integer :: &
    co, &                                                                                            !< constituent number
    so, ph, me


  do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
    ph = material_phaseAt(co,el)
    me = material_phaseMemberAt(co,ip,el)

    call mech_windForward(ph,me)

    do so = 1, phase_Nsources(material_phaseAt(co,el))
      damageState(ph)%p(so)%partitionedState0(:,me) = damageState(ph)%p(so)%state(:,me)
    enddo

  enddo

end subroutine constitutive_windForward


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations(co,ip,el)

  integer, intent(in) :: &
    co, &                                                                                            !< counter in integration point component loop
    ip, &                                                                                            !< counter in integration point loop
    el                                                                                               !< counter in element loop


  call crystallite_orientation(co,ip,el)%fromMatrix(transpose(math_rotationalPart(&
    mech_F_e(material_phaseAt(co,el),material_phaseMemberAt(co,ip,el)))))

  if (plasticState(material_phaseAt(1,el))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(crystallite_orientation, &
                                              phase_plasticityInstance(material_phaseAt(1,el)),ip,el)


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


  T = matmul(material_orientation0(co,ip,el)%asMatrix(),transpose(math_inv33(constitutive_mech_getF(co,ip,el)))) ! ToDo: initial orientation correct?

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
subroutine constitutive_restartWrite(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_addGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_addGroup(groupHandle(1),material_name_phase(ph))

    call mech_restartWrite(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine constitutive_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read data for restart
! ToDo: Merge data into one file for MPI
!--------------------------------------------------------------------------------------------------
subroutine constitutive_restartRead(fileHandle)

  integer(HID_T), intent(in) :: fileHandle

  integer(HID_T), dimension(2) :: groupHandle
  integer :: ph


  groupHandle(1) = HDF5_openGroup(fileHandle,'phase')

  do ph = 1, size(material_name_phase)

    groupHandle(2) = HDF5_openGroup(groupHandle(1),material_name_phase(ph))

    call mech_restartRead(groupHandle(2),ph)

    call HDF5_closeGroup(groupHandle(2))

  enddo

  call HDF5_closeGroup(groupHandle(1))

end subroutine constitutive_restartRead


end module constitutive
