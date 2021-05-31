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


  character(len=2), allocatable, dimension(:) :: phase_lattice
  real(pReal),      allocatable, dimension(:) :: phase_cOverA
  real(pReal),      allocatable, dimension(:) :: phase_rho

  type(tRotationContainer), dimension(:), allocatable :: &
    phase_orientation0, &
    phase_orientation

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

  logical, dimension(:), allocatable, public :: &                                                   ! ToDo: should be protected (bug in Intel Compiler)
    phase_localPlasticity                                                                           !< flags phases with local constitutive law

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tState),  allocatable, dimension(:), public :: &
    damageState


  interface

! == cleaned:begin =================================================================================
    module subroutine mechanical_init(materials,phases)
      class(tNode), pointer :: materials,phases
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

    module subroutine mechanical_forward()
    end subroutine mechanical_forward

    module subroutine damage_forward()
    end subroutine damage_forward

    module subroutine thermal_forward()
    end subroutine thermal_forward


    module subroutine mechanical_restore(ce,includeL)
      integer, intent(in) :: ce
      logical, intent(in) :: includeL
    end subroutine mechanical_restore


    module function phase_mechanical_dPdF(dt,co,ce) result(dPdF)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: &
        co, &                                                                                       !< counter in constituent loop
        ce
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


    module function mechanical_S(ph,en) result(S)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: S
    end function mechanical_S

    module function mechanical_L_p(ph,en) result(L_p)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: L_p
    end function mechanical_L_p

    module function mechanical_F_e(ph,en) result(F_e)
      integer, intent(in) :: ph,en
      real(pReal), dimension(3,3) :: F_e
    end function mechanical_F_e


    module function phase_F(co,ce) result(F)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: F
    end function phase_F

    module function phase_P(co,ce) result(P)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: P
    end function phase_P

    module function thermal_T(ph,en) result(T)
      integer, intent(in) :: ph,en
      real(pReal) :: T
    end function thermal_T

    module function thermal_dot_T(ph,en) result(dot_T)
      integer, intent(in) :: ph,en
      real(pReal) :: dot_T
    end function thermal_dot_T

    module function damage_phi(ph,me) result(phi)
      integer, intent(in) :: ph,me
      real(pReal) :: phi
    end function damage_phi


    module subroutine phase_set_F(F,co,ce)
      real(pReal), dimension(3,3), intent(in) :: F
      integer, intent(in) :: co, ce
    end subroutine phase_set_F

    module subroutine phase_thermal_setField(T,dot_T, co,ce)
      real(pReal), intent(in) :: T, dot_T
      integer, intent(in) :: co, ce
    end subroutine phase_thermal_setField

    module subroutine phase_set_phi(phi,co,ce)
      real(pReal), intent(in) :: phi
      integer, intent(in) :: co, ce
    end subroutine phase_set_phi


    module function phase_mu_phi(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pReal) :: mu
    end function phase_mu_phi

    module function phase_K_phi(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: K
    end function phase_K_phi


    module function phase_mu_T(co,ce) result(mu)
      integer, intent(in) :: co, ce
      real(pReal) :: mu
    end function phase_mu_T

    module function phase_K_T(co,ce) result(K)
      integer, intent(in) :: co, ce
      real(pReal), dimension(3,3) :: K
    end function phase_K_T

! == cleaned:end ===================================================================================

    module function thermal_stress(Delta_t,ph,en) result(converged_)

      real(pReal), intent(in) :: Delta_t
      integer, intent(in) :: ph, en
      logical :: converged_

    end function thermal_stress

    module function integrateDamageState(dt,co,ce) result(broken)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: &
        ce, &
        co
      logical :: broken
    end function integrateDamageState

    module function crystallite_stress(dt,co,ip,el) result(converged_)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: co, ip, el
      logical :: converged_
    end function crystallite_stress

    !ToDo: Try to merge the all stiffness functions
    module function phase_homogenizedC(ph,en) result(C)
      integer, intent(in) :: ph, en
      real(pReal), dimension(6,6) :: C
    end function phase_homogenizedC
    module function phase_damage_C(C_homogenized,ph,en) result(C)
      real(pReal), dimension(3,3,3,3), intent(in)  :: C_homogenized
      integer,                         intent(in)  :: ph,en
      real(pReal), dimension(3,3,3,3) :: C
    end function phase_damage_C

    module function phase_f_phi(phi,co,ce) result(f)
      integer, intent(in) :: ce,co
      real(pReal), intent(in) :: &
        phi                                                                                         !< damage parameter
      real(pReal) :: &
        f
    end function phase_f_phi

    module function phase_f_T(ph,en) result(f)
      integer, intent(in) :: ph, en
      real(pReal) :: f
    end function phase_f_T

    module subroutine plastic_nonlocal_updateCompatibility(orientation,ph,i,e)
      integer, intent(in) :: &
        ph, &
        i, &
        e
        type(tRotationContainer), dimension(:), intent(in) :: orientation
    end subroutine plastic_nonlocal_updateCompatibility

    module subroutine plastic_dependentState(co,ip,el)
      integer, intent(in) :: &
        co, &                                                                                       !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine plastic_dependentState


    module subroutine damage_anisobrittle_LiAndItsTangent(Ld, dLd_dTstar, S, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine damage_anisobrittle_LiAndItsTangent

    module subroutine damage_isoductile_LiAndItsTangent(Ld, dLd_dTstar, S, ph,me)
      integer, intent(in) :: ph, me
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine damage_isoductile_LiAndItsTangent

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
    phase_f_phi, &
    phase_f_T, &
    phase_K_phi, &
    phase_K_T, &
    phase_mu_phi, &
    phase_mu_T, &
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
    phase_set_phi, &
    phase_P, &
    phase_set_F, &
    phase_F

contains

!--------------------------------------------------------------------------------------------------
!> @brief Initialize constitutive models for individual physics
!--------------------------------------------------------------------------------------------------
subroutine phase_init

  integer :: &
    ph, ce, co, ma
  class (tNode), pointer :: &
    debug_constitutive, &
    materials, &
    phases, &
    phase


  print'(/,a)', ' <<<+-  phase init  -+>>>'; flush(IO_STDOUT)

  debug_constitutive => config_debug%get('phase', defaultVal=emptyList)
  debugConstitutive%basic     = debug_constitutive%contains('basic')
  debugConstitutive%extensive = debug_constitutive%contains('extensive')
  debugConstitutive%selective = debug_constitutive%contains('selective')
  debugConstitutive%element   = config_debug%get_asInt('element',         defaultVal = 1)
  debugConstitutive%ip        = config_debug%get_asInt('integrationpoint',defaultVal = 1)
  debugConstitutive%grain     = config_debug%get_asInt('constituent',     defaultVal = 1)


  materials => config_material%get('material')
  phases    => config_material%get('phase')

  allocate(phase_lattice(phases%length))
  allocate(phase_cOverA(phases%length),source=-1.0_pReal)
  allocate(phase_rho(phases%length))
  allocate(phase_orientation0(phases%length))

  do ph = 1,phases%length
    phase => phases%get(ph)
    phase_lattice(ph) = phase%get_asString('lattice')
    if (all(phase_lattice(ph) /= ['cF','cI','hP','tI'])) &
      call IO_error(130,ext_msg='phase_init: '//phase%get_asString('lattice'))
    if (any(phase_lattice(ph) == ['hP','tI'])) &
      phase_cOverA(ph) = phase%get_asFloat('c/a')
    phase_rho(ph) = phase%get_asFloat('rho',defaultVal=0.0_pReal)
    allocate(phase_orientation0(ph)%data(count(material_phaseID==ph)))
  enddo

  do ce = 1, size(material_phaseID,2)
    ma = discretization_materialAt((ce-1)/discretization_nIPs+1)
    do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
      ph = material_phaseID(co,ce)
      phase_orientation0(ph)%data(material_phaseEntry(co,ce)) = material_orientation0(ma)%data(co)
    enddo
  enddo

  allocate(phase_orientation(phases%length))
  do ph = 1,phases%length
    phase_orientation(ph)%data = phase_orientation0(ph)%data
  enddo

  call mechanical_init(materials,phases)
  call damage_init
  call thermal_init(phases)

end subroutine phase_init


!--------------------------------------------------------------------------------------------------
!> @brief Allocate the components of the state structure for a given phase
!--------------------------------------------------------------------------------------------------
subroutine phase_allocateState(state, &
                               NEntries,sizeState,sizeDotState,sizeDeltaState)

  class(tState), intent(out) :: &
    state
  integer, intent(in) :: &
    NEntries, &
    sizeState, &
    sizeDotState, &
    sizeDeltaState


  state%sizeState        = sizeState
  state%sizeDotState     = sizeDotState
  state%sizeDeltaState   = sizeDeltaState
  state%offsetDeltaState = sizeState-sizeDeltaState                                                 ! deltaState occupies latter part of state by definition

  allocate(state%atol             (sizeState),          source=0.0_pReal)
  allocate(state%state0           (sizeState,NEntries), source=0.0_pReal)
  allocate(state%state            (sizeState,NEntries), source=0.0_pReal)

  allocate(state%dotState      (sizeDotState,NEntries), source=0.0_pReal)

  allocate(state%deltaState  (sizeDeltaState,NEntries), source=0.0_pReal)


end subroutine phase_allocateState


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine phase_restore(ce,includeL)

  logical, intent(in) :: includeL
  integer, intent(in) :: ce

  integer :: &
    co


  do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
    if (damageState(material_phaseID(co,ce))%sizeState > 0) &
    damageState(material_phaseID(co,ce))%state( :,material_phaseEntry(co,ce)) = &
      damageState(material_phaseID(co,ce))%state0(:,material_phaseEntry(co,ce))
  enddo

  call mechanical_restore(ce,includeL)

end subroutine phase_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
!--------------------------------------------------------------------------------------------------
subroutine phase_forward()

  call mechanical_forward()
  call damage_forward()
  call thermal_forward()

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
    ce, &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el, &                                                                                           !< counter in element loop
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax                                                                                            !< maximum number of elements

  class(tNode), pointer :: &
    num_crystallite, &
    phases


  print'(/,a)', ' <<<+-  crystallite init  -+>>>'

  cMax = homogenization_maxNconstituents
  iMax = discretization_nIPs
  eMax = discretization_Nelems

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
    if (damageState(ph)%sizeState > 0) allocate(damageState(ph)%subState0,source=damageState(ph)%state0)  ! ToDo: hack
  enddo

  print'(a42,1x,i10)', '    # of elements:                       ', eMax
  print'(a42,1x,i10)', '    # of integration points/element:     ', iMax
  print'(a42,1x,i10)', 'max # of constituents/integration point: ', cMax
  flush(IO_STDOUT)


  !$OMP PARALLEL DO PRIVATE(ce)
  do el = 1, eMax
    do ip = 1, iMax
      ce = (el-1)*discretization_nIPs + ip
      do co = 1,homogenization_Nconstituents(material_homogenizationID(ce))
        call crystallite_orientations(co,ip,el)
        call plastic_dependentState(co,ip,el)                                                       ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations(co,ip,el)

  integer, intent(in) :: &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el                                                                                              !< counter in element loop

  integer :: ph, en


  ph = material_phaseID(co,(el-1)*discretization_nIPs + ip)
  en = material_phaseEntry(co,(el-1)*discretization_nIPs + ip)

  call phase_orientation(ph)%data(en)%fromMatrix(transpose(math_rotationalPart(mechanical_F_e(ph,en))))

  if (plasticState(material_phaseAt(1,el))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(phase_orientation,material_phaseAt(1,el),ip,el)


end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(co,ce, tensor33)

  real(pReal), dimension(3,3), intent(in) :: tensor33
  integer, intent(in):: &
    co, &
    ce
  real(pReal), dimension(3,3) :: crystallite_push33ToRef

  real(pReal), dimension(3,3)             :: T
  integer :: ph, en

  ph = material_phaseID(co,ce)
  en = material_phaseEntry(co,ce)
  T = matmul(phase_orientation0(ph)%data(en)%asMatrix(),transpose(math_inv33(phase_F(co,ce))))      ! ToDo: initial orientation correct?

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
!> @brief Write restart data to file.
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
!> @brief Read restart data from file.
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
