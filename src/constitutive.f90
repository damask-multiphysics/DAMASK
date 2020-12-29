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
  use DAMASK_interface
  use results

  implicit none
  private

  enum, bind(c); enumerator :: &
    PLASTICITY_UNDEFINED_ID, &
    PLASTICITY_NONE_ID, &
    PLASTICITY_ISOTROPIC_ID, &
    PLASTICITY_PHENOPOWERLAW_ID, &
    PLASTICITY_KINEHARDENING_ID, &
    PLASTICITY_DISLOTWIN_ID, &
    PLASTICITY_DISLOTUNGSTEN_ID, &
    PLASTICITY_NONLOCAL_ID, &
    SOURCE_UNDEFINED_ID ,&
    SOURCE_THERMAL_DISSIPATION_ID, &
    SOURCE_THERMAL_EXTERNALHEAT_ID, &
    SOURCE_DAMAGE_ISOBRITTLE_ID, &
    SOURCE_DAMAGE_ISODUCTILE_ID, &
    SOURCE_DAMAGE_ANISOBRITTLE_ID, &
    SOURCE_DAMAGE_ANISODUCTILE_ID, &
    KINEMATICS_UNDEFINED_ID ,&
    KINEMATICS_CLEAVAGE_OPENING_ID, &
    KINEMATICS_SLIPPLANE_OPENING_ID, &
    KINEMATICS_THERMAL_EXPANSION_ID
  end enum

  type(rotation),            dimension(:,:,:),        allocatable :: &
    crystallite_orientation                                                                         !< current orientation
  real(pReal),               dimension(:,:,:,:,:),    allocatable :: &
    crystallite_F0, &                                                                               !< def grad at start of FE inc
    crystallite_S0, &                                                                               !< 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partitionedS0                                                                       !< 2nd Piola-Kirchhoff stress vector at start of homog inc
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_P, &                                                                                !< 1st Piola-Kirchhoff stress per grain
    crystallite_S, &                                                                                !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_partitionedF0, &                                                                    !< def grad at start of homog inc
    crystallite_F                                                                                   !< def grad to be reached at end of homog inc

  type :: tTensorContainer
    real(pReal), dimension(:,:,:), allocatable :: data
  end type

  type(tTensorContainer), dimension(:), allocatable :: &
    constitutive_mech_Fe, &
    constitutive_mech_Fi, &
    constitutive_mech_Fp, &
    constitutive_mech_Li, &
    constitutive_mech_Lp, &
    constitutive_mech_Fi0, &
    constitutive_mech_Fp0, &
    constitutive_mech_Li0, &
    constitutive_mech_Lp0, &
    constitutive_mech_partitionedFi0, &
    constitutive_mech_partitionedFp0, &
    constitutive_mech_partitionedLi0, &
    constitutive_mech_partitionedLp0


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



  integer(kind(PLASTICITY_undefined_ID)), dimension(:),   allocatable, public :: &
    phase_plasticity                                                                                !< plasticity of each phase

  integer(kind(SOURCE_undefined_ID)),     dimension(:,:), allocatable :: &
    phase_source, &                                                                                 !< active sources mechanisms of each phase
    phase_kinematics                                                                                !< active kinematic mechanisms of each phase

  integer, dimension(:), allocatable, public :: &                                                   !< ToDo: should be protected (bug in Intel compiler)
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
    sourceState


  integer, public, protected :: &
    constitutive_plasticity_maxSizeDotState, &
    constitutive_source_maxSizeDotState

  interface

! == cleaned:begin =================================================================================
    module subroutine mech_init
    end subroutine mech_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init
    end subroutine thermal_init


    module subroutine mech_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine mech_results

    module subroutine damage_results(group,ph)
      character(len=*), intent(in) :: group
      integer,          intent(in) :: ph
    end subroutine damage_results


    module subroutine mech_restart_read(fileHandle)
      integer(HID_T), intent(in) :: fileHandle
    end subroutine mech_restart_read

    module subroutine mech_initializeRestorationPoints(ph,me)
      integer, intent(in) :: ph, me
    end subroutine mech_initializeRestorationPoints

    module subroutine constitutive_mech_windForward(ph,me)
      integer, intent(in) :: ph, me
    end subroutine constitutive_mech_windForward

    module subroutine constitutive_mech_forward
    end subroutine constitutive_mech_forward

    module subroutine mech_restore(ip,el,includeL)
      integer, intent(in) :: &
        ip, &
        el
      logical, intent(in) :: &
        includeL
    end subroutine mech_restore

! == cleaned:end ===================================================================================

    module function crystallite_stress(dt,co,ip,el) result(converged_)
      real(pReal), intent(in) :: dt
      integer, intent(in) :: co, ip, el
      logical :: converged_
    end function crystallite_stress

    module function constitutive_homogenizedC(co,ip,el) result(C)
      integer, intent(in) :: co, ip, el
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

    module subroutine constitutive_thermal_getRateAndItsTangents(TDot, dTDot_dT, T, ip, el)
      integer, intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal), intent(in) :: &
        T
      real(pReal), intent(inout) :: &
        TDot, &
        dTDot_dT
    end subroutine constitutive_thermal_getRateAndItsTangents



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


    module subroutine source_damage_isoBrittle_deltaState(C, Fe, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),  intent(in), dimension(3,3) :: &
        Fe
      real(pReal),  intent(in), dimension(6,6) :: &
        C
    end subroutine source_damage_isoBrittle_deltaState


    module subroutine constitutive_plastic_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                         S, Fi, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),   intent(in),  dimension(3,3) :: &
        S, &                                                                                        !< 2nd Piola-Kirchhoff stress
        Fi                                                                                          !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        Lp                                                                                          !< plastic velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLp_dS, &
        dLp_dFi                                                                                     !< derivative of Lp with respect to Fi
    end subroutine constitutive_plastic_LpAndItsTangents


    module subroutine constitutive_plastic_dependentState(F, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),   intent(in), dimension(3,3) :: &
        F                                                                                           !< elastic deformation gradient
    end subroutine constitutive_plastic_dependentState



    module subroutine constitutive_hooke_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, co, ip, el)
      integer, intent(in) :: &
        co, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),   intent(in),  dimension(3,3) :: &
        Fe, &                                                                                       !< elastic deformation gradient
        Fi                                                                                          !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        S                                                                                           !< 2nd Piola-Kirchhoff stress tensor
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dS_dFe, &                                                                                   !< derivative of 2nd P-K stress with respect to elastic deformation gradient
        dS_dFi                                                                                      !< derivative of 2nd P-K stress with respect to intermediate deformation gradient
    end subroutine constitutive_hooke_SandItsTangents

  end interface



  type(tDebugOptions) :: debugConstitutive

  public :: &
    constitutive_init, &
    constitutive_homogenizedC, &
    constitutive_LiAndItsTangents, &
    constitutive_damage_getRateAndItsTangents, &
    constitutive_thermal_getRateAndItsTangents, &
    constitutive_results, &
    constitutive_allocateState, &
    constitutive_forward, &
    constitutive_restore, &
    plastic_nonlocal_updateCompatibility, &
    source_active, &
    kinematics_active, &
    converged, &
    crystallite_init, &
    crystallite_stress, &
    crystallite_stressTangent, &
    crystallite_orientations, &
    crystallite_push33ToRef, &
    crystallite_restartWrite, &
    integrateSourceState, &
    constitutive_mech_getLp, &
    constitutive_mech_getS, &
    crystallite_restartRead, &
    constitutive_initializeRestorationPoints, &
    constitutive_windForward, &
    PLASTICITY_UNDEFINED_ID, &
    PLASTICITY_NONE_ID, &
    PLASTICITY_ISOTROPIC_ID, &
    PLASTICITY_PHENOPOWERLAW_ID, &
    PLASTICITY_KINEHARDENING_ID, &
    PLASTICITY_DISLOTWIN_ID, &
    PLASTICITY_DISLOTUNGSTEN_ID, &
    PLASTICITY_NONLOCAL_ID, &
    SOURCE_UNDEFINED_ID ,&
    SOURCE_THERMAL_DISSIPATION_ID, &
    SOURCE_THERMAL_EXTERNALHEAT_ID, &
    SOURCE_DAMAGE_ISOBRITTLE_ID, &
    SOURCE_DAMAGE_ISODUCTILE_ID, &
    SOURCE_DAMAGE_ANISOBRITTLE_ID, &
    SOURCE_DAMAGE_ANISODUCTILE_ID, &
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


  debug_constitutive => config_debug%get('constitutive', defaultVal=emptyList)
  debugConstitutive%basic      =  debug_constitutive%contains('basic')
  debugConstitutive%extensive  =  debug_constitutive%contains('extensive')
  debugConstitutive%selective  =  debug_constitutive%contains('selective')
  debugConstitutive%element    =  config_debug%get_asInt('element',defaultVal = 1)
  debugConstitutive%ip         =  config_debug%get_asInt('integrationpoint',defaultVal = 1)
  debugConstitutive%grain      =  config_debug%get_asInt('grain',defaultVal = 1)

!--------------------------------------------------------------------------------------------------
! initialize constitutive laws
  print'(/,a)', ' <<<+-  constitutive init  -+>>>'; flush(IO_STDOUT)
  call mech_init
  call damage_init
  call thermal_init


  phases => config_material%get('phase')
  constitutive_source_maxSizeDotState = 0
  PhaseLoop2:do ph = 1,phases%length
!--------------------------------------------------------------------------------------------------
! partition and initialize state
    plasticState(ph)%partitionedState0 = plasticState(ph)%state0
    plasticState(ph)%state             = plasticState(ph)%partitionedState0
    forall(so = 1:phase_Nsources(ph))
      sourceState(ph)%p(so)%partitionedState0 = sourceState(ph)%p(so)%state0
      sourceState(ph)%p(so)%state             = sourceState(ph)%p(so)%partitionedState0
    end forall

    constitutive_source_maxSizeDotState   = max(constitutive_source_maxSizeDotState, &
                                                maxval(sourceState(ph)%p%sizeDotState))
  enddo PhaseLoop2
  constitutive_plasticity_maxSizeDotState = maxval(plasticState%sizeDotState)

end subroutine constitutive_init


!--------------------------------------------------------------------------------------------------
!> @brief checks if a source mechanism is active or not
!--------------------------------------------------------------------------------------------------
function source_active(source_label,src_length)  result(active_source)

  character(len=*), intent(in)         :: source_label                                              !< name of source mechanism
  integer,          intent(in)         :: src_length                                                !< max. number of sources in system
  logical, dimension(:,:), allocatable :: active_source

  class(tNode), pointer :: &
    phases, &
    phase, &
    sources, &
    src
  integer :: p,s

  phases => config_material%get('phase')
  allocate(active_source(src_length,phases%length), source = .false. )
  do p = 1, phases%length
    phase => phases%get(p)
    sources => phase%get('source',defaultVal=emptyList)
    do s = 1, sources%length
      src => sources%get(s)
      if(src%get_asString('type') == source_label) active_source(s,p) = .true.
    enddo
  enddo


end function source_active


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
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                         S, Fi, co, ip, el)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal),   intent(in),  dimension(3,3) :: &
    S                                                                                               !< 2nd Piola-Kirchhoff stress
  real(pReal),   intent(in),  dimension(3,3) :: &
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),   intent(out), dimension(3,3) :: &
    Li                                                                                              !< intermediate velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLi_dS, &                                                                                       !< derivative of Li with respect to S
    dLi_dFi

  real(pReal), dimension(3,3) :: &
    my_Li, &                                                                                        !< intermediate velocity gradient
    FiInv, &
    temp_33
  real(pReal), dimension(3,3,3,3) :: &
    my_dLi_dS
  real(pReal) :: &
    detFi
  integer :: &
    k, i, j, &
    instance, of

  Li = 0.0_pReal
  dLi_dS  = 0.0_pReal
  dLi_dFi = 0.0_pReal

  plasticityType: select case (phase_plasticity(material_phaseAt(co,el)))
    case (PLASTICITY_isotropic_ID) plasticityType
      of = material_phasememberAt(co,ip,el)
      instance = phase_plasticityInstance(material_phaseAt(co,el))
      call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, S ,instance,of)
    case default plasticityType
      my_Li = 0.0_pReal
      my_dLi_dS = 0.0_pReal
  end select plasticityType

  Li = Li + my_Li
  dLi_dS = dLi_dS + my_dLi_dS

  KinematicsLoop: do k = 1, phase_Nkinematics(material_phaseAt(co,el))
    kinematicsType: select case (phase_kinematics(k,material_phaseAt(co,el)))
      case (KINEMATICS_cleavage_opening_ID) kinematicsType
        call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, co, ip, el)
      case (KINEMATICS_slipplane_opening_ID) kinematicsType
        call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, co, ip, el)
      case (KINEMATICS_thermal_expansion_ID) kinematicsType
        call kinematics_thermal_expansion_LiAndItsTangent(my_Li, my_dLi_dS, co, ip, el)
      case default kinematicsType
        my_Li = 0.0_pReal
        my_dLi_dS = 0.0_pReal
    end select kinematicsType
    Li = Li + my_Li
    dLi_dS = dLi_dS + my_dLi_dS
  enddo KinematicsLoop

  FiInv = math_inv33(Fi)
  detFi = math_det33(Fi)
  Li = matmul(matmul(Fi,Li),FiInv)*detFi                                                            !< push forward to intermediate configuration
  temp_33 = matmul(FiInv,Li)

  do i = 1,3; do j = 1,3
    dLi_dS(1:3,1:3,i,j)  = matmul(matmul(Fi,dLi_dS(1:3,1:3,i,j)),FiInv)*detFi
    dLi_dFi(1:3,1:3,i,j) = dLi_dFi(1:3,1:3,i,j) + Li*FiInv(j,i)
    dLi_dFi(1:3,i,1:3,j) = dLi_dFi(1:3,i,1:3,j) + math_I3*temp_33(j,i) + Li*FiInv(j,i)
  enddo; enddo

end subroutine constitutive_LiAndItsTangents


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_damage_collectDotState(S, co, ip, el,ph,of) result(broken)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    of
  real(pReal),  intent(in), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
  integer :: &
    so                                                                                               !< counter in source loop
  logical :: broken


  broken = .false.

  SourceLoop: do so = 1, phase_Nsources(ph)

    sourceType: select case (phase_source(so,ph))

      case (SOURCE_damage_anisoBrittle_ID) sourceType
        call source_damage_anisoBrittle_dotState(S, co, ip, el) ! correct stress?

      case (SOURCE_damage_isoDuctile_ID) sourceType
        call source_damage_isoDuctile_dotState(co, ip, el)

      case (SOURCE_damage_anisoDuctile_ID) sourceType
        call source_damage_anisoDuctile_dotState(co, ip, el)

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(sourceState(ph)%p(so)%dotState(:,of)))

  enddo SourceLoop

end function constitutive_damage_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_thermal_collectDotState(ph,me) result(broken)

  integer, intent(in) :: ph, me
  logical :: broken

  integer :: i


  broken = .false.

  SourceLoop: do i = 1, phase_Nsources(ph)

    if (phase_source(i,ph) == SOURCE_thermal_externalheat_ID) &
      call source_thermal_externalheat_dotState(ph,me)

    broken = broken .or. any(IEEE_is_NaN(sourceState(ph)%p(i)%dotState(:,me)))

  enddo SourceLoop

end function constitutive_thermal_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
function constitutive_damage_deltaState(Fe, co, ip, el, ph, of) result(broken)

  integer, intent(in) :: &
    co, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    ph, &
    of
  real(pReal),   intent(in), dimension(3,3) :: &
    Fe                                                                                              !< elastic deformation gradient
  integer :: &
    so, &
    myOffset, &
    mySize
  logical :: &
    broken


  broken = .false.

  sourceLoop: do so = 1, phase_Nsources(ph)

     sourceType: select case (phase_source(so,ph))

      case (SOURCE_damage_isoBrittle_ID) sourceType
        call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(co,ip,el), Fe, &
                                                   co, ip, el)
        broken = any(IEEE_is_NaN(sourceState(ph)%p(so)%deltaState(:,of)))
        if(.not. broken) then
          myOffset = sourceState(ph)%p(so)%offsetDeltaState
          mySize   = sourceState(ph)%p(so)%sizeDeltaState
          sourceState(ph)%p(so)%state(myOffset + 1: myOffset + mySize,of) = &
          sourceState(ph)%p(so)%state(myOffset + 1: myOffset + mySize,of) + sourceState(ph)%p(so)%deltaState(1:mySize,of)
        endif

    end select sourceType

  enddo SourceLoop

end function constitutive_damage_deltaState


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
      sourceState(material_phaseAt(co,el))%p(so)%state(          :,material_phasememberAt(co,ip,el)) = &
      sourceState(material_phaseAt(co,el))%p(so)%partitionedState0(:,material_phasememberAt(co,ip,el))
    enddo
  enddo

  call mech_restore(ip,el,includeL)

end subroutine constitutive_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_forward

  integer :: i, j

  crystallite_F0  = crystallite_F
  crystallite_S0  = crystallite_S

  call constitutive_mech_forward()

  do i = 1, size(sourceState)
    do j = 1,phase_Nsources(i)
      sourceState(i)%p(j)%state0 = sourceState(i)%p(j)%state
  enddo; enddo

end subroutine constitutive_forward


!--------------------------------------------------------------------------------------------------
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine constitutive_results

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
subroutine crystallite_init

  integer :: &
    Nconstituents, &
    ph, &
    me, &
    co, &                                                                                           !< counter in integration point component loop
    ip, &                                                                                           !< counter in integration point loop
    el, &                                                                                           !< counter in element loop
    cMax, &                                                                                         !< maximum number of  integration point components
    iMax, &                                                                                         !< maximum number of integration points
    eMax, &                                                                                         !< maximum number of elements
    so

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

  allocate(crystallite_F(3,3,cMax,iMax,eMax),source=0.0_pReal)

  allocate(crystallite_S0, &
           crystallite_F0, &
           crystallite_partitionedS0, &
           crystallite_partitionedF0,&
           crystallite_S,crystallite_P, &
           source = crystallite_F)

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

  allocate(constitutive_mech_Fe(phases%length))
  allocate(constitutive_mech_Fi(phases%length))
  allocate(constitutive_mech_Fi0(phases%length))
  allocate(constitutive_mech_partitionedFi0(phases%length))
  allocate(constitutive_mech_Fp(phases%length))
  allocate(constitutive_mech_Fp0(phases%length))
  allocate(constitutive_mech_partitionedFp0(phases%length))
  allocate(constitutive_mech_Li(phases%length))
  allocate(constitutive_mech_Li0(phases%length))
  allocate(constitutive_mech_partitionedLi0(phases%length))
  allocate(constitutive_mech_partitionedLp0(phases%length))
  allocate(constitutive_mech_Lp0(phases%length))
  allocate(constitutive_mech_Lp(phases%length))
  do ph = 1, phases%length
    Nconstituents = count(material_phaseAt == ph) * discretization_nIPs

    allocate(constitutive_mech_Fi(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Fe(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Fi0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partitionedFi0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Fp(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Fp0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partitionedFp0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Li(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Li0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partitionedLi0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partitionedLp0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Lp0(ph)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Lp(ph)%data(3,3,Nconstituents))
    do so = 1, phase_Nsources(ph)
      allocate(sourceState(ph)%p(so)%subState0,source=sourceState(ph)%p(so)%state0)                 ! ToDo: hack
    enddo
  enddo

  print'(a42,1x,i10)', '    # of elements:                       ', eMax
  print'(a42,1x,i10)', '    # of integration points/element:     ', iMax
  print'(a42,1x,i10)', 'max # of constituents/integration point: ', cMax
  flush(IO_STDOUT)

 !$OMP PARALLEL DO PRIVATE(ph,me)
  do el = 1, size(material_phaseMemberAt,3); do ip = 1, size(material_phaseMemberAt,2)
    do co = 1, homogenization_Nconstituents(material_homogenizationAt(el))

      ph = material_phaseAt(co,el)
      me = material_phaseMemberAt(co,ip,el)
      constitutive_mech_Fp0(ph)%data(1:3,1:3,me) = material_orientation0(co,ip,el)%asMatrix()                      ! Fp reflects initial orientation (see 10.1016/j.actamat.2006.01.005)
      constitutive_mech_Fp0(ph)%data(1:3,1:3,me) = constitutive_mech_Fp0(ph)%data(1:3,1:3,me) &
                                     / math_det33(constitutive_mech_Fp0(ph)%data(1:3,1:3,me))**(1.0_pReal/3.0_pReal)
      constitutive_mech_Fi0(ph)%data(1:3,1:3,me) = math_I3

      crystallite_F0(1:3,1:3,co,ip,el)  = math_I3

      constitutive_mech_Fe(ph)%data(1:3,1:3,me) = math_inv33(matmul(constitutive_mech_Fi0(ph)%data(1:3,1:3,me), &
                                                                    constitutive_mech_Fp0(ph)%data(1:3,1:3,me)))           ! assuming that euler angles are given in internal strain free configuration
      constitutive_mech_Fp(ph)%data(1:3,1:3,me) = constitutive_mech_Fp0(ph)%data(1:3,1:3,me)
      constitutive_mech_Fi(ph)%data(1:3,1:3,me) = constitutive_mech_Fi0(ph)%data(1:3,1:3,me)

      constitutive_mech_partitionedFi0(ph)%data(1:3,1:3,me) = constitutive_mech_Fi0(ph)%data(1:3,1:3,me)
      constitutive_mech_partitionedFp0(ph)%data(1:3,1:3,me) = constitutive_mech_Fp0(ph)%data(1:3,1:3,me)

    enddo
  enddo; enddo
  !$OMP END PARALLEL DO

  crystallite_partitionedF0  = crystallite_F0
  crystallite_F   = crystallite_F0


  !$OMP PARALLEL DO PRIVATE(ph,me)
  do el = 1, size(material_phaseMemberAt,3)
    do ip = 1, size(material_phaseMemberAt,2)
      do co = 1,homogenization_Nconstituents(material_homogenizationAt(el))
        ph = material_phaseAt(co,el)
        me = material_phaseMemberAt(co,ip,el)
        call crystallite_orientations(co,ip,el)
        call constitutive_plastic_dependentState(crystallite_partitionedF0(1:3,1:3,co,ip,el),co,ip,el)  ! update dependent state variables to be consistent with basic states
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
    crystallite_partitionedF0(1:3,1:3,co,ip,el)  = crystallite_F0(1:3,1:3,co,ip,el)
    crystallite_partitionedS0(1:3,1:3,co,ip,el)  = crystallite_S0(1:3,1:3,co,ip,el)

    call mech_initializeRestorationPoints(ph,me)

    do so = 1, phase_Nsources(material_phaseAt(co,el))
      sourceState(material_phaseAt(co,el))%p(so)%partitionedState0(:,material_phasememberAt(co,ip,el)) = &
      sourceState(material_phaseAt(co,el))%p(so)%state0(           :,material_phasememberAt(co,ip,el))
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
    crystallite_partitionedF0 (1:3,1:3,co,ip,el) = crystallite_F (1:3,1:3,co,ip,el)
    crystallite_partitionedS0 (1:3,1:3,co,ip,el) = crystallite_S (1:3,1:3,co,ip,el)

    call constitutive_mech_windForward(ph,me)
    do so = 1, phase_Nsources(material_phaseAt(co,el))
      sourceState(ph)%p(so)%partitionedState0(:,me) = sourceState(ph)%p(so)%state(:,me)
    enddo
  enddo

end subroutine constitutive_windForward


!--------------------------------------------------------------------------------------------------
!> @brief Calculate tangent (dPdF).
!--------------------------------------------------------------------------------------------------
function crystallite_stressTangent(dt,co,ip,el) result(dPdF)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
    co, &                                                                                            !< counter in constituent loop
    ip, &                                                                                            !< counter in integration point loop
    el                                                                                               !< counter in element loop
  real(pReal), dimension(3,3,3,3) :: dPdF

  integer :: &
    o, &
    p, ph, me
  real(pReal), dimension(3,3)     ::   devNull, &
                                       invSubFp0,invSubFi0,invFp,invFi, &
                                       temp_33_1, temp_33_2, temp_33_3
  real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                       dSdF, &
                                       dSdFi, &
                                       dLidS, &                                                     ! tangent in lattice configuration
                                       dLidFi, &
                                       dLpdS, &
                                       dLpdFi, &
                                       dFidS, &
                                       dFpinvdF, &
                                       rhs_3333, &
                                       lhs_3333, &
                                       temp_3333
  real(pReal), dimension(9,9)::        temp_99
  logical :: error


  ph = material_phaseAt(co,el)
  me = material_phaseMemberAt(co,ip,el)

  call constitutive_hooke_SandItsTangents(devNull,dSdFe,dSdFi, &
                                          constitutive_mech_Fp(ph)%data(1:3,1:3,me), &
                                          constitutive_mech_Fi(ph)%data(1:3,1:3,me),co,ip,el)
  call constitutive_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                     crystallite_S (1:3,1:3,co,ip,el), &
                                     constitutive_mech_Fi(ph)%data(1:3,1:3,me), &
                                     co,ip,el)

  invFp = math_inv33(constitutive_mech_Fp(ph)%data(1:3,1:3,me))
  invFi = math_inv33(constitutive_mech_Fi(ph)%data(1:3,1:3,me))
  invSubFp0 = math_inv33(constitutive_mech_partitionedFp0(ph)%data(1:3,1:3,me))
  invSubFi0 = math_inv33(constitutive_mech_partitionedFi0(ph)%data(1:3,1:3,me))

  if (sum(abs(dLidS)) < tol_math_check) then
    dFidS = 0.0_pReal
  else
    lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
    do o=1,3; do p=1,3
      lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                            + matmul(invSubFi0,dLidFi(1:3,1:3,o,p)) * dt
      lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                            + invFi*invFi(p,o)
      rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                            - matmul(invSubFi0,dLidS(1:3,1:3,o,p)) * dt
    enddo; enddo
    call math_invert(temp_99,error,math_3333to99(lhs_3333))
    if (error) then
      call IO_warning(warning_ID=600,el=el,ip=ip,g=co, &
                      ext_msg='inversion error in analytic tangent calculation')
      dFidS = 0.0_pReal
    else
      dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
    endif
    dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
  endif

  call constitutive_plastic_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                             crystallite_S (1:3,1:3,co,ip,el), &
                                             constitutive_mech_Fi(ph)%data(1:3,1:3,me),co,ip,el)
  dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
  temp_33_1 = transpose(matmul(invFp,invFi))
  temp_33_2 = matmul(crystallite_F(1:3,1:3,co,ip,el),invSubFp0)
  temp_33_3 = matmul(matmul(crystallite_F(1:3,1:3,co,ip,el),invFp), invSubFi0)

  do o=1,3; do p=1,3
    rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
    temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), invFi) &
                           + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
  enddo; enddo
  lhs_3333 = math_mul3333xx3333(dSdFe,temp_3333) * dt &
           + math_mul3333xx3333(dSdFi,dFidS)

  call math_invert(temp_99,error,math_eye(9)+math_3333to99(lhs_3333))
  if (error) then
    call IO_warning(warning_ID=600,el=el,ip=ip,g=co, &
                    ext_msg='inversion error in analytic tangent calculation')
    dSdF = rhs_3333
  else
    dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
  endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
  temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
  do o=1,3; do p=1,3
    dFpinvdF(1:3,1:3,p,o) = - matmul(invSubFp0, matmul(temp_3333(1:3,1:3,p,o),invFi)) * dt
  enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
  temp_33_1 = matmul(crystallite_S(1:3,1:3,co,ip,el),transpose(invFp))
  temp_33_2 = matmul(crystallite_F(1:3,1:3,co,ip,el),invFp)
  temp_33_3 = matmul(temp_33_2,crystallite_S(1:3,1:3,co,ip,el))

  dPdF = 0.0_pReal
  do p=1,3
    dPdF(p,1:3,p,1:3) = transpose(matmul(invFp,temp_33_1))
  enddo
  do o=1,3; do p=1,3
    dPdF(1:3,1:3,p,o) = dPdF(1:3,1:3,p,o) &
                      + matmul(matmul(crystallite_F(1:3,1:3,co,ip,el),dFpinvdF(1:3,1:3,p,o)),temp_33_1) &
                      + matmul(matmul(temp_33_2,dSdF(1:3,1:3,p,o)),transpose(invFp)) &
                      + matmul(temp_33_3,transpose(dFpinvdF(1:3,1:3,p,o)))
  enddo; enddo

end function crystallite_stressTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations(co,ip,el)

  integer, intent(in) :: &
    co, &                                                                                            !< counter in integration point component loop
    ip, &                                                                                            !< counter in integration point loop
    el                                                                                               !< counter in element loop


  call crystallite_orientation(co,ip,el)%fromMatrix(transpose(math_rotationalPart(&
    constitutive_mech_Fe(material_phaseAt(ip,el))%data(1:3,1:3,material_phaseMemberAt(co,ip,el)))))

  if (plasticState(material_phaseAt(1,el))%nonlocal) &
    call plastic_nonlocal_updateCompatibility(crystallite_orientation, &
                                              phase_plasticityInstance(material_phaseAt(1,el)),ip,el)


end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(co,ip,el, tensor33)

  real(pReal), dimension(3,3), intent(in) :: tensor33
  real(pReal), dimension(3,3)             :: T
  integer, intent(in):: &
    el, &
    ip, &
    co

  real(pReal), dimension(3,3) :: crystallite_push33ToRef


  T = matmul(material_orientation0(co,ip,el)%asMatrix(), &                                         ! ToDo: initial orientation correct?
             transpose(math_inv33(crystallite_F(1:3,1:3,co,ip,el))))
  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
function integrateSourceState(dt,co,ip,el) result(broken)

  real(pReal), intent(in) :: dt
  integer, intent(in) :: &
    el, &                                                                                            !< element index in element loop
    ip, &                                                                                            !< integration point index in ip loop
    co                                                                                               !< grain index in grain loop

  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    ph, &
    me, &
    so
  integer, dimension(maxval(phase_Nsources)) :: &
    size_so
  real(pReal) :: &
    zeta
  real(pReal), dimension(constitutive_source_maxSizeDotState) :: &
    r                                                                                               ! state residuum
  real(pReal), dimension(constitutive_source_maxSizeDotState,2,maxval(phase_Nsources)) :: source_dotState
  logical :: &
    broken, converged_


  ph = material_phaseAt(co,el)
  me = material_phaseMemberAt(co,ip,el)

  converged_ = .true.
  broken = constitutive_thermal_collectDotState(ph,me)
  broken = broken .or. constitutive_damage_collectDotState(crystallite_S(1:3,1:3,co,ip,el), co,ip,el,ph,me)
  if(broken) return

  do so = 1, phase_Nsources(ph)
    size_so(so) = sourceState(ph)%p(so)%sizeDotState
    sourceState(ph)%p(so)%state(1:size_so(so),me) = sourceState(ph)%p(so)%subState0(1:size_so(so),me) &
                                                  + sourceState(ph)%p(so)%dotState (1:size_so(so),me) * dt
    source_dotState(1:size_so(so),2,so) = 0.0_pReal
  enddo

  iteration: do NiterationState = 1, num%nState

    do so = 1, phase_Nsources(ph)
      if(nIterationState > 1) source_dotState(1:size_so(so),2,so) = source_dotState(1:size_so(so),1,so)
      source_dotState(1:size_so(so),1,so) = sourceState(ph)%p(so)%dotState(:,me)
    enddo

    broken = constitutive_thermal_collectDotState(ph,me)
    broken = broken .or. constitutive_damage_collectDotState(crystallite_S(1:3,1:3,co,ip,el), co,ip,el,ph,me)
    if(broken) exit iteration

    do so = 1, phase_Nsources(ph)
      zeta = damper(sourceState(ph)%p(so)%dotState(:,me), &
                    source_dotState(1:size_so(so),1,so),&
                    source_dotState(1:size_so(so),2,so))
      sourceState(ph)%p(so)%dotState(:,me) = sourceState(ph)%p(so)%dotState(:,me) * zeta &
                                        + source_dotState(1:size_so(so),1,so)* (1.0_pReal - zeta)
      r(1:size_so(so)) = sourceState(ph)%p(so)%state    (1:size_so(so),me)  &
                       - sourceState(ph)%p(so)%subState0(1:size_so(so),me)  &
                       - sourceState(ph)%p(so)%dotState (1:size_so(so),me) * dt
      sourceState(ph)%p(so)%state(1:size_so(so),me) = sourceState(ph)%p(so)%state(1:size_so(so),me) &
                                                - r(1:size_so(so))
      converged_ = converged_  .and. converged(r(1:size_so(so)), &
                                               sourceState(ph)%p(so)%state(1:size_so(so),me), &
                                               sourceState(ph)%p(so)%atol(1:size_so(so)))
    enddo

    if(converged_) then
      broken = constitutive_damage_deltaState(constitutive_mech_Fe(ph)%data(1:3,1:3,me),co,ip,el,ph,me)
      exit iteration
    endif

  enddo iteration

  broken = broken .or. .not. converged_


  contains

  !--------------------------------------------------------------------------------------------------
  !> @brief calculate the damping for correction of state and dot state
  !--------------------------------------------------------------------------------------------------
  real(pReal) pure function damper(current,previous,previous2)

  real(pReal), dimension(:), intent(in) ::&
    current, previous, previous2

  real(pReal) :: dot_prod12, dot_prod22

  dot_prod12 = dot_product(current  - previous,  previous - previous2)
  dot_prod22 = dot_product(previous - previous2, previous - previous2)
  if ((dot_product(current,previous) < 0.0_pReal .or. dot_prod12 < 0.0_pReal) .and. dot_prod22 > 0.0_pReal) then
    damper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
  else
    damper = 1.0_pReal
  endif

  end function damper

end function integrateSourceState


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
! ToDo: Merge data into one file for MPI, move state to constitutive and homogenization, respectively
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restartWrite

  integer :: ph
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print*, ' writing field and constitutive data required for restart to file';flush(IO_STDOUT)

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName,'a')

  call HDF5_write(fileHandle,crystallite_F,'F')
  call HDF5_write(fileHandle,crystallite_S,           'S')

  groupHandle = HDF5_addGroup(fileHandle,'phase')
  do ph = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') ph,'_omega'
    call HDF5_write(groupHandle,plasticState(ph)%state,datasetName)
    write(datasetName,'(i0,a)') ph,'_F_i'
    call HDF5_write(groupHandle,constitutive_mech_Fi(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_L_i'
    call HDF5_write(groupHandle,constitutive_mech_Li(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_L_p'
    call HDF5_write(groupHandle,constitutive_mech_Lp(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_F_p'
    call HDF5_write(groupHandle,constitutive_mech_Fp(ph)%data,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_addGroup(fileHandle,'homogenization')
  do ph = 1, size(material_name_homogenization)
    write(datasetName,'(i0,a)') ph,'_omega'
    call HDF5_write(groupHandle,homogState(ph)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read data for restart
! ToDo: Merge data into one file for MPI, move state to constitutive and homogenization, respectively
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restartRead

  integer :: ph
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print'(/,a,i0,a)', ' reading restart information of increment from file'

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName)

  call HDF5_read(fileHandle,crystallite_F0, 'F')
  call HDF5_read(fileHandle,crystallite_S0, 'S')

  groupHandle = HDF5_openGroup(fileHandle,'phase')
  do ph = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') ph,'_omega'
    call HDF5_read(groupHandle,plasticState(ph)%state0,datasetName)
    write(datasetName,'(i0,a)') ph,'_F_i'
    call HDF5_read(groupHandle,constitutive_mech_Fi0(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_L_i'
    call HDF5_read(groupHandle,constitutive_mech_Li0(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_L_p'
    call HDF5_read(groupHandle,constitutive_mech_Lp0(ph)%data,datasetName)
    write(datasetName,'(i0,a)') ph,'_F_p'
    call HDF5_read(groupHandle,constitutive_mech_Fp0(ph)%data,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_openGroup(fileHandle,'homogenization')
  do ph = 1,size(material_name_homogenization)
    write(datasetName,'(i0,a)') ph,'_omega'
    call HDF5_read(groupHandle,homogState(ph)%state0,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartRead


! getter for non-mech (e.g. thermal)
function constitutive_mech_getS(co,ip,el) result(S)

  integer, intent(in) :: co, ip, el
  real(pReal), dimension(3,3) :: S

  S = crystallite_S(1:3,1:3,co,ip,el)

end function constitutive_mech_getS

! getter for non-mech (e.g. thermal)
function constitutive_mech_getLp(co,ip,el) result(Lp)

  integer, intent(in) :: co, ip, el
  real(pReal), dimension(3,3) :: Lp

  Lp = constitutive_mech_Lp(material_phaseAt(co,el))%data(1:3,1:3,material_phaseMemberAt(co,ip,el))

end function constitutive_mech_getLp

end module constitutive
