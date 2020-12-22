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
  use FEsolving
  use results

  implicit none
  private
  real(pReal),               dimension(:,:,:),        allocatable, public :: &
    crystallite_dt                                                                                  !< requested time increment of each grain
  real(pReal),               dimension(:,:,:),        allocatable :: &
    crystallite_subdt, &                                                                            !< substepped time increment of each grain
    crystallite_subStep                                                                             !< size of next integration step
  type(rotation),            dimension(:,:,:),        allocatable :: &
    crystallite_orientation                                                                         !< current orientation
  real(pReal),               dimension(:,:,:,:,:),    allocatable :: &
    crystallite_F0, &                                                                               !< def grad at start of FE inc
    crystallite_subF,  &                                                                            !< def grad to be reached at end of crystallite inc
    crystallite_subF0, &                                                                            !< def grad at start of crystallite inc
    !
    crystallite_Fe, &                                                                               !< current "elastic" def grad (end of converged time step)
    !
    crystallite_Fp, &                                                                               !< current plastic def grad (end of converged time step)
    crystallite_Fp0, &                                                                              !< plastic def grad at start of FE inc
    crystallite_partitionedFp0,&                                                                    !< plastic def grad at start of homog inc
    crystallite_subFp0,&                                                                            !< plastic def grad at start of crystallite inc
    !
    crystallite_subFi0,&                                                                            !< intermediate def grad at start of crystallite inc
    !
    crystallite_Lp0, &                                                                              !< plastic velocitiy grad at start of FE inc
    crystallite_partitionedLp0, &                                                                   !< plastic velocity grad at start of homog inc
    !
    crystallite_S0, &                                                                               !< 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partitionedS0                                                                       !< 2nd Piola-Kirchhoff stress vector at start of homog inc
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_P, &                                                                                !< 1st Piola-Kirchhoff stress per grain
    crystallite_Lp, &                                                                               !< current plastic velocitiy grad (end of converged time step)
    crystallite_S, &                                                                                !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_partitionedF0                                                                       !< def grad at start of homog inc
  real(pReal),               dimension(:,:,:,:,:),    allocatable, public :: &
    crystallite_partitionedF                                                                        !< def grad to be reached at end of homog inc

  logical,                    dimension(:,:,:),       allocatable, public :: &
    crystallite_requested                                                                           !< used by upper level (homogenization) to request crystallite calculation
  logical,                    dimension(:,:,:),       allocatable :: &
    crystallite_converged                                                                           !< convergence flag



  type :: tTensorContainer
    real(pReal), dimension(:,:,:), allocatable :: data
  end type

  type(tTensorContainer), dimension(:), allocatable :: &
    constitutive_mech_Fi, &
    constitutive_mech_Fi0, &
    constitutive_mech_partionedFi0, &
    constitutive_mech_Li, &
    constitutive_mech_Li0, &
    constitutive_mech_partionedLi0

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

  procedure(integrateStateFPI), pointer :: integrateState

  integer(kind(PLASTICITY_undefined_ID)), dimension(:),   allocatable :: &
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

! == cleaned:end ===================================================================================


module function constitutive_deltaState(S, Fi, ipc, ip, el, phase, of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),   intent(in), dimension(3,3) :: &
    S, &                                                                                            !< 2nd Piola Kirchhoff stress
    Fi                                                                                              !< intermediate deformation gradient
  logical :: &
    broken


end function constitutive_deltaState


    module subroutine source_damage_anisoBrittle_dotState(S, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),  intent(in), dimension(3,3) :: &
        S
    end subroutine source_damage_anisoBrittle_dotState

    module subroutine source_damage_anisoDuctile_dotState(ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end subroutine source_damage_anisoDuctile_dotState

    module subroutine source_damage_isoDuctile_dotState(ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
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

    module subroutine constitutive_thermal_getRateAndItsTangents(TDot, dTDot_dT, T, S, Lp, ip, el)
      integer, intent(in) :: &
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal), intent(in) :: &
        T
      real(pReal),  intent(in), dimension(:,:,:,:,:) :: &
        S, &                                                                                        !< current 2nd Piola Kitchoff stress vector
        Lp                                                                                          !< plastic velocity gradient
      real(pReal), intent(inout) :: &
        TDot, &
        dTDot_dT
    end subroutine constitutive_thermal_getRateAndItsTangents

    module function plastic_dislotwin_homogenizedC(ipc,ip,el) result(homogenizedC)
      real(pReal), dimension(6,6) :: &
        homogenizedC
      integer,     intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
    end function plastic_dislotwin_homogenizedC

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

    module subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_cleavage_opening_LiAndItsTangent

    module subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(in),  dimension(3,3) :: &
        S
      real(pReal),   intent(out), dimension(3,3) :: &
        Ld                                                                                          !< damage velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLd_dTstar                                                                                  !< derivative of Ld with respect to Tstar (4th-order tensor)
    end subroutine kinematics_slipplane_opening_LiAndItsTangent

    module subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< grain number
        ip, &                                                                                       !< integration point number
        el                                                                                          !< element number
      real(pReal),   intent(out), dimension(3,3) :: &
        Li                                                                                          !< thermal velocity gradient
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dLi_dTstar                                                                                  !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)
    end subroutine kinematics_thermal_expansion_LiAndItsTangent


    module subroutine source_damage_isoBrittle_deltaState(C, Fe, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),  intent(in), dimension(3,3) :: &
        Fe
      real(pReal),  intent(in), dimension(6,6) :: &
        C
    end subroutine source_damage_isoBrittle_deltaState


    module subroutine constitutive_plastic_LpAndItsTangents(Lp, dLp_dS, dLp_dFi, &
                                         S, Fi, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
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


    module subroutine constitutive_plastic_dependentState(F, Fp, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),   intent(in), dimension(3,3) :: &
        F, &                                                                                        !< elastic deformation gradient
        Fp                                                                                          !< plastic deformation gradient
    end subroutine constitutive_plastic_dependentState



    module subroutine constitutive_hooke_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
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

    module subroutine integrateStateFPI(g,i,e)
      integer, intent(in) :: e, i, g
    end subroutine integrateStateFPI

    module subroutine integrateStateEuler(g,i,e)
      integer, intent(in) :: e, i, g
    end subroutine integrateStateEuler

    module subroutine integrateStateAdaptiveEuler(g,i,e)
      integer, intent(in) :: e, i, g
    end subroutine integrateStateAdaptiveEuler

    module subroutine integrateStateRK4(g,i,e)
      integer, intent(in) :: e, i, g
    end subroutine integrateStateRK4

    module subroutine integrateStateRKCK45(g,i,e)
      integer, intent(in) :: e, i, g
    end subroutine integrateStateRKCK45

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
    crystallite_restartRead, &
    crystallite_forward, &
    constitutive_initializeRestorationPoints, &
    constitutive_windForward, &
    crystallite_restore

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialze constitutive models for individual physics
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init

  integer :: &
    p, &                                                                                            !< counter in phase loop
    s                                                                                               !< counter in source loop
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
  call mech_init
  call damage_init
  call thermal_init

  print'(/,a)', ' <<<+-  constitutive init  -+>>>'; flush(IO_STDOUT)

  phases => config_material%get('phase')
  constitutive_source_maxSizeDotState = 0
  PhaseLoop2:do p = 1,phases%length
!--------------------------------------------------------------------------------------------------
! partition and initialize state
    plasticState(p)%partitionedState0 = plasticState(p)%state0
    plasticState(p)%state             = plasticState(p)%partitionedState0
    forall(s = 1:phase_Nsources(p))
      sourceState(p)%p(s)%partitionedState0 = sourceState(p)%p(s)%state0
      sourceState(p)%p(s)%state             = sourceState(p)%p(s)%partitionedState0
    end forall

    constitutive_source_maxSizeDotState   = max(constitutive_source_maxSizeDotState, &
                                                maxval(sourceState(p)%p%sizeDotState))
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
!> @brief returns the homogenize elasticity matrix
!> ToDo: homogenizedC66 would be more consistent
!--------------------------------------------------------------------------------------------------
function constitutive_homogenizedC(ipc,ip,el)

  real(pReal), dimension(6,6) :: &
    constitutive_homogenizedC
  integer,      intent(in)     :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element

  plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
    case (PLASTICITY_DISLOTWIN_ID) plasticityType
     constitutive_homogenizedC = plastic_dislotwin_homogenizedC(ipc,ip,el)
    case default plasticityType
     constitutive_homogenizedC = lattice_C66(1:6,1:6,material_phaseAt(ipc,el))
  end select plasticityType

end function constitutive_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
! ToDo: MD: S is Mi?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LiAndItsTangents(Li, dLi_dS, dLi_dFi, &
                                         S, Fi, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
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

  plasticityType: select case (phase_plasticity(material_phaseAt(ipc,el)))
    case (PLASTICITY_isotropic_ID) plasticityType
      of = material_phasememberAt(ipc,ip,el)
      instance = phase_plasticityInstance(material_phaseAt(ipc,el))
      call plastic_isotropic_LiAndItsTangent(my_Li, my_dLi_dS, S ,instance,of)
    case default plasticityType
      my_Li = 0.0_pReal
      my_dLi_dS = 0.0_pReal
  end select plasticityType

  Li = Li + my_Li
  dLi_dS = dLi_dS + my_dLi_dS

  KinematicsLoop: do k = 1, phase_Nkinematics(material_phaseAt(ipc,el))
    kinematicsType: select case (phase_kinematics(k,material_phaseAt(ipc,el)))
      case (KINEMATICS_cleavage_opening_ID) kinematicsType
        call kinematics_cleavage_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ipc, ip, el)
      case (KINEMATICS_slipplane_opening_ID) kinematicsType
        call kinematics_slipplane_opening_LiAndItsTangent(my_Li, my_dLi_dS, S, ipc, ip, el)
      case (KINEMATICS_thermal_expansion_ID) kinematicsType
        call kinematics_thermal_expansion_LiAndItsTangent(my_Li, my_dLi_dS, ipc, ip, el)
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
function constitutive_damage_collectDotState(S, ipc, ip, el,phase,of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),  intent(in), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
  integer :: &
    i                                                                                               !< counter in source loop
  logical :: broken


  broken = .false.

  SourceLoop: do i = 1, phase_Nsources(phase)

    sourceType: select case (phase_source(i,phase))

      case (SOURCE_damage_anisoBrittle_ID) sourceType
        call source_damage_anisoBrittle_dotState(S, ipc, ip, el) ! correct stress?

      case (SOURCE_damage_isoDuctile_ID) sourceType
        call source_damage_isoDuctile_dotState(ipc, ip, el)

      case (SOURCE_damage_anisoDuctile_ID) sourceType
        call source_damage_anisoDuctile_dotState(ipc, ip, el)

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(sourceState(phase)%p(i)%dotState(:,of)))

  enddo SourceLoop

end function constitutive_damage_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_thermal_collectDotState(S, ipc, ip, el,phase,of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),  intent(in), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
  integer :: &
    i                                                                                               !< counter in source loop
  logical :: broken


  broken = .false.

  SourceLoop: do i = 1, phase_Nsources(phase)

    sourceType: select case (phase_source(i,phase))

      case (SOURCE_thermal_externalheat_ID) sourceType
        call source_thermal_externalheat_dotState(phase,of)

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(sourceState(phase)%p(i)%dotState(:,of)))

  enddo SourceLoop

end function constitutive_thermal_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
function constitutive_deltaState_source(Fe, ipc, ip, el, phase, of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),   intent(in), dimension(3,3) :: &
    Fe                                                                                              !< elastic deformation gradient
  integer :: &
    i, &
    myOffset, &
    mySize
  logical :: &
    broken


  broken = .false.

  sourceLoop: do i = 1, phase_Nsources(phase)

     sourceType: select case (phase_source(i,phase))

      case (SOURCE_damage_isoBrittle_ID) sourceType
        call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(ipc,ip,el), Fe, &
                                                   ipc, ip, el)
        broken = any(IEEE_is_NaN(sourceState(phase)%p(i)%deltaState(:,of)))
        if(.not. broken) then
          myOffset = sourceState(phase)%p(i)%offsetDeltaState
          mySize   = sourceState(phase)%p(i)%sizeDeltaState
          sourceState(phase)%p(i)%state(myOffset + 1: myOffset + mySize,of) = &
          sourceState(phase)%p(i)%state(myOffset + 1: myOffset + mySize,of) + sourceState(phase)%p(i)%deltaState(1:mySize,of)
        endif

    end select sourceType

  enddo SourceLoop

end function constitutive_deltaState_source


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

  allocate(state%atol           (sizeState),             source=0.0_pReal)
  allocate(state%state0         (sizeState,Nconstituents), source=0.0_pReal)
  allocate(state%partitionedState0(sizeState,Nconstituents), source=0.0_pReal)
  allocate(state%subState0      (sizeState,Nconstituents), source=0.0_pReal)
  allocate(state%state          (sizeState,Nconstituents), source=0.0_pReal)

  allocate(state%dotState    (sizeDotState,Nconstituents), source=0.0_pReal)

  allocate(state%deltaState(sizeDeltaState,Nconstituents), source=0.0_pReal)


end subroutine constitutive_allocateState


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_restore(i,e)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  integer :: &
    c, &                                                                                            !< constituent number
    s

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%state(          :,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine constitutive_restore


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine constitutive_forward

  integer :: i, j

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
    p, &
    m, &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
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

  allocate(crystallite_partitionedF(3,3,cMax,iMax,eMax),source=0.0_pReal)

  allocate(crystallite_S0, &
           crystallite_F0,crystallite_Fp0,crystallite_Lp0, &
           crystallite_partitionedS0, &
           crystallite_partitionedF0,crystallite_partitionedFp0,&
                                   crystallite_partitionedLp0, &
           crystallite_S,crystallite_P, &
           crystallite_Fe,crystallite_Fp,crystallite_Lp, &
           crystallite_subF,crystallite_subF0, &
           crystallite_subFp0,crystallite_subFi0, &
           source = crystallite_partitionedF)

  allocate(crystallite_dt(cMax,iMax,eMax),source=0.0_pReal)
  allocate(crystallite_subdt,crystallite_subStep, &
           source = crystallite_dt)

  allocate(crystallite_orientation(cMax,iMax,eMax))

  allocate(crystallite_requested(cMax,iMax,eMax),             source=.false.)
  allocate(crystallite_converged(cMax,iMax,eMax),             source=.true.)

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

  select case(num_crystallite%get_asString('integrator',defaultVal='FPI'))
    case('FPI')
      integrateState => integrateStateFPI
    case('Euler')
      integrateState => integrateStateEuler
    case('AdaptiveEuler')
      integrateState => integrateStateAdaptiveEuler
    case('RK4')
      integrateState => integrateStateRK4
    case('RKCK45')
      integrateState => integrateStateRKCK45
    case default
     call IO_error(301,ext_msg='integrator')
  end select

  phases => config_material%get('phase')

  allocate(constitutive_mech_Fi(phases%length))
  allocate(constitutive_mech_Fi0(phases%length))
  allocate(constitutive_mech_partionedFi0(phases%length))
  allocate(constitutive_mech_Li(phases%length))
  allocate(constitutive_mech_Li0(phases%length))
  allocate(constitutive_mech_partionedLi0(phases%length))
  do p = 1, phases%length
    Nconstituents = count(material_phaseAt == p) * discretization_nIPs

    allocate(constitutive_mech_Fi(p)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Fi0(p)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partionedFi0(p)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Li(p)%data(3,3,Nconstituents))
    allocate(constitutive_mech_Li0(p)%data(3,3,Nconstituents))
    allocate(constitutive_mech_partionedLi0(p)%data(3,3,Nconstituents))
  enddo

  print'(a42,1x,i10)', '    # of elements:                       ', eMax
  print'(a42,1x,i10)', '    # of integration points/element:     ', iMax
  print'(a42,1x,i10)', 'max # of constituents/integration point: ', cMax
  flush(IO_STDOUT)

 !$OMP PARALLEL DO PRIVATE(p,m)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1), FEsolving_execIP(2); do c = 1, homogenization_Nconstituents(material_homogenizationAt(e))

      p = material_phaseAt(c,e)
      m = material_phaseMemberAt(c,i,e)
      crystallite_Fp0(1:3,1:3,c,i,e) = material_orientation0(c,i,e)%asMatrix()                      ! Fp reflects initial orientation (see 10.1016/j.actamat.2006.01.005)
      crystallite_Fp0(1:3,1:3,c,i,e) = crystallite_Fp0(1:3,1:3,c,i,e) &
                                     / math_det33(crystallite_Fp0(1:3,1:3,c,i,e))**(1.0_pReal/3.0_pReal)
      constitutive_mech_Fi0(p)%data(1:3,1:3,m) = math_I3

      crystallite_F0(1:3,1:3,c,i,e)  = math_I3

      crystallite_Fe(1:3,1:3,c,i,e)  = math_inv33(matmul(constitutive_mech_Fi0(p)%data(1:3,1:3,m), &
                                                         crystallite_Fp0(1:3,1:3,c,i,e)))           ! assuming that euler angles are given in internal strain free configuration
      crystallite_Fp(1:3,1:3,c,i,e)  = crystallite_Fp0(1:3,1:3,c,i,e)
      constitutive_mech_Fi(p)%data(1:3,1:3,m) = constitutive_mech_Fi0(p)%data(1:3,1:3,m)

      constitutive_mech_partionedFi0(p)%data(1:3,1:3,m) = constitutive_mech_Fi0(p)%data(1:3,1:3,m)


      crystallite_requested(c,i,e) = .true.
    enddo; enddo
  enddo
  !$OMP END PARALLEL DO


  crystallite_partitionedFp0 = crystallite_Fp0
  crystallite_partitionedF0  = crystallite_F0
  crystallite_partitionedF   = crystallite_F0

  call crystallite_orientations()

  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
        call constitutive_plastic_dependentState(crystallite_partitionedF0(1:3,1:3,c,i,e), &
                                         crystallite_partitionedFp0(1:3,1:3,c,i,e), &
                                         c,i,e)                                                     ! update dependent state variables to be consistent with basic states
     enddo
    enddo
  enddo
  !$OMP END PARALLEL DO


end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P)
!--------------------------------------------------------------------------------------------------
function crystallite_stress()

  logical, dimension(discretization_nIPs,discretization_Nelems) :: crystallite_stress
  real(pReal) :: &
    formerSubStep
  integer :: &
    NiterationCrystallite, &                                                                        ! number of iterations in crystallite loop
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e, &                                                                                            !< counter in element loop
    s, p, m
  logical, dimension(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems) :: todo     !ToDo: need to set some values to false for different Ngrains
  real(pReal), dimension(homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems) :: subFrac     !ToDo: need to set some values to false for different Ngrains
  real(pReal),               dimension(:,:,:,:,:),    allocatable :: &
    subLp0,&                                                                                        !< plastic velocity grad at start of crystallite inc
    subLi0                                                                                          !< intermediate velocity grad at start of crystallite inc

  todo = .false.

  allocate(subLi0(3,3,homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems))
  subLp0 = crystallite_partitionedLp0

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
  crystallite_subStep = 0.0_pReal
  !$OMP PARALLEL DO PRIVATE(p,m)
  elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2); do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
            p = material_phaseAt(c,e)
            m = material_phaseMemberAt(c,i,e)
      subLi0(1:3,1:3,c,i,e) = constitutive_mech_partionedLi0(p)%data(1:3,1:3,m)
      homogenizationRequestsCalculation: if (crystallite_requested(c,i,e)) then
        plasticState    (material_phaseAt(c,e))%subState0(      :,material_phaseMemberAt(c,i,e)) = &
        plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phaseMemberAt(c,i,e))

        do s = 1, phase_Nsources(material_phaseAt(c,e))
          sourceState(material_phaseAt(c,e))%p(s)%subState0(      :,material_phaseMemberAt(c,i,e)) = &
          sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phaseMemberAt(c,i,e))
        enddo
        crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_partitionedFp0(1:3,1:3,c,i,e)
        crystallite_subFi0(1:3,1:3,c,i,e) = constitutive_mech_partionedFi0(p)%data(1:3,1:3,m)
        crystallite_subF0(1:3,1:3,c,i,e)  = crystallite_partitionedF0(1:3,1:3,c,i,e)
        subFrac(c,i,e) = 0.0_pReal
        crystallite_subStep(c,i,e) = 1.0_pReal/num%subStepSizeCryst
        todo(c,i,e) = .true.
        crystallite_converged(c,i,e) = .false.                                                      ! pretend failed step of 1/subStepSizeCryst
      endif homogenizationRequestsCalculation
    enddo; enddo
  enddo elementLooping1
  !$OMP END PARALLEL DO

  NiterationCrystallite = 0
  cutbackLooping: do while (any(todo(:,FEsolving_execIP(1):FEsolving_execIP(2),FEsolving_execELem(1):FEsolving_execElem(2))))
    NiterationCrystallite = NiterationCrystallite + 1

#ifdef DEBUG
    if (debugCrystallite%extensive) &
      print'(a,i6)', '<< CRYST stress >> crystallite iteration ',NiterationCrystallite
#endif
    !$OMP PARALLEL DO PRIVATE(formerSubStep,p,m)
    elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
      do i = FEsolving_execIP(1),FEsolving_execIP(2)
        do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
            p = material_phaseAt(c,e)
            m = material_phaseMemberAt(c,i,e)
!--------------------------------------------------------------------------------------------------
!  wind forward
          if (crystallite_converged(c,i,e)) then
            formerSubStep = crystallite_subStep(c,i,e)
            subFrac(c,i,e) = subFrac(c,i,e) + crystallite_subStep(c,i,e)
            crystallite_subStep(c,i,e) = min(1.0_pReal - subFrac(c,i,e), &
                                             num%stepIncreaseCryst * crystallite_subStep(c,i,e))

            todo(c,i,e) = crystallite_subStep(c,i,e) > 0.0_pReal                        ! still time left to integrate on?

            if (todo(c,i,e)) then
              crystallite_subF0 (1:3,1:3,c,i,e) = crystallite_subF(1:3,1:3,c,i,e)
              subLp0(1:3,1:3,c,i,e) = crystallite_Lp  (1:3,1:3,c,i,e)
              subLi0(1:3,1:3,c,i,e) = constitutive_mech_Li(p)%data(1:3,1:3,m)
              crystallite_subFp0(1:3,1:3,c,i,e) = crystallite_Fp  (1:3,1:3,c,i,e)
              crystallite_subFi0(1:3,1:3,c,i,e) = constitutive_mech_Fi(p)%data(1:3,1:3,m)
              plasticState(    material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e)) &
                = plasticState(material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e))
              do s = 1, phase_Nsources(material_phaseAt(c,e))
                sourceState(    material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e)) &
                  = sourceState(material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e))
              enddo
            endif

!--------------------------------------------------------------------------------------------------
!  cut back (reduced time and restore)
          else
            crystallite_subStep(c,i,e)       = num%subStepSizeCryst * crystallite_subStep(c,i,e)
            crystallite_Fp   (1:3,1:3,c,i,e) =            crystallite_subFp0(1:3,1:3,c,i,e)
            constitutive_mech_Fi(p)%data(1:3,1:3,m) =            crystallite_subFi0(1:3,1:3,c,i,e)
            crystallite_S    (1:3,1:3,c,i,e) =            crystallite_S0    (1:3,1:3,c,i,e)
            if (crystallite_subStep(c,i,e) < 1.0_pReal) then                                        ! actual (not initial) cutback
              crystallite_Lp (1:3,1:3,c,i,e) =            subLp0(1:3,1:3,c,i,e)
              constitutive_mech_Li(p)%data(1:3,1:3,m) =            subLi0(1:3,1:3,c,i,e)
            endif
            plasticState    (material_phaseAt(c,e))%state(    :,material_phaseMemberAt(c,i,e)) &
              = plasticState(material_phaseAt(c,e))%subState0(:,material_phaseMemberAt(c,i,e))
            do s = 1, phase_Nsources(material_phaseAt(c,e))
              sourceState(    material_phaseAt(c,e))%p(s)%state(    :,material_phaseMemberAt(c,i,e)) &
                = sourceState(material_phaseAt(c,e))%p(s)%subState0(:,material_phaseMemberAt(c,i,e))
            enddo

                                                                                                    ! cant restore dotState here, since not yet calculated in first cutback after initialization
            todo(c,i,e) = crystallite_subStep(c,i,e) > num%subStepMinCryst                          ! still on track or already done (beyond repair)
          endif

!--------------------------------------------------------------------------------------------------
!  prepare for integration
          if (todo(c,i,e)) then
            crystallite_subF(1:3,1:3,c,i,e) = crystallite_subF0(1:3,1:3,c,i,e) &
                                            + crystallite_subStep(c,i,e) *( crystallite_partitionedF (1:3,1:3,c,i,e) &
                                                                           -crystallite_partitionedF0(1:3,1:3,c,i,e))
            crystallite_Fe(1:3,1:3,c,i,e) = matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                                   math_inv33(matmul(constitutive_mech_Fi(p)%data(1:3,1:3,m), &
                                                                     crystallite_Fp(1:3,1:3,c,i,e))))
            crystallite_subdt(c,i,e) = crystallite_subStep(c,i,e) * crystallite_dt(c,i,e)
            crystallite_converged(c,i,e) = .false.
            call integrateState(c,i,e)
            call integrateSourceState(c,i,e)
          endif

        enddo
      enddo
    enddo elementLooping3
    !$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------
!  integrate --- requires fully defined state array (basic + dependent state)
    where(.not. crystallite_converged .and. crystallite_subStep > num%subStepMinCryst) &            ! do not try non-converged but fully cutbacked any further
      todo = .true.                                                                                 ! TODO: again unroll this into proper elementloop to avoid N^2 for single point evaluation
  enddo cutbackLooping

! return whether converged or not
  crystallite_stress = .false.
  elementLooping5: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      crystallite_stress(i,e) = all(crystallite_converged(:,i,e))
    enddo
  enddo elementLooping5

end function crystallite_stress


!--------------------------------------------------------------------------------------------------
!> @brief Backup data for homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_initializeRestorationPoints(i,e)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  integer :: &
    c, &                                                                                            !< constituent number
    s,ph, me

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    ph = material_phaseAt(c,e)
    me = material_phaseMemberAt(c,i,e)
    crystallite_partitionedFp0(1:3,1:3,c,i,e) = crystallite_Fp0(1:3,1:3,c,i,e)
    crystallite_partitionedLp0(1:3,1:3,c,i,e) = crystallite_Lp0(1:3,1:3,c,i,e)
    crystallite_partitionedF0(1:3,1:3,c,i,e)  = crystallite_F0(1:3,1:3,c,i,e)
    crystallite_partitionedS0(1:3,1:3,c,i,e)  = crystallite_S0(1:3,1:3,c,i,e)

    call mech_initializeRestorationPoints(ph,me)
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%state0(           :,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine constitutive_initializeRestorationPoints


!--------------------------------------------------------------------------------------------------
!> @brief Wind homog inc forward.
!--------------------------------------------------------------------------------------------------
subroutine constitutive_windForward(i,e)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  integer :: &
    c, &                                                                                            !< constituent number
    s, p, m
  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
    p = material_phaseAt(c,e)
    m = material_phaseMemberAt(c,i,e)
    crystallite_partitionedF0 (1:3,1:3,c,i,e) = crystallite_partitionedF(1:3,1:3,c,i,e)
    crystallite_partitionedFp0(1:3,1:3,c,i,e) = crystallite_Fp        (1:3,1:3,c,i,e)
    crystallite_partitionedLp0(1:3,1:3,c,i,e) = crystallite_Lp        (1:3,1:3,c,i,e)
    constitutive_mech_partionedFi0(p)%data(1:3,1:3,m) = constitutive_mech_Fi(p)%data(1:3,1:3,m)
    constitutive_mech_partionedLi0(p)%data(1:3,1:3,m) = constitutive_mech_Li(p)%data(1:3,1:3,m)
    crystallite_partitionedS0 (1:3,1:3,c,i,e) = crystallite_S         (1:3,1:3,c,i,e)

    plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phasememberAt(c,i,e)) = &
    plasticState    (material_phaseAt(c,e))%state          (:,material_phasememberAt(c,i,e))
    do s = 1, phase_Nsources(material_phaseAt(c,e))
      sourceState(material_phaseAt(c,e))%p(s)%partitionedState0(:,material_phasememberAt(c,i,e)) = &
      sourceState(material_phaseAt(c,e))%p(s)%state          (:,material_phasememberAt(c,i,e))
    enddo
  enddo

end subroutine constitutive_windForward


!--------------------------------------------------------------------------------------------------
!> @brief Restore data after homog cutback.
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restore(i,e,includeL)

  integer, intent(in) :: &
    i, &                                                                                            !< integration point number
    e                                                                                               !< element number
  logical, intent(in) :: &
    includeL                                                                                        !< protect agains fake cutback
  integer :: &
    c, p, m                                                                                            !< constituent number

  do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
  p = material_phaseAt(c,e)
  m = material_phaseMemberAt(c,i,e)
    if (includeL) then
      crystallite_Lp(1:3,1:3,c,i,e) = crystallite_partitionedLp0(1:3,1:3,c,i,e)
      constitutive_mech_Li(p)%data(1:3,1:3,m) = constitutive_mech_partionedLi0(p)%data(1:3,1:3,m)
    endif                                                                                           ! maybe protecting everything from overwriting makes more sense

    crystallite_Fp(1:3,1:3,c,i,e)   = crystallite_partitionedFp0(1:3,1:3,c,i,e)
    constitutive_mech_Fi(p)%data(1:3,1:3,m)   = constitutive_mech_partionedFi0(p)%data(1:3,1:3,m)
    crystallite_S (1:3,1:3,c,i,e)   = crystallite_partitionedS0 (1:3,1:3,c,i,e)

    plasticState    (material_phaseAt(c,e))%state(          :,material_phasememberAt(c,i,e)) = &
    plasticState    (material_phaseAt(c,e))%partitionedState0(:,material_phasememberAt(c,i,e))
  enddo

end subroutine crystallite_restore


!--------------------------------------------------------------------------------------------------
!> @brief Calculate tangent (dPdF).
!--------------------------------------------------------------------------------------------------
function crystallite_stressTangent(c,i,e) result(dPdF)

  real(pReal), dimension(3,3,3,3) :: dPdF
  integer, intent(in) :: &
    c, &                                                                                            !< counter in constituent loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop
  integer :: &
    o, &
    p, pp, m

  real(pReal), dimension(3,3)     ::   devNull, &
                                       invSubFp0,invSubFi0,invFp,invFi, &
                                       temp_33_1, temp_33_2, temp_33_3, temp_33_4
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

  pp = material_phaseAt(c,e)
  m = material_phaseMemberAt(c,i,e)

        call constitutive_hooke_SandItsTangents(devNull,dSdFe,dSdFi, &
                                         crystallite_Fe(1:3,1:3,c,i,e), &
                                         constitutive_mech_Fi(pp)%data(1:3,1:3,m),c,i,e)
        call constitutive_LiAndItsTangents(devNull,dLidS,dLidFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           constitutive_mech_Fi(pp)%data(1:3,1:3,m), &
                                           c,i,e)

        invFp = math_inv33(crystallite_Fp(1:3,1:3,c,i,e))
        invFi = math_inv33(constitutive_mech_Fi(pp)%data(1:3,1:3,m))
        invSubFp0 = math_inv33(crystallite_subFp0(1:3,1:3,c,i,e))
        invSubFi0 = math_inv33(crystallite_subFi0(1:3,1:3,c,i,e))

        if (sum(abs(dLidS)) < tol_math_check) then
          dFidS = 0.0_pReal
        else
          lhs_3333 = 0.0_pReal; rhs_3333 = 0.0_pReal
          do o=1,3; do p=1,3
            lhs_3333(1:3,1:3,o,p) = lhs_3333(1:3,1:3,o,p) &
                                  + crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidFi(1:3,1:3,o,p))
            lhs_3333(1:3,o,1:3,p) = lhs_3333(1:3,o,1:3,p) &
                                  + invFi*invFi(p,o)
            rhs_3333(1:3,1:3,o,p) = rhs_3333(1:3,1:3,o,p) &
                                  - crystallite_subdt(c,i,e)*matmul(invSubFi0,dLidS(1:3,1:3,o,p))
          enddo; enddo
          call math_invert(temp_99,error,math_3333to99(lhs_3333))
          if (error) then
            call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                            ext_msg='inversion error in analytic tangent calculation')
            dFidS = 0.0_pReal
          else
            dFidS = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
          endif
          dLidS = math_mul3333xx3333(dLidFi,dFidS) + dLidS
        endif

        call constitutive_plastic_LpAndItsTangents(devNull,dLpdS,dLpdFi, &
                                           crystallite_S (1:3,1:3,c,i,e), &
                                           constitutive_mech_Fi(pp)%data(1:3,1:3,m),c,i,e)
        dLpdS = math_mul3333xx3333(dLpdFi,dFidS) + dLpdS

!--------------------------------------------------------------------------------------------------
! calculate dSdF
        temp_33_1 = transpose(matmul(invFp,invFi))
        temp_33_2 = matmul(crystallite_subF(1:3,1:3,c,i,e),invSubFp0)
        temp_33_3 = matmul(matmul(crystallite_subF(1:3,1:3,c,i,e),invFp), invSubFi0)

        do o=1,3; do p=1,3
          rhs_3333(p,o,1:3,1:3)  = matmul(dSdFe(p,o,1:3,1:3),temp_33_1)
          temp_3333(1:3,1:3,p,o) = matmul(matmul(temp_33_2,dLpdS(1:3,1:3,p,o)), invFi) &
                                 + matmul(temp_33_3,dLidS(1:3,1:3,p,o))
        enddo; enddo
        lhs_3333 = crystallite_subdt(c,i,e)*math_mul3333xx3333(dSdFe,temp_3333) &
                 + math_mul3333xx3333(dSdFi,dFidS)

        call math_invert(temp_99,error,math_eye(9)+math_3333to99(lhs_3333))
        if (error) then
          call IO_warning(warning_ID=600,el=e,ip=i,g=c, &
                          ext_msg='inversion error in analytic tangent calculation')
          dSdF = rhs_3333
        else
          dSdF = math_mul3333xx3333(math_99to3333(temp_99),rhs_3333)
        endif

!--------------------------------------------------------------------------------------------------
! calculate dFpinvdF
        temp_3333 = math_mul3333xx3333(dLpdS,dSdF)
        do o=1,3; do p=1,3
          dFpinvdF(1:3,1:3,p,o) = -crystallite_subdt(c,i,e) &
                                * matmul(invSubFp0, matmul(temp_3333(1:3,1:3,p,o),invFi))
        enddo; enddo

!--------------------------------------------------------------------------------------------------
! assemble dPdF
        temp_33_1 = matmul(crystallite_S(1:3,1:3,c,i,e),transpose(invFp))
        temp_33_2 = matmul(invFp,temp_33_1)
        temp_33_3 = matmul(crystallite_subF(1:3,1:3,c,i,e),invFp)
        temp_33_4 = matmul(temp_33_3,crystallite_S(1:3,1:3,c,i,e))

        dPdF = 0.0_pReal
        do p=1,3
          dPdF(p,1:3,p,1:3) = transpose(temp_33_2)
        enddo
        do o=1,3; do p=1,3
          dPdF(1:3,1:3,p,o) = dPdF(1:3,1:3,p,o) &
                            + matmul(matmul(crystallite_subF(1:3,1:3,c,i,e), &
                                     dFpinvdF(1:3,1:3,p,o)),temp_33_1) &
                            + matmul(matmul(temp_33_3,dSdF(1:3,1:3,p,o)), &
                                     transpose(invFp)) &
                            + matmul(temp_33_4,transpose(dFpinvdF(1:3,1:3,p,o)))
        enddo; enddo

end function crystallite_stressTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations

  integer &
    c, &                                                                                            !< counter in integration point component loop
    i, &                                                                                            !< counter in integration point loop
    e                                                                                               !< counter in element loop

  !$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1),FEsolving_execIP(2)
      do c = 1,homogenization_Nconstituents(material_homogenizationAt(e))
        call crystallite_orientation(c,i,e)%fromMatrix(transpose(math_rotationalPart(crystallite_Fe(1:3,1:3,c,i,e))))
  enddo; enddo; enddo
  !$OMP END PARALLEL DO

  nonlocalPresent: if (any(plasticState%nonlocal)) then
    !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)
      if (plasticState(material_phaseAt(1,e))%nonlocal) then
        do i = FEsolving_execIP(1),FEsolving_execIP(2)
          call plastic_nonlocal_updateCompatibility(crystallite_orientation, &
                                                    phase_plasticityInstance(material_phaseAt(1,e)),i,e)
        enddo
      endif
    enddo
    !$OMP END PARALLEL DO
  endif nonlocalPresent

end subroutine crystallite_orientations


!--------------------------------------------------------------------------------------------------
!> @brief Map 2nd order tensor to reference config
!--------------------------------------------------------------------------------------------------
function crystallite_push33ToRef(ipc,ip,el, tensor33)

  real(pReal), dimension(3,3) :: crystallite_push33ToRef
  real(pReal), dimension(3,3), intent(in) :: tensor33
  real(pReal), dimension(3,3)             :: T
  integer, intent(in):: &
    el, &
    ip, &
    ipc

  T = matmul(material_orientation0(ipc,ip,el)%asMatrix(), &                                         ! ToDo: initial orientation correct?
             transpose(math_inv33(crystallite_subF(1:3,1:3,ipc,ip,el))))
  crystallite_push33ToRef = matmul(transpose(T),matmul(tensor33,T))

end function crystallite_push33ToRef


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state with adaptive 1st order explicit Euler method
!> using Fixed Point Iteration to adapt the stepsize
!--------------------------------------------------------------------------------------------------
subroutine integrateSourceState(g,i,e)

  integer, intent(in) :: &
    e, &                                                                                            !< element index in element loop
    i, &                                                                                            !< integration point index in ip loop
    g                                                                                               !< grain index in grain loop
  integer :: &
    NiterationState, &                                                                              !< number of iterations in state loop
    p, &
    c, &
    s, &
    size_pl
  integer, dimension(maxval(phase_Nsources)) :: &
    size_so
  real(pReal) :: &
    zeta
  real(pReal), dimension(max(constitutive_plasticity_maxSizeDotState,constitutive_source_maxSizeDotState)) :: &
    r                                                                                               ! state residuum
  real(pReal), dimension(constitutive_source_maxSizeDotState,2,maxval(phase_Nsources)) :: source_dotState
  logical :: &
    broken

  p = material_phaseAt(g,e)
  c = material_phaseMemberAt(g,i,e)

  broken = constitutive_thermal_collectDotState(crystallite_S(1:3,1:3,g,i,e), g,i,e,p,c)
  broken = broken .or. constitutive_damage_collectDotState(crystallite_S(1:3,1:3,g,i,e), g,i,e,p,c)
  if(broken) return

  do s = 1, phase_Nsources(p)
    size_so(s) = sourceState(p)%p(s)%sizeDotState
    sourceState(p)%p(s)%state(1:size_so(s),c) = sourceState(p)%p(s)%subState0(1:size_so(s),c) &
                                              + sourceState(p)%p(s)%dotState (1:size_so(s),c) &
                                              * crystallite_subdt(g,i,e)
    source_dotState(1:size_so(s),2,s) = 0.0_pReal
  enddo

  iteration: do NiterationState = 1, num%nState

    do s = 1, phase_Nsources(p)
      if(nIterationState > 1) source_dotState(1:size_so(s),2,s) = source_dotState(1:size_so(s),1,s)
      source_dotState(1:size_so(s),1,s) = sourceState(p)%p(s)%dotState(:,c)
    enddo

    broken = constitutive_thermal_collectDotState(crystallite_S(1:3,1:3,g,i,e), g,i,e,p,c)
    broken = broken .or. constitutive_damage_collectDotState(crystallite_S(1:3,1:3,g,i,e), g,i,e,p,c)
    if(broken) exit iteration

    do s = 1, phase_Nsources(p)
      zeta = damper(sourceState(p)%p(s)%dotState(:,c), &
                    source_dotState(1:size_so(s),1,s),&
                    source_dotState(1:size_so(s),2,s))
      sourceState(p)%p(s)%dotState(:,c) = sourceState(p)%p(s)%dotState(:,c) * zeta &
                                        + source_dotState(1:size_so(s),1,s)* (1.0_pReal - zeta)
      r(1:size_so(s)) = sourceState(p)%p(s)%state    (1:size_so(s),c)  &
                      - sourceState(p)%p(s)%subState0(1:size_so(s),c)  &
                      - sourceState(p)%p(s)%dotState (1:size_so(s),c) * crystallite_subdt(g,i,e)
      sourceState(p)%p(s)%state(1:size_so(s),c) = sourceState(p)%p(s)%state(1:size_so(s),c) &
                                                - r(1:size_so(s))
      crystallite_converged(g,i,e) = &
      crystallite_converged(g,i,e) .and. converged(r(1:size_so(s)), &
                                                   sourceState(p)%p(s)%state(1:size_so(s),c), &
                                                   sourceState(p)%p(s)%atol(1:size_so(s)))
    enddo

    if(crystallite_converged(g,i,e)) then
      broken = constitutive_deltaState_source(crystallite_Fe(1:3,1:3,g,i,e),g,i,e,p,c)
      exit iteration
    endif

  enddo iteration


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

end subroutine integrateSourceState


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

  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print*, ' writing field and constitutive data required for restart to file';flush(IO_STDOUT)

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName,'a')

  call HDF5_write(fileHandle,crystallite_partitionedF,'F')
  call HDF5_write(fileHandle,crystallite_Fp,        'F_p')
  call HDF5_write(fileHandle,crystallite_Lp,        'L_p')
  call HDF5_write(fileHandle,crystallite_S,           'S')

  groupHandle = HDF5_addGroup(fileHandle,'phase')
  do i = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') i,'_omega'
    call HDF5_write(groupHandle,plasticState(i)%state,datasetName)
    write(datasetName,'(i0,a)') i,'_F_i'
    call HDF5_write(groupHandle,constitutive_mech_Fi(i)%data,datasetName)
    write(datasetName,'(i0,a)') i,'_L_i'
    call HDF5_write(groupHandle,constitutive_mech_Li(i)%data,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_addGroup(fileHandle,'homogenization')
  do i = 1, size(material_name_homogenization)
    write(datasetName,'(i0,a)') i,'_omega'
    call HDF5_write(groupHandle,homogState(i)%state,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartWrite


!--------------------------------------------------------------------------------------------------
!> @brief Read data for restart
! ToDo: Merge data into one file for MPI, move state to constitutive and homogenization, respectively
!--------------------------------------------------------------------------------------------------
subroutine crystallite_restartRead

  integer :: i
  integer(HID_T) :: fileHandle, groupHandle
  character(len=pStringLen) :: fileName, datasetName

  print'(/,a,i0,a)', ' reading restart information of increment from file'

  write(fileName,'(a,i0,a)') trim(getSolverJobName())//'_',worldrank,'.hdf5'
  fileHandle = HDF5_openFile(fileName)

  call HDF5_read(fileHandle,crystallite_F0, 'F')
  call HDF5_read(fileHandle,crystallite_Fp0,'F_p')
  call HDF5_read(fileHandle,crystallite_Lp0,'L_p')
  call HDF5_read(fileHandle,crystallite_S0, 'S')

  groupHandle = HDF5_openGroup(fileHandle,'phase')
  do i = 1,size(material_name_phase)
    write(datasetName,'(i0,a)') i,'_omega'
    call HDF5_read(groupHandle,plasticState(i)%state0,datasetName)
    write(datasetName,'(i0,a)') i,'_F_i'
    call HDF5_read(groupHandle,constitutive_mech_Fi0(i)%data,datasetName)
    write(datasetName,'(i0,a)') i,'_L_i'
    call HDF5_read(groupHandle,constitutive_mech_Li0(i)%data,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  groupHandle = HDF5_openGroup(fileHandle,'homogenization')
  do i = 1,size(material_name_homogenization)
    write(datasetName,'(i0,a)') i,'_omega'
    call HDF5_read(groupHandle,homogState(i)%state0,datasetName)
  enddo
  call HDF5_closeGroup(groupHandle)

  call HDF5_closeFile(fileHandle)

end subroutine crystallite_restartRead


!--------------------------------------------------------------------------------------------------
!> @brief Forward data after successful increment.
! ToDo: Any guessing for the current states possible?
!--------------------------------------------------------------------------------------------------
subroutine crystallite_forward

  integer :: i, j

  crystallite_F0  = crystallite_partitionedF
  crystallite_Fp0 = crystallite_Fp
  crystallite_Lp0 = crystallite_Lp
  crystallite_S0  = crystallite_S

  do i = 1, size(plasticState)
    plasticState(i)%state0 = plasticState(i)%state
    constitutive_mech_Fi0(i) = constitutive_mech_Fi(i)
    constitutive_mech_Li0(i) = constitutive_mech_Li(i)
  enddo
  do i = 1,size(material_name_homogenization)
    homogState  (i)%state0 = homogState  (i)%state
    damageState (i)%state0 = damageState (i)%state
  enddo

end subroutine crystallite_forward

end module constitutive
