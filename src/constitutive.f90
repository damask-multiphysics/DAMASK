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
  use geometry_plastic_nonlocal, only: &
    geometry_plastic_nonlocal_disable

  implicit none
  private

  integer(kind(ELASTICITY_undefined_ID)),     dimension(:),   allocatable :: &
    phase_elasticity                                                                                !< elasticity of each phase

  integer(kind(PLASTICITY_undefined_ID)),     dimension(:),   allocatable :: &                      !ToDo: old intel compiler complains about protected
    phase_plasticity                                                                                !< plasticity of each phase

  integer(kind(SOURCE_undefined_ID)), dimension(:,:), allocatable :: &                              ! ToDo: old intel compiler complains about protected
    phase_source, &                                                                                 !< active sources mechanisms of each phase
    phase_kinematics, &                                                                             !< active kinematic mechanisms of each phase
    phase_stiffnessDegradation                                                                      !< active stiffness degradation mechanisms of each phase

  integer, dimension(:), allocatable,  public :: &                                                  ! ToDo: old intel compiler complains about protected
    phase_Nsources, &                                                                               !< number of source mechanisms active in each phase
    phase_Nkinematics, &                                                                            !< number of kinematic mechanisms active in each phase
    phase_NstiffnessDegradations, &                                                                 !< number of stiffness degradation mechanisms active in each phase
    phase_plasticityInstance, &                                                                     !< instance of particular plasticity of each phase
    phase_elasticityInstance                                                                        !< instance of particular elasticity of each phase

  logical, dimension(:), allocatable, public :: &                                                   ! ToDo: old intel compiler complains about protected
    phase_localPlasticity                                                                           !< flags phases with local constitutive law

  type(tPlasticState), allocatable, dimension(:), public :: &
    plasticState
  type(tSourceState),  allocatable, dimension(:), public :: &
    sourceState


  integer, public, protected :: &
    constitutive_plasticity_maxSizeDotState, &
    constitutive_source_maxSizeDotState

  interface

    module subroutine mech_init
    end subroutine mech_init

    module subroutine damage_init
    end subroutine damage_init

    module subroutine thermal_init
    end subroutine thermal_init


    module function plastic_active(plastic_label)  result(active_plastic)
      character(len=*), intent(in)        :: plastic_label
      logical, dimension(:), allocatable  :: active_plastic
    end function plastic_active

    module function source_active(source_label,src_length)  result(active_source)
      character(len=*), intent(in)          :: source_label
      integer,          intent(in)          :: src_length
      logical, dimension(:,:), allocatable  :: active_source
    end function source_active

    module function kinematics_active(kinematics_label,kinematics_length)  result(active_kinematics)
      character(len=*), intent(in)          :: kinematics_label
      integer,          intent(in)          :: kinematics_length
      logical, dimension(:,:), allocatable  :: active_kinematics
    end function kinematics_active

    module subroutine plastic_isotropic_dotState(Mp,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_isotropic_dotState

    module subroutine plastic_phenopowerlaw_dotState(Mp,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_phenopowerlaw_dotState

    module subroutine plastic_kinehardening_dotState(Mp,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_kinehardening_dotState

    module subroutine plastic_dislotwin_dotState(Mp,T,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_dislotwin_dotState

    module subroutine plastic_disloTungsten_dotState(Mp,T,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      real(pReal),                  intent(in) :: &
        T
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_disloTungsten_dotState

    module subroutine plastic_nonlocal_dotState(Mp, F, Fp, Temperature,timestep, &
                                                        instance,of,ip,el)
      real(pReal), dimension(3,3), intent(in) :: &
        Mp                                                                                          !< MandelStress
      real(pReal), dimension(3,3,homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems), intent(in) :: &
        F, &                                                                                        !< deformation gradient
        Fp                                                                                          !< plastic deformation gradient
      real(pReal), intent(in) :: &
        Temperature, &                                                                              !< temperature
        timestep                                                                                    !< substepped crystallite time increment
      integer, intent(in) :: &
        instance, &
        of, &
        ip, &                                                                                       !< current integration point
        el                                                                                          !< current element number
    end subroutine plastic_nonlocal_dotState


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

    pure module function kinematics_thermal_expansion_initialStrain(homog,phase,offset) result(initialStrain)
     integer, intent(in) :: &
       phase, &
       homog, &
       offset
     real(pReal), dimension(3,3) :: &
       initialStrain
    end function kinematics_thermal_expansion_initialStrain

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


    module subroutine plastic_kinehardening_deltaState(Mp,instance,of)
      real(pReal), dimension(3,3),  intent(in) :: &
        Mp                                                                                          !< Mandel stress
      integer,                      intent(in) :: &
        instance, &
        of
    end subroutine plastic_kinehardening_deltaState

    module subroutine plastic_nonlocal_deltaState(Mp,instance,of,ip,el)
      real(pReal), dimension(3,3), intent(in) :: &
        Mp
      integer, intent(in) :: &
        instance, &
        of, &
        ip, &
        el
    end subroutine plastic_nonlocal_deltaState

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

    module subroutine constitutive_SandItsTangents(S, dS_dFe, dS_dFi, Fe, Fi, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                          !< component-ID of integration point
        ip, &                                                                                           !< integration point
        el                                                                                              !< element
      real(pReal),   intent(in),  dimension(3,3) :: &
        Fe, &                                                                                           !< elastic deformation gradient
        Fi                                                                                              !< intermediate deformation gradient
      real(pReal),   intent(out), dimension(3,3) :: &
        S                                                                                               !< 2nd Piola-Kirchhoff stress tensor
      real(pReal),   intent(out), dimension(3,3,3,3) :: &
        dS_dFe, &                                                                                       !< derivative of 2nd P-K stress with respect to elastic deformation gradient
        dS_dFi                                                                                          !< derivative of 2nd P-K stress with respect to intermediate deformation gradient
    end subroutine constitutive_SandItsTangents

    module subroutine plastic_results
    end subroutine plastic_results

    module subroutine damage_results
    end subroutine damage_results

  end interface

  interface constitutive_LpAndItsTangents

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

  end interface constitutive_LpAndItsTangents

  interface constitutive_dependentState

    module subroutine constitutive_plastic_dependentState(F, Fp, ipc, ip, el)
      integer, intent(in) :: &
        ipc, &                                                                                      !< component-ID of integration point
        ip, &                                                                                       !< integration point
        el                                                                                          !< element
      real(pReal),   intent(in), dimension(3,3) :: &
        F, &                                                                                        !< elastic deformation gradient
        Fp                                                                                          !< plastic deformation gradient
    end subroutine constitutive_plastic_dependentState

  end interface constitutive_dependentState


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

  type(tDebugOptions) :: debugConstitutive

  public :: &
    constitutive_init, &
    constitutive_homogenizedC, &
    constitutive_LpAndItsTangents, &
    constitutive_dependentState, &
    constitutive_LiAndItsTangents, &
    constitutive_initialFi, &
    constitutive_SandItsTangents, &
    constitutive_collectDotState, &
    constitutive_deltaState, &
    constitutive_damage_getRateAndItsTangents, &
    constitutive_thermal_getRateAndItsTangents, &
    constitutive_results, &
    constitutive_allocateState, &
    plastic_nonlocal_updateCompatibility, &
    plastic_active, &
    source_active, &
    kinematics_active
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
module function source_active(source_label,src_length)  result(active_source)

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

module function kinematics_active(kinematics_label,kinematics_length)  result(active_kinematics)

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
!> @brief  collects initial intermediate deformation gradient
!--------------------------------------------------------------------------------------------------
pure function constitutive_initialFi(ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el                                                                                              !< element
  real(pReal), dimension(3,3) :: &
    constitutive_initialFi                                                                          !< composite initial intermediate deformation gradient
  integer :: &
    k                                                                                               !< counter in kinematics loop
  integer :: &
    phase, &
    homog, offset

  constitutive_initialFi = math_I3
  phase = material_phaseAt(ipc,el)

  KinematicsLoop: do k = 1, phase_Nkinematics(phase)                                                !< Warning: small initial strain assumption
    kinematicsType: select case (phase_kinematics(k,phase))
      case (KINEMATICS_thermal_expansion_ID) kinematicsType
        homog = material_homogenizationAt(el)
        offset = thermalMapping(homog)%p(ip,el)
        constitutive_initialFi = &
          constitutive_initialFi + kinematics_thermal_expansion_initialStrain(homog,phase,offset)
    end select kinematicsType
  enddo KinematicsLoop

end function constitutive_initialFi


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_collectDotState(S, FArray, Fi, FpArray, subdt, ipc, ip, el,phase,of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),  intent(in) :: &
    subdt                                                                                           !< timestep
  real(pReal),  intent(in), dimension(3,3,homogenization_maxNconstituents,discretization_nIPs,discretization_Nelems) :: &
    FArray, &                                                                                       !< elastic deformation gradient
    FpArray                                                                                         !< plastic deformation gradient
  real(pReal),  intent(in), dimension(3,3) :: &
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),  intent(in), dimension(3,3) :: &
    S                                                                                               !< 2nd Piola Kirchhoff stress (vector notation)
  real(pReal),              dimension(3,3) :: &
    Mp
  integer :: &
    ho, &                                                                                           !< homogenization
    tme, &                                                                                          !< thermal member position
    i, &                                                                                            !< counter in source loop
    instance
  logical :: broken

  ho = material_homogenizationAt(el)
  tme = thermalMapping(ho)%p(ip,el)
  instance = phase_plasticityInstance(phase)

  Mp = matmul(matmul(transpose(Fi),Fi),S)

  plasticityType: select case (phase_plasticity(phase))

    case (PLASTICITY_ISOTROPIC_ID) plasticityType
      call plastic_isotropic_dotState(Mp,instance,of)

    case (PLASTICITY_PHENOPOWERLAW_ID) plasticityType
      call plastic_phenopowerlaw_dotState(Mp,instance,of)

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call plastic_kinehardening_dotState(Mp,instance,of)

    case (PLASTICITY_DISLOTWIN_ID) plasticityType
      call plastic_dislotwin_dotState(Mp,temperature(ho)%p(tme),instance,of)

    case (PLASTICITY_DISLOTUNGSTEN_ID) plasticityType
      call plastic_disloTungsten_dotState(Mp,temperature(ho)%p(tme),instance,of)

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call plastic_nonlocal_dotState(Mp,FArray,FpArray,temperature(ho)%p(tme),subdt, &
                                          instance,of,ip,el)
  end select plasticityType
  broken = any(IEEE_is_NaN(plasticState(phase)%dotState(:,of)))

  SourceLoop: do i = 1, phase_Nsources(phase)

    sourceType: select case (phase_source(i,phase))

      case (SOURCE_damage_anisoBrittle_ID) sourceType
        call source_damage_anisoBrittle_dotState(S, ipc, ip, el) ! correct stress?

      case (SOURCE_damage_isoDuctile_ID) sourceType
        call source_damage_isoDuctile_dotState(ipc, ip, el)

      case (SOURCE_damage_anisoDuctile_ID) sourceType
        call source_damage_anisoDuctile_dotState(ipc, ip, el)

      case (SOURCE_thermal_externalheat_ID) sourceType
        call source_thermal_externalheat_dotState(phase,of)

    end select sourceType

    broken = broken .or. any(IEEE_is_NaN(sourceState(phase)%p(i)%dotState(:,of)))

  enddo SourceLoop

end function constitutive_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
function constitutive_deltaState(S, Fe, Fi, ipc, ip, el, phase, of) result(broken)

  integer, intent(in) :: &
    ipc, &                                                                                          !< component-ID of integration point
    ip, &                                                                                           !< integration point
    el, &                                                                                           !< element
    phase, &
    of
  real(pReal),   intent(in), dimension(3,3) :: &
    S, &                                                                                            !< 2nd Piola Kirchhoff stress
    Fe, &                                                                                           !< elastic deformation gradient
    Fi                                                                                              !< intermediate deformation gradient
  real(pReal),               dimension(3,3) :: &
    Mp
  integer :: &
    i, &
    instance, &
    myOffset, &
    mySize
  logical :: &
    broken

  Mp  = matmul(matmul(transpose(Fi),Fi),S)
  instance = phase_plasticityInstance(phase)

  plasticityType: select case (phase_plasticity(phase))

    case (PLASTICITY_KINEHARDENING_ID) plasticityType
      call plastic_kinehardening_deltaState(Mp,instance,of)
      broken = any(IEEE_is_NaN(plasticState(phase)%deltaState(:,of)))

    case (PLASTICITY_NONLOCAL_ID) plasticityType
      call plastic_nonlocal_deltaState(Mp,instance,of,ip,el)
      broken = any(IEEE_is_NaN(plasticState(phase)%deltaState(:,of)))

    case default
      broken = .false.

  end select plasticityType

  if(.not. broken) then
    select case(phase_plasticity(phase))
      case (PLASTICITY_NONLOCAL_ID,PLASTICITY_KINEHARDENING_ID)

        myOffset = plasticState(phase)%offsetDeltaState
        mySize   = plasticState(phase)%sizeDeltaState
        plasticState(phase)%state(myOffset + 1:myOffset + mySize,of) = &
        plasticState(phase)%state(myOffset + 1:myOffset + mySize,of) + plasticState(phase)%deltaState(1:mySize,of)
    end select
  endif


  sourceLoop: do i = 1, phase_Nsources(phase)

     sourceType: select case (phase_source(i,phase))

      case (SOURCE_damage_isoBrittle_ID) sourceType
        call source_damage_isoBrittle_deltaState  (constitutive_homogenizedC(ipc,ip,el), Fe, &
                                                   ipc, ip, el)
        broken = broken .or. any(IEEE_is_NaN(sourceState(phase)%p(i)%deltaState(:,of)))
        if(.not. broken) then
          myOffset = sourceState(phase)%p(i)%offsetDeltaState
          mySize   = sourceState(phase)%p(i)%sizeDeltaState
          sourceState(phase)%p(i)%state(myOffset + 1: myOffset + mySize,of) = &
          sourceState(phase)%p(i)%state(myOffset + 1: myOffset + mySize,of) + sourceState(phase)%p(i)%deltaState(1:mySize,of)
        endif

    end select sourceType

  enddo SourceLoop

end function constitutive_deltaState


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
!> @brief writes constitutive results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine constitutive_results

  call plastic_results
  call damage_results

end subroutine constitutive_results


end module constitutive
