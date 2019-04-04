!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystal plasticity model for bcc metals, especially Tungsten
!--------------------------------------------------------------------------------------------------
module plastic_disloUCLA
  use prec, only: &
    pReal
 
  implicit none
  private
  integer,                       dimension(:,:),   allocatable, target, public :: &
    plastic_disloUCLA_sizePostResult                                                                !< size of each post result output
  character(len=64),             dimension(:,:),   allocatable, target, public :: &
    plastic_disloUCLA_output                                                                        !< name of each post result output
 
  real(pReal),                                                 parameter,           private :: &
    kB = 1.38e-23_pReal                                                                             !< Boltzmann constant in J/Kelvin
 
  enum, bind(c)
    enumerator :: &
      undefined_ID, &
      rho_mob_ID, &
      rho_dip_ID, &
      dot_gamma_sl_ID, &
      gamma_sl_ID, &
      Lambda_sl_ID, &
      thresholdstress_ID
  end enum
 
  type, private :: tParameters
    real(pReal) :: &
      aTol_rho, &
      D, &                                                                                          !< grain size
      mu, &
      D_0, &                                                                                        !< prefactor for self-diffusion coefficient
      Q_cl                                                                                          !< activation energy for dislocation climb
    real(pReal),                 dimension(:),  allocatable :: &
      rho_mob_0, &                                                                                  !< initial dislocation density
      rho_dip_0, &                                                                                  !< initial dipole density
      b_sl, &                                                                                       !< magnitude of burgers vector [m]
      nonSchmidCoeff, &
      D_a, &
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations
      atomicVolume, &
      tau_0, &
      !* mobility law parameters
      delta_F, &                                                                                    !< activation energy for glide [J]
      v0, &                                                                                         !< dislocation velocity prefactor [m/s]
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      B, &                                                                                          !< friction coefficient
      kink_height, &                                                                                !< height of the kink pair
      w, &                                                                                          !< width of the kink pair
      omega                                                                                         !< attempt frequency for kink pair nucleation
    real(pReal),                  dimension(:,:),  allocatable :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      forestProjectionEdge
    real(pReal),                  dimension(:,:,:),  allocatable :: &
      Schmid, &
      nonSchmid_pos, &
      nonSchmid_neg
    integer :: &
      sum_N_sl                                                                                      !< total number of active slip system
    integer,                      dimension(:), allocatable :: &
      N_sl                                                                                          !< number of active slip systems for each family
    integer(kind(undefined_ID)),  dimension(:),allocatable :: &
      outputID                                                                                      !< ID of each post result output
    logical :: &
      dipoleFormation                                                                               !< flag indicating consideration of dipole formation
  end type                                                                                          !< container type for internal constitutive parameters
 
  type, private :: tDisloUCLAState
    real(pReal), dimension(:,:), pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tDisloUCLAState
 
  type, private :: tDisloUCLAdependentState
    real(pReal), dimension(:,:), allocatable :: &
      Lambda_sl, &
      threshold_stress
  end type tDisloUCLAdependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),              allocatable, dimension(:), private :: param
  type(tDisloUCLAState),          allocatable, dimension(:), private :: &
    dotState, &
    state
  type(tDisloUCLAdependentState), allocatable, dimension(:), private :: dependentState
 
  public :: &
    plastic_disloUCLA_init, &
    plastic_disloUCLA_dependentState, &
    plastic_disloUCLA_LpAndItsTangent, &
    plastic_disloUCLA_dotState, &
    plastic_disloUCLA_postResults, &
    plastic_disloUCLA_results
  private :: &
    kinetics

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_init()
  use prec, only: &
    pStringLen
  use debug, only: &
    debug_level,&
    debug_constitutive,&
    debug_levelBasic
  use math, only: &
    math_expand
  use IO, only: &
    IO_error
  use material, only: &
    phase_plasticity, &
    phase_plasticityInstance, &
    phase_Noutput, &
    material_allocatePlasticState, &
    PLASTICITY_DISLOUCLA_label, &
    PLASTICITY_DISLOUCLA_ID, &
    material_phase, &
    plasticState
  use config, only: &
    config_phase
  use lattice
 
  implicit none
  integer :: &
    Ninstance, &
    p, i, &
    NipcMyPhase, &
    sizeState, sizeDotState, &
    startIndex, endIndex
 
  integer,              dimension(0), parameter :: emptyIntArray    = [integer::]
  real(pReal),          dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
  character(len=65536), dimension(0), parameter :: emptyStringArray = [character(len=65536)::]
 
  integer(kind(undefined_ID)) :: &
    outputID
 
  character(len=pStringLen) :: &
    extmsg = ''
  character(len=65536), dimension(:), allocatable :: &
    outputs
 
  write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_DISLOUCLA_label//' init  -+>>>'
 
  write(6,'(/,a)')   ' Cereceda et al., International Journal of Plasticity 78:242–256, 2016'
  write(6,'(a)')     ' https://dx.doi.org/10.1016/j.ijplas.2015.09.002'
 
  Ninstance = count(phase_plasticity == PLASTICITY_DISLOUCLA_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance
 
  allocate(plastic_disloUCLA_sizePostResult(maxval(phase_Noutput),Ninstance),source=0)
  allocate(plastic_disloUCLA_output(maxval(phase_Noutput),Ninstance))
           plastic_disloUCLA_output = ''
 
  allocate(param(Ninstance))
  allocate(state(Ninstance))
  allocate(dotState(Ninstance))
  allocate(dependentState(Ninstance))
 
 
  do p = 1, size(phase_plasticity)
    if (phase_plasticity(p) /= PLASTICITY_DISLOUCLA_ID) cycle
    associate(prm => param(phase_plasticityInstance(p)), &
              dot => dotState(phase_plasticityInstance(p)), &
              stt => state(phase_plasticityInstance(p)), &
              dst => dependentState(phase_plasticityInstance(p)), &
              config => config_phase(p))

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
    prm%mu = lattice_mu(p)
 
    prm%aTol_rho = config%getFloat('atol_rho')
 
    ! sanity checks
    if (prm%aTol_rho <= 0.0_pReal) extmsg = trim(extmsg)//' atol_rho'

!--------------------------------------------------------------------------------------------------
! slip related parameters
    prm%N_sl      = config%getInts('nslip',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(prm%N_sl)
    slipActive: if (prm%sum_N_sl > 0) then
      prm%Schmid = lattice_SchmidMatrix_slip(prm%N_sl,config%getString('lattice_structure'),&
                                             config%getFloat('c/a',defaultVal=0.0_pReal))
 
      if(trim(config%getString('lattice_structure')) == 'bcc') then
        prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                                  defaultVal = emptyRealArray)
        prm%nonSchmid_pos  = lattice_nonSchmidMatrix(prm%N_sl,prm%nonSchmidCoeff,+1)
        prm%nonSchmid_neg  = lattice_nonSchmidMatrix(prm%N_sl,prm%nonSchmidCoeff,-1)
      else
        prm%nonSchmid_pos  = prm%Schmid
        prm%nonSchmid_neg  = prm%Schmid
      endif
 
      prm%h_sl_sl     = transpose(lattice_interaction_SlipBySlip(prm%N_sl, &
                                                       config%getFloats('interaction_slipslip'), &
                                                       config%getString('lattice_structure')))
      prm%forestProjectionEdge = lattice_forestProjection(prm%N_sl,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))
 
      prm%rho_mob_0   = config%getFloats('rhoedge0',       requiredSize=size(prm%N_sl))
      prm%rho_dip_0   = config%getFloats('rhoedgedip0',    requiredSize=size(prm%N_sl))
      prm%v0          = config%getFloats('v0',             requiredSize=size(prm%N_sl))
      prm%b_sl        = config%getFloats('slipburgers',    requiredSize=size(prm%N_sl))
      prm%delta_F     = config%getFloats('qedge',          requiredSize=size(prm%N_sl))
 
      prm%i_sl        = config%getFloats('clambdaslip',    requiredSize=size(prm%N_sl))
      prm%tau_0       = config%getFloats('tau_peierls',    requiredSize=size(prm%N_sl))
      prm%p           = config%getFloats('p_slip',         requiredSize=size(prm%N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(prm%N_sl))])
      prm%q           = config%getFloats('q_slip',         requiredSize=size(prm%N_sl), &
                                         defaultVal=[(1.0_pReal,i=1,size(prm%N_sl))])
      prm%kink_height = config%getFloats('kink_height',    requiredSize=size(prm%N_sl))
      prm%w           = config%getFloats('kink_width',     requiredSize=size(prm%N_sl))
      prm%omega       = config%getFloats('omega',          requiredSize=size(prm%N_sl))
      prm%B           = config%getFloats('friction_coeff', requiredSize=size(prm%N_sl))
 
      prm%D                   = config%getFloat('grainsize')
      prm%D_0                 = config%getFloat('d0')
      prm%Q_cl                = config%getFloat('qsd')
      prm%atomicVolume        = config%getFloat('catomicvolume')       * prm%b_sl**3.0_pReal
      prm%D_a                 = config%getFloat('cedgedipmindistance') * prm%b_sl
      prm%dipoleformation     = config%getFloat('dipoleformationfactor') > 0.0_pReal !should be on by default, ToDo: change to /key/-type key
 
      ! expand: family => system
      prm%rho_mob_0      = math_expand(prm%rho_mob_0,      prm%N_sl)
      prm%rho_dip_0      = math_expand(prm%rho_dip_0,      prm%N_sl)
      prm%q              = math_expand(prm%q,              prm%N_sl)
      prm%p              = math_expand(prm%p,              prm%N_sl)
      prm%delta_F        = math_expand(prm%delta_F,        prm%N_sl)
      prm%b_sl           = math_expand(prm%b_sl,           prm%N_sl)
      prm%kink_height    = math_expand(prm%kink_height,    prm%N_sl)
      prm%w              = math_expand(prm%w,              prm%N_sl)
      prm%omega          = math_expand(prm%omega,          prm%N_sl)
      prm%tau_0          = math_expand(prm%tau_0,          prm%N_sl)
      prm%v0             = math_expand(prm%v0,             prm%N_sl)
      prm%B              = math_expand(prm%B,              prm%N_sl)
      prm%i_sl           = math_expand(prm%i_sl,           prm%N_sl)
      prm%atomicVolume   = math_expand(prm%atomicVolume,   prm%N_sl)
      prm%D_a            = math_expand(prm%D_a,            prm%N_sl)
 
 
      ! sanity checks
      if (    prm%D_0            <= 0.0_pReal)  extmsg = trim(extmsg)//' D_0'
      if (    prm%Q_cl           <= 0.0_pReal)  extmsg = trim(extmsg)//' Q_cl'
      if (any(prm%rho_mob_0      <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedge0'
      if (any(prm%rho_dip_0      <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedgedip0'
      if (any(prm%v0             <  0.0_pReal)) extmsg = trim(extmsg)//' v0'
      if (any(prm%b_sl           <= 0.0_pReal)) extmsg = trim(extmsg)//' slipb_sl'
      if (any(prm%delta_F        <= 0.0_pReal)) extmsg = trim(extmsg)//' qedge'
      if (any(prm%tau_0          <  0.0_pReal)) extmsg = trim(extmsg)//' tau_0'
      if (any(prm%D_a            <= 0.0_pReal)) extmsg = trim(extmsg)//' cedgedipmindistance or slipb_sl'
      if (any(prm%atomicVolume   <= 0.0_pReal)) extmsg = trim(extmsg)//' catomicvolume or slipb_sl'
 
    else slipActive
      allocate(prm%rho_mob_0(0))
      allocate(prm%rho_dip_0(0))
    endif slipActive

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') &
      call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_DISLOUCLA_label//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
    outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
    allocate(prm%outputID(0))
    do i=1, size(outputs)
      outputID = undefined_ID
      select case(trim(outputs(i)))
 
        case ('edge_density')
          outputID  = merge(rho_mob_ID,undefined_ID,prm%sum_N_sl>0)
        case ('dipole_density')
          outputID = merge(rho_dip_ID,undefined_ID,prm%sum_N_sl>0)
        case ('shear_rate','shearrate','shear_rate_slip','shearrate_slip')
          outputID = merge(dot_gamma_sl_ID,undefined_ID,prm%sum_N_sl>0)
        case ('accumulated_shear','accumulatedshear','accumulated_shear_slip')
          outputID = merge(gamma_sl_ID,undefined_ID,prm%sum_N_sl>0)
        case ('mfp','mfp_slip')
          outputID = merge(Lambda_sl_ID,undefined_ID,prm%sum_N_sl>0)
        case ('threshold_stress','threshold_stress_slip')
          outputID = merge(thresholdstress_ID,undefined_ID,prm%sum_N_sl>0)
 
      end select
 
      if (outputID /= undefined_ID) then
        plastic_disloUCLA_output(i,phase_plasticityInstance(p)) = outputs(i)
        plastic_disloUCLA_sizePostResult(i,phase_plasticityInstance(p)) = prm%sum_N_sl
        prm%outputID = [prm%outputID, outputID]
      endif
 
    enddo

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    NipcMyPhase = count(material_phase == p)
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl
    sizeState = sizeDotState
 
    call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0, &
                                       prm%sum_N_sl,0,0)
    plasticState(p)%sizePostResults = sum(plastic_disloUCLA_sizePostResult(:,phase_plasticityInstance(p)))

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
    startIndex = 1
    endIndex   = prm%sum_N_sl
    stt%rho_mob=>plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_mob= spread(prm%rho_mob_0,2,NipcMyPhase)
    dot%rho_mob=>plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho
 
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%rho_dip=>plasticState(p)%state(startIndex:endIndex,:)
    stt%rho_dip= spread(prm%rho_dip_0,2,NipcMyPhase)
    dot%rho_dip=>plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%aTolState(startIndex:endIndex) = prm%aTol_rho
 
    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    stt%gamma_sl=>plasticState(p)%state(startIndex:endIndex,:)
    dot%gamma_sl=>plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal                                    ! Don't use for convergence check
    ! global alias
    plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)
    plasticState(p)%accumulatedSlip => plasticState(p)%state(startIndex:endIndex,:)
 
    allocate(dst%Lambda_sl(prm%sum_N_sl,NipcMyPhase),         source=0.0_pReal)
    allocate(dst%threshold_stress(prm%sum_N_sl,NipcMyPhase),  source=0.0_pReal)
 
    plasticState(p)%state0 = plasticState(p)%state                                                  ! ToDo: this could be done centrally
 
    end associate
 
  enddo

end subroutine plastic_disloUCLA_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
pure subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp, &
                                                  Mp,T,instance,of)
  implicit none
  real(pReal), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pReal), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress
 
  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                 intent(in) :: &
    T                                                                                               !< temperature
  integer,                     intent(in) :: &
    instance, &
    of
 
  integer :: &
    i,k,l,m,n
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    dot_gamma_pos,dot_gamma_neg, &
    ddot_gamma_dtau_pos,ddot_gamma_dtau_neg
 
  Lp = 0.0_pReal
  dLp_dMp = 0.0_pReal
 
  associate(prm => param(instance))
 
  call kinetics(Mp,T,instance,of,dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg)
  do i = 1, prm%sum_N_sl
    Lp = Lp + (dot_gamma_pos(i)+dot_gamma_neg(i))*prm%Schmid(1:3,1:3,i)
    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                       + ddot_gamma_dtau_pos(i) * prm%Schmid(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                       + ddot_gamma_dtau_neg(i) * prm%Schmid(k,l,i) * prm%nonSchmid_neg(m,n,i)
  enddo
 
  end associate

end subroutine plastic_disloUCLA_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dotState(Mp,T,instance,of)
  use prec, only: &
    tol_math_check, &
    dEq0
  use math, only: &
    PI, &
    math_clip
 
  implicit none
  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                               !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                                !< temperature
  integer,                      intent(in) :: &
    instance, &
    of
 
  real(pReal) :: &
    VacancyDiffusion
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos, gdot_neg,&
    tau_pos,&
    tau_neg, &
    v_cl, &
    dot_rho_dip_formation, &
    dot_rho_dip_climb, &
    dip_distance
 
  associate(prm => param(instance), stt => state(instance),dot => dotState(instance), dst => dependentState(instance))
 
  call kinetics(Mp,T,instance,of,&
                gdot_pos,gdot_neg, &
                tau_pos_out = tau_pos,tau_neg_out = tau_neg)
 
  dot%gamma_sl(:,of) = (gdot_pos+gdot_neg)                                                          ! ToDo: needs to be abs
  VacancyDiffusion = prm%D_0*exp(-prm%Q_cl/(kB*T))
 
  where(dEq0(tau_pos))                                                                              ! ToDo: use avg of pos and neg
    dot_rho_dip_formation = 0.0_pReal
    dot_rho_dip_climb     = 0.0_pReal
  else where
    dip_distance = math_clip(3.0_pReal*prm%mu*prm%b_sl/(16.0_pReal*PI*abs(tau_pos)), &
                             prm%D_a, &                                                          ! lower limit
                             dst%Lambda_sl(:,of))                                                ! upper limit
    dot_rho_dip_formation = merge(2.0_pReal*dip_distance* stt%rho_mob(:,of)*abs(dot%gamma_sl(:,of))/prm%b_sl, & ! ToDo: ignore region of spontaneous annihilation
                                  0.0_pReal, &
                                  prm%dipoleformation)
    v_cl = (3.0_pReal*prm%mu*VacancyDiffusion*prm%atomicVolume/(2.0_pReal*pi*kB*T)) &
                  * (1.0_pReal/(dip_distance+prm%D_a))
    dot_rho_dip_climb = (4.0_pReal*v_cl*stt%rho_dip(:,of))/(dip_distance-prm%D_a)               ! ToDo: Discuss with Franz: Stress dependency?
  end where
 
  dot%rho_mob(:,of) = abs(dot%gamma_sl(:,of))/(prm%b_sl*dst%Lambda_sl(:,of)) &                      ! multiplication
                    - dot_rho_dip_formation &
                    - (2.0_pReal*prm%D_a)/prm%b_sl*stt%rho_mob(:,of)*abs(dot%gamma_sl(:,of))        ! Spontaneous annihilation of 2 single edge dislocations
  dot%rho_dip(:,of) = dot_rho_dip_formation &
                    - (2.0_pReal*prm%D_a)/prm%b_sl*stt%rho_dip(:,of)*abs(dot%gamma_sl(:,of)) &      ! Spontaneous annihilation of a single edge dislocation with a dipole constituent
                    - dot_rho_dip_climb
 
  end associate

end subroutine plastic_disloUCLA_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dependentState(instance,of)
 
  implicit none
  integer,      intent(in) :: &
    instance, &
    of
 
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    dislocationSpacing
  integer :: &
    i
 
  associate(prm => param(instance), stt => state(instance),dst => dependentState(instance))
 
  forall (i = 1:prm%sum_N_sl) &
    dislocationSpacing(i) = sqrt(dot_product(stt%rho_mob(:,of)+stt%rho_dip(:,of), &
                                        prm%forestProjectionEdge(:,i)))
  dst%threshold_stress(:,of) = prm%mu*prm%b_sl &
                             * sqrt(matmul(prm%h_sl_sl,stt%rho_mob(:,of)+stt%rho_dip(:,of)))
 
  dst%Lambda_sl(:,of) = prm%D/(1.0_pReal+prm%D*dislocationSpacing/prm%i_sl)
 
  end associate

end subroutine plastic_disloUCLA_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_disloUCLA_postResults(Mp,T,instance,of) result(postResults)
  use prec, only: &
    dEq, dNeq0
  use math, only: &
    PI, &
    math_mul33xx33
 
  implicit none
  real(pReal), dimension(3,3), intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                 intent(in) :: &
    T                                                                                               !< temperature
  integer,                     intent(in) :: &
    instance, &
    of
 
  real(pReal), dimension(sum(plastic_disloUCLA_sizePostResult(:,instance))) :: &
   postResults
 
  integer :: &
    o,c
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    gdot_pos,gdot_neg
 
  c = 0
 
  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
 
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
 
      case (rho_mob_ID)
        postResults(c+1:c+prm%sum_N_sl) = stt%rho_mob(1:prm%sum_N_sl,of)
      case (rho_dip_ID)
        postResults(c+1:c+prm%sum_N_sl) = stt%rho_dip(1:prm%sum_N_sl,of)
      case (dot_gamma_sl_ID)
        call kinetics(Mp,T,instance,of,gdot_pos,gdot_neg)
        postResults(c+1:c+prm%sum_N_sl) = gdot_pos + gdot_neg
      case (gamma_sl_ID)
        postResults(c+1:c+prm%sum_N_sl) = stt%gamma_sl(1:prm%sum_N_sl, of)
      case (Lambda_sl_ID)
        postResults(c+1:c+prm%sum_N_sl) = dst%Lambda_sl(1:prm%sum_N_sl, of)
      case (thresholdstress_ID)
        postResults(c+1:c+prm%sum_N_sl) = dst%threshold_stress(1:prm%sum_N_sl,of)
 
    end select
 
    c = c + prm%sum_N_sl
 
  enddo outputsLoop

  end associate

end function plastic_disloUCLA_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_results(instance,group)
#if defined(PETSc) || defined(DAMASK_HDF5)
  use results

  implicit none
  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
  
  integer :: o

  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))
      case (rho_mob_ID)
        call results_writeDataset(group,stt%rho_mob,'rho_mob',&
                                 'mobile dislocation density','1/m^2')
      case (rho_dip_ID)
        call results_writeDataset(group,stt%rho_dip,'rho_dip',&
                                  'dislocation dipole density''1/m^2')
      case (dot_gamma_sl_ID)
        call results_writeDataset(group,stt%gamma_sl,'dot_gamma_sl',&
                                  'plastic slip','1')
      case (Lambda_sl_ID)
        call results_writeDataset(group,dst%Lambda_sl,'Lambda_sl',&
                                  'mean free path for slip','m')
      case (thresholdstress_ID)
        call results_writeDataset(group,dst%threshold_stress,'threshold_stress',&
                                  'threshold stress for slip','Pa')
    end select
  enddo outputsLoop
  end associate
  
#else
  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
#endif

end subroutine plastic_disloUCLA_results


!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on slip systems, their derivatives with respect to resolved stress and the
!  resolved stresss
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,T,instance,of, &
                 dot_gamma_pos,dot_gamma_neg,ddot_gamma_dtau_pos,ddot_gamma_dtau_neg,tau_pos_out,tau_neg_out)
  use prec, only: &
    tol_math_check, &
    dEq, dNeq0
  use math, only: &
    PI, &
    math_mul33xx33
 
  implicit none
  real(pReal), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pReal),                  intent(in) :: &
    T                                                                                               !< temperature
  integer,                      intent(in) :: &
    instance, &
    of
 
  real(pReal),                  intent(out), dimension(param(instance)%sum_N_sl) :: &
    dot_gamma_pos, &
    dot_gamma_neg
  real(pReal),                  intent(out), optional, dimension(param(instance)%sum_N_sl) :: &
    ddot_gamma_dtau_pos, &
    ddot_gamma_dtau_neg, &
    tau_pos_out, &
    tau_neg_out
  real(pReal), dimension(param(instance)%sum_N_sl) :: &
    StressRatio, &
    StressRatio_p,StressRatio_pminus1, &
    dvel, vel, &
    tau_pos,tau_neg, &
    t_n, t_k, dtk,dtn, &
    needsGoodName                                                                                   ! ToDo: @Karo: any idea?
  integer :: j
 
  associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))
 
  do j = 1, prm%sum_N_sl
    tau_pos(j) = math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,j))
    tau_neg(j) = math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,j))
  enddo
 
 
  if (present(tau_pos_out)) tau_pos_out = tau_pos
  if (present(tau_neg_out)) tau_neg_out = tau_neg
 
  associate(BoltzmannRatio  => prm%delta_F/(kB*T), &
            dot_gamma_0     => stt%rho_mob(:,of)*prm%b_sl*prm%v0, &
            effectiveLength => dst%Lambda_sl(:,of) - prm%w)
 
  significantPositiveTau: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
    StressRatio = (abs(tau_pos)-dst%threshold_stress(:,of))/prm%tau_0
    StressRatio_p       = StressRatio** prm%p
    StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
    needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)
    
    t_n = prm%b_sl/(needsGoodName*prm%omega*effectiveLength)
    t_k = effectiveLength * prm%B /(2.0_pReal*prm%b_sl*tau_pos)     ! our definition of tk is different with the one in dislotwin
 
    vel = prm%kink_height/(t_n + t_k)
 
    dot_gamma_pos = dot_gamma_0 * sign(vel,tau_pos) * 0.5_pReal
  else where significantPositiveTau
    dot_gamma_pos = 0.0_pReal
  end where significantPositiveTau
 
  if (present(ddot_gamma_dtau_pos)) then
  significantPositiveTau2: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
    dtn = -1.0_pReal * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pReal-StressRatio_p)**(prm%q - 1.0_pReal) &
        * (StressRatio)**(prm%p - 1.0_pReal) / prm%tau_0
    dtk = -1.0_pReal * t_k / tau_pos
   
    dvel = -1.0_pReal * prm%kink_height * (dtk + dtn) / (t_n + t_k)**2.0_pReal
 
    ddot_gamma_dtau_pos = dot_gamma_0 * dvel* 0.5_pReal
  else where significantPositiveTau2
    ddot_gamma_dtau_pos = 0.0_pReal
  end where significantPositiveTau2
  endif
 
  significantNegativeTau: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
    StressRatio = (abs(tau_neg)-dst%threshold_stress(:,of))/prm%tau_0
    StressRatio_p       = StressRatio** prm%p
    StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
    needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)
 
    t_n = prm%b_sl/(needsGoodName*prm%omega*effectiveLength)
    t_k = effectiveLength * prm%B /(2.0_pReal*prm%b_sl*tau_pos)     ! our definition of tk is different with the one in dislotwin
 
    vel = prm%kink_height/(t_n + t_k)
 
    dot_gamma_neg = dot_gamma_0 * sign(vel,tau_neg) * 0.5_pReal
  else where significantNegativeTau
    dot_gamma_neg = 0.0_pReal
  end where significantNegativeTau
 
  if (present(ddot_gamma_dtau_neg)) then
  significantNegativeTau2: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
    dtn = -1.0_pReal * t_n * BoltzmannRatio * prm%p * prm%q * (1.0_pReal-StressRatio_p)**(prm%q - 1.0_pReal) &
        * (StressRatio)**(prm%p - 1.0_pReal) / prm%tau_0
    dtk = -1.0_pReal * t_k / tau_neg
   
    dvel = -1.0_pReal * prm%kink_height * (dtk + dtn) / (t_n + t_k)**2.0_pReal
 
    ddot_gamma_dtau_neg = dot_gamma_0 * dvel * 0.5_pReal
  else where significantNegativeTau2
    ddot_gamma_dtau_neg = 0.0_pReal
  end where significantNegativeTau2
  end if
  
  end associate
  end associate

end subroutine kinetics

end module plastic_disloUCLA
