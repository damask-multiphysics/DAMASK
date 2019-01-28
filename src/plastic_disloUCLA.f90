!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystal plasticity model for bcc metals, especially Tungsten
!--------------------------------------------------------------------------------------------------
module plastic_disloUCLA
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:,:),         allocatable, target, public :: &
   plastic_disloUCLA_sizePostResult                                                                 !< size of each post result output
 character(len=64),                   dimension(:,:),         allocatable, target, public :: &
   plastic_disloUCLA_output                                                                         !< name of each post result output

 real(pReal),                                                 parameter,           private :: &
   kB = 1.38e-23_pReal                                                                              !< Boltzmann constant in J/Kelvin

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     rho_ID, &
     rhoDip_ID, &
     shearrate_ID, &
     accumulatedshear_ID, &
     mfp_ID, &
     thresholdstress_ID
 end enum

 type, private :: tParameters
   real(pReal) :: &
     aTolRho, &
     grainSize, &
     SolidSolutionStrength, &                                                                       !< Strength due to elements in solid solution
     mu, &
     D0, &                                                                                          !< prefactor for self-diffusion coefficient
     Qsd                                                                                            !< activation energy for dislocation climb
   real(pReal),                 allocatable, dimension(:) :: &
     rho0, &                                                                                        !< initial edge dislocation density
     rhoDip0, &                                                                                     !< initial edge dipole density
     burgers, &                                                                                     !< absolute length of burgers vector [m]
     nonSchmidCoeff, &
     minDipDistance, &
     CLambda, &                                                                                     !< Adj. parameter for distance between 2 forest dislocations
     atomicVolume, &
     tau_Peierls, &
     tau0, &
     !* mobility law parameters
     H0kp, &                                                                                        !< activation energy for glide [J]
     v0, &                                                                                          !< dislocation velocity prefactor [m/s]
     p, &                                                                                           !< p-exponent in glide velocity
     q, &                                                                                           !< q-exponent in glide velocity
     B, &                                                                                           !< friction coefficient
     kink_height, &                                                                                 !< height of the kink pair
     w, &                                                                                           !< width of the kink pair
     omega                                                                                          !< attempt frequency for kink pair nucleation
   real(pReal),                 allocatable, dimension(:,:) :: &
     interaction_SlipSlip, &                                                                        !< slip resistance from slip activity
     forestProjectionEdge
   real(pReal),                 allocatable, dimension(:,:,:) :: &
     Schmid, &
     nonSchmid_pos, &
     nonSchmid_neg
   integer(pInt) :: &
     totalNslip                                                                                     !< total number of active slip system
   integer(pInt),               allocatable, dimension(:) :: &
     Nslip                                                                                          !< number of active slip systems for each family
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
   logical :: &
     dipoleFormation                                                                                !< flag indicating consideration of dipole formation
 end type                                                                                           !< container type for internal constitutive parameters

 type, private :: tDisloUCLAState
   real(pReal), pointer, dimension(:,:) :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear
 end type tDisloUCLAState

 type, private :: tDisloUCLAdependentState
   real(pReal), allocatable, dimension(:,:) :: &
     mfp, &
     dislocationSpacing, &
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
   plastic_disloUCLA_postResults
 private :: &
   kinetics

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   pStringLen
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_expand
 use IO, only: &
   IO_error, &
   IO_timeStamp
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
   MATERIAL_partPhase, &
   config_phase
 use lattice

 implicit none
 integer(pInt) :: &
   Ninstance, &
   p, i, &
   NipcMyPhase, &
   sizeState, sizeDotState, &
   startIndex, endIndex

 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_DISLOUCLA_label//' init  -+>>>'
 write(6,'(/,a)')   ' Cereceda et al., International Journal of Plasticity 78, 2016, 242-256'
 write(6,'(/,a)')   ' http://dx.doi.org/10.1016/j.ijplas.2015.09.002'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 Ninstance = int(count(phase_plasticity == PLASTICITY_DISLOUCLA_ID),pInt)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

 allocate(plastic_disloUCLA_sizePostResult(maxval(phase_Noutput),Ninstance),source=0_pInt)
 allocate(plastic_disloUCLA_output(maxval(phase_Noutput),Ninstance))
          plastic_disloUCLA_output = ''

 allocate(param(Ninstance))
 allocate(state(Ninstance))
 allocate(dotState(Ninstance))
 allocate(dependentState(Ninstance))


 do p = 1_pInt, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_DISLOUCLA_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             dst => dependentState(phase_plasticityInstance(p)), &
             config => config_phase(p))

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
   prm%mu = lattice_mu(p)

   prm%aTolRho = config%getFloat('atol_rho')

   ! sanity checks
   if (prm%aTolRho <= 0.0_pReal) extmsg = trim(extmsg)//' atol_rho'

!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%Nslip      = config%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%Nslip)
   slipActive: if (prm%totalNslip > 0_pInt) then
     prm%Schmid = lattice_SchmidMatrix_slip(prm%Nslip,config%getString('lattice_structure'),&
                                            config%getFloat('c/a',defaultVal=0.0_pReal))

     if(trim(config%getString('lattice_structure')) == 'bcc') then
       prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                                 defaultVal = emptyRealArray)
       prm%nonSchmid_pos  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1_pInt)
       prm%nonSchmid_neg  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1_pInt)
     else
       prm%nonSchmid_pos  = prm%Schmid
       prm%nonSchmid_neg  = prm%Schmid
     endif

     prm%interaction_SlipSlip = lattice_interaction_SlipSlip(prm%Nslip, &
                                                             config%getFloats('interaction_slipslip'), &
                                                             config%getString('lattice_structure'))
     prm%forestProjectionEdge = lattice_forestProjection(prm%Nslip,config%getString('lattice_structure'),&
                                                         config%getFloat('c/a',defaultVal=0.0_pReal))

     prm%rho0        = config%getFloats('rhoedge0',       requiredSize=size(prm%Nslip))
     prm%rhoDip0     = config%getFloats('rhoedgedip0',    requiredSize=size(prm%Nslip))
     prm%v0          = config%getFloats('v0',             requiredSize=size(prm%Nslip))
     prm%burgers     = config%getFloats('slipburgers',    requiredSize=size(prm%Nslip))
     prm%H0kp        = config%getFloats('qedge',          requiredSize=size(prm%Nslip))

     prm%clambda     = config%getFloats('clambdaslip',    requiredSize=size(prm%Nslip))
     prm%tau_Peierls = config%getFloats('tau_peierls',    requiredSize=size(prm%Nslip))  ! ToDo: Deprecated
     prm%p           = config%getFloats('p_slip',         requiredSize=size(prm%Nslip), &
                                        defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%q           = config%getFloats('q_slip',         requiredSize=size(prm%Nslip), &
                                        defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%kink_height = config%getFloats('kink_height',    requiredSize=size(prm%Nslip))
     prm%w           = config%getFloats('kink_width',     requiredSize=size(prm%Nslip))
     prm%omega       = config%getFloats('omega',          requiredSize=size(prm%Nslip))
     prm%B           = config%getFloats('friction_coeff', requiredSize=size(prm%Nslip))

     prm%SolidSolutionStrength  = config%getFloat('solidsolutionstrength')                 ! ToDo: Deprecated
     prm%grainSize              = config%getFloat('grainsize')
     prm%D0                     = config%getFloat('d0')
     prm%Qsd                    = config%getFloat('qsd')
     prm%atomicVolume           = config%getFloat('catomicvolume')       * prm%burgers**3.0_pReal
     prm%minDipDistance         = config%getFloat('cedgedipmindistance') * prm%burgers
     prm%dipoleformation        = config%getFloat('dipoleformationfactor') > 0.0_pReal !should be on by default, ToDo: change to /key/-type key

     ! expand: family => system
     prm%rho0           = math_expand(prm%rho0,           prm%Nslip)
     prm%rhoDip0        = math_expand(prm%rhoDip0,        prm%Nslip)
     prm%q              = math_expand(prm%q,              prm%Nslip)
     prm%p              = math_expand(prm%p,              prm%Nslip)
     prm%H0kp           = math_expand(prm%H0kp,           prm%Nslip)
     prm%burgers        = math_expand(prm%burgers,        prm%Nslip)
     prm%kink_height    = math_expand(prm%kink_height,    prm%Nslip)
     prm%w              = math_expand(prm%w,              prm%Nslip)
     prm%omega          = math_expand(prm%omega,          prm%Nslip)
     prm%tau_Peierls    = math_expand(prm%tau_Peierls,    prm%Nslip)
     prm%v0             = math_expand(prm%v0,             prm%Nslip)
     prm%B              = math_expand(prm%B,              prm%Nslip)
     prm%clambda        = math_expand(prm%clambda,        prm%Nslip)
     prm%atomicVolume   = math_expand(prm%atomicVolume,   prm%Nslip)
     prm%minDipDistance = math_expand(prm%minDipDistance, prm%Nslip)

     prm%tau0 = prm%tau_peierls + prm%SolidSolutionStrength

     ! sanity checks
     if (    prm%D0             <= 0.0_pReal)  extmsg = trim(extmsg)//' d0'
     if (    prm%Qsd            <= 0.0_pReal)  extmsg = trim(extmsg)//' qsd'
     if (any(prm%rho0           <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedge0'
     if (any(prm%rhoDip0        <  0.0_pReal)) extmsg = trim(extmsg)//' rhoedgedip0'
     if (any(prm%v0             <  0.0_pReal)) extmsg = trim(extmsg)//' v0'
     if (any(prm%burgers        <= 0.0_pReal)) extmsg = trim(extmsg)//' slipburgers'
     if (any(prm%H0kp           <= 0.0_pReal)) extmsg = trim(extmsg)//' qedge'
     if (any(prm%tau_peierls    <  0.0_pReal)) extmsg = trim(extmsg)//' tau_peierls'
     if (any(prm%minDipDistance <= 0.0_pReal)) extmsg = trim(extmsg)//' cedgedipmindistance or slipburgers'
     if (any(prm%atomicVolume   <= 0.0_pReal)) extmsg = trim(extmsg)//' catomicvolume or slipburgers'

   else slipActive
     allocate(prm%rho0(0))
     allocate(prm%rhoDip0(0))
   endif slipActive

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211_pInt,ext_msg=trim(extmsg)//'('//PLASTICITY_DISLOUCLA_label//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(trim(outputs(i)))

       case ('edge_density')
         outputID  = merge(rho_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('dipole_density')
         outputID = merge(rhoDip_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('shear_rate','shearrate','shear_rate_slip','shearrate_slip')
         outputID = merge(shearrate_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('accumulated_shear','accumulatedshear','accumulated_shear_slip')
         outputID = merge(accumulatedshear_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('mfp','mfp_slip')
         outputID = merge(mfp_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('threshold_stress','threshold_stress_slip')
         outputID = merge(thresholdstress_ID,undefined_ID,prm%totalNslip>0_pInt)

     end select

     if (outputID /= undefined_ID) then
       plastic_disloUCLA_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_disloUCLA_sizePostResult(i,phase_plasticityInstance(p)) = prm%totalNslip
       prm%outputID = [prm%outputID, outputID]
     endif

   enddo

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phase == p)
   sizeDotState = int(size(['rhoEdge     ','rhoEdgeDip  ','accshearslip']),pInt) * prm%totalNslip
   sizeState = sizeDotState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0_pInt, &
                                      prm%totalNslip,0_pInt,0_pInt)
   plasticState(p)%sizePostResults = sum(plastic_disloUCLA_sizePostResult(:,phase_plasticityInstance(p)))

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1_pInt
   endIndex   = prm%totalNslip
   stt%rhoEdge=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdge= spread(prm%rho0,2,NipcMyPhase)
   dot%rhoEdge=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + prm%totalNslip
   stt%rhoEdgeDip=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdgeDip= spread(prm%rhoDip0,2,NipcMyPhase)
   dot%rhoEdgeDip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + prm%totalNslip
   stt%accshear=>plasticState(p)%state(startIndex:endIndex,:)
   dot%accshear=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1.0e6_pReal  !ToDo: better make optional parameter
   ! global alias
   plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%accumulatedSlip => plasticState(p)%state(startIndex:endIndex,:)

   allocate(dst%mfp(prm%totalNslip,NipcMyPhase),               source=0.0_pReal)
   allocate(dst%dislocationSpacing(prm%totalNslip,NipcMyPhase),source=0.0_pReal)
   allocate(dst%threshold_stress(prm%totalNslip,NipcMyPhase),  source=0.0_pReal)

   plasticState(p)%state0 = plasticState(p)%state                                                   ! ToDo: this could be done centrally

   end associate

 enddo

end subroutine plastic_disloUCLA_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
pure subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp,Mp,Temperature,instance,of)

 implicit none
 real(pReal), dimension(3,3),     intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to the Mandel stress

 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature
 integer(pInt),               intent(in) :: &
   instance, &
   of

 integer(pInt) :: &
   i,k,l,m,n
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg, &
   dgdot_dtau_pos,dgdot_dtau_neg

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal

 associate(prm => param(instance))

 call kinetics(Mp,Temperature,instance,of,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)
 do i = 1_pInt, prm%totalNslip
   Lp = Lp + (gdot_pos(i)+gdot_neg(i))*prm%Schmid(1:3,1:3,i)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_pos(i) * prm%Schmid(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                      + dgdot_dtau_neg(i) * prm%Schmid(k,l,i) * prm%nonSchmid_neg(m,n,i)
 enddo

 end associate

end subroutine plastic_disloUCLA_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dotState(Mp,Temperature,instance,of)
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
   temperature                                                                                      !< temperature
 integer(pInt),                intent(in) :: &
   instance, &
   of

 real(pReal) :: &
   VacancyDiffusion
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos, gdot_neg,&
   tau_pos,&
   tau_neg, &
   DotRhoDipFormation, ClimbVelocity, EdgeDipDistance, &
   DotRhoEdgeDipClimb

 associate(prm => param(instance), stt => state(instance),dot => dotState(instance), dst => dependentState(instance))

 call kinetics(Mp,Temperature,instance,of,&
               gdot_pos,gdot_neg, &
               tau_pos1 = tau_pos,tau_neg1 = tau_neg)

 dot%accshear(:,of) = (gdot_pos+gdot_neg)                                                      ! ToDo: needs to be abs
 VacancyDiffusion = prm%D0*exp(-prm%Qsd/(kB*Temperature))

 where(dEq0(tau_pos))                                                                          ! ToDo: use avg of pos and neg
   DotRhoDipFormation = 0.0_pReal
   DotRhoEdgeDipClimb = 0.0_pReal
 else where
   EdgeDipDistance = math_clip((3.0_pReal*prm%mu*prm%burgers)/(16.0_pReal*PI*abs(tau_pos)), &
                               prm%minDipDistance, &                                                ! lower limit
                               dst%mfp(:,of))                                                       ! upper limit
   DotRhoDipFormation = merge(((2.0_pReal*EdgeDipDistance)/prm%burgers)* stt%rhoEdge(:,of)*abs(dot%accshear(:,of)), & ! ToDo: ignore region of spontaneous annihilation
                              0.0_pReal, &
                              prm%dipoleformation)
   ClimbVelocity = (3.0_pReal*prm%mu*VacancyDiffusion*prm%atomicVolume/(2.0_pReal*pi*kB*Temperature)) &
                 * (1.0_pReal/(EdgeDipDistance+prm%minDipDistance))
   DotRhoEdgeDipClimb = (4.0_pReal*ClimbVelocity*stt%rhoEdgeDip(:,of))/(EdgeDipDistance-prm%minDipDistance) ! ToDo: Discuss with Franz: Stress dependency?
 end where

 dot%rhoEdge(:,of) = abs(dot%accshear(:,of))/(prm%burgers*dst%mfp(:,of)) &                     ! multiplication
                   - DotRhoDipFormation &
                   - (2.0_pReal*prm%minDipDistance)/prm%burgers*stt%rhoEdge(:,of)*abs(dot%accshear(:,of)) !* Spontaneous annihilation of 2 single edge dislocations
 dot%rhoEdgeDip(:,of) = DotRhoDipFormation &
                      - (2.0_pReal*prm%minDipDistance)/prm%burgers*stt%rhoEdgeDip(:,of)*abs(dot%accshear(:,of)) & !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
                      - DotRhoEdgeDipClimb

 end associate

end subroutine plastic_disloUCLA_dotState


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dependentState(instance,of)

 implicit none
 integer(pInt),                intent(in) :: &
   instance, &
   of

 integer(pInt) :: &
   i

 associate(prm => param(instance), stt => state(instance),dst => dependentState(instance))

 forall (i = 1_pInt:prm%totalNslip)
   dst%dislocationSpacing(i,of) = sqrt(dot_product(stt%rhoEdge(:,of)+stt%rhoEdgeDip(:,of), &
                                       prm%forestProjectionEdge(:,i)))
   dst%threshold_stress(i,of) = prm%mu*prm%burgers(i) &
                              * sqrt(dot_product(stt%rhoEdge(:,of)+stt%rhoEdgeDip(:,of), &
                                                 prm%interaction_SlipSlip(i,:)))
 end forall

 dst%mfp(:,of) = prm%grainSize/(1.0_pReal+prm%grainSize*dst%dislocationSpacing(:,of)/prm%Clambda)
 dst%dislocationSpacing(:,of) = dst%mfp(:,of)                                                       ! ToDo: Hack to recover wrong behavior for the moment

 end associate

end subroutine plastic_disloUCLA_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_disloUCLA_postResults(Mp,Temperature,instance,of) result(postResults)
 use prec, only: &
   dEq, dNeq0
 use math, only: &
   PI, &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                 intent(in) :: &
   temperature                                                                                      !< temperature
 integer(pInt),               intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_disloUCLA_sizePostResult(:,instance))) :: &
  postResults

 integer(pInt) :: &
   o,c,i
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg

 c = 0_pInt

 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))

     case (rho_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdge(1_pInt:prm%totalNslip,of)
     case (rhoDip_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdgeDip(1_pInt:prm%totalNslip,of)
     case (shearrate_ID)
       call kinetics(Mp,Temperature,instance,of,gdot_pos,gdot_neg)
       postResults(c+1:c+prm%totalNslip) = gdot_pos + gdot_neg
     case (accumulatedshear_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%accshear(1_pInt:prm%totalNslip, of)
     case (mfp_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = dst%mfp(1_pInt:prm%totalNslip, of)
     case (thresholdstress_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = dst%threshold_stress(1_pInt:prm%totalNslip,of)

   end select

   c = c + prm%totalNslip

 enddo outputsLoop

 end associate

end function plastic_disloUCLA_postResults


!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on slip systems, their derivatives with respect to resolved stress and the
!  resolved stresss
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,Temperature,instance,of, &
                 gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg,tau_pos1,tau_neg1)
 use prec, only: &
   tol_math_check, &
   dEq, dNeq0
 use math, only: &
   PI, &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   temperature                                                                                      !< temperature
 integer(pInt),                intent(in) :: &
   instance, &
   of

 real(pReal),                  intent(out), dimension(param(instance)%totalNslip) :: &
   gdot_pos, &
   gdot_neg
 real(pReal),                  intent(out), optional, dimension(param(instance)%totalNslip) :: &
   dgdot_dtau_pos, &
   dgdot_dtau_neg, &
   tau_pos1, &
   tau_neg1
 real(pReal), dimension(param(instance)%totalNslip) :: &
   StressRatio, &
   StressRatio_p,StressRatio_pminus1, &
   dvel, vel, &
   tau_pos,tau_neg, &
   needsGoodName                                                                                    ! ToDo: @Karo: any idea?
 integer(pInt) :: j

 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 do j = 1_pInt, prm%totalNslip
   tau_pos(j) = math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,j))
   tau_neg(j) = math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,j))
 enddo


 if (present(tau_pos1)) tau_pos1 = tau_pos
 if (present(tau_neg1)) tau_neg1 = tau_neg

 associate(BoltzmannRatio  => prm%H0kp/(kB*Temperature), &
           DotGamma0       => stt%rhoEdge(:,of)*prm%burgers*prm%v0, &
           effectiveLength => dst%mfp(:,of) - prm%w)

 significantPositiveTau: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
   StressRatio = (abs(tau_pos)-dst%threshold_stress(:,of))/prm%tau0
   StressRatio_p       = StressRatio** prm%p
   StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
   needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)

   vel = 2.0_pReal*prm%burgers * prm%kink_height * prm%omega  &
            * effectiveLength * tau_pos * needsGoodName  &
          / (  2.0_pReal*(prm%burgers**2.0_pReal)*tau_pos &
             + prm%omega * prm%B * effectiveLength**2.0_pReal* needsGoodName  &
            )

   gdot_pos = DotGamma0 * sign(vel,tau_pos) * 0.5_pReal
 else where significantPositiveTau
   gdot_pos = 0.0_pReal
 end where significantPositiveTau

 if (present(dgdot_dtau_pos)) then
 significantPositiveTau2: where(abs(tau_pos)-dst%threshold_stress(:,of) > tol_math_check)
   dvel = 2.0_pReal*prm%burgers * prm%kink_height * prm%omega* effectiveLength &
             * ( &
                 (needsGoodName + tau_pos * abs(needsGoodName)*BoltzmannRatio*prm%p &
                                * prm%q/prm%tau0 &
                                * StressRatio_pminus1*(1-StressRatio_p)**(prm%q-1.0_pReal)   &
                 ) &
               * (  2.0_pReal*(prm%burgers**2.0_pReal)*tau_pos &
                  +  prm%omega * prm%B* effectiveLength **2.0_pReal* needsGoodName &
                 ) &
               - tau_pos * needsGoodName * (2.0_pReal*prm%burgers**2.0_pReal &
                                                 +  prm%omega * prm%B *effectiveLength **2.0_pReal&
                                                  * (abs(needsGoodName)*BoltzmannRatio*prm%p *prm%q/prm%tau0 &
                                                *StressRatio_pminus1*(1-StressRatio_p)**(prm%q-1.0_pReal)  )&
                                               ) &
              )  &
            /(2.0_pReal*prm%burgers**2.0_pReal*tau_pos &
              + prm%omega * prm%B* effectiveLength**2.0_pReal* needsGoodName )**2.0_pReal

   dgdot_dtau_pos = DotGamma0 * dvel* 0.5_pReal
 else where significantPositiveTau2
   dgdot_dtau_pos = 0.0_pReal
 end where significantPositiveTau2
 endif

 significantNegativeTau: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
   StressRatio = (abs(tau_neg)-dst%threshold_stress(:,of))/prm%tau0
   StressRatio_p       = StressRatio** prm%p
   StressRatio_pminus1 = StressRatio**(prm%p-1.0_pReal)
   needsGoodName       = exp(-BoltzmannRatio*(1-StressRatio_p) ** prm%q)

   vel = 2.0_pReal*prm%burgers * prm%kink_height * prm%omega  &
            * effectiveLength * tau_neg * needsGoodName  &
          / (  2.0_pReal*(prm%burgers**2.0_pReal)*tau_neg &
             + prm%omega * prm%B * effectiveLength**2.0_pReal* needsGoodName  &
            )

   gdot_neg = DotGamma0 * sign(vel,tau_neg) * 0.5_pReal
 else where significantNegativeTau
   gdot_neg = 0.0_pReal
 end where significantNegativeTau

 if (present(dgdot_dtau_neg)) then
 significantNegativeTau2: where(abs(tau_neg)-dst%threshold_stress(:,of) > tol_math_check)
   dvel = 2.0_pReal*prm%burgers * prm%kink_height * prm%omega* effectiveLength &
             * ( &
                 (needsGoodName + tau_neg * abs(needsGoodName)*BoltzmannRatio*prm%p &
                                * prm%q/prm%tau0 &
                                * StressRatio_pminus1*(1-StressRatio_p)**(prm%q-1.0_pReal)   &
                 ) &
               * (  2.0_pReal*(prm%burgers**2.0_pReal)*tau_neg &
                  +  prm%omega * prm%B* effectiveLength **2.0_pReal* needsGoodName &
                 ) &
               - tau_neg * needsGoodName * (2.0_pReal*prm%burgers**2.0_pReal &
                                                 +  prm%omega * prm%B *effectiveLength **2.0_pReal&
                                                  * (abs(needsGoodName)*BoltzmannRatio*prm%p *prm%q/prm%tau0 &
                                                *StressRatio_pminus1*(1-StressRatio_p)**(prm%q-1.0_pReal)  )&
                                               ) &
              )  &
            /(2.0_pReal*prm%burgers**2.0_pReal*tau_neg &
              + prm%omega * prm%B* effectiveLength**2.0_pReal* needsGoodName )**2.0_pReal

   dgdot_dtau_neg = DotGamma0 * dvel * 0.5_pReal
 else where significantNegativeTau2
   dgdot_dtau_neg = 0.0_pReal
 end where significantNegativeTau2
 end if
 end associate
 end associate

end subroutine kinetics

end module plastic_disloUCLA
