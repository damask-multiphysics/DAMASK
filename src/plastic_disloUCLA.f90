!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incoprorating dislocation and twinning physics
!> @details to be done
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

 integer(pInt),                       dimension(:),           allocatable,         private :: &
   plastic_disloUCLA_totalNslip                                                                     !< total number of active slip systems for each instance


 enum, bind(c)
   enumerator :: undefined_ID, &
     rho_ID, &
     rhoDip_ID, &
     shearrate_ID, &
     accumulatedshear_ID, &
     mfp_ID, &
     resolvedstress_ID, &
     thresholdstress_ID, &
     dipoledistance_ID, &
     stressexponent_ID
 end enum

 type, private :: tParameters
   real(pReal) :: &
     aTolRho, &
     grainSize, &
     SolidSolutionStrength, &                                                       !< Strength due to elements in solid solution
     mu, &
     D0, &                                                                          !< prefactor for self-diffusion coefficient
     Qsd                                                                            !< activation energy for dislocation climb
   real(pReal),                 allocatable, dimension(:) :: &
     rho0, &                                                                        !< initial edge dislocation density per slip system for each family and instance
     rhoDip0, &                                                                     !< initial edge dipole density per slip system for each family and instance
     burgers, &                                                                     !< absolute length of burgers vector [m] for each slip system and instance
     nonSchmidCoeff, &
     minDipDistance, &
     CLambda, &                                                                     !< Adj. parameter for distance between 2 forest dislocations for each slip system and instance
     atomicVolume, &
     !* mobility law parameters
     H0kp, &                                                                        !< activation energy for glide [J] for each slip system and instance
     v0, &                                                                          !< dislocation velocity prefactor [m/s] for each family and instance
     p, &                                                                           !< p-exponent in glide velocity
     q, &                                                                           !< q-exponent in glide velocity
     B, &                                                                            !< friction coeff. B (kMC)
     kink_height, &                                                                  !< height of the kink pair
     kink_width, &                                                                   !< width of the kink pair
     omega, &                                                                        !< attempt frequency for kink pair nucleation
     tau_Peierls 
   real(pReal),                 allocatable, dimension(:,:) :: &
     interaction_SlipSlip, &                                                         !< slip resistance from slip activity
     forestProjectionEdge
   real(pReal),                 allocatable, dimension(:,:,:) :: &
     Schmid_slip, &
     Schmid_twin, &
     nonSchmid_pos, &
     nonSchmid_neg
   integer(pInt) :: &
     totalNslip                                                                                     !< total number of active slip system
   integer(pInt),               allocatable, dimension(:) :: &
     Nslip                                                                                          !< number of active slip systems for each family
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
   logical :: &
     dipoleformation
 end type                                                                                           !< container type for internal constitutive parameters

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 type, private :: tDisloUCLAState
   real(pReal), pointer, dimension(:,:) :: &
     rhoEdge, &
     rhoEdgeDip, &
     accshear_slip, &
     whole
 end type

 type, private :: tDisloUCLAdependentState
   real(pReal), allocatable,     dimension(:,:) :: &
     mfp, &
     threshold_stress
 end type tDisloUCLAdependentState

 type(tDisloUCLAState ), allocatable, dimension(:), private :: &
   state, &
   dotState

 type(tDisloUCLAdependentState), allocatable, dimension(:), private :: &
   dependentState

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
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_mul3x3, &
   math_expand
 use IO, only: &
   IO_error, &
   IO_timeStamp
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_DISLOUCLA_label, &
   PLASTICITY_DISLOUCLA_ID, &
   material_phase, &
   plasticState, &
   material_allocatePlasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice

 implicit none
 integer(pInt) :: maxNinstance,&
                  f,j,k,o, i, &
                  outputSize, &
                  offset_slip, index_myFamily, index_otherFamily, &
                  startIndex, endIndex, p, &
                  sizeState, sizeDotState, &
                  NofMyPhase
 character(len=65536) :: &
   structure = ''
 character(len=65536), dimension(:), allocatable :: outputs
 integer(kind(undefined_ID))  :: outputID
 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_DISLOUCLA_label//' init  -+>>>'
 write(6,'(/,a)')   ' Cereceda et al., International Journal of Plasticity 78, 2016, 242-256'
 write(6,'(/,a)')   ' http://dx.doi.org/10.1016/j.ijplas.2015.09.002'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == PLASTICITY_DISLOUCLA_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_disloUCLA_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(plastic_disloUCLA_output(maxval(phase_Noutput),maxNinstance))
          plastic_disloUCLA_output = ''


 allocate(plastic_disloUCLA_totalNslip(maxNinstance),                          source=0_pInt)

 allocate(param(maxNinstance))
 allocate(state(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(dependentState(maxNinstance))


 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_DISLOUCLA_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             dst => dependentState(phase_plasticityInstance(p)))

   structure          = config_phase(p)%getString('lattice_structure')
   prm%mu = lattice_mu(p)

   prm%aTolRho = config_phase(p)%getFloat('atol_rho')
!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%Nslip      = config_phase(p)%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%Nslip)
   slipActive: if (prm%totalNslip > 0_pInt) then
     prm%Schmid_slip          = lattice_SchmidMatrix_slip(prm%Nslip,structure(1:3),&
                                                          config_phase(p)%getFloat('c/a',defaultVal=0.0_pReal))
     if(structure=='bcc') then
       prm%nonSchmidCoeff     = config_phase(p)%getFloats('nonschmid_coefficients',&
                                                          defaultVal = emptyRealArray)
       prm%nonSchmid_pos      = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1_pInt)
       prm%nonSchmid_neg      = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1_pInt)
     else
       prm%nonSchmid_pos      = prm%Schmid_slip
       prm%nonSchmid_neg      = prm%Schmid_slip
     endif
     prm%interaction_SlipSlip = lattice_interaction_SlipSlip(prm%Nslip, &
                                                             config_phase(p)%getFloats('interaction_slipslip'), &
                                                             structure(1:3))
     prm%rho0        = config_phase(p)%getFloats('rhoedge0',       requiredShape=shape(prm%Nslip)) 
     prm%rhoDip0     = config_phase(p)%getFloats('rhoedgedip0',    requiredShape=shape(prm%Nslip)) 
     prm%burgers     = config_phase(p)%getFloats('slipburgers',    requiredShape=shape(prm%Nslip)) 
     prm%H0kp        = config_phase(p)%getFloats('qedge',          requiredShape=shape(prm%Nslip)) 
     prm%v0          = config_phase(p)%getFloats('v0',             requiredShape=shape(prm%Nslip)) 
     prm%clambda     = config_phase(p)%getFloats('clambdaslip',    requiredShape=shape(prm%Nslip)) 
     prm%tau_Peierls = config_phase(p)%getFloats('tau_peierls',    requiredShape=shape(prm%Nslip))
     prm%p           = config_phase(p)%getFloats('p_slip',         requiredShape=shape(prm%Nslip), &
                                                 defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%q           = config_phase(p)%getFloats('q_slip',         requiredShape=shape(prm%Nslip), &
                                                 defaultVal=[(1.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%kink_height = config_phase(p)%getFloats('kink_height',    requiredShape=shape(prm%Nslip))
     prm%kink_width  = config_phase(p)%getFloats('kink_width',     requiredShape=shape(prm%Nslip))
     prm%omega       = config_phase(p)%getFloats('omega',          requiredShape=shape(prm%Nslip))
     prm%B           = config_phase(p)%getFloats('friction_coeff', requiredShape=shape(prm%Nslip))

     prm%SolidSolutionStrength  = config_phase(p)%getFloat('solidsolutionstrength')
     prm%grainSize              = config_phase(p)%getFloat('grainsize')
     prm%D0                     = config_phase(p)%getFloat('d0')
     prm%Qsd                    = config_phase(p)%getFloat('qsd')
     prm%atomicVolume           = config_phase(p)%getFloat('catomicvolume')       * prm%burgers**3.0_pReal
     prm%minDipDistance         = config_phase(p)%getFloat('cedgedipmindistance') * prm%burgers
     prm%dipoleformation        = config_phase(p)%getFloat('dipoleformationfactor') > 0.0_pReal !should be on by default


     ! expand: family => system
     prm%rho0   = math_expand(prm%rho0,  prm%Nslip)
     prm%rhoDip0   = math_expand(prm%rhoDip0,  prm%Nslip)
     prm%q   = math_expand(prm%q,  prm%Nslip)
     prm%p   = math_expand(prm%p,  prm%Nslip)
     prm%H0kp   = math_expand(prm%H0kp,  prm%Nslip)
     prm%burgers   = math_expand(prm%burgers,  prm%Nslip)
     prm%kink_height   = math_expand(prm%kink_height,  prm%Nslip)
     prm%kink_width   = math_expand(prm%kink_width,  prm%Nslip)
     prm%omega   = math_expand(prm%omega,  prm%Nslip)
     prm%tau_Peierls   = math_expand(prm%tau_Peierls,  prm%Nslip)
     prm%v0   = math_expand(prm%v0,  prm%Nslip)
     prm%B   = math_expand(prm%B,  prm%Nslip)
     prm%clambda   = math_expand(prm%clambda,  prm%Nslip)
     prm%atomicVolume   = math_expand(prm%atomicVolume,  prm%Nslip)
     prm%minDipDistance   = math_expand(prm%minDipDistance,  prm%Nslip)

     plastic_disloUCLA_totalNslip(phase_plasticityInstance(p)) = prm%totalNslip
      !if (plastic_disloUCLA_CAtomicVolume(instance) <= 0.0_pReal) &
      !  call IO_error(211_pInt,el=instance,ext_msg='cAtomicVolume ('//PLASTICITY_DISLOUCLA_label//')')
     ! if (prm%D0 <= 0.0_pReal) &
     !   call IO_error(211_pInt,el=instance,ext_msg='D0 ('//PLASTICITY_DISLOUCLA_label//')')
     ! if (plastic_disloUCLA_Qsd(instance) <= 0.0_pReal) &
     !   call IO_error(211_pInt,el=instance,ext_msg='Qsd ('//PLASTICITY_DISLOUCLA_label//')')
     ! if (plastic_disloUCLA_aTolRho(instance) <= 0.0_pReal) &
     !   call IO_error(211_pInt,el=instance,ext_msg='aTolRho ('//PLASTICITY_DISLOUCLA_label//')')
     !if (plastic_disloUCLA_rhoEdge0(f,instance) < 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='rhoEdge0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_rhoEdgeDip0(f,instance) < 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='rhoEdgeDip0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_burgersPerSlipFamily(f,instance) <= 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='slipBurgers ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_v0PerSlipFamily(f,instance) <= 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='v0 ('//PLASTICITY_DISLOUCLA_label//')')
          !if (plastic_disloUCLA_tau_peierlsPerSlipFamily(f,instance) < 0.0_pReal) &
          !  call IO_error(211_pInt,el=instance,ext_msg='tau_peierls ('//PLASTICITY_DISLOUCLA_label//')')


   else slipActive
     allocate(prm%rho0(0))
     allocate(prm%rhoDip0(0))
   endif slipActive


#if defined(__GFORTRAN__)
   outputs = ['GfortranBug86277']
   outputs = config_phase(p)%getStrings('(output)',defaultVal=outputs)
   if (outputs(1) == 'GfortranBug86277') outputs = emptyStringArray
#else
   outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
#endif
   allocate(prm%outputID(0))

   do i = 1_pInt, size(outputs)
     outputID = undefined_ID
     outputSize = prm%totalNslip
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
       case ('resolved_stress','resolved_stress_slip')
         outputID = merge(resolvedstress_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('threshold_stress','threshold_stress_slip')
         outputID = merge(thresholdstress_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('edge_dipole_distance')
         outputID = merge(dipoleDistance_ID,undefined_ID,prm%totalNslip>0_pInt)
       case ('stress_exponent')
         outputID = merge(stressexponent_ID,undefined_ID,prm%totalNslip>0_pInt)
     end select

     if (outputID /= undefined_ID) then
       plastic_disloUCLA_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_disloUCLA_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
       prm%outputID = [prm%outputID, outputID]
     endif

   enddo

   NofMyPhase=count(material_phase==p)

!--------------------------------------------------------------------------------------------------
! allocate state arrays

   sizeDotState = int(size(['rhoEdge     ','rhoEdgeDip  ','accshearslip']),pInt) * prm%totalNslip
   sizeState    = sizeDotState

   call material_allocatePlasticState(p,NofMyPhase,sizeState,sizeDotState,0_pInt, &
                                      prm%totalNslip,0_pInt,0_pInt)

   plasticState(p)%sizePostResults = sum(plastic_disloUCLA_sizePostResult(:,phase_plasticityInstance(p)))

   allocate(prm%forestProjectionEdge(prm%totalNslip,prm%totalNslip),source = 0.0_pReal)
   
   i = 0_pInt
   mySlipFamilies: do f = 1_pInt,size(prm%Nslip,1)
     index_myFamily = sum(prm%Nslip(1:f-1_pInt))

     slipSystemsLoop: do j = 1_pInt,prm%Nslip(f)
       i = i + 1_pInt
       do o = 1_pInt, size(prm%Nslip,1)
         index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
         do k = 1_pInt,prm%Nslip(o)                                   ! loop over (active) systems in other family (slip)
           prm%forestProjectionEdge(index_myFamily+j,index_otherFamily+k) = &
             abs(math_mul3x3(lattice_sn(:,sum(lattice_NslipSystem(1:f-1,p))+j,p), &
                             lattice_st(:,sum(lattice_NslipSystem(1:o-1,p))+k,p)))
       enddo; enddo
     enddo slipSystemsLoop
   enddo mySlipFamilies 

   offset_slip = 2_pInt*plasticState(p)%nSlip
   plasticState(p)%slipRate => &
     plasticState(p)%dotState(offset_slip+1:offset_slip+plasticState(p)%nSlip,1:NofMyPhase)
   plasticState(p)%accumulatedSlip => &
     plasticState(p)%state   (offset_slip+1:offset_slip+plasticState(p)%nSlip,1:NofMyPhase)

   startIndex=1_pInt
   endIndex=prm%totalNslip
   stt%rhoEdge=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdge= spread(prm%rho0,2,NofMyPhase)
   dot%rhoEdge=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

   startIndex=endIndex+1_pInt
   endIndex=endIndex+prm%totalNslip
   stt%rhoEdgeDip=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rhoEdgeDip= spread(prm%rhoDip0,2,NofMyPhase)
   dot%rhoEdgeDip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolRho

   startIndex=endIndex+1_pInt
   endIndex=endIndex+prm%totalNslip
   stt%accshear_slip=>plasticState(p)%state(startIndex:endIndex,:)
   dot%accshear_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = 1e6_pReal

   dot%whole => plasticState(p)%dotState


   allocate(dst%mfp(prm%totalNslip,NofMyPhase),source=0.0_pReal)
   allocate(dst%threshold_stress(prm%totalNslip,NofMyPhase),source=0.0_pReal)


   plasticState(p)%state0 = plasticState(p)%state                                                 ! ToDo: this could be done centrally
   end associate
 enddo

end subroutine plastic_disloUCLA_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_dependentState(instance,of)

 implicit none
 integer(pInt), intent(in)                  :: instance, of

 integer(pInt) :: &
   i
 real(pReal), dimension(param(instance)%totalNslip) :: &
  invLambdaSlip ! 1/mean free distance between 2 forest dislocations seen by a moving dislocation

 associate(prm => param(instance), stt => state(instance),dst => dependentState(instance))

 forall (i = 1_pInt:prm%totalNslip)
   invLambdaSlip(i) = sqrt(dot_product(stt%rhoEdge(:,of)+stt%rhoEdgeDip(:,of), &
                                       prm%forestProjectionEdge(:,i))) &
                    / prm%Clambda(i)
   dst%threshold_stress(i,of) = prm%mu*prm%burgers(i) &
                              * sqrt(dot_product(stt%rhoEdge(:,of)+stt%rhoEdgeDip(:,of), &
                                                 prm%interaction_SlipSlip(i,:)))
 end forall

 dst%mfp(:,of) = prm%grainSize/(1.0_pReal+prm%grainSize*invLambdaSlip)
 end associate


end subroutine plastic_disloUCLA_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_disloUCLA_LpAndItsTangent(Lp,dLp_dMp,Mp,Temperature,instance,of)

 implicit none
 integer(pInt), intent(in)                  :: instance, of
 real(pReal), intent(in)                    :: Temperature
 real(pReal), dimension(3,3),   intent(in)    :: Mp
 real(pReal), dimension(3,3), intent(out)   :: Lp
 real(pReal), dimension(3,3,3,3), intent(out)   :: dLp_dMp

 integer(pInt) :: i,k,l,m,n

 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_slip_pos,gdot_slip_neg,tau_slip_pos,tau_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg

 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal

 call kinetics(prm,stt,dst,Mp,Temperature,of, &
               gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)
 slipSystems: do i = 1_pInt, prm%totalNslip
   Lp = Lp + (gdot_slip_pos(i)+gdot_slip_neg(i))*prm%Schmid_slip(1:3,1:3,i)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtauslip_pos(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                      + dgdot_dtauslip_neg(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_neg(m,n,i)
 enddo slipSystems
 end associate

 Lp = 0.5_pReal * Lp
 dLp_dMp = 0.5_pReal * dLp_dMp

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
 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   instance, of

 real(pReal) :: &
   VacancyDiffusion
 real(pReal), dimension(plastic_disloUCLA_totalNslip(instance)) :: &
   gdot_slip_pos, gdot_slip_neg,&
   tau_slip_pos,&
   tau_slip_neg, &
   dgdot_dtauslip_neg,dgdot_dtauslip_pos,DotRhoDipFormation, ClimbVelocity, EdgeDipDistance, &
   DotRhoEdgeDipClimb

 associate(prm => param(instance), stt => state(instance),dot => dotState(instance), dst => dependentState(instance))

 call kinetics(prm,stt,dst,Mp,Temperature,of, &
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)

 dot%whole(:,of) = 0.0_pReal
 dot%accshear_slip(:,of) = (gdot_slip_pos+gdot_slip_neg)*0.5_pReal

 VacancyDiffusion = prm%D0*exp(-prm%Qsd/(kB*Temperature))

 where(dEq0(tau_slip_pos))
   EdgeDipDistance = dst%mfp(:,of)        !ToDo MD@FR: correct? was not handled properly before
   DotRhoDipFormation = 0.0_pReal
   DotRhoEdgeDipClimb = 0.0_pReal
 else where
   EdgeDipDistance = math_clip((3.0_pReal*prm%mu*prm%burgers)/(16.0_pReal*PI*abs(tau_slip_pos)), &
                               prm%minDipDistance, &        ! lower limit
                               dst%mfp(:,of))               ! upper limit
   DotRhoDipFormation = merge(((2.0_pReal*EdgeDipDistance)/prm%burgers)* stt%rhoEdge(:,of)*abs(dot%accshear_slip(:,of)), &
                              0.0_pReal, &
                              prm%dipoleformation)
   ClimbVelocity = (3.0_pReal*prm%mu*VacancyDiffusion*prm%atomicVolume/(2.0_pReal*pi*kB*Temperature)) &
                 * (1.0_pReal/(EdgeDipDistance+prm%minDipDistance))
   DotRhoEdgeDipClimb = (4.0_pReal*ClimbVelocity*stt%rhoEdgeDip(:,of))/(EdgeDipDistance-prm%minDipDistance)
 end where

 dot%rhoEdge(:,of) = abs(dot%accshear_slip(:,of))/(prm%burgers*dst%mfp(:,of)) & ! multiplication
                   - DotRhoDipFormation &
                   - (2.0_pReal*prm%minDipDistance)/prm%burgers*stt%rhoEdge(:,of)*abs(dot%accshear_slip(:,of)) !* Spontaneous annihilation of 2 single edge dislocations

 dot%rhoEdgeDip(:,of) = DotRhoDipFormation &
                      - (2.0_pReal*prm%minDipDistance)/prm%burgers* stt%rhoEdgeDip(:,of)*abs(dot%accshear_slip(:,of)) & !* Spontaneous annihilation of a single edge dislocation with a dipole constituent
                      - DotRhoEdgeDipClimb

end associate

end subroutine plastic_disloUCLA_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_disloUCLA_postResults(Mp,Temperature,instance,of) result(postResults)
 use prec, only: &
   dEq, dNeq0
 use math, only: &
   pi, &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                 intent(in) :: &
   Temperature                                                                                      !< Mandel stress
 integer(pInt),               intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_disloUCLA_sizePostResult(:,instance))) :: &
  postResults

 integer(pInt) :: &
   o,c,i
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos, &
   gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg

 associate(prm => param(instance), stt => state(instance), dst => dependentState(instance))

 postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))

     case (rho_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdge(1_pInt:prm%totalNslip,of)
     case (rhoDip_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%rhoEdgeDip(1_pInt:prm%totalNslip,of)
     case (shearrate_ID,stressexponent_ID)
       call kinetics(prm,stt,dst,Mp,Temperature,of, &
                     gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)

       if    (prm%outputID(o) == shearrate_ID) then
         postResults(c+1:c+prm%totalNslip) = (gdot_slip_pos + gdot_slip_neg)*0.5_pReal
       elseif(prm%outputID(o) == stressexponent_ID) then
         where (dNeq0(gdot_slip_pos+gdot_slip_neg))
           postResults(c+1_pInt:c + prm%totalNslip) = (tau_slip_pos+tau_slip_neg) * 0.5_pReal &
                                             / (gdot_slip_pos+gdot_slip_neg) &
                                             * (dgdot_dtauslip_pos+dgdot_dtauslip_neg)
         else where
           postResults(c+1_pInt:c + prm%totalNslip) = 0.0_pReal
         end where
       endif
     case (accumulatedshear_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%accshear_slip(1_pInt:prm%totalNslip, of)
     case (mfp_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = dst%mfp(1_pInt:prm%totalNslip, of)
     case (resolvedstress_ID)
       do i = 1_pInt, prm%totalNslip
         postResults(c+i) =math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))
       enddo
     case (thresholdstress_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = dst%threshold_stress(1_pInt:prm%totalNslip,of)
     case (dipoleDistance_ID)
       do i = 1_pInt, prm%totalNslip
         if (dNeq0(abs(math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))))) then
           postResults(c+i) = (3.0_pReal*prm%mu*prm%burgers(i)) &
                            / (16.0_pReal*pi*abs(math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))))
         else
           postResults(c+i) = huge(1.0_pReal)
         endif
         postResults(c+i)=min(postResults(c+i),dst%mfp(i,of))
       enddo
   end select

   c = c + prm%totalNslip
 enddo outputsLoop
 end associate

end function plastic_disloUCLA_postResults


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
subroutine kinetics(prm,stt,dst,Mp,Temperature,of, &
                 gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg)
 use prec, only: &
   tol_math_check, &
   dEq, dNeq0
 use math, only: &
   pi, &
math_mul33xx33

 implicit none
 type(tParameters), intent(in) :: &
   prm
 type(tDisloUCLAState), intent(in) :: &
   stt
 type(tDisloUCLAdependentState), intent(in) :: &
   dst
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),                intent(in) :: &
   temperature                                                                                      !< temperature at integration point
 integer(pInt),              intent(in) :: &
   of

 integer(pInt) :: &
   j
 real(pReal) :: dvel_slip, vel_slip
 real(pReal), intent(out), dimension(prm%totalNslip) :: &
   gdot_slip_pos,dgdot_dtauslip_pos,tau_slip_pos,gdot_slip_neg,dgdot_dtauslip_neg,tau_slip_neg
 real(pReal), dimension(prm%totalNslip) :: &
   StressRatio, BoltzmannRatio, &
   StressRatio_p,StressRatio_pminus1, &
   DotGamma0

 gdot_slip_pos = 0.0_pReal
 gdot_slip_neg = 0.0_pReal
 dgdot_dtauslip_pos = 0.0_pReal
 dgdot_dtauslip_neg = 0.0_pReal

 do j = 1_pInt, prm%totalNslip
   tau_slip_pos(j) = math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,j))
   tau_slip_neg(j) = math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,j))
 enddo

 BoltzmannRatio = prm%H0kp/(kB*Temperature)
 DotGamma0 = stt%rhoEdge(:,of)*prm%burgers*prm%v0

 do j = 1_pInt, prm%totalNslip

   significantPositiveTau: if((abs(tau_slip_pos(j))-dst%threshold_stress(j, of)) > tol_math_check) then

     StressRatio(j) = ((abs(tau_slip_pos(j))-dst%threshold_stress(j, of)) &
                 / (prm%solidSolutionStrength+prm%tau_Peierls(j)))
     StressRatio_p(j)       = StressRatio(j)** prm%p(j)
     StressRatio_pminus1(j) = StressRatio(j)**(prm%p(j)-1.0_pReal)

     vel_slip = 2.0_pReal*prm%burgers(j) * prm%kink_height(j) * prm%omega(j)  &
              * ( dst%mfp(j,of) - prm%kink_width(j) ) &
            * (tau_slip_pos(j)  &
            * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) ) &
            / ( &
            2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
            + prm%omega(j) * prm%B(j) &
            *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
            * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
            )

     gdot_slip_pos(j) = DotGamma0(j) * sign(vel_slip,tau_slip_pos(j))

     dvel_slip = 2.0_pReal*prm%burgers(j) * prm%kink_height(j) * prm%omega(j)  &
               * ( dst%mfp(j,of) - prm%kink_width(j) ) &
          * ( &
          (exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) &
          + tau_slip_pos(j) &
          * (abs(exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)))&
          *BoltzmannRatio(j)*prm%p(j)&
          *prm%q(j)/&
          (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
          StressRatio_pminus1(j)*(1-StressRatio_p(j))**(prm%q(j)-1.0_pReal)  ) &
          ) &
          *  (2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
          +  prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
          ) &
          -  (tau_slip_pos(j) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) )  &
          *  (2.0_pReal*(prm%burgers(j)**2.0_pReal) &
          +  prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * (abs(exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)))&
          *BoltzmannRatio(j)*prm%p(j)&
          *prm%q(j)/&
          (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
          StressRatio_pminus1(j)*(1-StressRatio_p(j))**(prm%q(j)-1.0_pReal)  )&
          ) &
          )  &
          / (  &
          ( &
          2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_pos(j) &
          + prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
          )**2.0_pReal &
          )

     dgdot_dtauslip_pos(j) = DotGamma0(j) * dvel_slip

   endif significantPositiveTau


   significantNegativeTau: if((abs(tau_slip_neg(j))-dst%threshold_stress(j, of)) > tol_math_check) then

     StressRatio(j) = ((abs(tau_slip_neg(j))-dst%threshold_stress(j, of)) &
                 / (prm%solidSolutionStrength+prm%tau_Peierls(j)))
     StressRatio_p(j)       = StressRatio(j)** prm%p(j)
     StressRatio_pminus1(j) = StressRatio(j)**(prm%p(j)-1.0_pReal)

     vel_slip = 2.0_pReal*prm%burgers(j) * prm%kink_height(j) * prm%omega(j)  &
              * ( dst%mfp(j,of) - prm%kink_width(j) ) &
            * (tau_slip_neg(j)  &
            * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) ) &
            / ( &
            2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
            + prm%omega(j) * prm%B(j) &
            *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
            * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
            )

     gdot_slip_neg(j) = DotGamma0(j) * sign(vel_slip,tau_slip_neg(j))

     dvel_slip = 2.0_pReal*prm%burgers(j) * prm%kink_height(j) * prm%omega(j)  &
               * ( dst%mfp(j,of) - prm%kink_width(j) ) &
          * ( &
          (exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) &
          + tau_slip_neg(j) &
          * (abs(exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)))&
          *BoltzmannRatio(j)*prm%p(j)&
          *prm%q(j)/&
          (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
          StressRatio_pminus1(j)*(1-StressRatio_p(j))**(prm%q(j)-1.0_pReal)  ) &
          ) &
          *  (2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
          +  prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
          ) &
          -  (tau_slip_neg(j) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)) )  &
          *  (2.0_pReal*(prm%burgers(j)**2.0_pReal) &
          +  prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * (abs(exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j)))&
          *BoltzmannRatio(j)*prm%p(j)&
          *prm%q(j)/&
          (prm%solidSolutionStrength+prm%tau_Peierls(j))*&
          StressRatio_pminus1(j)*(1-StressRatio_p(j))**(prm%q(j)-1.0_pReal)  )&
          ) &
          )  &
          / (  &
          ( &
          2.0_pReal*(prm%burgers(j)**2.0_pReal)*tau_slip_neg(j) &
          + prm%omega(j) * prm%B(j) &
          *(( dst%mfp(j,of) - prm%kink_width(j) )**2.0_pReal) &
          * exp(-BoltzmannRatio(j)*(1-StressRatio_p(j)) ** prm%q(j))  &
          )**2.0_pReal &
          )


     dgdot_dtauslip_neg(j) = DotGamma0(j) * dvel_slip
   endif significantNegativeTau
 enddo

end subroutine kinetics

end module plastic_disloUCLA
