!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Michigan State University
!> @author Zhuowen Zhao, Michigan State University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  Phenomenological crystal plasticity using a power law formulation for the shear rates
!!         and a Voce-type kinematic hardening rule
!--------------------------------------------------------------------------------------------------
module plastic_kinehardening
 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_sizePostResult                                                             !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_kinehardening_output                                                                     !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_kinehardening_Noutput                                                                    !< number of outputs per instance
 
  

 integer(pInt),                       dimension(:,:),   allocatable,         private :: &
   plastic_kinehardening_Nslip                                                                      !< active number of slip systems per family (input parameter, per family)
   

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     crss_ID, &                                                                         !< critical resolved stress
     crss_back_ID, &                                                                    !< critical resolved back stress
     sense_ID, &                                                                        !< sense of acting shear stress (-1 or +1)
     chi0_ID, &                                                                         !< backstress at last switch of stress sense (positive?)
     gamma0_ID, &                                                                       !< accumulated shear at last switch of stress sense (at current switch?)
     accshear_ID, &
     shearrate_ID, &
     resolvedstress_ID
 end enum


 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   real(pReal) :: &
     gdot0, &                                                                                       !< reference shear strain rate for slip (input parameter)
     n_slip, &                                                                                      !< stress exponent for slip (input parameter)
     aTolResistance, &
     aTolShear
   real(pReal),                         dimension(:),   allocatable,          private :: &
     crss0, &                                                                                        !< initial critical shear stress for slip (input parameter, per family)
     theta0, &                                                                                      !< initial hardening rate of forward stress for each slip
     theta1, &                                                                                      !< asymptotic hardening rate of forward stress for each slip >
     theta0_b, &                                                                                    !< initial hardening rate of back stress for each slip > 
     theta1_b, &                                                                                    !< asymptotic hardening rate of back stress for each slip >
     tau1, &
     tau1_b, &
     nonSchmidCoeff
   real(pReal),                         dimension(:,:),   allocatable,          private :: &
     interaction_slipslip                                                                        !< latent hardening matrix
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
 end type

 type, private :: tKinehardeningState
   real(pReal), pointer, dimension(:,:) :: &                                                        !< vectors along NipcMyInstance
     crss, &                                                                                        !< critical resolved stress
     crss_back, &                                                                                   !< critical resolved back stress
     sense, &                                                                                       !< sense of acting shear stress (-1 or +1)
     chi0, &                                                                                        !< backstress at last switch of stress sense
     gamma0, &                                                                                      !< accumulated shear at last switch of stress sense
     accshear                                                                                       !< accumulated (absolute) shear

 end type

 type(tParameters), dimension(:), allocatable, private :: &
   param                                                                                            !< containers of constitutive parameters (len Ninstance)

 type(tKinehardeningState), allocatable, dimension(:), private :: &
   dotState, &
   deltaState, &
   state

  
 public :: &
   plastic_kinehardening_init, &
   plastic_kinehardening_LpAndItsTangent, &
   plastic_kinehardening_dotState, &
   plastic_kinehardening_deltaState, &
   plastic_kinehardening_postResults
 private :: &
   kinetics


contains



!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   dEq0
 use debug, only: &
   debug_level, &
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
   PLASTICITY_kinehardening_label, &
   PLASTICITY_kinehardening_ID, &
   material_phase, &
   plasticState
 use config, only: &
   config_phase, &
   MATERIAL_partPhase
 use lattice

 implicit none

 integer(kind(undefined_ID)) :: &
   output_ID
 integer(pInt) :: &
   o, i, p,  &
   instance, &
   Ninstance, &
   NipcMyPhase, &
   outputSize, &
   offset_slip, &
   startIndex, endIndex, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState

 integer(pInt),          dimension(0), parameter :: emptyIntArray    = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 integer(kind(undefined_ID)) :: &
   outputID                                                                                         !< ID of each post result output
   
 character(len=65536), dimension(:), allocatable :: &
  outputs
 character(len=65536) :: &
   extmsg    = '', &
   structure    = ''

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_KINEHARDENING_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 Ninstance = int(count(phase_plasticity == PLASTICITY_KINEHARDENING_ID),pInt)
 if (Ninstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a,1x,i5,/)') '# instances:',Ninstance                                      
   
 allocate(plastic_kinehardening_sizePostResult(maxval(phase_Noutput),Ninstance), &
                                                                             source=0_pInt)
 allocate(plastic_kinehardening_output(maxval(phase_Noutput),Ninstance))
          plastic_kinehardening_output                                             = ''
 allocate(plastic_kinehardening_Noutput(Ninstance),                       source=0_pInt)
 allocate(plastic_kinehardening_Nslip(lattice_maxNslipFamily,Ninstance),  source=0_pInt)
 
 allocate(param(Ninstance))                                                                      ! one container of parameters per instance
 allocate(state(Ninstance))
 allocate(dotState(Ninstance))
 allocate(deltaState(Ninstance))
 
 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_KINEHARDENING_ID) cycle
     instance = phase_plasticityInstance(p)                                                     ! which instance of my phase
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             delta => deltaState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)))

   structure          = config_phase(p)%getString('lattice_structure')

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
   prm%aTolResistance = config_phase(p)%getFloat('atol_resistance',defaultVal=1.0_pReal)
   prm%aTolShear      = config_phase(p)%getFloat('atol_shear',     defaultVal=1.0e-6_pReal)

   ! sanity checks
   if (prm%aTolResistance <= 0.0_pReal) extmsg = trim(extmsg)//'aTolresistance '
   if (prm%aTolShear      <= 0.0_pReal) extmsg = trim(extmsg)//'aTolShear '

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

     prm%crss0                = config_phase(p)%getFloats('crss0',   requiredShape=shape(prm%Nslip))
     prm%tau1                 = config_phase(p)%getFloats('tau1', requiredShape=shape(prm%Nslip))
     prm%tau1_b               = config_phase(p)%getFloats('tau1_b', requiredShape=shape(prm%Nslip))
     prm%theta0               = config_phase(p)%getFloats('theta0', requiredShape=shape(prm%Nslip))
     prm%theta1               = config_phase(p)%getFloats('theta1', requiredShape=shape(prm%Nslip))
     prm%theta0_b             = config_phase(p)%getFloats('theta0_b', requiredShape=shape(prm%Nslip))
     prm%theta1_b             = config_phase(p)%getFloats('theta1_b', requiredShape=shape(prm%Nslip))

     ! expand: family => system
     prm%crss0   = math_expand(prm%crss0,  prm%Nslip)
     prm%tau1 = math_expand(prm%tau1,prm%Nslip)
     prm%tau1_b       = math_expand(prm%tau1_b,      prm%Nslip)
     prm%theta0 = math_expand(prm%theta0,prm%Nslip)
     prm%theta1 = math_expand(prm%theta1,prm%Nslip)
     prm%theta0_b = math_expand(prm%theta0_b,prm%Nslip)
     prm%theta1_b = math_expand(prm%theta1_b,prm%Nslip)

     prm%gdot0           = config_phase(p)%getFloat('gdot0')
     prm%n_slip          = config_phase(p)%getFloat('n_slip')

   endif slipActive

   
!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
        case ('resistance')
          outputID = merge(crss_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('accumulatedshear')
          outputID = merge(accshear_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('shearrate')
          outputID = merge(shearrate_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('resolvedstress')
          outputID = merge(resolvedstress_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('backstress')
          outputID = merge(crss_back_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('sense')
          outputID = merge(sense_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('chi0')
          outputID = merge(chi0_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip
        case ('gamma0')
          outputID = merge(gamma0_ID,undefined_ID,prm%totalNslip>0_pInt)
          outputSize = prm%totalNslip

      end select

      if (outputID /= undefined_ID) then
        plastic_kinehardening_Noutput(instance) = plastic_kinehardening_Noutput(instance) + 1_pInt
        plastic_kinehardening_output(i,phase_plasticityInstance(p)) = outputs(i)
        plastic_kinehardening_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
        prm%outputID = [prm%outputID , outputID]
      endif

   end do

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phase == p)                                                         ! number of constituents with my phase
   sizeDotState   = size(['crss     ','crss_back', 'accshear ']) * prm%TotalNslip
   sizeDeltaState = size(['sense ',   'chi0  ',    'gamma0'   ]) * prm%TotalNslip
   sizeState = sizeDotState + sizeDeltaState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,sizeDeltaState, &
                                      prm%totalNslip,0_pInt,0_pInt)
                                      
   plasticState(p)%sizePostResults = sum(plastic_kinehardening_sizePostResult(:,phase_plasticityInstance(p)))
   plasticState(p)%offsetDeltaState = sizeDotState


     startIndex = 1_pInt
     endIndex   =  prm%totalNslip
     stt%crss          => plasticState(p)%state    (startIndex:endIndex,1:NipcMyPhase)
     dot%crss          => plasticState(p)%dotState (startIndex:endIndex,1:NipcMyPhase)
     stt%crss = spread(prm%crss0, 2, NipcMyPhase)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex +  prm%totalNslip 
     stt%crss_back          => plasticState(p)%state    (startIndex:endIndex,1:NipcMyPhase)
     dot%crss_back          => plasticState(p)%dotState (startIndex:endIndex,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex +  prm%totalNslip
     stt%accshear          => plasticState(p)%state    (startIndex:endIndex,1:NipcMyPhase)
     dot%accshear          => plasticState(p)%dotState (startIndex:endIndex,1:NipcMyPhase)
     plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear
     
!----------------------------------------------------------------------------------------------
!locally define deltaState alias
     o = endIndex
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex +  prm%totalNslip
     stt%sense          => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     delta%sense        => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)
     
     startIndex = endIndex + 1_pInt
     endIndex   = endIndex +  prm%totalNslip
     stt%chi0           => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)
     delta%chi0         => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex +  prm%totalNslip
     stt%gamma0         => plasticState(p)%state     (startIndex  :endIndex  ,1:NipcMyPhase)         
     delta%gamma0       => plasticState(p)%deltaState(startIndex-o:endIndex-o,1:NipcMyPhase)
     
     plasticState(p)%state0 = plasticState(p)%state                                                 ! ToDo: this could be done centrally


     offset_slip = plasticState(p)%nSlip
     plasticState(p)%slipRate => &
       plasticState(p)%dotState(offset_slip+1:offset_slip+plasticState(p)%nSlip,1:NipcMyPhase)
     plasticState(p)%accumulatedSlip => &
       plasticState(p)%state(offset_slip+1:offset_slip+plasticState(p)%nSlip,1:NipcMyPhase) 
       
   end associate
 end do
 end subroutine plastic_kinehardening_init
     
!--------------------------------------------------------------------------------------------------
!  sanity checks
 
   !  if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
   !          .and. param(instance)%crss0 (1:nSlipFamilies)          < 0.0_pReal)) extmsg = trim(extmsg)//' crss0'
   !  if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
   !          .and. param(instance)%tau1  (1:nSlipFamilies)         <= 0.0_pReal)) extmsg = trim(extmsg)//' tau1'
   !  if (any(plastic_kinehardening_Nslip (1:nSlipFamilies,instance) > 0_pInt &
   !          .and. param(instance)%tau1_b(1:nSlipFamilies)         < 0.0_pReal)) extmsg = trim(extmsg)//' tau1_b'
   !  if (param(instance)%gdot0            <= 0.0_pReal) extmsg = trim(extmsg)//' gdot0'
   !  if (param(instance)%n_slip           <= 0.0_pReal) extmsg = trim(extmsg)//' n_slip'               



!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)
 use prec, only: &
   dNeq0
   
 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to the Mandel stress

 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer(pInt),               intent(in) :: &
   instance, &
   of

 integer(pInt) :: &
   j,k,l,m,n

   
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg, &
   dgdot_dtau_pos,dgdot_dtau_neg

 associate(prm => param(instance), stt => state(instance))
 Lp = 0.0_pReal 
 dLp_dMp = 0.0_pReal

 call kinetics(Mp,instance,of,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)

 do j = 1_pInt, prm%totalNslip
   Lp = Lp + (gdot_pos(j)+gdot_neg(j))*prm%Schmid_slip(1:3,1:3,j)
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtau_pos(j)*prm%Schmid_slip(k,l,j)*prm%nonSchmid_pos(m,n,j) &
                      + dgdot_dtau_neg(j)*prm%Schmid_slip(k,l,j)*prm%nonSchmid_neg(m,n,j)
 enddo
end associate

end subroutine plastic_kinehardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_deltaState(Mp,instance,of)
 use prec, only: &
   dNeq, &
   dEq0
   
 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer(pInt),              intent(in) :: &
   instance, &
   of

 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg, &
   sense

 associate( prm => param(instance), stt => state(instance), del => deltaState(instance))
 
 call kinetics(Mp,instance,of,gdot_pos,gdot_neg)
 sense = merge(state(instance)%sense(:,of), &                                                       ! keep existing...
               sign(1.0_pReal,gdot_pos+gdot_neg), &                                                 ! ...or have a defined 
               dEq0(gdot_pos+gdot_neg,1e-10_pReal))                                                 ! current sense of shear direction

#ifdef DEBUG
!         if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &               ! ToDo: We need an inverse mapping of ->el, ip, co
!            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
!                   .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
!           write(6,'(a)') '======= kinehardening delta state ======='
!         endif
!         if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
!            .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
!                   .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
!           write(6,'(i2,1x,f7.4,1x,f7.4)') j,sense(j),state(instance)%sense(j,of)
!         endif
#endif

!--------------------------------------------------------------------------------------------------
! switch in sense of shear?
 where(dNeq(sense,stt%sense(:,of),0.1_pReal))
   del%sense (:,of) = sense - stt%sense(:,of)                                                       ! switch sense
   del%chi0  (:,of) = abs(stt%crss_back(:,of)) - stt%chi0(:,of)                                     ! remember current backstress magnitude
   del%gamma0(:,of) = stt%accshear(:,of) - stt%gamma0(:,of)                                         ! remember current accumulated shear
 else where
   del%sense (:,of) = 0.0_pReal
   del%chi0  (:,of) = 0.0_pReal
   del%gamma0(:,of) = 0.0_pReal
 end where
 
 end associate

end subroutine plastic_kinehardening_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_kinehardening_dotState(Mp,instance,of)

 implicit none
 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer(pInt),              intent(in) :: &
   instance, &
   of

 integer(pInt) :: &
   j
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg
 real(pReal) :: &
   sumGamma
 

 associate( prm => param(instance), stt => state(instance), dot => dotState(instance))

 call kinetics(Mp,instance,of,gdot_pos,gdot_neg)
 dot%accshear(:,of) = abs(gdot_pos+gdot_neg)
 sumGamma = sum(stt%accshear(:,of))       
                 
 do j = 1_pInt, prm%totalNslip
   dot%crss(j,of) = dot_product(prm%interaction_SlipSlip(j,:),dot%accshear(:,of)) &
                  * (  prm%theta1(j) &
                     + (prm%theta0(j) - prm%theta1(j) + prm%theta0(j)*prm%theta1(j)*sumGamma/prm%tau1(j)) &
                                                     *exp(-sumGamma*prm%theta0(j)/prm%tau1(j)) &
                    )
 enddo
 dot%crss_back(:,of) = stt%sense(:,of)*dot%accshear(:,of) * &
          ( prm%theta1_b + &
            (prm%theta0_b - prm%theta1_b &
              + prm%theta0_b*prm%theta1_b/(prm%tau1_b+stt%chi0(:,of))*(stt%accshear(:,of)-stt%gamma0(:,of))&
            ) *exp(-(stt%accshear(:,of)-stt%gamma0(:,of)) *prm%theta0_b/(prm%tau1_b+stt%chi0(:,of))) &
          ) 
    
 end associate

end subroutine plastic_kinehardening_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_kinehardening_postResults(Mp,instance,of) result(postResults)
 use math, only: &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer(pInt),             intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_kinehardening_sizePostResult(:,instance))) :: &
   postResults
 integer(pInt) :: &
   o,c,j
   
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_pos,gdot_neg

 postResults = 0.0_pReal
 c = 0_pInt


 associate( prm => param(instance), stt => state(instance))
 
 call kinetics(Mp,instance,of,gdot_pos,gdot_neg)
 
 outputsLoop: do o = 1_pInt,plastic_kinehardening_Noutput(instance)
   select case(prm%outputID(o))
     case (crss_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%crss(:,of)
       c = c + prm%totalNslip
       
     case(crss_back_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%crss_back(:,of)
       c = c + prm%totalNslip
       
     case (sense_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%sense(:,of)
       c = c + prm%totalNslip
                                                                        
     case (chi0_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%chi0(:,of)
       c = c + prm%totalNslip
       
     case (gamma0_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%gamma0(:,of)
       c = c + prm%totalNslip
     
     case (accshear_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = stt%accshear(:,of)
       c = c + prm%totalNslip
       
     case (shearrate_ID)
       postResults(c+1_pInt:c+prm%totalNslip) = gdot_pos+gdot_neg
       c = c + prm%totalNslip

     case (resolvedstress_ID)
       do j = 1_pInt, prm%totalNslip
         postResults(c+j) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,j))
       enddo
       c = c + prm%totalNslip
       
   end select
 enddo outputsLoop
 end associate

end function plastic_kinehardening_postResults


!--------------------------------------------------------------------------------------------------
!> @brief calculates shear rates on slip systems and derivatives with respect to resolved stress
!> @details: Shear rates are calculated only optionally. NOTE: Against the common convention, the
!> result (i.e. intent(out)) variables are the last to have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,instance,of,gdot_pos,gdot_neg,dgdot_dtau_pos,dgdot_dtau_neg)
 use prec, only: &
  dNeq0
 use math, only: &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp
 integer(pInt),     intent(in) :: &
   instance, & 
   of
 real(pReal), dimension(param(instance)%totalNslip), intent(out) :: &
   gdot_pos, &
   gdot_neg
 real(pReal), dimension(param(instance)%totalNslip), optional, intent(out) :: &
   dgdot_dtau_pos, &
   dgdot_dtau_neg


 real(pReal), dimension(param(instance)%totalNslip) :: &
   tau_pos, &
   tau_neg
 integer(pInt) :: i
 logical       :: nonSchmidActive

 associate( prm => param(instance), stt => state(instance))
 
 nonSchmidActive = size(prm%nonSchmidCoeff) > 0_pInt

 do i = 1_pInt, prm%totalNslip
   tau_pos(i) =       math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i)) - stt%crss_back(i,of)
   tau_neg(i) = merge(math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,i)) - stt%crss_back(i,of), &
                           0.0_pReal, nonSchmidActive)
 enddo

 where(dNeq0(tau_pos))
   gdot_pos = prm%gdot0 * merge(0.5_pReal,1.0_pReal, nonSchmidActive) &                             ! 1/2 if non-Schmid active
            * sign(abs(tau_pos/stt%crss(:,of))**prm%n_slip,  tau_pos)
 else where
   gdot_pos = 0.0_pReal
 end where

 where(dNeq0(tau_neg))
   gdot_neg = prm%gdot0 * 0.5_pReal &                                                               ! only used if non-Schmid active, always 1/2
            * sign(abs(tau_neg/stt%crss(:,of))**prm%n_slip,  tau_neg)
 else where
   gdot_neg = 0.0_pReal
 end where

 if (present(dgdot_dtau_pos)) then
   where(dNeq0(gdot_pos))
     dgdot_dtau_pos = gdot_pos*prm%n_slip/tau_pos
   else where
     dgdot_dtau_pos = 0.0_pReal
   end where
 endif
 if (present(dgdot_dtau_neg)) then
   where(dNeq0(gdot_neg))
     dgdot_dtau_neg = gdot_neg*prm%n_slip/tau_neg
   else where
     dgdot_dtau_neg = 0.0_pReal
   end where
 endif
 end associate

end subroutine kinetics

end module plastic_kinehardening
