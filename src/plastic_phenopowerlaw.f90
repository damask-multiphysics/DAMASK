!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief  phenomenological crystal plasticity formulation using a powerlaw fitting
!--------------------------------------------------------------------------------------------------
module plastic_phenopowerlaw
 use prec
 use debug
 use math
 use IO
 use material
 use config
 use lattice
 use discretization
 use results

 implicit none
 private
 
 integer,          dimension(:,:),   allocatable, target, public :: &
   plastic_phenopowerlaw_sizePostResult                                                             !< size of each post result output
 character(len=64), dimension(:,:),   allocatable, target, public :: &
   plastic_phenopowerlaw_output                                                                     !< name of each post result output

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     resistance_slip_ID, &
     accumulatedshear_slip_ID, &
     shearrate_slip_ID, &
     resolvedstress_slip_ID, &
     resistance_twin_ID, &
     accumulatedshear_twin_ID, &
     shearrate_twin_ID, &
     resolvedstress_twin_ID
 end enum

 type :: tParameters
   real(pReal) :: &
     gdot0_slip, &                                                                                  !< reference shear strain rate for slip
     gdot0_twin, &                                                                                  !< reference shear strain rate for twin
     n_slip, &                                                                                      !< stress exponent for slip
     n_twin, &                                                                                      !< stress exponent for twin
     spr, &                                                                                         !< push-up factor for slip saturation due to twinning
     twinB, &
     twinC, &
     twinD, &
     twinE, &
     h0_SlipSlip, &                                                                                 !< reference hardening slip - slip
     h0_TwinSlip, &                                                                                 !< reference hardening twin - slip
     h0_TwinTwin, &                                                                                 !< reference hardening twin - twin
     a_slip, &
     aTolResistance, &                                                                              !< absolute tolerance for integration of xi
     aTolShear, &                                                                                   !< absolute tolerance for integration of gamma
     aTolTwinfrac                                                                                   !< absolute tolerance for integration of f
   real(pReal), allocatable, dimension(:) :: &
     xi_slip_0, &                                                                                   !< initial critical shear stress for slip
     xi_twin_0, &                                                                                   !< initial critical shear stress for twin
     xi_slip_sat, &                                                                                 !< maximum critical shear stress for slip
     nonSchmidCoeff, &
     H_int, &                                                                                       !< per family hardening activity (optional)
     gamma_twin_char                                                                                !< characteristic shear for twins
   real(pReal), allocatable, dimension(:,:) :: &
     interaction_SlipSlip, &                                                                        !< slip resistance from slip activity
     interaction_SlipTwin, &                                                                        !< slip resistance from twin activity
     interaction_TwinSlip, &                                                                        !< twin resistance from slip activity
     interaction_TwinTwin                                                                           !< twin resistance from twin activity
   real(pReal),  allocatable, dimension(:,:,:) :: &
     Schmid_slip, &
     Schmid_twin, &
     nonSchmid_pos, &
     nonSchmid_neg
   integer :: &
     totalNslip, &                                                                                  !< total number of active slip system
     totalNtwin                                                                                     !< total number of active twin systems
   integer,      allocatable, dimension(:) :: &
     Nslip, &                                                                                       !< number of active slip systems for each family
     Ntwin                                                                                          !< number of active twin systems for each family
   integer(kind(undefined_ID)), allocatable, dimension(:) :: &
     outputID                                                                                       !< ID of each post result output
 end type tParameters

 type :: tPhenopowerlawState
   real(pReal), pointer, dimension(:,:) :: &
     xi_slip, &
     xi_twin, &
     gamma_slip, &
     gamma_twin
 end type tPhenopowerlawState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
 type(tParameters),         allocatable, dimension(:) :: param
 type(tPhenopowerlawState), allocatable, dimension(:) :: &
   dotState, &
   state

 public :: &
   plastic_phenopowerlaw_init, &
   plastic_phenopowerlaw_LpAndItsTangent, &
   plastic_phenopowerlaw_dotState, &
   plastic_phenopowerlaw_postResults, &
   plastic_phenopowerlaw_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_init

 integer :: &
   Ninstance, &
   p, i, &
   NipcMyPhase, outputSize, &
   sizeState, sizeDotState, &
   startIndex, endIndex

 integer,                dimension(0), parameter :: emptyIntArray    = [integer::]
 real(pReal),            dimension(0), parameter :: emptyRealArray   = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 integer(kind(undefined_ID)) :: &
   outputID

 character(len=pStringLen) :: &
   extmsg = ''
 character(len=65536), dimension(:), allocatable :: &
   outputs

 write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_PHENOPOWERLAW_label//' init  -+>>>'

 Ninstance = count(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

 allocate(plastic_phenopowerlaw_sizePostResult(maxval(phase_Noutput),Ninstance),source=0)
 allocate(plastic_phenopowerlaw_output(maxval(phase_Noutput),Ninstance))
          plastic_phenopowerlaw_output = ''

 allocate(param(Ninstance))
 allocate(state(Ninstance))
 allocate(dotState(Ninstance))

 do p = 1, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_PHENOPOWERLAW_ID) cycle
   associate(prm => param(phase_plasticityInstance(p)), &
             dot => dotState(phase_plasticityInstance(p)), &
             stt => state(phase_plasticityInstance(p)), &
             config => config_phase(p))

!--------------------------------------------------------------------------------------------------
!  optional parameters that need to be defined
   prm%twinB          = config%getFloat('twin_b',defaultVal=1.0_pReal)
   prm%twinC          = config%getFloat('twin_c',defaultVal=0.0_pReal)
   prm%twinD          = config%getFloat('twin_d',defaultVal=0.0_pReal)
   prm%twinE          = config%getFloat('twin_e',defaultVal=0.0_pReal)

   prm%aTolResistance = config%getFloat('atol_resistance',defaultVal=1.0_pReal)
   prm%aTolShear      = config%getFloat('atol_shear',     defaultVal=1.0e-6_pReal)
   prm%aTolTwinfrac   = config%getFloat('atol_twinfrac',  defaultVal=1.0e-6_pReal)

   ! sanity checks
   if (prm%aTolResistance <= 0.0_pReal) extmsg = trim(extmsg)//' aTolresistance'
   if (prm%aTolShear      <= 0.0_pReal) extmsg = trim(extmsg)//' aTolShear'
   if (prm%aTolTwinfrac   <= 0.0_pReal) extmsg = trim(extmsg)//' atoltwinfrac'

!--------------------------------------------------------------------------------------------------
! slip related parameters
   prm%Nslip      = config%getInts('nslip',defaultVal=emptyIntArray)
   prm%totalNslip = sum(prm%Nslip)
   slipActive: if (prm%totalNslip > 0) then
     prm%Schmid_slip = lattice_SchmidMatrix_slip(prm%Nslip,config%getString('lattice_structure'),&
                                                 config%getFloat('c/a',defaultVal=0.0_pReal))

     if(trim(config%getString('lattice_structure')) == 'bcc') then
       prm%nonSchmidCoeff = config%getFloats('nonschmid_coefficients',&
                                                 defaultVal = emptyRealArray)
       prm%nonSchmid_pos  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,+1)
       prm%nonSchmid_neg  = lattice_nonSchmidMatrix(prm%Nslip,prm%nonSchmidCoeff,-1)
     else
       prm%nonSchmid_pos  = prm%Schmid_slip
       prm%nonSchmid_neg  = prm%Schmid_slip
     endif
     prm%interaction_SlipSlip = lattice_interaction_SlipBySlip(prm%Nslip, &
                                                               config%getFloats('interaction_slipslip'), &
                                                               config%getString('lattice_structure'))

     prm%xi_slip_0   = config%getFloats('tau0_slip',   requiredSize=size(prm%Nslip))
     prm%xi_slip_sat = config%getFloats('tausat_slip', requiredSize=size(prm%Nslip))
     prm%H_int       = config%getFloats('h_int',       requiredSize=size(prm%Nslip), &
                                        defaultVal=[(0.0_pReal,i=1,size(prm%Nslip))])

     prm%gdot0_slip  = config%getFloat('gdot0_slip')
     prm%n_slip      = config%getFloat('n_slip')
     prm%a_slip      = config%getFloat('a_slip')
     prm%h0_SlipSlip = config%getFloat('h0_slipslip')

     ! expand: family => system
     prm%xi_slip_0   = math_expand(prm%xi_slip_0,  prm%Nslip)
     prm%xi_slip_sat = math_expand(prm%xi_slip_sat,prm%Nslip)
     prm%H_int       = math_expand(prm%H_int,      prm%Nslip)

     ! sanity checks
     if (    prm%gdot0_slip  <= 0.0_pReal)      extmsg = trim(extmsg)//' gdot0_slip'
     if (    prm%a_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' a_slip'
     if (    prm%n_slip      <= 0.0_pReal)      extmsg = trim(extmsg)//' n_slip'
     if (any(prm%xi_slip_0   <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_slip_0'
     if (any(prm%xi_slip_sat <= 0.0_pReal))     extmsg = trim(extmsg)//' xi_slip_sat'
   else slipActive
     allocate(prm%interaction_SlipSlip(0,0))
     allocate(prm%xi_slip_0(0))
   endif slipActive

!--------------------------------------------------------------------------------------------------
! twin related parameters
   prm%Ntwin      = config%getInts('ntwin', defaultVal=emptyIntArray)
   prm%totalNtwin = sum(prm%Ntwin)
   twinActive: if (prm%totalNtwin > 0) then
     prm%Schmid_twin          = lattice_SchmidMatrix_twin(prm%Ntwin,config%getString('lattice_structure'),&
                                                          config%getFloat('c/a',defaultVal=0.0_pReal))
     prm%interaction_TwinTwin = lattice_interaction_TwinByTwin(prm%Ntwin,&
                                                               config%getFloats('interaction_twintwin'), &
                                                               config%getString('lattice_structure'))
     prm%gamma_twin_char      = lattice_characteristicShear_twin(prm%Ntwin,config%getString('lattice_structure'),&
                                                            config%getFloat('c/a'))

     prm%xi_twin_0            = config%getFloats('tau0_twin',requiredSize=size(prm%Ntwin))

     prm%gdot0_twin           = config%getFloat('gdot0_twin')
     prm%n_twin               = config%getFloat('n_twin')
     prm%spr                  = config%getFloat('s_pr')
     prm%h0_TwinTwin          = config%getFloat('h0_twintwin')

     ! expand: family => system
     prm%xi_twin_0   = math_expand(prm%xi_twin_0,  prm%Ntwin)

     ! sanity checks
     if (prm%gdot0_twin <= 0.0_pReal)  extmsg = trim(extmsg)//' gdot0_twin'
     if (prm%n_twin     <= 0.0_pReal)  extmsg = trim(extmsg)//' n_twin'
   else twinActive
     allocate(prm%interaction_TwinTwin(0,0))
     allocate(prm%xi_twin_0(0))
     allocate(prm%gamma_twin_char(0))
   endif twinActive

!--------------------------------------------------------------------------------------------------
! slip-twin related parameters
   slipAndTwinActive: if (prm%totalNslip > 0 .and. prm%totalNtwin > 0) then
     prm%interaction_SlipTwin = lattice_interaction_SlipByTwin(prm%Nslip,prm%Ntwin,&
                                                               config%getFloats('interaction_sliptwin'), &
                                                               config%getString('lattice_structure'))
     prm%interaction_TwinSlip = lattice_interaction_TwinBySlip(prm%Ntwin,prm%Nslip,&
                                                               config%getFloats('interaction_twinslip'), &
                                                               config%getString('lattice_structure'))
   else slipAndTwinActive
     allocate(prm%interaction_SlipTwin(prm%TotalNslip,prm%TotalNtwin))                              ! at least one dimension is 0
     allocate(prm%interaction_TwinSlip(prm%TotalNtwin,prm%TotalNslip))                              ! at least one dimension is 0
     prm%h0_TwinSlip = 0.0_pReal
   endif slipAndTwinActive

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') &
     call IO_error(211,ext_msg=trim(extmsg)//'('//PLASTICITY_PHENOPOWERLAW_label//')')

!--------------------------------------------------------------------------------------------------
!  output pararameters
   outputs = config%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))

       case ('resistance_slip')
         outputID = merge(resistance_slip_ID,undefined_ID,prm%totalNslip>0)
         outputSize = prm%totalNslip
       case ('accumulatedshear_slip')
         outputID = merge(accumulatedshear_slip_ID,undefined_ID,prm%totalNslip>0)
         outputSize = prm%totalNslip
       case ('shearrate_slip')
         outputID = merge(shearrate_slip_ID,undefined_ID,prm%totalNslip>0)
         outputSize = prm%totalNslip
       case ('resolvedstress_slip')
         outputID = merge(resolvedstress_slip_ID,undefined_ID,prm%totalNslip>0)
         outputSize = prm%totalNslip

       case ('resistance_twin')
         outputID = merge(resistance_twin_ID,undefined_ID,prm%totalNtwin>0)
         outputSize = prm%totalNtwin
       case ('accumulatedshear_twin')
         outputID = merge(accumulatedshear_twin_ID,undefined_ID,prm%totalNtwin>0)
         outputSize = prm%totalNtwin
       case ('shearrate_twin')
         outputID = merge(shearrate_twin_ID,undefined_ID,prm%totalNtwin>0)
         outputSize = prm%totalNtwin
       case ('resolvedstress_twin')
         outputID = merge(resolvedstress_twin_ID,undefined_ID,prm%totalNtwin>0)
         outputSize = prm%totalNtwin

     end select

     if (outputID /= undefined_ID) then
       plastic_phenopowerlaw_output(i,phase_plasticityInstance(p)) = outputs(i)
       plastic_phenopowerlaw_sizePostResult(i,phase_plasticityInstance(p)) = outputSize
       prm%outputID = [prm%outputID, outputID]
     endif

   enddo

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
   sizeDotState = size(['tau_slip  ','gamma_slip']) * prm%totalNslip &
                + size(['tau_twin  ','gamma_twin']) * prm%totalNtwin
   sizeState = sizeDotState

   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0, &
                                      prm%totalNslip,prm%totalNtwin,0)
   plasticState(p)%sizePostResults = sum(plastic_phenopowerlaw_sizePostResult(:,phase_plasticityInstance(p)))

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1
   endIndex   = prm%totalNslip
   stt%xi_slip => plasticState(p)%state   (startIndex:endIndex,:)
   stt%xi_slip = spread(prm%xi_slip_0, 2, NipcMyPhase)
   dot%xi_slip => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNtwin
   stt%xi_twin => plasticState(p)%state   (startIndex:endIndex,:)
   stt%xi_twin = spread(prm%xi_twin_0, 2, NipcMyPhase)
   dot%xi_twin => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNslip
   stt%gamma_slip => plasticState(p)%state   (startIndex:endIndex,:)
   dot%gamma_slip => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear
   ! global alias
   plasticState(p)%slipRate        => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%accumulatedSlip => plasticState(p)%state(startIndex:endIndex,:)

   startIndex = endIndex + 1
   endIndex   = endIndex + prm%totalNtwin
   stt%gamma_twin => plasticState(p)%state   (startIndex:endIndex,:)
   dot%gamma_twin => plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear

   plasticState(p)%state0 = plasticState(p)%state                                                   ! ToDo: this could be done centrally

   end associate

 enddo

end subroutine plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!> @details asummes that deformation by dislocation glide affects twinned and untwinned volume
!  equally (Taylor assumption). Twinning happens only in untwinned volume
!--------------------------------------------------------------------------------------------------
pure subroutine plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dMp,Mp,instance,of)

 real(pReal), dimension(3,3),     intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to the Mandel stress

 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer,               intent(in) :: &
   instance, &
   of

 integer :: &
   i,k,l,m,n
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_slip_pos,gdot_slip_neg, &
   dgdot_dtauslip_pos,dgdot_dtauslip_neg
 real(pReal), dimension(param(instance)%totalNtwin) :: &
   gdot_twin,dgdot_dtautwin

 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal

 associate(prm => param(instance))

 call kinetics_slip(Mp,instance,of,gdot_slip_pos,gdot_slip_neg,dgdot_dtauslip_pos,dgdot_dtauslip_neg)
 slipSystems: do i = 1, prm%totalNslip
   Lp = Lp + (gdot_slip_pos(i)+gdot_slip_neg(i))*prm%Schmid_slip(1:3,1:3,i)
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtauslip_pos(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_pos(m,n,i) &
                      + dgdot_dtauslip_neg(i) * prm%Schmid_slip(k,l,i) * prm%nonSchmid_neg(m,n,i)
 enddo slipSystems

 call kinetics_twin(Mp,instance,of,gdot_twin,dgdot_dtautwin)
 twinSystems: do i = 1, prm%totalNtwin
   Lp = Lp + gdot_twin(i)*prm%Schmid_twin(1:3,1:3,i)
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + dgdot_dtautwin(i)*prm%Schmid_twin(k,l,i)*prm%Schmid_twin(m,n,i)
 enddo twinSystems

 end associate

end subroutine plastic_phenopowerlaw_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_dotState(Mp,instance,of)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer,                      intent(in) :: &
   instance, &
   of

 real(pReal) :: &
   c_SlipSlip,c_TwinSlip,c_TwinTwin, &
   xi_slip_sat_offset,&
   sumGamma,sumF
 real(pReal), dimension(param(instance)%totalNslip) :: &
   left_SlipSlip,right_SlipSlip, &
   gdot_slip_pos,gdot_slip_neg

 associate(prm => param(instance), stt => state(instance), dot => dotState(instance))

 sumGamma = sum(stt%gamma_slip(:,of))
 sumF     = sum(stt%gamma_twin(:,of)/prm%gamma_twin_char)

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = prm%h0_slipslip * (1.0_pReal + prm%twinC*sumF** prm%twinB)
 c_TwinSlip = prm%h0_TwinSlip * sumGamma**prm%twinE
 c_TwinTwin = prm%h0_TwinTwin * sumF**prm%twinD

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors
 left_SlipSlip  = 1.0_pReal + prm%H_int
 xi_slip_sat_offset = prm%spr*sqrt(sumF)
 right_SlipSlip = abs(1.0_pReal-stt%xi_slip(:,of) / (prm%xi_slip_sat+xi_slip_sat_offset)) **prm%a_slip &
                * sign(1.0_pReal,1.0_pReal-stt%xi_slip(:,of) / (prm%xi_slip_sat+xi_slip_sat_offset))

!--------------------------------------------------------------------------------------------------
! shear rates
 call kinetics_slip(Mp,instance,of,gdot_slip_pos,gdot_slip_neg)
 dot%gamma_slip(:,of) = abs(gdot_slip_pos+gdot_slip_neg)
 call kinetics_twin(Mp,instance,of,dot%gamma_twin(:,of))

!--------------------------------------------------------------------------------------------------
! hardening
 dot%xi_slip(:,of) = c_SlipSlip * left_SlipSlip * &
                     matmul(prm%interaction_SlipSlip,dot%gamma_slip(:,of)*right_SlipSlip) &
                   + matmul(prm%interaction_SlipTwin,dot%gamma_twin(:,of))

 dot%xi_twin(:,of) = c_TwinSlip * matmul(prm%interaction_TwinSlip,dot%gamma_slip(:,of)) &
                   + c_TwinTwin * matmul(prm%interaction_TwinTwin,dot%gamma_twin(:,of))
 end associate

end subroutine plastic_phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_phenopowerlaw_postResults(Mp,instance,of) result(postResults)

 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer,                     intent(in) :: &
   instance, &
   of

 real(pReal), dimension(sum(plastic_phenopowerlaw_sizePostResult(:,instance))) :: &
   postResults

 integer :: &
   o,c,i
 real(pReal), dimension(param(instance)%totalNslip) :: &
   gdot_slip_pos,gdot_slip_neg

 c = 0

 associate(prm => param(instance), stt => state(instance))

 outputsLoop: do o = 1,size(prm%outputID)
   select case(prm%outputID(o))

     case (resistance_slip_ID)
       postResults(c+1:c+prm%totalNslip) = stt%xi_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (accumulatedshear_slip_ID)
       postResults(c+1:c+prm%totalNslip) = stt%gamma_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip
     case (shearrate_slip_ID)
       call kinetics_slip(Mp,instance,of,gdot_slip_pos,gdot_slip_neg)
       postResults(c+1:c+prm%totalNslip) = gdot_slip_pos+gdot_slip_neg
       c = c + prm%totalNslip
     case (resolvedstress_slip_ID)
       do i = 1, prm%totalNslip
         postResults(c+i) = math_mul33xx33(Mp,prm%Schmid_slip(1:3,1:3,i))
       enddo
       c = c + prm%totalNslip

     case (resistance_twin_ID)
       postResults(c+1:c+prm%totalNtwin) = stt%xi_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (accumulatedshear_twin_ID)
       postResults(c+1:c+prm%totalNtwin) = stt%gamma_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (shearrate_twin_ID)
       call kinetics_twin(Mp,instance,of,postResults(c+1:c+prm%totalNtwin))
       c = c + prm%totalNtwin
     case (resolvedstress_twin_ID)
       do i = 1, prm%totalNtwin
         postResults(c+i) = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,i))
       enddo
       c = c + prm%totalNtwin

   end select
 enddo outputsLoop

 end associate

end function plastic_phenopowerlaw_postResults


!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_results(instance,group)
#if defined(PETSc) || defined(DAMASK_HDF5)

  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
  
  integer :: o

  associate(prm => param(instance), stt => state(instance))
  outputsLoop: do o = 1,size(prm%outputID)
    select case(prm%outputID(o))

      case (resistance_slip_ID)
        call results_writeDataset(group,stt%xi_slip,   'xi_sl', &
                                  'resistance against plastic slip','Pa')
      case (accumulatedshear_slip_ID)
        call results_writeDataset(group,stt%gamma_slip,'gamma_sl', &
                                  'plastic shear','1')
                                  
      case (resistance_twin_ID)
        call results_writeDataset(group,stt%xi_twin,   'xi_tw', &
                                  'resistance against twinning','Pa')
      case (accumulatedshear_twin_ID)
        call results_writeDataset(group,stt%gamma_twin,'gamma_tw', &
                                  'twinning shear','1')
        
    end select
  enddo outputsLoop
  end associate
  
#else
  integer,          intent(in) :: instance
  character(len=*), intent(in) :: group
#endif

end subroutine plastic_phenopowerlaw_results


!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on slip systems and their derivatives with respect to resolved stress
!> @details Derivatives are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_slip(Mp,instance,of, &
                              gdot_slip_pos,gdot_slip_neg,dgdot_dtau_slip_pos,dgdot_dtau_slip_neg)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer,                      intent(in) :: &
   instance, &
   of

 real(pReal),                  intent(out), dimension(param(instance)%totalNslip) :: &
   gdot_slip_pos, &
   gdot_slip_neg
 real(pReal),                  intent(out), optional, dimension(param(instance)%totalNslip) :: &
   dgdot_dtau_slip_pos, &
   dgdot_dtau_slip_neg

 real(pReal), dimension(param(instance)%totalNslip) :: &
   tau_slip_pos, &
   tau_slip_neg
 integer :: i
 logical :: nonSchmidActive

 associate(prm => param(instance), stt => state(instance))

 nonSchmidActive = size(prm%nonSchmidCoeff) > 0

 do i = 1, prm%totalNslip
   tau_slip_pos(i) =       math_mul33xx33(Mp,prm%nonSchmid_pos(1:3,1:3,i))
   tau_slip_neg(i) = merge(math_mul33xx33(Mp,prm%nonSchmid_neg(1:3,1:3,i)), &
                           0.0_pReal, nonSchmidActive)
 enddo

 where(dNeq0(tau_slip_pos))
   gdot_slip_pos = prm%gdot0_slip * merge(0.5_pReal,1.0_pReal, nonSchmidActive) &                   ! 1/2 if non-Schmid active
                 * sign(abs(tau_slip_pos/stt%xi_slip(:,of))**prm%n_slip,  tau_slip_pos)
 else where
   gdot_slip_pos = 0.0_pReal
 end where

 where(dNeq0(tau_slip_neg))
   gdot_slip_neg = prm%gdot0_slip * 0.5_pReal &                                                     ! only used if non-Schmid active, always 1/2
                 * sign(abs(tau_slip_neg/stt%xi_slip(:,of))**prm%n_slip,  tau_slip_neg)
 else where
   gdot_slip_neg = 0.0_pReal
 end where

 if (present(dgdot_dtau_slip_pos)) then
   where(dNeq0(gdot_slip_pos))
     dgdot_dtau_slip_pos = gdot_slip_pos*prm%n_slip/tau_slip_pos
   else where
     dgdot_dtau_slip_pos = 0.0_pReal
   end where
 endif
 if (present(dgdot_dtau_slip_neg)) then
   where(dNeq0(gdot_slip_neg))
     dgdot_dtau_slip_neg = gdot_slip_neg*prm%n_slip/tau_slip_neg
   else where
     dgdot_dtau_slip_neg = 0.0_pReal
   end where
 endif
 end associate

end subroutine kinetics_slip


!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on twin systems and their derivatives with respect to resolved stress.
!  twinning is assumed to take place only in untwinned volume.
!> @details Derivates are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_twin(Mp,instance,of,&
                              gdot_twin,dgdot_dtau_twin)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 integer,                      intent(in) :: &
   instance, &
   of

 real(pReal), dimension(param(instance)%totalNtwin), intent(out) :: &
   gdot_twin
 real(pReal), dimension(param(instance)%totalNtwin), intent(out), optional :: &
   dgdot_dtau_twin

 real(pReal), dimension(param(instance)%totalNtwin) :: &
   tau_twin
 integer :: i

 associate(prm => param(instance), stt => state(instance))

 do i = 1, prm%totalNtwin
   tau_twin(i)  = math_mul33xx33(Mp,prm%Schmid_twin(1:3,1:3,i))
 enddo

 where(tau_twin > 0.0_pReal)
   gdot_twin = (1.0_pReal-sum(stt%gamma_twin(:,of)/prm%gamma_twin_char)) &                          ! only twin in untwinned volume fraction
             * prm%gdot0_twin*(abs(tau_twin)/stt%xi_twin(:,of))**prm%n_twin
 else where
   gdot_twin = 0.0_pReal
 end where

 if (present(dgdot_dtau_twin)) then
   where(dNeq0(gdot_twin))
     dgdot_dtau_twin = gdot_twin*prm%n_twin/tau_twin
   else where
     dgdot_dtau_twin = 0.0_pReal
   end where
 endif

 end associate

end subroutine kinetics_twin

end module plastic_phenopowerlaw
