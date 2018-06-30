!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for phenomenological crystal plasticity formulation using a powerlaw
!! fitting
!--------------------------------------------------------------------------------------------------
module plastic_phenopowerlaw
 use prec, only: &
   pReal,&
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_phenopowerlaw_sizePostResult                                                             !< size of each post result output
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_phenopowerlaw_output                                                                     !< name of each post result output

 enum, bind(c)
   enumerator :: &
     undefined_ID, &
     resistance_slip_ID, &
     accumulatedshear_slip_ID, &
     shearrate_slip_ID, &
     resolvedstress_slip_ID, &
     totalshear_ID, &
     resistance_twin_ID, &
     accumulatedshear_twin_ID, &
     shearrate_twin_ID, &
     resolvedstress_twin_ID, &
     totalvolfrac_twin_ID
 end enum

 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(pInt) :: &
     totalNslip, &
     totalNtwin
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
     aTolResistance, &                                                                              ! default absolute tolerance 1 Pa
     aTolShear, &                                                                                   ! default absolute tolerance 1e-6
     aTolTwinfrac                                                                                   ! default absolute tolerance 1e-6
   integer(pInt), dimension(:),   allocatable :: &
     Nslip, &                                                                                       !< active number of slip systems per family
     Ntwin                                                                                          !< active number of twin systems per family
   real(pReal),   dimension(:),   allocatable :: &
     tau0_slip, &                                                                                   !< initial critical shear stress for slip
     tau0_twin, &                                                                                   !< initial critical shear stress for twin
     tausat_slip, &                                                                                 !< maximum critical shear stress for slip
     nonSchmidCoeff, &
     H_int                                                                                         !< per family hardening activity (optional)
   real(pReal),   dimension(:,:),   allocatable :: &
     interaction_SlipSlip, &                                                                        !< slip resistance from slip activity
     interaction_SlipTwin, &                                                                        !< slip resistance from twin activity
     interaction_TwinSlip, &                                                                        !< twin resistance from slip activity
     interaction_TwinTwin                                                                           !< twin resistance from twin activity

 integer(kind(undefined_ID)),         dimension(:),   allocatable          :: &
   outputID                                                                                     !< ID of each post result output
 end type

 type(tParameters), dimension(:), allocatable, target, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 type, private :: tPhenopowerlawState
   real(pReal), pointer,     dimension(:,:) :: &
     s_slip, &
     s_twin, &
     accshear_slip, &
     accshear_twin, &
     whole
   real(pReal), pointer,     dimension(:) :: &
     sumGamma, &
     sumF
 end type

 type(tPhenopowerlawState), allocatable, dimension(:), target,  private :: &
   dotState, &
   state

 public :: &
   plastic_phenopowerlaw_init, &
   plastic_phenopowerlaw_LpAndItsTangent, &
   plastic_phenopowerlaw_dotState, &
   plastic_phenopowerlaw_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_init
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use prec, only: &
   dEq0
 use debug, only: &
   debug_level, &
   debug_constitutive,&
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333, &
   math_expand
 use IO, only: &
   IO_warning, &
   IO_error, &
   IO_timeStamp
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   material_phase, &
   plasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice
 use numerics,only: &
   numerics_integrator

 implicit none

 integer(pInt) :: &
   maxNinstance, &
   instance,p,j,k, f,o, i,&
   NipcMyPhase, outputSize, &
   index_myFamily, index_otherFamily, &
   sizeState,sizeDotState, &
   startIndex, endIndex

 real(pReal), dimension(:,:), allocatable :: temp1, temp2

 integer(pInt),          dimension(0), parameter :: emptyIntArray = [integer(pInt)::]
 real(pReal),            dimension(0), parameter :: emptyRealArray = [real(pReal)::]
 character(len=65536),   dimension(0), parameter :: emptyStringArray = [character(len=65536)::]

 type(tParameters), pointer :: prm

 integer(kind(undefined_ID))                 :: &
   outputID                                                                                     !< ID of each post result output

 character(len=512) :: &
   extmsg    = ''
 character(len=65536), dimension(:), allocatable :: outputs

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_PHENOPOWERLAW_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_plasticity == PLASTICITY_PHENOPOWERLAW_ID),pInt)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_phenopowerlaw_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(plastic_phenopowerlaw_output(maxval(phase_Noutput),maxNinstance))
          plastic_phenopowerlaw_output = ''

 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
 allocate(state(maxNinstance))
 allocate(dotState(maxNinstance))

 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_PHENOPOWERLAW_ID) cycle
   instance = phase_plasticityInstance(p)
   prm => param(instance)

   prm%Nslip =  config_phase(p)%getInts('nslip',defaultVal=emptyIntArray)
   if (size(prm%Nslip) > count(lattice_NslipSystem(:,p) > 0_pInt))       call IO_error(150_pInt,ext_msg='Nslip')
   if (any(lattice_NslipSystem(1:size(prm%Nslip),p)-prm%Nslip < 0_pInt)) call IO_error(150_pInt,ext_msg='Nslip')
   prm%totalNslip = sum(prm%Nslip)

   if (prm%totalNslip > 0_pInt) then
     prm%tau0_slip            = config_phase(p)%getFloats('tau0_slip')
     prm%tausat_slip          = config_phase(p)%getFloats('tausat_slip')
     prm%interaction_SlipSlip = spread(config_phase(p)%getFloats('interaction_slipslip'),2,1)
     prm%H_int                = config_phase(p)%getFloats('h_int',&
                              defaultVal=[(0.0_pReal,i=1_pInt,size(prm%Nslip))])
     prm%nonSchmidCoeff       = config_phase(p)%getFloats('nonschmid_coefficients',&
                              defaultVal = emptyRealArray )

     prm%gdot0_slip           = config_phase(p)%getFloat('gdot0_slip')
     prm%n_slip               = config_phase(p)%getFloat('n_slip')
     prm%a_slip               = config_phase(p)%getFloat('a_slip')
     prm%h0_SlipSlip          = config_phase(p)%getFloat('h0_slipslip')
   endif

   prm%Ntwin            =  config_phase(p)%getInts('ntwin', defaultVal=emptyIntArray)
   if (size(prm%Ntwin) > count(lattice_NtwinSystem(:,p) > 0_pInt))       call IO_error(150_pInt,ext_msg='Ntwin')
   if (any(lattice_NtwinSystem(1:size(prm%Ntwin),p)-prm%Ntwin < 0_pInt)) call IO_error(150_pInt,ext_msg='Ntwin')
   prm%totalNtwin = sum(prm%Ntwin)

   if (prm%totalNtwin > 0_pInt) then
     prm%tau0_twin       =  config_phase(p)%getFloats('tau0_twin')
     prm%interaction_TwinTwin = spread(config_phase(p)%getFloats('interaction_twintwin'),2,1)

     prm%gdot0_twin  = config_phase(p)%getFloat('gdot0_twin')
     prm%n_twin      = config_phase(p)%getFloat('n_twin')
     prm%spr         = config_phase(p)%getFloat('s_pr')
     prm%twinB       = config_phase(p)%getFloat('twin_b')
     prm%twinC       = config_phase(p)%getFloat('twin_c')
     prm%twinD       = config_phase(p)%getFloat('twin_d')
     prm%twinE       = config_phase(p)%getFloat('twin_e')
     prm%h0_TwinTwin = config_phase(p)%getFloat('h0_twintwin')
   endif

   if (prm%totalNslip > 0_pInt .and. prm%totalNtwin > 0_pInt) then
     prm%interaction_SlipTwin = spread(config_phase(p)%getFloats('interaction_sliptwin'),2,1)
     prm%interaction_TwinSlip = spread(config_phase(p)%getFloats('interaction_twinslip'),2,1)
     prm%h0_TwinSlip = config_phase(p)%getFloat('h0_twinslip')
   endif


   prm%aTolResistance = config_phase(p)%getFloat('atol_resistance',defaultVal=1.0_pReal)
   prm%aTolShear      = config_phase(p)%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
   prm%aTolTwinfrac   = config_phase(p)%getFloat('atol_twinfrac',defaultVal=1.0e-6_pReal)

   outputs = config_phase(p)%getStrings('(output)',defaultVal=emptyStringArray)
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
        case ('resistance_slip')
          outputID = resistance_slip_ID
          outputSize = sum(prm%Nslip)
        case ('accumulatedshear_slip')
          outputID = accumulatedshear_slip_ID
          outputSize = sum(prm%Nslip)
        case ('shearrate_slip')
          outputID = shearrate_slip_ID
          outputSize = sum(prm%Nslip)
        case ('resolvedstress_slip')
          outputID = resolvedstress_slip_ID
          outputSize = sum(prm%Nslip)

        case ('resistance_twin')
          outputID = resistance_twin_ID
          outputSize = sum(prm%Ntwin)
        case ('accumulatedshear_twin')
          outputID = accumulatedshear_twin_ID
          outputSize = sum(prm%Ntwin)
        case ('shearrate_twin')
          outputID = shearrate_twin_ID
          outputSize = sum(prm%Ntwin)
        case ('resolvedstress_twin')
          outputID = resolvedstress_twin_ID
          outputSize = sum(prm%Ntwin)

        case ('totalvolfrac_twin')
          outputID = totalvolfrac_twin_ID
          outputSize = 1_pInt
        case ('totalshear')
          outputID = totalshear_ID
          outputSize = 1_pInt
      end select

      if (outputID /= undefined_ID) then
        plastic_phenopowerlaw_output(i,instance) = outputs(i)
        plastic_phenopowerlaw_sizePostResult(i,instance) = outputSize
        prm%outputID = [prm%outputID , outputID]
      endif

   end do

   extmsg = ''
   if (sum(prm%Nslip) > 0_pInt) then
     if (size(prm%tau0_slip) /= size(prm%Nslip))   call IO_error(211_pInt,ip=instance, &
            ext_msg='shape(tau0_slip) ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (size(prm%tausat_slip) /= size(prm%Nslip)) call IO_error(211_pInt,ip=instance, &
            ext_msg='shape(tausat_slip) ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (size(prm%H_int) /= size(prm%Nslip))       call IO_error(211_pInt,ip=instance, &
            ext_msg='shape(H_int) ('//PLASTICITY_PHENOPOWERLAW_label//')')

     if (any(prm%tau0_slip < 0.0_pReal .and. prm%Nslip > 0_pInt)) &
       extmsg = trim(extmsg)//"tau0_slip "
     if (any(prm%tausat_slip < prm%tau0_slip .and. prm%Nslip > 0_pInt)) &
       extmsg = trim(extmsg)//"tausat_slip "

     if (prm%gdot0_slip <= 0.0_pReal) extmsg = trim(extmsg)//" gdot0_slip "
     if (dEq0(prm%a_slip))            extmsg = trim(extmsg)//" a_slip " ! ToDo: negative values ok?
     if (dEq0(prm%n_slip))            extmsg = trim(extmsg)//" n_slip " ! ToDo: negative values ok?
   endif 

   if (sum(prm%Ntwin) > 0_pInt) then
     if (size(prm%tau0_twin) /= size(prm%ntwin)) call IO_error(211_pInt,ip=instance,&
            ext_msg='shape(tau0_twin) ('//PLASTICITY_PHENOPOWERLAW_label//')')

     if (any(prm%tau0_twin < 0.0_pReal .and. prm%Ntwin > 0_pInt)) &
       extmsg = trim(extmsg)//"tau0_twin "

     if (prm%gdot0_twin <= 0.0_pReal) extmsg = trim(extmsg)//"gdot0_twin "
     if (dEq0(prm%n_twin))            extmsg = trim(extmsg)//"n_twin " ! ToDo: negative values ok?
   endif

   if (prm%aTolResistance <= 0.0_pReal) extmsg = trim(extmsg)//"aTolresistance "
   if (prm%aTolShear      <= 0.0_pReal) extmsg = trim(extmsg)//"aTolShear "
   if (prm%aTolTwinfrac   <= 0.0_pReal) extmsg = trim(extmsg)//"atoltwinfrac "

   if (extmsg /= '') call IO_error(211_pInt,ip=instance,&
                                 ext_msg=trim(extmsg)//'('//PLASTICITY_PHENOPOWERLAW_label//')')

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phase == p)                                                   ! number of IPCs containing my phase
   sizeState = size(['tau_slip     ','accshear_slip']) * prm%TotalNslip &
             + size(['tau_twin     ','accshear_twin']) * prm%TotalNtwin &
             + size(['sum(gamma)', 'sum(f)    '])

   sizeDotState = sizeState
   plasticState(p)%sizeState = sizeState
   plasticState(p)%sizeDotState = sizeDotState
   plasticState(p)%sizePostResults = sum(plastic_phenopowerlaw_sizePostResult(:,instance))
   plasticState(p)%nSlip = sum(prm%Nslip)
   plasticState(p)%nTwin = sum(prm%Ntwin)
   allocate(plasticState(p)%aTolState          (   sizeState),             source=0.0_pReal)
   allocate(plasticState(p)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
   allocate(plasticState(p)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
   allocate(plasticState(p)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
   allocate(plasticState(p)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)

   allocate(plasticState(p)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
   allocate(plasticState(p)%deltaState               (0_pInt,NipcMyPhase), source=0.0_pReal)
   if (any(numerics_integrator == 1_pInt)) then
     allocate(plasticState(p)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(p)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
   endif
   if (any(numerics_integrator == 4_pInt)) &
     allocate(plasticState(p)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
   if (any(numerics_integrator == 5_pInt)) &
     allocate(plasticState(p)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)


!--------------------------------------------------------------------------------------------------
! calculate hardening matrices
   allocate(temp1(sum(prm%Nslip),sum(prm%Nslip)),source =0.0_pReal)
   allocate(temp2(sum(prm%Nslip),sum(prm%Ntwin)),source =0.0_pReal)
   mySlipFamilies: do f = 1_pInt,size(prm%Nslip,1)                                    ! >>> interaction slip -- X
     index_myFamily = sum(prm%Nslip(1:f-1_pInt))

     mySlipSystems: do j = 1_pInt,prm%Nslip(f)
       otherSlipFamilies: do o = 1_pInt,size(prm%Nslip,1)
         index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
         otherSlipSystems: do k = 1_pInt,prm%Nslip(o)
           temp1(index_myFamily+j,index_otherFamily+k) = &
               prm%interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                 sum(lattice_NslipSystem(1:f-1,p))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1,p))+k, &
                                                                 p),1)
       enddo otherSlipSystems; enddo otherSlipFamilies

       twinFamilies: do o = 1_pInt,size(prm%Ntwin,1)
         index_otherFamily = sum(prm%Ntwin(1:o-1_pInt))
         twinSystems: do k = 1_pInt,prm%Ntwin(o)
           temp2(index_myFamily+j,index_otherFamily+k) = &
               prm%interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                 sum(lattice_NslipSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,p))+k, &
                                                                 p),1)
       enddo twinSystems; enddo twinFamilies
     enddo mySlipSystems
   enddo mySlipFamilies
   prm%interaction_SlipSlip = temp1; deallocate(temp1)
   prm%interaction_SlipTwin = temp2; deallocate(temp2)
   

   allocate(temp1(sum(prm%Ntwin),sum(prm%Nslip)),source =0.0_pReal)
   allocate(temp2(sum(prm%Ntwin),sum(prm%Ntwin)),source =0.0_pReal)
   myTwinFamilies: do f = 1_pInt,size(prm%Ntwin,1)                                    ! >>> interaction twin -- X
     index_myFamily = sum(prm%Ntwin(1:f-1_pInt))
     myTwinSystems: do j = 1_pInt,prm%Ntwin(f)
       slipFamilies: do o = 1_pInt,size(prm%Nslip,1)
         index_otherFamily = sum(prm%Nslip(1:o-1_pInt))
         slipSystems: do k = 1_pInt,prm%Nslip(o)
           temp1(index_myFamily+j,index_otherFamily+k) = &
               prm%interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NslipSystem(1:o-1_pInt,p))+k, &
                                                                 p),1)
       enddo slipSystems; enddo slipFamilies

       otherTwinFamilies: do o = 1_pInt,size(prm%Ntwin,1)
         index_otherFamily = sum(prm%Ntwin(1:o-1_pInt))
         otherTwinSystems: do k = 1_pInt,prm%Ntwin(o)
           temp2(index_myFamily+j,index_otherFamily+k) = &
               prm%interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                 sum(lattice_NtwinSystem(1:f-1_pInt,p))+j, &
                                                                 sum(lattice_NtwinSystem(1:o-1_pInt,p))+k, &
                                                                 p),1)
       enddo otherTwinSystems; enddo otherTwinFamilies
     enddo myTwinSystems
   enddo myTwinFamilies
   prm%interaction_TwinSlip = temp1; deallocate(temp1)
   prm%interaction_TwinTwin = temp2; deallocate(temp2)

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1_pInt
   endIndex   = plasticState(p)%nSlip
   state   (instance)%s_slip=>plasticState(p)%state   (startIndex:endIndex,:)
   dotState(instance)%s_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%state0(startIndex:endIndex,:) = &
     spread(math_expand(prm%tau0_slip, prm%Nslip), 2, NipcMyPhase)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + plasticState(p)%nTwin
   state   (instance)%s_twin=>plasticState(p)%state   (startIndex:endIndex,:)
   dotState(instance)%s_twin=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%state0(startIndex:endIndex,:) = &
     spread(math_expand(prm%tau0_twin, prm%Ntwin), 2, NipcMyPhase)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolResistance

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + 1_pInt
   state   (instance)%sumGamma=>plasticState(p)%state   (startIndex,:)
   dotState(instance)%sumGamma=>plasticState(p)%dotState(startIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + 1_pInt
   state   (instance)%sumF=>plasticState(p)%state   (startIndex,:)
   dotState(instance)%sumF=>plasticState(p)%dotState(startIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolTwinFrac

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + plasticState(p)%nSlip
   state   (instance)%accshear_slip=>plasticState(p)%state   (startIndex:endIndex,:)
   dotState(instance)%accshear_slip=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear
   ! global alias
   plasticState(p)%slipRate =>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%accumulatedSlip =>plasticState(p)%state(startIndex:endIndex,:)

   startIndex = endIndex + 1_pInt
   endIndex   = endIndex + plasticState(p)%nTwin
   state   (instance)%accshear_twin=>plasticState(p)%state   (startIndex:endIndex,:)
   dotState(instance)%accshear_twin=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(p)%aTolState(startIndex:endIndex) = prm%aTolShear
   
   dotState(instance)%whole        =>plasticState(p)%dotState

 enddo

end subroutine plastic_phenopowerlaw_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,ipc,ip,el)
 use prec, only: &
   dNeq0
 use math, only: &
   math_Plain3333to99, &
   math_Mandel6to33
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_Stwin, &
   lattice_Stwin_v, &
   lattice_NslipSystem, &
   lattice_NtwinSystem
 use material, only: &
   phasememberAt, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation

 integer(pInt) :: &
   index_myFamily, &
   f,i,j,k,l,m,n, &
   of, &
   ph
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg, &
   gdot_slip_pos,gdot_slip_neg, &
   dgdot_dtauslip_pos,dgdot_dtauslip_neg, &
   gdot_twin,dgdot_dtautwin,tau_twin
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333                                                                                   !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal), dimension(3,3,2) :: &
   nonSchmid_tensor
 type(tParameters),         pointer :: prm
 type(tPhenopowerlawState), pointer :: stt

 of = phasememberAt(ipc,ip,el)
 ph = material_phase(ipc,ip,el)

 prm     => param(phase_plasticityInstance(ph))
 stt     => state(phase_plasticityInstance(ph))

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! Slip part
 j = 0_pInt
 slipFamilies: do f = 1_pInt,size(prm%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                          ! at which index starts my family
   slipSystems: do i = 1_pInt,prm%Nslip(f)
     j = j+1_pInt

     ! Calculation of Lp
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     do k = 1,size(prm%nonSchmidCoeff)
       tau_slip_pos = tau_slip_pos + prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + prm%nonSchmidCoeff(k)*&
                                           lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + prm%nonSchmidCoeff(k)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo
     gdot_slip_pos = 0.5_pReal*prm%gdot0_slip* &
                    ((abs(tau_slip_pos)/(stt%s_slip(j,of)))**prm%n_slip)*sign(1.0_pReal,tau_slip_pos)

     gdot_slip_neg = 0.5_pReal*prm%gdot0_slip* &
                    ((abs(tau_slip_neg)/(stt%s_slip(j,of)))**prm%n_slip)*sign(1.0_pReal,tau_slip_neg)

     Lp = Lp + (1.0_pReal-stt%sumF(of))*& 
               (gdot_slip_pos+gdot_slip_neg)*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (dNeq0(gdot_slip_pos)) then  !@ Philip: Needed? No division
       dgdot_dtauslip_pos = gdot_slip_pos*prm%n_slip/tau_slip_pos
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_pos*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                     nonSchmid_tensor(m,n,1)
     endif

     if (dNeq0(gdot_slip_neg)) then !@ Philip: Needed? No division
       dgdot_dtauslip_neg = gdot_slip_neg*prm%n_slip/tau_slip_neg
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_neg*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                     nonSchmid_tensor(m,n,2)
     endif
   enddo slipSystems
 enddo slipFamilies

!--------------------------------------------------------------------------------------------------
! Twinning part
 j = 0_pInt
 twinFamilies: do f = 1_pInt,size(prm%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                      ! at which index starts my family
   twinSystems: do i = 1_pInt,prm%Ntwin(f)
     j = j+1_pInt

     ! Calculation of Lp
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
     gdot_twin = (1.0_pReal-stt%sumF(of))*prm%gdot0_twin*&
                    (abs(tau_twin)/stt%s_twin(j,of))**&
                    prm%n_twin*max(0.0_pReal,sign(1.0_pReal,tau_twin))
     Lp = Lp + gdot_twin*lattice_Stwin(1:3,1:3,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (dNeq0(gdot_twin)) then !@ Philip: Needed? No division
       dgdot_dtautwin = gdot_twin*prm%n_twin/tau_twin
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtautwin*lattice_Stwin(k,l,index_myFamily+i,ph)* &
                                                  lattice_Stwin(m,n,index_myFamily+i,ph)
     endif
   enddo twinSystems
 enddo twinFamilies

 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)


end subroutine plastic_phenopowerlaw_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_phenopowerlaw_dotState(Tstar_v,ipc,ip,el)
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin
 use material, only: &
   material_phase, &
   phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                    !< microstructure state

 integer(pInt) :: &
   ph, &
   f,i,j,k, &
   index_myFamily, &
   of
 real(pReal) :: &
   c_SlipSlip,c_TwinSlip,c_TwinTwin, &
   ssat_offset, &
   tau_slip_pos,tau_slip_neg,tau_twin

 real(pReal), dimension(param(phase_plasticityInstance(material_phase(ipc,ip,el)))%totalNslip) :: &
   gdot_slip,left_SlipSlip,right_SlipSlip
 real(pReal), dimension(param(phase_plasticityInstance(material_phase(ipc,ip,el)))%totalNtwin) :: &
   gdot_twin

 type(tParameters),         pointer :: prm
 type(tPhenopowerlawState), pointer :: dst,stt

 of = phasememberAt(ipc,ip,el)
 ph = material_phase(ipc,ip,el)

 prm => param(phase_plasticityInstance(ph))
 stt => state(phase_plasticityInstance(ph))
 dst => dotState(phase_plasticityInstance(ph))

 dst%whole(:,of) = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = prm%h0_slipslip * (1.0_pReal + prm%twinC*stt%sumF(of)** prm%twinB)
 c_TwinSlip = prm%h0_TwinSlip * stt%sumGamma(of)**prm%twinE
 c_TwinTwin = prm%h0_TwinTwin * stt%sumF(of)**prm%twinD

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors and calculate dot gammas
 ssat_offset = prm%spr*sqrt(stt%sumF(of))
 j = 0_pInt
 slipFamilies1: do f =1_pInt,size(prm%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems1: do i = 1_pInt,prm%Nslip(f)
     j = j+1_pInt
     left_SlipSlip(j) = 1.0_pReal + prm%H_int(f)                         ! modified no system-dependent left part
     right_SlipSlip(j) = abs(1.0_pReal-stt%s_slip(j,of) / (prm%tausat_slip(f)+ssat_offset)) **prm%a_slip &
                       * sign(1.0_pReal,1.0_pReal-stt%s_slip(j,of) / (prm%tausat_slip(f)+ssat_offset))

!--------------------------------------------------------------------------------------------------
! Calculation of dot gamma
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmidSystems: do k = 1,size(prm%nonSchmidCoeff)
       tau_slip_pos = tau_slip_pos + prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg +prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems
     gdot_slip(j) = prm%gdot0_slip*0.5_pReal* &
                  ( (abs(tau_slip_pos)/(stt%s_slip(j,of)))**prm%n_slip*sign(1.0_pReal,tau_slip_pos) &
                   +(abs(tau_slip_neg)/(stt%s_slip(j,of)))**prm%n_slip*sign(1.0_pReal,tau_slip_neg))
   enddo slipSystems1
 enddo slipFamilies1

 j = 0_pInt
 twinFamilies1: do f = 1_pInt,size(prm%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems1: do i = 1_pInt,prm%Ntwin(f)
     j = j+1_pInt

!--------------------------------------------------------------------------------------------------
! Calculation of dot vol frac
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
     gdot_twin(j) = (1.0_pReal-stt%sumF(of))*&                                       ! 1-F
                    prm%gdot0_twin*&
                    (abs(tau_twin)/stt%s_twin(j,of))**&
                    prm%n_twin*max(0.0_pReal,sign(1.0_pReal,tau_twin))
    enddo twinSystems1
  enddo twinFamilies1

!--------------------------------------------------------------------------------------------------
! calculate the overall hardening based on above
 do j = 1_pInt,prm%totalNslip
   dst%s_slip(j,of) = c_SlipSlip * left_SlipSlip(j) * &                                         ! evolution of slip resistance j
     dot_product(prm%interaction_SlipSlip(j,1:prm%totalNslip),right_SlipSlip*abs(gdot_slip)) + &    ! dot gamma_slip modulated by right-side slip factor
     dot_product(prm%interaction_SlipTwin(j,1:prm%totalNtwin),gdot_twin)                            ! dot gamma_twin modulated by right-side twin factor
 enddo
 dst%sumGamma(of) = dst%sumGamma(of) + sum(abs(gdot_slip))
 dst%accshear_slip(1:prm%totalNslip,of) = abs(gdot_slip)

 j = 0_pInt
 twinFamilies2: do f = 1_pInt,size(prm%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems2: do i = 1_pInt,prm%Ntwin(f)
     j = j+1_pInt
     dst%s_twin(j,of) = &                                                                       ! evolution of twin resistance j
       c_TwinSlip * dot_product(prm%interaction_TwinSlip(j,1:prm%totalNslip),abs(gdot_slip)) + &    ! dot gamma_slip modulated by right-side slip factor
       c_TwinTwin * dot_product(prm%interaction_TwinTwin(j,1:prm%totalNtwin),gdot_twin)             ! dot gamma_twin modulated by right-side twin factor
     if (stt%sumF(of) < 0.98_pReal) &                                                               ! ensure twin volume fractions stays below 1.0
       dst%sumF(of) = dst%sumF(of) + gdot_twin(j)/lattice_shearTwin(index_myFamily+i,ph)
      dst%accshear_twin(j,of) = abs(gdot_twin(j))
   enddo twinSystems2
 enddo twinFamilies2

end subroutine plastic_phenopowerlaw_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_phenopowerlaw_postResults(Tstar_v,ipc,ip,el)
 use material, only: &
   material_phase, &
   plasticState, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_Sslip_v, &
   lattice_Stwin_v, &
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid

 implicit none
 real(pReal), dimension(6), intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                        !< microstructure state

 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
   plastic_phenopowerlaw_postResults

 integer(pInt) :: &
   instance,ph, of, &
   o,f,i,c,j,k, &
   index_myFamily
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg,tau
 type(tParameters), pointer :: prm

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)
 prm => param(instance)


 plastic_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
     case (resistance_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+prm%totalNslip) = state(instance)%s_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip

     case (accumulatedshear_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+prm%totalNslip) = state(instance)%accshear_slip(1:prm%totalNslip,of)
       c = c + prm%totalNslip

     case (shearrate_slip_ID)
       j = 0_pInt
       slipFamilies1: do f = 1_pInt,size(prm%Nslip,1)
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems1: do i = 1_pInt,prm%Nslip(f)
           j = j + 1_pInt
           tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
           tau_slip_neg  = tau_slip_pos
           do k = 1,lattice_NnonSchmid(ph)
             tau_slip_pos = tau_slip_pos +prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
             tau_slip_neg = tau_slip_neg +prm%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
           enddo
           plastic_phenopowerlaw_postResults(c+j) = prm%gdot0_slip*0.5_pReal* &
                    ((abs(tau_slip_pos)/state(instance)%s_slip(j,of))**prm%n_slip &
                    *sign(1.0_pReal,tau_slip_pos) &
                    +(abs(tau_slip_neg)/(state(instance)%s_slip(j,of)))**prm%n_slip &
                    *sign(1.0_pReal,tau_slip_neg))
         enddo slipSystems1
       enddo slipFamilies1
       c = c + prm%totalNslip

     case (resolvedstress_slip_ID)
       j = 0_pInt
       slipFamilies2: do f = 1_pInt,size(prm%Nslip,1)
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems2: do i = 1_pInt,prm%Nslip(f)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
         enddo slipSystems2
       enddo slipFamilies2
       c = c + prm%totalNslip

     case (totalshear_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = &
                             state(instance)%sumGamma(of)
       c = c + 1_pInt

     case (resistance_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+prm%totalNtwin) = &
                            state(instance)%s_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin

     case (accumulatedshear_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+prm%totalNtwin) = &
                             state(instance)%accshear_twin(1:prm%totalNtwin,of)
       c = c + prm%totalNtwin
     case (shearrate_twin_ID)
       j = 0_pInt
       twinFamilies1: do f = 1_pInt,size(prm%Ntwin,1)
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems1: do i = 1_pInt,prm%Ntwin(f)
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
           plastic_phenopowerlaw_postResults(c+j) = (1.0_pReal-state(instance)%sumF(of))*&  ! 1-F
                                                         prm%gdot0_twin*&
                                                         (abs(tau)/state(instance)%s_twin(j,of))**&
                                           prm%n_twin*max(0.0_pReal,sign(1.0_pReal,tau))
         enddo twinSystems1
       enddo twinFamilies1
       c = c + prm%totalNtwin

     case (resolvedstress_twin_ID)
       j = 0_pInt
       twinFamilies2: do f = 1_pInt,size(prm%Ntwin,1)
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems2: do i = 1_pInt,prm%Ntwin(f)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
         enddo twinSystems2
       enddo twinFamilies2
       c = c + prm%totalNtwin

     case (totalvolfrac_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = state(instance)%sumF(of)
       c = c + 1_pInt

   end select
 enddo outputsLoop

end function plastic_phenopowerlaw_postResults

end module plastic_phenopowerlaw
