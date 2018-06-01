!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
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
     H_int, &                                                                                       !< per family hardening activity (optional)
     interaction_SlipSlip, &                                                                        !< slip resistance from slip activity
     interaction_SlipTwin, &                                                                        !< slip resistance from twin activity
     interaction_TwinSlip, &                                                                        !< twin resistance from slip activity
     interaction_TwinTwin                                                                           !< twin resistance from twin activity
   real(pReal),   dimension(:,:),   allocatable :: &
     matrix_SlipSlip, &                                                                        !< slip resistance from slip activity
     matrix_SlipTwin, &                                                                        !< slip resistance from twin activity
     matrix_TwinSlip, &                                                                        !< twin resistance from slip activity
     matrix_TwinTwin                                                                           !< twin resistance from twin activity

 integer(kind(undefined_ID)),         dimension(:),   allocatable          :: &
   outputID                                                                                     !< ID of each post result output
 end type

 type(tParameters), dimension(:), allocatable, target, private :: param                                     !< containers of constitutive parameters (len Ninstance)

 type, private :: tPhenopowerlawState
   real(pReal), pointer,     dimension(:,:) :: &
     s_slip, &
     s_twin, &
     accshear_slip, &
     accshear_twin
   real(pReal), pointer,     dimension(:) :: &
     sumGamma, &
     sumF
 end type

 type(tPhenopowerlawState), allocatable, dimension(:), private :: &
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
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_warning, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_PHENOPOWERLAW_label, &
   PLASTICITY_PHENOPOWERLAW_ID, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase, &
   phaseConfig

 use lattice
 use numerics,only: &
   numerics_integrator

 implicit none

 integer(pInt) :: &
   maxNinstance, &
   instance,phase,j,k, f,o, i,&
   NipcMyPhase, outputSize, &
   offset_slip, index_myFamily, index_otherFamily, &
   sizeState,sizeDotState, sizeDeltaState, &
   startIndex, endIndex
 integer(pInt), dimension(0), parameter :: emptyInt = [integer(pInt)::]
 real(pReal), dimension(0), parameter :: emptyReal = [real(pReal)::]

 type(tParameters), pointer :: p

 integer(kind(undefined_ID))                 :: &
   outputID                                                                                     !< ID of each post result output

 character(len=65536) :: &
   extmsg    = ''
 character(len=64), dimension(:), allocatable :: outputs

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

 do phase = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(phase) == PLASTICITY_PHENOPOWERLAW_ID) then
     instance = phase_plasticityInstance(phase)
     p => param(instance)

     p%Nslip            =  phaseConfig(phase)%getIntArray('nslip',defaultVal=emptyInt)
     !if (size > Nchunks_SlipFamilies + 1_pInt) call IO_error(150_pInt,ext_msg=extmsg)
     if (sum(p%Nslip) > 0_pInt) then
       p%tau0_slip   =  phaseConfig(phase)%getFloatArray('tau0_slip')
       p%tausat_slip =  phaseConfig(phase)%getFloatArray('tausat_slip')
       p%H_int       = phaseConfig(phase)%getFloatArray('h_int',defaultVal=[(0.0_pReal,i=1_pInt,size(p%Nslip))])
       print*, (shape(p%H_int))
       print*, (shape(p%Nslip))
       p%interaction_SlipSlip = phaseConfig(phase)%getFloatArray('interaction_slipslip')
       p%nonSchmidCoeff       = phaseConfig(phase)%getFloatArray('nonschmid_coefficients',&
                                 defaultVal = [real(pReal)::1] )
       p%gdot0_slip = phaseConfig(phase)%getFloat('gdot0_slip')
       p%n_slip = phaseConfig(phase)%getFloat('n_slip')
       p%a_slip = phaseConfig(phase)%getFloat('a_slip')
        p%h0_SlipSlip = phaseConfig(phase)%getFloat('h0_slipslip')
     endif

     p%Ntwin            =  phaseConfig(phase)%getIntArray('ntwin', defaultVal=emptyInt)
     !if (size > Nchunks_SlipFamilies + 1_pInt) call IO_error(150_pInt,ext_msg=extmsg)
     if (sum(p%Ntwin) > 0_pInt) then
       p%tau0_twin       =  phaseConfig(phase)%getFloatArray('tau0_twin')
       p%interaction_TwinTwin = phaseConfig(phase)%getFloatArray('interaction_twintwin')
         p%gdot0_twin = phaseConfig(phase)%getFloat('gdot0_twin')
         p%n_twin = phaseConfig(phase)%getFloat('n_twin')
         p%spr = phaseConfig(phase)%getFloat('s_pr')
         p%twinB = phaseConfig(phase)%getFloat('twin_b')
         p%twinC = phaseConfig(phase)%getFloat('twin_c')
         p%twinD = phaseConfig(phase)%getFloat('twin_d')
         p%twinE = phaseConfig(phase)%getFloat('twin_e')
         p%h0_TwinTwin = phaseConfig(phase)%getFloat('h0_twintwin')
      endif
     if (sum(p%Nslip) > 0_pInt .and. sum(p%Ntwin) > 0_pInt) then
       p%interaction_SlipTwin = phaseConfig(phase)%getFloatArray('interaction_sliptwin')
       p%interaction_TwinSlip = phaseConfig(phase)%getFloatArray('interaction_twinslip')
       p%h0_TwinSlip = phaseConfig(phase)%getFloat('h0_twinslip')
     endif

     allocate(p%matrix_SlipSlip(sum(p%Nslip),sum(p%Nslip)),source =0.0_pReal)
     allocate(p%matrix_SlipTwin(sum(p%Nslip),sum(p%Ntwin)),source =0.0_pReal)
     allocate(p%matrix_TwinSlip(sum(p%Ntwin),sum(p%Nslip)),source =0.0_pReal)
     allocate(p%matrix_TwinTwin(sum(p%Ntwin),sum(p%Ntwin)),source =0.0_pReal)
         p%aTolResistance = phaseConfig(phase)%getFloat('atol_resistance',defaultVal=1.0_pReal)
         p%aTolShear      = phaseConfig(phase)%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
         p%aTolTwinfrac   = phaseConfig(phase)%getFloat('atol_twinfrac',defaultVal=1.0e-6_pReal)
     outputs = phaseConfig(phase)%getStrings('(output)')
     allocate(p%outputID(0))
     do i=1_pInt, size(outputs)
       outputID = undefined_ID
       select case(outputs(i))
          case ('resistance_slip')
            outputID = resistance_slip_ID
            outputSize = sum(p%Nslip)
          case ('acumulatedshear_slip','accumulated_shear_slip')
            outputID = accumulatedshear_slip_ID
            outputSize = sum(p%Nslip)
          case ('shearrate_slip')
            outputID = shearrate_slip_ID
            outputSize = sum(p%Nslip)
          case ('resolvedstress_slip')
            outputID = resolvedstress_slip_ID
            outputSize = sum(p%Nslip)

          case ('resistance_twin')
            outputID = resistance_twin_ID
            outputSize = sum(p%Ntwin)
          case ('accumulatedshear_twin','accumulated_shear_twin')
            outputID = accumulatedshear_twin_ID
            outputSize = sum(p%Ntwin)
          case ('shearrate_twin')
            outputID = shearrate_twin_ID
            outputSize = sum(p%Ntwin)
          case ('resolvedstress_twin')
            outputID = resolvedstress_twin_ID
            outputSize = sum(p%Ntwin)

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
          p%outputID = [p%outputID , outputID]
        endif

      end do

!--------------------------------------------------------------------------------------------------
! parameters independent of number of slip/twin systems
extmsg = ''
if (size(p%tau0_slip) /= size(p%nslip)) extmsg = trim(extmsg)//" shape(tau0_slip) "
if (size(p%tausat_slip) /= size(p%nslip)) extmsg = trim(extmsg)//" shape(tausat_slip) "
if (size(p%H_int) /= size(p%nslip)) extmsg = trim(extmsg)//" shape(h_int) "
if (size(p%tau0_twin) /= size(p%ntwin)) extmsg = trim(extmsg)//" shape(tau0_twin) "
     if (extmsg /= '') call IO_error(211_pInt,ip=instance,&
                            ext_msg=trim(extmsg)//'('//PLASTICITY_PHENOPOWERLAW_label//')')

if (any(p%tau0_slip < 0.0_pReal .and. p%Nslip > 0_pInt)) &
  extmsg = trim(extmsg)//" 'tau0_slip' "
if (any(p%tau0_slip < p%tausat_slip .and. p%Nslip > 0_pInt)) &
  extmsg = trim(extmsg)//" 'tausat_slip' "
if (any(p%gdot0_slip <= 0.0_pReal .and. p%Nslip > 0_pInt)) &
  extmsg = trim(extmsg)//" 'tausat_slip' "
if (p%n_slip <= 0.0_pReal) extmsg = trim(extmsg)//" 'n_slip' "

       !if (any(dEq0(p%a_slip) .and. sum(p%Nslip) > 0)) &
       !  call IO_error(211_pInt,el=instance,ext_msg='a_slip ('//PLASTICITY_PHENOPOWERLAW_label//')')

     !    if (any(p%tau0_twin < 0.0_pReal .and. &
    !           p%Ntwin(:) > 0)) &
    !     call IO_error(211_pInt,el=instance,ext_msg='tau0_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')
    !   if (    p%gdot0_twin <= 0.0_pReal .and. &
    !       any(p%Ntwin(:) > 0)) &
    !     call IO_error(211_pInt,el=instance,ext_msg='gdot0_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')
    !   if (    p%n_twin <= 0.0_pReal .and. &
    !      any(p%Ntwin(:) > 0)) &
    !     call IO_error(211_pInt,el=instance,ext_msg='n_twin ('//PLASTICITY_PHENOPOWERLAW_label//')')

     if (p%aTolResistance <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='aTolResistance ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (p%aTolShear      <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='aTolShear ('//PLASTICITY_PHENOPOWERLAW_label//')')
     if (p%aTolTwinfrac   <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='aTolTwinfrac ('//PLASTICITY_PHENOPOWERLAW_label//')')




     NipcMyPhase = count(material_phase == phase)                                                   ! number of IPCs containing my phase

!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = size(['tau_slip     ','accshear_slip']) * sum(p%nslip) &
               + size(['tau_twin     ','accshear_twin']) * sum(p%ntwin) &
               + size(['sum(gamma)', 'sum(f)    '])

     sizeDotState = sizeState
     sizeDeltaState = 0_pInt
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%nSlip = sum(p%Nslip)
     plasticState(phase)%nTwin = sum(p%Ntwin)
     allocate(plasticState(phase)%aTolState          (   sizeState),             source=0.0_pReal)
     allocate(plasticState(phase)%state0             (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%partionedState0    (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NipcMyPhase), source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NipcMyPhase), source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase), source=0.0_pReal)

     offset_slip = plasticState(phase)%nSlip+plasticState(phase)%nTwin+2_pInt
     plasticState(phase)%slipRate => &
       plasticState(phase)%dotState(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip => &
       plasticState(phase)%state(offset_slip+1:offset_slip+plasticState(phase)%nSlip,1:NipcMyPhase)

!--------------------------------------------------------------------------------------------------
! calculate hardening matrices
     mySlipFamilies: do f = 1_pInt,size(p%Nslip,1)                                    ! >>> interaction slip -- X
       index_myFamily = sum(p%Nslip(1:f-1_pInt))

       mySlipSystems: do j = 1_pInt,p%Nslip(f)
         otherSlipFamilies: do o = 1_pInt,size(p%Nslip,1)
           index_otherFamily = sum(p%Nslip(1:o-1_pInt))
           otherSlipSystems: do k = 1_pInt,p%Nslip(o)
             p%matrix_SlipSlip(index_myFamily+j,index_otherFamily+k) = &
                 p%interaction_SlipSlip(lattice_interactionSlipSlip( &
                                                                   sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                                                   sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                                                   phase))
         enddo otherSlipSystems; enddo otherSlipFamilies

         twinFamilies: do o = 1_pInt,size(p%Ntwin,1)
           index_otherFamily = sum(p%Ntwin(1:o-1_pInt))
           twinSystems: do k = 1_pInt,p%Ntwin(o)
             p%matrix_SlipTwin(index_myFamily+j,index_otherFamily+k) = &
                 p%interaction_SlipTwin(lattice_interactionSlipTwin( &
                                                                   sum(lattice_NslipSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                   phase))
         enddo twinSystems; enddo twinFamilies
       enddo mySlipSystems
     enddo mySlipFamilies

     myTwinFamilies: do f = 1_pInt,size(p%Ntwin,1)                                    ! >>> interaction twin -- X
       index_myFamily = sum(p%Ntwin(1:f-1_pInt))
       myTwinSystems: do j = 1_pInt,p%Ntwin(f)
         slipFamilies: do o = 1_pInt,size(p%Nslip,1)
           index_otherFamily = sum(p%Nslip(1:o-1_pInt))
           slipSystems: do k = 1_pInt,p%Nslip(o)
             p%matrix_TwinSlip(index_myFamily+j,index_otherFamily+k) = &
                 p%interaction_TwinSlip(lattice_interactionTwinSlip( &
                                                                   sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NslipSystem(1:o-1_pInt,phase))+k, &
                                                                   phase))
         enddo slipSystems; enddo slipFamilies

         otherTwinFamilies: do o = 1_pInt,size(p%Ntwin,1)
           index_otherFamily = sum(p%Ntwin(1:o-1_pInt))
           otherTwinSystems: do k = 1_pInt,p%Ntwin(o)
             p%matrix_TwinTwin(index_myFamily+j,index_otherFamily+k) = &
                 p%interaction_TwinTwin(lattice_interactionTwinTwin( &
                                                                   sum(lattice_NtwinSystem(1:f-1_pInt,phase))+j, &
                                                                   sum(lattice_NtwinSystem(1:o-1_pInt,phase))+k, &
                                                                   phase))
         enddo otherTwinSystems; enddo otherTwinFamilies
       enddo myTwinSystems
     enddo myTwinFamilies

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
     startIndex = 1_pInt
     endIndex   = plasticState(phase)%nSlip
     state   (instance)%s_slip=>plasticState(phase)%state   (startIndex:endIndex,:)
     dotState(instance)%s_slip=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(phase)%state0(startIndex:endIndex,:) = &
       spread(math_expand(p%tau0_slip, p%Nslip), 2, NipcMyPhase)

     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolResistance

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + plasticState(phase)%nTwin
     state   (instance)%s_twin=>plasticState(phase)%state   (startIndex:endIndex,:)
     dotState(instance)%s_twin=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(phase)%state0(startIndex:endIndex,:) = &
       spread(math_expand(p%tau0_twin, p%Ntwin), 2, NipcMyPhase)
     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolResistance

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + 1_pInt
     state   (instance)%sumGamma=>plasticState(phase)%state   (startIndex,:)
     dotState(instance)%sumGamma=>plasticState(phase)%dotState(startIndex,:)
     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolShear

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + 1_pInt
     state   (instance)%sumF=>plasticState(phase)%state   (startIndex,:)
     dotState(instance)%sumF=>plasticState(phase)%dotState(startIndex,:)
     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolTwinFrac

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + plasticState(phase)%nSlip
     state   (instance)%accshear_slip=>plasticState(phase)%state   (startIndex:endIndex,:)
     dotState(instance)%accshear_slip=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolShear
     ! global alias
     plasticState(phase)%slipRate =>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(phase)%accumulatedSlip =>plasticState(phase)%state(startIndex:endIndex,:)

     startIndex = endIndex + 1_pInt
     endIndex   = endIndex + plasticState(phase)%nTwin
     state   (instance)%accshear_twin=>plasticState(phase)%state   (startIndex:endIndex,:)
     dotState(instance)%accshear_twin=>plasticState(phase)%dotState(startIndex:endIndex,:)
     plasticState(phase)%aTolState(startIndex:endIndex) = p%aTolShear

   endif
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
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_NnonSchmid
 use material, only: &
   phaseAt, phasememberAt, &
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
   instance, &
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

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! Slip part
 j = 0_pInt
 slipFamilies: do f = 1_pInt,size(param(instance)%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                          ! at which index starts my family
   slipSystems: do i = 1_pInt,param(instance)%Nslip(f)
     j = j+1_pInt

     ! Calculation of Lp
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmid_tensor(1:3,1:3,1) = lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,1)
     do k = 1,size(param(instance)%nonSchmidCoeff)
       tau_slip_pos = tau_slip_pos + param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg + param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
       nonSchmid_tensor(1:3,1:3,1) = nonSchmid_tensor(1:3,1:3,1) + param(instance)%nonSchmidCoeff(k)*&
                                           lattice_Sslip(1:3,1:3,2*k,index_myFamily+i,ph)
       nonSchmid_tensor(1:3,1:3,2) = nonSchmid_tensor(1:3,1:3,2) + param(instance)%nonSchmidCoeff(k)*&
                                           lattice_Sslip(1:3,1:3,2*k+1,index_myFamily+i,ph)
     enddo
     gdot_slip_pos = 0.5_pReal*param(instance)%gdot0_slip* &
                    ((abs(tau_slip_pos)/(state(instance)%s_slip(j,of))) &
                    **param(instance)%n_slip)*sign(1.0_pReal,tau_slip_pos)

     gdot_slip_neg = 0.5_pReal*param(instance)%gdot0_slip* &
                    ((abs(tau_slip_neg)/(state(instance)%s_slip(j,of))) &
                    **param(instance)%n_slip)*sign(1.0_pReal,tau_slip_neg)

     Lp = Lp + (1.0_pReal-state(instance)%sumF(of))*&                                             ! 1-F
               (gdot_slip_pos+gdot_slip_neg)*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (dNeq0(gdot_slip_pos)) then
       dgdot_dtauslip_pos = gdot_slip_pos*param(instance)%n_slip/tau_slip_pos
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip_pos*lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                     nonSchmid_tensor(m,n,1)
     endif

     if (dNeq0(gdot_slip_neg)) then
       dgdot_dtauslip_neg = gdot_slip_neg*param(instance)%n_slip/tau_slip_neg
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
 twinFamilies: do f = 1_pInt,size(param(instance)%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                      ! at which index starts my family
   twinSystems: do i = 1_pInt,param(instance)%Ntwin(f)
     j = j+1_pInt

     ! Calculation of Lp
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
     gdot_twin = (1.0_pReal-state(instance)%sumF(of))*&                                          ! 1-F
                    param(instance)%gdot0_twin*&
                    (abs(tau_twin)/state(instance)%s_twin(j,of))**&
                    param(instance)%n_twin*max(0.0_pReal,sign(1.0_pReal,tau_twin))
     Lp = Lp + gdot_twin*lattice_Stwin(1:3,1:3,index_myFamily+i,ph)

     ! Calculation of the tangent of Lp
     if (dNeq0(gdot_twin)) then
       dgdot_dtautwin = gdot_twin*param(instance)%n_twin/tau_twin
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
   lattice_maxNslipFamily, &
   lattice_maxNtwinFamily, &
   lattice_NslipSystem, &
   lattice_NtwinSystem, &
   lattice_shearTwin, &
   lattice_NnonSchmid
 use material, only: &
   material_phase, &
   phaseAt, phasememberAt, &
   plasticState, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element                                                                                    !< microstructure state

 integer(pInt) :: &
   instance,ph, &
   f,i,j,k, &
   index_myFamily, nslip,ntwin,&
   of
 real(pReal) :: &
   c_SlipSlip,c_TwinSlip,c_TwinTwin, &
   ssat_offset, &
   tau_slip_pos,tau_slip_neg,tau_twin

 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Nslip) :: &
   gdot_slip,left_SlipSlip,left_SlipTwin,right_SlipSlip,right_TwinSlip
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%Ntwin) :: &
   gdot_twin,left_TwinSlip,left_TwinTwin,right_SlipTwin,right_TwinTwin

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 nSlip= sum(param(instance)%nslip)
 nTwin= sum(param(instance)%nTwin)

 plasticState(ph)%dotState(:,of) = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! system-independent (nonlinear) prefactors to M_Xx (X influenced by x) matrices
 c_SlipSlip = param(instance)%h0_slipslip*&
              (1.0_pReal + param(instance)%twinC*state(instance)%sumF(of)**&
                                                           param(instance)%twinB)
 c_TwinSlip = param(instance)%h0_TwinSlip*&
              state(instance)%sumGamma(of)**param(instance)%twinE
 c_TwinTwin = param(instance)%h0_TwinTwin*&
              state(instance)%sumF(of)**param(instance)%twinD

!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors and calculate dot gammas
 ssat_offset = param(instance)%spr*sqrt(state(instance)%sumF(of))
 j = 0_pInt
 slipFamilies1: do f =1_pInt,size(param(instance)%Nslip,1)
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems1: do i = 1_pInt,param(instance)%Nslip(f)
     j = j+1_pInt
     left_SlipSlip(j) = 1.0_pReal + param(instance)%H_int(f)                         ! modified no system-dependent left part
     left_SlipTwin(j) = 1.0_pReal                                                                   ! no system-dependent left part
     right_SlipSlip(j) = abs(1.0_pReal-state(instance)%s_slip(j,of) / &
                                    (param(instance)%tausat_slip(f)+ssat_offset)) &
                         **param(instance)%a_slip&
                         *sign(1.0_pReal,1.0_pReal-state(instance)%s_slip(j,of) / &
                                    (param(instance)%tausat_slip(f)+ssat_offset))
     right_TwinSlip(j) = 1.0_pReal                                                                  ! no system-dependent part

!--------------------------------------------------------------------------------------------------
! Calculation of dot gamma
     tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
     tau_slip_neg  = tau_slip_pos
     nonSchmidSystems: do k = 1,lattice_NnonSchmid(ph)
       tau_slip_pos = tau_slip_pos + param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,  index_myFamily+i,ph))
       tau_slip_neg = tau_slip_neg +param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
     enddo nonSchmidSystems
     gdot_slip(j) = param(instance)%gdot0_slip*0.5_pReal* &
                  ((abs(tau_slip_pos)/(state(instance)%s_slip(j,of)))**param(instance)%n_slip &
                  *sign(1.0_pReal,tau_slip_pos) &
                  +(abs(tau_slip_neg)/(state(instance)%s_slip(j,of)))**param(instance)%n_slip &
                  *sign(1.0_pReal,tau_slip_neg))
   enddo slipSystems1
 enddo slipFamilies1



 j = 0_pInt
 twinFamilies1: do f = 1_pInt,size(param(instance)%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems1: do i = 1_pInt,param(instance)%Ntwin(f)
     j = j+1_pInt
     left_TwinSlip(j)  = 1.0_pReal                                                                  ! no system-dependent left part
     left_TwinTwin(j)  = 1.0_pReal                                                                  ! no system-dependent left part
     right_SlipTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part
     right_TwinTwin(j) = 1.0_pReal                                                                  ! no system-dependent right part

!--------------------------------------------------------------------------------------------------
! Calculation of dot vol frac
     tau_twin  = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
     gdot_twin(j) = (1.0_pReal-state(instance)%sumF(of))*&                                       ! 1-F
                    param(instance)%gdot0_twin*&
                    (abs(tau_twin)/state(instance)%s_twin(j,of))**&
                    param(instance)%n_twin*max(0.0_pReal,sign(1.0_pReal,tau_twin))
    enddo twinSystems1
  enddo twinFamilies1

!--------------------------------------------------------------------------------------------------
! calculate the overall hardening based on above
 j = 0_pInt
 slipFamilies2: do f = 1_pInt,size(param(instance)%Nslip,1)
   slipSystems2: do i = 1_pInt,param(instance)%Nslip(f)
     j = j+1_pInt
     dotState(instance)%s_slip(j,of) = &                                                            ! evolution of slip resistance j
       c_SlipSlip * left_SlipSlip(j) * &
       dot_product(param(instance)%matrix_SlipSlip(j,1:nslip), &
                   right_SlipSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       dot_product(param(instance)%matrix_SlipTwin(j,1:ntwin), &
                   right_SlipTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     dotState(instance)%sumGamma(of) = dotState(instance)%sumGamma(of) + &
                                                        abs(gdot_slip(j))
     dotState(instance)%accshear_slip(j,of) = abs(gdot_slip(j))
   enddo slipSystems2
 enddo slipFamilies2

 j = 0_pInt
 twinFamilies2: do f = 1_pInt,size(param(instance)%Ntwin,1)
   index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   twinSystems2: do i = 1_pInt,param(instance)%Ntwin(f)
     j = j+1_pInt
     dotState(instance)%s_twin(j,of) = &                                                      ! evolution of twin resistance j
       c_TwinSlip * left_TwinSlip(j) * &
       dot_product(param(instance)%matrix_TwinSlip(j,1:nslip), &
                   right_TwinSlip*abs(gdot_slip)) + &                                               ! dot gamma_slip modulated by right-side slip factor
       c_TwinTwin * left_TwinTwin(j) * &
       dot_product(param(instance)%matrix_TwinTwin(j,1:ntwin), &
                   right_TwinTwin*gdot_twin)                                                        ! dot gamma_twin modulated by right-side twin factor
     if (state(instance)%sumF(of) < 0.98_pReal) &                                                   ! ensure twin volume fractions stays below 1.0
       dotState(instance)%sumF(of) = dotState(instance)%sumF(of) + &
                                                      gdot_twin(j)/lattice_shearTwin(index_myFamily+i,ph)
      dotState(instance)%accshear_twin(j,of) = abs(gdot_twin(j))
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
   nSlip,nTwin, &
   o,f,i,c,j,k, &
   index_myFamily
 real(pReal) :: &
   tau_slip_pos,tau_slip_neg,tau

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(ph)

 nSlip= sum(param(instance)%nslip)
 nTwin= sum(param(instance)%nTwin)

 plastic_phenopowerlaw_postResults = 0.0_pReal
 c = 0_pInt

 outputsLoop: do o = 1_pInt,size(param(instance)%outputID)
   select case(param(instance)%outputID(o))
     case (resistance_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(instance)%s_slip(1:nSlip,of)
       c = c + nSlip

     case (accumulatedshear_slip_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nSlip) = state(instance)%accshear_slip(1:nSlip,of)
       c = c + nSlip

     case (shearrate_slip_ID)
       j = 0_pInt
       slipFamilies1: do f = 1_pInt,size(param(instance)%Nslip,1)
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems1: do i = 1_pInt,param(instance)%Nslip(f)
           j = j + 1_pInt
           tau_slip_pos  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
           tau_slip_neg  = tau_slip_pos
           do k = 1,lattice_NnonSchmid(ph)
             tau_slip_pos = tau_slip_pos +param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k,index_myFamily+i,ph))
             tau_slip_neg = tau_slip_neg +param(instance)%nonSchmidCoeff(k)* &
                                   dot_product(Tstar_v,lattice_Sslip_v(1:6,2*k+1,index_myFamily+i,ph))
           enddo
           plastic_phenopowerlaw_postResults(c+j) = param(instance)%gdot0_slip*0.5_pReal* &
                    ((abs(tau_slip_pos)/state(instance)%s_slip(j,of))**param(instance)%n_slip &
                    *sign(1.0_pReal,tau_slip_pos) &
                    +(abs(tau_slip_neg)/(state(instance)%s_slip(j,of)))**param(instance)%n_slip &
                    *sign(1.0_pReal,tau_slip_neg))
         enddo slipSystems1
       enddo slipFamilies1
       c = c + nSlip

     case (resolvedstress_slip_ID)
       j = 0_pInt
       slipFamilies2: do f = 1_pInt,size(param(instance)%Nslip,1)
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems2: do i = 1_pInt,param(instance)%Nslip(f)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
         enddo slipSystems2
       enddo slipFamilies2
       c = c + nSlip

     case (totalshear_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = &
                             state(instance)%sumGamma(of)
       c = c + 1_pInt

     case (resistance_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                            state(instance)%s_twin(1:nTwin,of)
       c = c + nTwin

     case (accumulatedshear_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt:c+nTwin) = &
                             state(instance)%accshear_twin(1:nTwin,of)
       c = c + nTwin
     case (shearrate_twin_ID)
       j = 0_pInt
       twinFamilies1: do f = 1_pInt,size(param(instance)%Ntwin,1)
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems1: do i = 1_pInt,param(instance)%Ntwin(f)
           j = j + 1_pInt
           tau = dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
           plastic_phenopowerlaw_postResults(c+j) = (1.0_pReal-state(instance)%sumF(of))*&  ! 1-F
                                                         param(instance)%gdot0_twin*&
                                                         (abs(tau)/state(instance)%s_twin(j,of))**&
                                           param(instance)%n_twin*max(0.0_pReal,sign(1.0_pReal,tau))
         enddo twinSystems1
       enddo twinFamilies1
       c = c + nTwin

     case (resolvedstress_twin_ID)
       j = 0_pInt
       twinFamilies2: do f = 1_pInt,size(param(instance)%Ntwin,1)
         index_myFamily = sum(lattice_NtwinSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         twinSystems2: do i = 1_pInt,param(instance)%Ntwin(f)
           j = j + 1_pInt
           plastic_phenopowerlaw_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Stwin_v(1:6,index_myFamily+i,ph))
         enddo twinSystems2
       enddo twinFamilies2
       c = c + nTwin

     case (totalvolfrac_twin_ID)
       plastic_phenopowerlaw_postResults(c+1_pInt) = state(instance)%sumF(of)
       c = c + 1_pInt

   end select
 enddo outputsLoop

end function plastic_phenopowerlaw_postResults

end module plastic_phenopowerlaw
