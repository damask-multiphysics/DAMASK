!--------------------------------------------------------------------------------------------------
!> @author Tias Maiti, Michigan State University
!> @author Philip Eisenlohr, Michigan State University
!> @brief material subroutine implementing the diagonal hardening concept outlined by
!> Raul Radovitzky in the paper
!--------------------------------------------------------------------------------------------------

module plastic_MITdiag
 use prec, only: &
   pReal,&
   pInt
 use lattice, only: &
   lattice_maxNinteraction, &
   lattice_maxNslipFamily
   
 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_MITdiag_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_MITdiag_sizePostResult                                                                   !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_MITdiag_output                                                                           !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_MITdiag_Noutput                                                                          !< number of outputs per instance

 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_MITdiag_totalNslip                                                            !< no. of slip system used in simulation

 real(pReal),                         dimension(:,:,:), allocatable,         private :: &
   plastic_MITdiag_hardeningMatrix

 enum, bind(c) 
   enumerator :: undefined_ID, &
                 resistance_ID, &
                 accumulatedshear_ID, &
                 shearrate_ID, &
                 resolvedstress_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), allocatable, dimension(:) :: & 
     outputID
   real(pReal), dimension(lattice_maxNinteraction) :: & 
     interaction
   real(pReal), dimension(:,:), allocatable :: &
     hardeningMatrix
   integer(pInt), dimension(lattice_maxNslipFamily) :: & 
     n_slip          = 0_pInt
   real(pReal),   dimension(lattice_maxNslipFamily) :: & 
     b, &
     g0, &
     rho0, &
     rhosat, &
     gammasat
   real(pReal) :: &
     aTolResistance  = 1.0_pReal, &
     aTolShear       = 1.0e-6_pReal, &
     gdot0, &
     n
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                      !< containers of constitutive parameters (len Ninstance)
 
 type, private :: tMITdiagState                                                                      !< internal state aliases
   real(pReal), pointer,     dimension(:,:) :: &                                                     ! scalars along NipcMyInstance
     resistance, &
     accumulatedShear
 end type
 type, private :: tMITdiagAbsTol                                                                     !< internal alias for abs tolerance in state
   real(pReal), pointer :: &                                                                         ! scalars along NipcMyInstance
     resistance, &
     accumulatedShear
 end type
 type(tMITdiagState), allocatable, dimension(:), private :: &                                        !< state aliases per instance
   state, &
   state0, &
   dotState
 type(tMITdiagAbsTol), allocatable, dimension(:), private :: &                                       !< state aliases per instance
   stateAbsTol

 public  :: &
   plastic_MITdiag_init, &
   plastic_MITdiag_LpAndItsTangent, &
   plastic_MITdiag_dotState, &
   plastic_MITdiag_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_MITdiag_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   numerics_integrator
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_intValue, &
   IO_floatValue, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF, &
   IO_warning
 use lattice
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_MITdiag_label, &
   PLASTICITY_MITdiag_ID, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase
   
 use lattice  

 implicit none
 integer(pInt), intent(in) :: fileUnit
 
 
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   f,i,j,o,k, &
   phase, & 
   instance, &
   maxNinstance, &
   mySize = 0_pInt, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState, &
   totalNslip, &
   index_myFamily, &
   index_otherFamily, &
   Nchunks_SlipSlip = 0_pInt, &
   Nchunks_SlipFamilies = 0_pInt
 character(len=65536) :: &
   tag       = '', &
   line      = '', &
   extmsg    = ''
 real(pReal), dimension(:), allocatable :: tempPerSlip
 character(len=64) :: &
   outputtag = ''
  integer(pInt) :: NipcMyPhase

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_MITdiag_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_MITdiag_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_MITdiag_sizePostResults(maxNinstance),                      source=0_pInt)
 allocate(plastic_MITdiag_sizePostResult(maxval(phase_Noutput), maxNinstance),source=0_pInt)
 allocate(plastic_MITdiag_output(maxval(phase_Noutput), maxNinstance))
          plastic_MITdiag_output = ''
 allocate(plastic_MITdiag_Noutput(maxNinstance),                              source=0_pInt)
 allocate(plastic_MITdiag_totalNslip(maxNinstance),                           source=0_pInt)
 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     phase = phase + 1_pInt                                                                         ! advance section counter
     if (phase_plasticity(phase) == PLASTICITY_MITdiag_ID) then
       instance = phase_plasticityInstance(phase)                                                   ! count instances of my constitutive law
       allocate(param(instance)%outputID(phase_Noutput(phase)))                                     ! allocate space for IDs of every requested output
       Nchunks_SlipFamilies = count(lattice_NslipSystem(:,phase) > 0_pInt)                          ! maximum number of slip families according to lattice type of current phase
       Nchunks_SlipSlip =     maxval(lattice_interactionSlipSlip(:,:,phase))
       if(allocated(tempPerSlip)) deallocate(tempPerSlip)
       allocate(tempPerSlip(Nchunks_SlipFamilies))
     endif
     cycle                                                                                          ! skip to next line
   endif
   
   if (phase > 0_pInt) then; if (phase_plasticity(phase) == PLASTICITY_MITdiag_ID) then           ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     chunkPos = IO_stringPos(line) 
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     extmsg = trim(tag)//' ('//PLASTICITY_MITdiag_label//')'                                      ! prepare error message identifier

     select case(tag)
       case ('(output)')
         outputtag = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
         select case(outputtag)
           case ('resistance')
             plastic_MITdiag_Noutput(instance) = plastic_MITdiag_Noutput(instance) + 1_pInt
             param(instance)%outputID(plastic_MITdiag_Noutput(instance)) = resistance_ID
             plastic_MITdiag_output(plastic_MITdiag_Noutput(instance),instance) = outputtag
           case ('accumulatedshear')
             plastic_MITdiag_Noutput(instance) = plastic_MITdiag_Noutput(instance) + 1_pInt
             param(instance)%outputID(plastic_MITdiag_Noutput(instance)) = accumulatedshear_ID
             plastic_MITdiag_output(plastic_MITdiag_Noutput(instance),instance) = outputtag
           case ('shearrate')
             plastic_MITdiag_Noutput(instance) = plastic_MITdiag_Noutput(instance) + 1_pInt
             param(instance)%outputID(plastic_MITdiag_Noutput(instance)) = shearrate_ID
             plastic_MITdiag_output(plastic_MITdiag_Noutput(instance),instance) = outputtag
           case ('resolvedstress')
             plastic_MITdiag_Noutput(instance) = plastic_MITdiag_Noutput(instance) + 1_pInt
             param(instance)%outputID(plastic_MITdiag_Noutput(instance)) = resolvedstress_ID
             plastic_MITdiag_output(plastic_MITdiag_Noutput(instance),instance) = outputtag

         end select

!--------------------------------------------------------------------------------------------------
! parameters depending on number of slip families
       case ('nslip')
         if (chunkPos(1) < Nchunks_SlipFamilies + 1_pInt) &
           call IO_warning(50_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_MITdiag_label//')')
         if (chunkPos(1) > Nchunks_SlipFamilies + 1_pInt) &
           call IO_error(150_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_MITdiag_label//')')
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt                                                 ! user specified number of (possibly) active slip families (e.g. 6 0 6 --> 3)
         do j = 1_pInt, Nchunks_SlipFamilies
           param(instance)%n_slip(j) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo
       case ('b','g0','rho0','rhosat','gammasat')
         tempPerSlip = 0.0_pReal
         do j = 1_pInt, Nchunks_SlipFamilies
           if (param(instance)%n_slip(j) > 0_pInt) &
             tempPerSlip(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         select case(tag)
           case ('b')
             param(instance)%b(1:Nchunks_SlipFamilies)        = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('g0')
            param(instance)%g0(1:Nchunks_SlipFamilies)        = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('rho0')
             param(instance)%rho0(1:Nchunks_SlipFamilies)     = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('rhosat')
             param(instance)%rhosat(1:Nchunks_SlipFamilies)   = tempPerSlip(1:Nchunks_SlipFamilies)
           case ('gammasat')
             param(instance)%gammasat(1:Nchunks_SlipFamilies) = tempPerSlip(1:Nchunks_SlipFamilies)
         end select

       case ('interaction_slipslip')
         if (chunkPos(1) < 1_pInt + Nchunks_SlipSlip) &
           call IO_warning(52_pInt,ext_msg=trim(tag)//' ('//PLASTICITY_MITDIAG_label//')')
         do j = 1_pInt, Nchunks_SlipSlip
           param(instance)%interaction(j) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo

       case ('gdot0')
         param(instance)%gdot0           = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%gdot0 <= 0.0_pReal)                call IO_error(211_pInt,ext_msg=extmsg)

       case ('n')
         param(instance)%n               = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%n <= 0.0_pReal)                    call IO_error(211_pInt,ext_msg=extmsg)

       case ('atol_resistance')
         param(instance)%aTolResistance  = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%aTolResistance <= 0.0_pReal)       call IO_error(211_pInt,ext_msg=extmsg)

       case ('atol_shear')
         param(instance)%aTolShear       = IO_floatValue(line,chunkPos,2_pInt)

       case default

     end select
   endif; endif
 enddo parsingFile

 allocate(state(maxNinstance))                                                                      ! internal state aliases
 allocate(state0(maxNinstance))
 allocate(dotState(maxNinstance))
 allocate(stateAbsTol(maxNinstance)) 

 initializeInstances: do phase = 1_pInt, size(phase_plasticity)                                     ! loop over every plasticity
   myPhase: if (phase_plasticity(phase) == PLASTICITY_MITdiag_ID) then                              ! isolate instances of own constitutive description
     NipcMyPhase = count(material_phase == phase)                                                   ! number of own material points (including point components ipc)
     instance = phase_plasticityInstance(phase)
!--------------------------------------------------------------------------------------------------
!  sanity checks
     if (param(instance)%aTolShear <= 0.0_pReal) &
       param(instance)%aTolShear = 1.0e-6_pReal                                                     ! default absolute tolerance 1e-6
     param(instance)%n_slip(1:lattice_maxNslipFamily) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase),&                                    ! limit active slip systems per family to min of available and requested
           param(instance)%n_slip(1:lattice_maxNslipFamily))
     plastic_MITdiag_totalNslip(instance)  = sum(param(instance)%n_slip(1:lattice_maxNslipFamily))  ! how many slip systems altogether
     totalNslip = plastic_MITdiag_totalNslip(instance)

     allocate(plastic_MITdiag_hardeningMatrix(maxval(plastic_MITdiag_totalNslip),&                  ! slip resistance from slip activity
                                              maxval(plastic_MITdiag_totalNslip),&
                                              maxNinstance), source=0.0_pReal)
                                              
!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,plastic_MITdiag_Noutput(instance)
       select case(param(instance)%outputID(o))
         case(resistance_ID, &
              accumulatedshear_ID, &
              shearrate_ID, &
              resolvedstress_ID &
              )
           mySize = plastic_MITdiag_totalNslip(instance)
         case default
       end select
  
       outputFound: if (mySize > 0_pInt) then
         plastic_MITdiag_sizePostResult(o,instance) = mySize
         plastic_MITdiag_sizePostResults(instance) = &
         plastic_MITdiag_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState = totalNslip &                                                                        ! g (resistance)
               + totalNslip                                                                          ! accshear
     sizeDotState = sizeState                                                                        ! both evolve
     sizeDeltaState = 0_pInt                                                                         ! no sudden jumps in state
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_MITdiag_sizePostResults(instance)
     plasticState(phase)%nSlip = totalNslip
     plasticState(phase)%nTwin = 0
     plasticState(phase)%nTrans= 0
     allocate(plasticState(phase)%aTolState          (   sizeState))

     allocate(plasticState(phase)%state0             (   sizeState,NipcMyPhase),source=0.0_pReal)

     allocate(plasticState(phase)%partionedState0    (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NipcMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NipcMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase),source=0.0_pReal)

     do f = 1_pInt,lattice_maxNslipFamily                                                                    ! >>> interaction slip -- X
       index_myFamily = sum(param(instance)%n_slip(1:f-1_pInt))
       do j = 1_pInt,param(instance)%n_slip(f)                                            ! loop over (active) systems in my family (slip)
         do o = 1_pInt,lattice_maxNslipFamily
           index_otherFamily = sum(param(instance)%n_slip(1:o-1_pInt))
           do k = 1_pInt,param(instance)%n_slip(o)                                        ! loop over (active) systems in other family (slip)
             plastic_MITdiag_hardeningMatrix(index_myFamily+j, index_otherFamily+k, instance) = &
                 param(instance)%interaction(lattice_interactionSlipSlip( &
                                             sum(lattice_NslipSystem(1:f-1,phase))+j, &
                                             sum(lattice_NslipSystem(1:o-1,phase))+k, &
                                             phase))
           enddo
         enddo
       enddo
     enddo

!--------------------------------------------------------------------------------------------------
! globally required state aliases
     plasticState(phase)%slipRate           => plasticState(phase)%dotState(totalNslip+1:2*totalNslip,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip    => plasticState(phase)%state   (totalNslip+1:2*totalNslip,1:NipcMyPhase)

!--------------------------------------------------------------------------------------------------
! locally defined state aliases
     state(instance)%resistance             => plasticState(phase)%state    (1:totalNslip,1:NipcMyPhase)
     state0(instance)%resistance            => plasticState(phase)%state0   (1:totalNslip,1:NipcMyPhase)
     dotState(instance)%resistance          => plasticState(phase)%dotState (1:totalNslip,1:NipcMyPhase)
     stateAbsTol(instance)%resistance       => plasticState(phase)%aTolState(1)

     state(instance)%accumulatedShear       => plasticState(phase)%state    (totalNslip+1:2*totalNslip,1:NipcMyPhase)
     state0(instance)%accumulatedShear      => plasticState(phase)%state0   (totalNslip+1:2*totalNslip,1:NipcMyPhase)
     dotState(instance)%accumulatedShear    => plasticState(phase)%dotState (totalNslip+1:2*totalNslip,1:NipcMyPhase)
     stateAbsTol(instance)%accumulatedShear => plasticState(phase)%aTolState(2)

!--------------------------------------------------------------------------------------------------
! init state
     state0(instance)%resistance            = param(instance)%g0(1)
     state0(instance)%accumulatedShear      = 1.0e-150_pReal

!--------------------------------------------------------------------------------------------------
! init absolute state tolerances
     stateAbsTol(instance)%resistance       = param(instance)%aTolResistance
     stateAbsTol(instance)%accumulatedShear = param(instance)%aTolShear

   endif myPhase
 enddo initializeInstances

end subroutine plastic_MITdiag_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_MITdiag_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,ipc,ip,el)
 use prec, only: &
   dNeq0
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_g
 use math, only: &
   math_mul6x6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_deviatoric33, &
   math_mul33xx33, &
   math_transpose33
 use lattice, only: &
   lattice_Sslip, &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem
 use material, only: &
   phaseAt, phasememberAt, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) :: &
   dLp_dTstar99                                                                                     !< derivative of Lp with respect to 2nd Piola Kirchhoff stress

 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation

 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: &
   instance, &
   index_myFamily, &
   f,i,j,k,l,m,n, &
   of, &
   ph

 real(pReal) :: &
   tau_slip, &
   gdot_slip, &
   dgdot_dtauslip

 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar3333                                                                                   !< derivative of Lp with respect to Tstar as 4th order tensor

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))

 gdot_slip = 0.0_pReal
 Lp = 0.0_pReal
 dLp_dTstar3333 = 0.0_pReal
 dLp_dTstar99 = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! Slip part
 j = 0_pInt
 slipFamilies: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                          ! at which index starts my family
   slipSystems: do i = 1_pInt,param(instance)%n_slip(f)
     j = j+1_pInt
     
     ! Calculation of Lp
     tau_slip  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))

     if (dNeq0(abs(tau_slip))) then
       gdot_slip = param(instance)%gdot0 * &
                   ((abs(tau_slip)/(state(instance)%resistance(j,of))) &
                    **param(instance)%n)*sign(1.0_pReal,tau_slip)

       Lp = Lp + gdot_slip*lattice_Sslip(1:3,1:3,1,index_myFamily+i,ph)
     end if

     ! Calculation of the tangent of Lp
     if (dNeq0(abs(tau_slip))) then
       dgdot_dtauslip = gdot_slip * param(instance)%n / tau_slip
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLp_dTstar3333(k,l,m,n) = dLp_dTstar3333(k,l,m,n) + &
                                   dgdot_dtauslip * lattice_Sslip(k,l,1,index_myFamily+i,ph)* &
                                                    lattice_Sslip(m,n,1,index_myFamily+i,ph)
     endif

   enddo slipSystems
 enddo slipFamilies
 
 dLp_dTstar99 = math_Plain3333to99(dLp_dTstar3333)

 if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
     .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
          .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
   write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Lp', &
                                 math_transpose33(Lp(1:3,1:3))
   write(6,'(/,a,/,f12.5)') '<< CONST phenopowerlaw >> gdot', gdot_slip
 end if

end subroutine plastic_MITdiag_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_MITdiag_dotState(Tstar_v,ipc,ip,el)
 use prec, only: &
   dNeq0
 use math, only: &
   pi
 use material, only: &
   phaseAt, phasememberAt, &
   material_phase, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_mu, &
   lattice_Sslip_v, &
   lattice_maxNslipFamily, &
   lattice_NslipSystem
 
 implicit none
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 integer(pInt) :: &
   instance,ph, &
   nSlip, &
   f,i,j, &
   index_myFamily, &
   of
 real(pReal) :: &
   tau_slip, &
   forresthardening

 real(pReal), dimension(plastic_MITdiag_totalNslip(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
   gdot_slip, &
   tau_c, &
   gamma_c, &
   h, rho

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 ph = phaseAt(ipc,ip,el)
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 nSlip = plastic_MITdiag_totalNslip(instance)
 
!--------------------------------------------------------------------------------------------------
!  calculate left and right vectors and calculate dot gammas
 gdot_slip = 0.0_pReal
 tau_slip  = 0.0_pReal
 rho       = 0.0_pReal
 
 j = 0_pInt
 slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems1: do i = 1_pInt,param(instance)%n_slip(f)
     j = j+1_pInt

     tau_slip  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))

     if (dNeq0(abs(tau_slip))) then
        gdot_slip(j) = param(instance)%gdot0 * &
                       ((abs(tau_slip)/(state(instance)%resistance(j,of))) &
                       **param(instance)%n * sign(1.0_pReal,tau_slip))
     
          
     endif
     rho(j) = param(instance)%rhosat(f)* &
          (1.0_pReal - (1.0_pReal - (param(instance)%rho0(f) / param(instance)%rhosat(f))) * &
          exp(-state(instance)%accumulatedShear(j,of) / param(instance)%gammasat(f)))

     dotState(ph)%accumulatedShear(j,of) = abs(gdot_slip(j))

   enddo slipSystems1
 enddo slipFamilies1

 tau_c     = 0.0_pReal
 gamma_c   = 0.0_pReal
 h         = 0.0_pReal

  
 j = 0_pInt
 slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                         ! at which index starts my family
   slipSystems2: do i = 1_pInt,param(instance)%n_slip(f)
     j = j+1_pInt
     
     forresthardening = &
                        dot_product(plastic_MITdiag_hardeningMatrix(j,1:nSlip,instance), rho)     
     tau_c(j) = 0.3_pReal * lattice_mu(ph) * param(instance)%b(f) * &
                sqrt(pi * forresthardening)                         
     gamma_c(j) = param(instance)%b(f) * rho(j) / &
                  (2.0_pReal * sqrt(forresthardening))     
     h(j)  = (tau_c(j) / gamma_c(j)) * (state(instance)%resistance(j,of) / tau_c(j))**3.0_pReal * &
                  (cosh((tau_c(j) / state(instance)%resistance(j,of))**2.0_pReal) - 1.0_pReal)
     dotState(ph)%resistance(j,of) = h(j) * abs(gdot_slip(j))
     
   enddo slipSystems2
 enddo slipFamilies2
   
end subroutine plastic_MITdiag_dotState

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_MITdiag_postResults(Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6
 use material, only: &
   material_phase, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance
 use lattice, only: &
   lattice_NslipSystem, &
   lattice_maxNslipFamily, &
   lattice_Sslip_v
   
 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plastic_MITdiag_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           plastic_MITdiag_postResults

 integer(pInt) :: &
   instance,ph, of, &
   nSlip, &
   o,f,i,c,j, &
   index_myFamily
 real(pReal) :: &
   tau_slip

 of = phasememberAt(ipc,ip,el)
 ph = phaseAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 
 nSlip = plastic_MITdiag_totalNslip(instance) 
 
 c = 0_pInt
 plastic_MITdiag_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,plastic_MITdiag_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (resistance_ID)
       plastic_MITdiag_postResults(c+1_pInt:c+nSlip) = state(instance)%resistance(:,of)
       c = c + nSlip

     case (accumulatedshear_ID)
       plastic_MITdiag_postResults(c+1_pInt:c+nSlip) = state(instance)%accumulatedShear(:,of)
       c = c + nSlip

     case (shearrate_ID)
       j = 0_pInt
       slipFamilies1: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems1: do i = 1_pInt,param(instance)%n_slip(f)
           j = j + 1_pInt
           tau_slip  = dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
           plastic_MITdiag_postResults(c+j) = param(instance)%gdot0 * &
                    ((abs(tau_slip) / state(instance)%resistance(j,of)) ** param(instance)%n &
                    *sign(1.0_pReal,tau_slip))
         enddo slipSystems1
       enddo slipFamilies1
       c = c + nSlip

     case (resolvedstress_ID)
       j = 0_pInt
       slipFamilies2: do f = 1_pInt,lattice_maxNslipFamily
         index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,ph))                                ! at which index starts my family
         slipSystems2: do i = 1_pInt,param(instance)%n_slip(f)
           j = j + 1_pInt
           plastic_MITdiag_postResults(c+j) = &
                             dot_product(Tstar_v,lattice_Sslip_v(1:6,1,index_myFamily+i,ph))
         enddo slipSystems2
       enddo slipFamilies2
       c = c + nSlip

   end select
 enddo outputsLoop

end function plastic_MITdiag_postResults


end module plastic_MITdiag
