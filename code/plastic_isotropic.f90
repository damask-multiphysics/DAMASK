!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic (ISOTROPIC) plasticity
!> @details Isotropic (ISOTROPIC) Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
module plastic_isotropic
#ifdef HDF
 use hdf5, only: &
   HID_T
#endif

 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
 integer(pInt),                       dimension(:),     allocatable,         public, protected :: &
   plastic_isotropic_sizePostResults                                                                  !< cumulative size of post results
   
 integer(pInt),                       dimension(:,:),   allocatable, target, public :: &
   plastic_isotropic_sizePostResult                                                                   !< size of each post result output
   
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_isotropic_output                                                                           !< name of each post result output
 
 integer(pInt),                       dimension(:),     allocatable, target, public :: &
   plastic_isotropic_Noutput                                                                          !< number of outputs per instance
 
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 flowstress_ID, &
                 strainrate_ID
 end enum

 type, private :: tParameters                                                                         !< container type for internal constitutive parameters
   integer(kind(undefined_ID)), allocatable, dimension(:) :: & 
     outputID
  real(pReal) :: &
     fTaylor, &
     tau0, &
     gdot0, &
     n, &
     h0, &
     h0_slopeLnRate, &
     tausat, &
     a, &
     aTolFlowstress, &
     aTolShear     , &
     tausat_SinhFitA, &
     tausat_SinhFitB, &
     tausat_SinhFitC, &
     tausat_SinhFitD
  logical :: &
     dilatation
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                       !< containers of constitutive parameters (len Ninstance)
 
 type, private :: tIsotropicState                                                                     !< internal state aliases
   real(pReal), pointer,     dimension(:) :: &                                                        ! scalars along NipcMyInstance
     flowstress, &
     accumulatedShear
 end type
 type, private :: tIsotropicAbsTol                                                                  !< internal alias for abs tolerance in state
   real(pReal), pointer :: &                                                                        ! scalars along NipcMyInstance
     flowstress, &
     accumulatedShear
 end type
 type(tIsotropicState), allocatable, dimension(:), private :: &                                       !< state aliases per instance
   state, &
   state0, &
   dotState
 type(tIsotropicAbsTol), allocatable, dimension(:), private :: &                                       !< state aliases per instance
   stateAbsTol

 public  :: &
   plastic_isotropic_init, &
   plastic_isotropic_LpAndItsTangent, &
   plastic_isotropic_LiAndItsTangent, &
   plastic_isotropic_dotState, &
   plastic_isotropic_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   analyticJaco, &
   worldrank, &
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
   IO_floatValue, &
   IO_error, &
   IO_timeStamp, &
   IO_EOF
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_ISOTROPIC_ID, &
   material_phase, &
   plasticState, &
   MATERIAL_partPhase
   
 use lattice  

 implicit none
 integer(pInt), intent(in) :: fileUnit
 
 
 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: &
   o, &
   phase, & 
   instance, &
   maxNinstance, &
   mySize, &
   sizeDotState, &
   sizeState, &
   sizeDeltaState
 character(len=65536) :: &
   tag       = '', &
   outputtag = '', &
   line      = '', &
   extmsg    = ''
  integer(pInt) :: NipcMyPhase

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_ISOTROPIC_label//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_ISOTROPIC_ID),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

 allocate(plastic_isotropic_sizePostResults(maxNinstance),                      source=0_pInt)
 allocate(plastic_isotropic_sizePostResult(maxval(phase_Noutput), maxNinstance),source=0_pInt)
 allocate(plastic_isotropic_output(maxval(phase_Noutput), maxNinstance))
          plastic_isotropic_output = ''
 allocate(plastic_isotropic_Noutput(maxNinstance),                              source=0_pInt)

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
     if (phase_plasticity(phase) == PLASTICITY_ISOTROPIC_ID) then
       instance = phase_plasticityInstance(phase)

     endif
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt) then; if (phase_plasticity(phase) == PLASTICITY_ISOTROPIC_ID) then           ! one of my phases. Do not short-circuit here (.and. between if-statements), it's not safe in Fortran
     instance = phase_plasticityInstance(phase)                                                     ! which instance of my plasticity is present phase
     allocate(param(instance)%outputID(phase_Noutput(phase)))                                       ! allocate space for IDs of every requested output
     chunkPos = IO_stringPos(line) 
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     extmsg = trim(tag)//' ('//PLASTICITY_ISOTROPIC_label//')'                                      ! prepare error message identifier

     select case(tag)
       case ('(output)')
         outputtag = IO_lc(IO_stringValue(line,chunkPos,2_pInt))
         select case(outputtag)
           case ('flowstress')
             plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
             param(instance)%outputID (plastic_isotropic_Noutput(instance)) = flowstress_ID
             plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = outputtag
           case ('strainrate')
             plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
             param(instance)%outputID (plastic_isotropic_Noutput(instance)) = strainrate_ID
             plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = outputtag

         end select

       case ('/dilatation/')
         param(instance)%dilatation      = .true.

       case ('tau0')
         param(instance)%tau0            = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%tau0 < 0.0_pReal)                  call IO_error(211_pInt,ext_msg=extmsg)

       case ('gdot0')
         param(instance)%gdot0           = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%gdot0 <= 0.0_pReal)                call IO_error(211_pInt,ext_msg=extmsg)

       case ('n')
         param(instance)%n               = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%n <= 0.0_pReal)                    call IO_error(211_pInt,ext_msg=extmsg)

       case ('h0')
         param(instance)%h0              = IO_floatValue(line,chunkPos,2_pInt)

       case ('h0_slope','slopelnrate')
         param(instance)%h0_slopeLnRate  = IO_floatValue(line,chunkPos,2_pInt)

       case ('tausat')
         param(instance)%tausat          = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%tausat <= 0.0_pReal)               call IO_error(211_pInt,ext_msg=extmsg)

       case ('tausat_sinhfita')
         param(instance)%tausat_SinhFitA = IO_floatValue(line,chunkPos,2_pInt)

       case ('tausat_sinhfitb')
         param(instance)%tausat_SinhFitB = IO_floatValue(line,chunkPos,2_pInt)

       case ('tausat_sinhfitc')
         param(instance)%tausat_SinhFitC = IO_floatValue(line,chunkPos,2_pInt)

       case ('tausat_sinhfitd')
         param(instance)%tausat_SinhFitD = IO_floatValue(line,chunkPos,2_pInt)

       case ('a', 'w0')
         param(instance)%a               = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%a <= 0.0_pReal)                    call IO_error(211_pInt,ext_msg=extmsg)

       case ('taylorfactor')
         param(instance)%fTaylor         = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%fTaylor <= 0.0_pReal)              call IO_error(211_pInt,ext_msg=extmsg)

       case ('atol_flowstress')
         param(instance)%aTolFlowstress  = IO_floatValue(line,chunkPos,2_pInt)
         if (param(instance)%aTolFlowstress <= 0.0_pReal)       call IO_error(211_pInt,ext_msg=extmsg)

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
   myPhase: if (phase_plasticity(phase) == PLASTICITY_isotropic_ID) then                            ! isolate instances of own constitutive description
     NipcMyPhase = count(material_phase == phase)                                                   ! number of own material points (including point components ipc)
     instance = phase_plasticityInstance(phase)
!--------------------------------------------------------------------------------------------------
!  sanity checks
     if (param(instance)%aTolShear <= 0.0_pReal) &
       param(instance)%aTolShear = 1.0e-6_pReal                                         ! default absolute tolerance 1e-6

!--------------------------------------------------------------------------------------------------
!  Determine size of postResults array
     outputsLoop: do o = 1_pInt,plastic_isotropic_Noutput(instance)
       select case(param(instance)%outputID(o))
         case(flowstress_ID,strainrate_ID)
           mySize = 1_pInt
         case default
       end select
  
       outputFound: if (mySize > 0_pInt) then
         plastic_isotropic_sizePostResult(o,instance) = mySize
         plastic_isotropic_sizePostResults(instance) = &
         plastic_isotropic_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

!--------------------------------------------------------------------------------------------------
! allocate state arrays
     sizeState    = 2_pInt                                                                           ! flowstress, accumulated_shear
     sizeDotState = sizeState                                                                        ! both evolve
     sizeDeltaState = 0_pInt                                                                         ! no sudden jumps in state
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%sizePostResults = plastic_isotropic_sizePostResults(instance)
     plasticState(phase)%nSlip = 1
     plasticState(phase)%nTwin = 0
     plasticState(phase)%nTrans= 0
     allocate(plasticState(phase)%aTolState          (   sizeState))

     allocate(plasticState(phase)%state0             (   sizeState,NipcMyPhase),source=0.0_pReal)

     allocate(plasticState(phase)%partionedState0    (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%subState0          (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%state              (   sizeState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%dotState           (sizeDotState,NipcMyPhase),source=0.0_pReal)
     allocate(plasticState(phase)%deltaState       (sizeDeltaState,NipcMyPhase),source=0.0_pReal)
     if (.not. analyticJaco) then
       allocate(plasticState(phase)%state_backup     (   sizeState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%dotState_backup  (sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 1_pInt)) then
       allocate(plasticState(phase)%previousDotState (sizeDotState,NipcMyPhase),source=0.0_pReal)
       allocate(plasticState(phase)%previousDotState2(sizeDotState,NipcMyPhase),source=0.0_pReal)
     endif
     if (any(numerics_integrator == 4_pInt)) &
       allocate(plasticState(phase)%RK4dotState      (sizeDotState,NipcMyPhase),source=0.0_pReal)
     if (any(numerics_integrator == 5_pInt)) &
       allocate(plasticState(phase)%RKCK45dotState (6,sizeDotState,NipcMyPhase),source=0.0_pReal)

!--------------------------------------------------------------------------------------------------
! globally required state aliases
     plasticState(phase)%slipRate           => plasticState(phase)%dotState(2:2,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip    => plasticState(phase)%state   (2:2,1:NipcMyPhase)

!--------------------------------------------------------------------------------------------------
! locally defined state aliases
     state(instance)%flowstress             => plasticState(phase)%state    (1,1:NipcMyPhase)
     state0(instance)%flowstress            => plasticState(phase)%state0   (1,1:NipcMyPhase)
     dotState(instance)%flowstress          => plasticState(phase)%dotState (1,1:NipcMyPhase)
     stateAbsTol(instance)%flowstress       => plasticState(phase)%aTolState(1)

     state(instance)%accumulatedShear       => plasticState(phase)%state    (2,1:NipcMyPhase)
     state0(instance)%accumulatedShear      => plasticState(phase)%state0   (2,1:NipcMyPhase)
     dotState(instance)%accumulatedShear    => plasticState(phase)%dotState (2,1:NipcMyPhase)
     stateAbsTol(instance)%accumulatedShear => plasticState(phase)%aTolState(2)

!--------------------------------------------------------------------------------------------------
! init state
     state0(instance)%flowstress       = param(instance)%tau0
     state0(instance)%accumulatedShear = 0.0_pReal

!--------------------------------------------------------------------------------------------------
! init absolute state tolerances
     stateAbsTol(instance)%flowstress       = param(instance)%aTolFlowstress
     stateAbsTol(instance)%accumulatedShear = param(instance)%aTolShear

   endif myPhase
 enddo initializeInstances

end subroutine plastic_isotropic_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LpAndItsTangent(Lp,dLp_dTstar99,Tstar_v,ipc,ip,el)
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
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
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

 real(pReal), dimension(3,3) :: &
   Tstar_dev_33                                                                                     !< deviatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal), dimension(3,3,3,3) :: &
   dLp_dTstar_3333                                                                                  !< derivative of Lp with respect to Tstar as 4th order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_dev, &                                                                                !< euclidean norm of Tstar_dev
   squarenorm_Tstar_dev                                                                             !< square of the euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, of, &
   k, l, m, n

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(phaseAt(ipc,ip,el))                                            ! "phaseAt" equivalent to "material_phase" !!

 Tstar_dev_33 = math_deviatoric33(math_Mandel6to33(Tstar_v))                                        ! deviatoric part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_dev = math_mul33xx33(Tstar_dev_33,Tstar_dev_33)
 norm_Tstar_dev = sqrt(squarenorm_Tstar_dev) 

 if (norm_Tstar_dev <= 0.0_pReal) then                                                              ! Tstar == 0 --> both Lp and dLp_dTstar are zero
   Lp = 0.0_pReal
   dLp_dTstar99 = 0.0_pReal
 else
   gamma_dot = param(instance)%gdot0 &
             * ( sqrt(1.5_pReal) * norm_Tstar_dev / param(instance)%fTaylor / state(instance)%flowstress(of) ) &
             **param(instance)%n

   Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/param(instance)%fTaylor

   if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
       .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CONST isotropic >> at el ip g ',el,ip,ipc
     write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Tstar (dev) / MPa', &
                                      math_transpose33(Tstar_dev_33(1:3,1:3))*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> norm Tstar / MPa', norm_Tstar_dev*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> gdot', gamma_dot
   end if
!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,m,n) = (param(instance)%n-1.0_pReal) * &
                                      Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
   forall (k=1_pInt:3_pInt,m=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,k,m,m) = dLp_dTstar_3333(k,k,m,m) - 1.0_pReal/3.0_pReal
   dLp_dTstar99 = math_Plain3333to99(gamma_dot / param(instance)%fTaylor * &
                                      dLp_dTstar_3333 / norm_Tstar_dev)
 end if
end subroutine plastic_isotropic_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dTstar_3333,Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6, &
   math_Mandel6to33, &
   math_Plain3333to99, &
   math_spherical33, &
   math_mul33xx33
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Li                                                                                               !< plastic velocity gradient

 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 real(pReal), dimension(3,3) :: &
   Tstar_sph_33                                                                                     !< sphiatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal), dimension(3,3,3,3), intent(out)  :: &
   dLi_dTstar_3333                                                                                  !< derivative of Li with respect to Tstar as 4th order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_sph, &                                                                                !< euclidean norm of Tstar_sph
   squarenorm_Tstar_sph                                                                             !< square of the euclidean norm of Tstar_sph
 integer(pInt) :: &
   instance, of, &
   k, l, m, n

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(phaseAt(ipc,ip,el))                                            ! "phaseAt" equivalent to "material_phase" !!

 Tstar_sph_33 = math_spherical33(math_Mandel6to33(Tstar_v))                                         ! spherical part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_sph = math_mul33xx33(Tstar_sph_33,Tstar_sph_33)
 norm_Tstar_sph = sqrt(squarenorm_Tstar_sph) 

 if (param(instance)%dilatation) then
     if (norm_Tstar_sph <= 0.0_pReal) then                                                          ! Tstar == 0 --> both Li and dLi_dTstar are zero
       Li = 0.0_pReal
       dLi_dTstar_3333 = 0.0_pReal
     else
       gamma_dot = param(instance)%gdot0 &
                   * (sqrt(1.5_pReal) * norm_Tstar_sph / param(instance)%fTaylor / state(instance)%flowstress(of) ) &
                   **param(instance)%n

       Li = Tstar_sph_33/norm_Tstar_sph * gamma_dot/param(instance)%fTaylor

       !--------------------------------------------------------------------------------------------------
       ! Calculation of the tangent of Li
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLi_dTstar_3333(k,l,m,n) = (param(instance)%n-1.0_pReal) * &
                                          Tstar_sph_33(k,l)*Tstar_sph_33(m,n) / squarenorm_Tstar_sph
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
         dLi_dTstar_3333(k,l,k,l) = dLi_dTstar_3333(k,l,k,l) + 1.0_pReal

       dLi_dTstar_3333 = gamma_dot / param(instance)%fTaylor * &
                                          dLi_dTstar_3333 / norm_Tstar_sph
     endif
 endif
 
end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_dotState(Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6
 use material, only: &
   phaseAt, phasememberAt, &
   plasticState, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   hardening, &                                                                                     !< hardening coefficient
   saturation, &                                                                                    !< saturation flowstress
   norm_Tstar_v                                                                                     !< euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of                                                                                               !< shortcut notation for offset position in state array

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(phaseAt(ipc,ip,el))                                            ! "phaseAt" equivalent to "material_phase" !!

!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (param(instance)%dilatation) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
!--------------------------------------------------------------------------------------------------
! strain rate 
 gamma_dot = param(instance)%gdot0 * ( sqrt(1.5_pReal) * norm_Tstar_v & 
            / &!-----------------------------------------------------------------------------------
           (param(instance)%fTaylor*state(instance)%flowstress(of) ))**param(instance)%n
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (abs(param(instance)%tausat_SinhFitA) <= tiny(0.0_pReal)) then
     saturation = param(instance)%tausat
   else
     saturation = (  param(instance)%tausat &
                   + ( log(  ( gamma_dot / param(instance)%tausat_SinhFitA&
                               )**(1.0_pReal / param(instance)%tausat_SinhFitD)&
                            + sqrt(  ( gamma_dot / param(instance)%tausat_SinhFitA &
                                      )**(2.0_pReal / param(instance)%tausat_SinhFitD) &
                                   + 1.0_pReal ) &
                            ) & ! asinh(K) = ln(K + sqrt(K^2 +1))
                       )**(1.0_pReal / param(instance)%tausat_SinhFitC) &
                   / (  param(instance)%tausat_SinhFitB &
                      * (gamma_dot / param(instance)%gdot0)**(1.0_pReal / param(instance)%n) &
                      ) &
                   )
   endif
   hardening = ( param(instance)%h0 + param(instance)%h0_slopeLnRate * log(gamma_dot) ) &
               * abs( 1.0_pReal - state(instance)%flowstress(of)/saturation )**param(instance)%a &
               * sign(1.0_pReal, 1.0_pReal - state(instance)%flowstress(of)/saturation)
 else
   hardening = 0.0_pReal
 endif

 dotState(instance)%flowstress      (of) = hardening * gamma_dot
 dotState(instance)%accumulatedShear(of) =             gamma_dot

end subroutine plastic_isotropic_dotState

!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_isotropic_postResults(Tstar_v,ipc,ip,el)
 use math, only: &
   math_mul6x6
 use material, only: &
   material_phase, &
   plasticState, &
   phaseAt, phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal), dimension(plastic_isotropic_sizePostResults(phase_plasticityInstance(material_phase(ipc,ip,el)))) :: &
                                           plastic_isotropic_postResults

 real(pReal), dimension(6) :: &
   Tstar_dev_v                                                                                      !< deviatoric 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal) :: &
   norm_Tstar_v                                                                                     ! euclidean norm of Tstar_dev
 integer(pInt) :: &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of, &                                                                                            !< shortcut notation for offset position in state array
   c, &
   o

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(phaseAt(ipc,ip,el))                                            ! "phaseAt" equivalent to "material_phase" !!
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (param(instance)%dilatation) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
 
 c = 0_pInt
 plastic_isotropic_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,plastic_isotropic_Noutput(instance)
   select case(param(instance)%outputID(o))
     case (flowstress_ID)
       plastic_isotropic_postResults(c+1_pInt) = state(instance)%flowstress(of)
       c = c + 1_pInt
     case (strainrate_ID)
       plastic_isotropic_postResults(c+1_pInt) = &
                param(instance)%gdot0 * (            sqrt(1.5_pReal) * norm_Tstar_v & 
             / &!----------------------------------------------------------------------------------
              (param(instance)%fTaylor * state(instance)%flowstress(of)) ) ** param(instance)%n
       c = c + 1_pInt
   end select
 enddo outputsLoop

end function plastic_isotropic_postResults


end module plastic_isotropic
