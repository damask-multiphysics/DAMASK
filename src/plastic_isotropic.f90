!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic (ISOTROPIC) plasticity
!> @details Isotropic (ISOTROPIC) Plasticity which resembles the phenopowerlaw plasticity without
!! resolving the stress on the slip systems. Will give the response of phenopowerlaw for an
!! untextured polycrystal
!--------------------------------------------------------------------------------------------------
module plastic_isotropic
 use prec, only: &
   pReal,&
   pInt
 
 implicit none
 private
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
     aTolShear, &
     tausat_SinhFitA, &
     tausat_SinhFitB, &
     tausat_SinhFitC, &
     tausat_SinhFitD
  logical :: &
     dilatation
 end type

 type(tParameters), dimension(:), allocatable, target, private :: param                               !< containers of constitutive parameters (len Ninstance)
 
 type, private :: tIsotropicState                                                                     !< internal state aliases
   real(pReal), pointer,     dimension(:) :: &                                                        ! scalars along NipcMyInstance
     flowstress, &
     accumulatedShear
 end type

 type(tIsotropicState), allocatable, dimension(:), private :: &                                       !< state aliases per instance
   state, &
   dotState

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
subroutine plastic_isotropic_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
use IO
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   numerics_integrator
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_ISOTROPIC_ID, &
   material_phase, &
   plasticState
 use config_material, only: &
   MATERIAL_partPhase, &
   phaseConfig
   
 use lattice  

 implicit none
 
 type(tParameters), pointer :: p
 
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
   extmsg    = ''
 integer(pInt) :: NipcMyPhase,i
 character(len=64), dimension(:), allocatable :: outputs

 write(6,'(/,a)')   ' <<<+-  constitutive_'//PLASTICITY_ISOTROPIC_label//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == PLASTICITY_ISOTROPIC_ID),pInt)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance

! public variables
 allocate(plastic_isotropic_sizePostResult(maxval(phase_Noutput), maxNinstance),source=0_pInt)
 allocate(plastic_isotropic_output(maxval(phase_Noutput), maxNinstance))
          plastic_isotropic_output = ''
 allocate(plastic_isotropic_Noutput(maxNinstance),                              source=0_pInt)

! inernal variable 
 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
 allocate(state(maxNinstance))                                                                      ! internal state aliases
 allocate(dotState(maxNinstance))

 do phase = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(phase) == PLASTICITY_ISOTROPIC_ID) then
     instance = phase_plasticityInstance(phase)
     p => param(instance)                                                                           ! shorthand pointer to parameter object of my constitutive law
     p%tau0            =  phaseConfig(phase)%getFloat('tau0')
     p%tausat          =  phaseConfig(phase)%getFloat('tausat')
     p%gdot0           =  phaseConfig(phase)%getFloat('gdot0')
     p%n               =  phaseConfig(phase)%getFloat('n')
     p%h0              =  phaseConfig(phase)%getFloat('h0')
     p%fTaylor         =  phaseConfig(phase)%getFloat('m')
     p%h0_slopeLnRate  =  phaseConfig(phase)%getFloat('h0_slopelnrate', defaultVal=0.0_pReal)          ! ToDo: alias allowed?
     p%tausat_SinhFitA =  phaseConfig(phase)%getFloat('tausat_sinhfita',defaultVal=0.0_pReal)
     p%tausat_SinhFitB =  phaseConfig(phase)%getFloat('tausat_sinhfitb',defaultVal=0.0_pReal)
     p%tausat_SinhFitC =  phaseConfig(phase)%getFloat('tausat_sinhfitc',defaultVal=0.0_pReal)
     p%tausat_SinhFitD =  phaseConfig(phase)%getFloat('tausat_sinhfitd',defaultVal=0.0_pReal)
     p%a               =  phaseConfig(phase)%getFloat('a')                                          ! ToDo: alias
     p%aTolFlowStress =  phaseConfig(phase)%getFloat('atol_flowstress',defaultVal=1.0_pReal)
     p%aTolShear =  phaseConfig(phase)%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
    
     p%dilatation      = phaseConfig(phase)%keyExists('/dilatation/')

     outputs = phaseConfig(phase)%getStrings('(output)')
     allocate(p%outputID(0))
     do i=1_pInt, size(outputs)
       select case(outputs(i))
         case ('flowstress')
           plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
           plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = outputs(i)
           plasticState(phase)%sizePostResults =  plasticState(phase)%sizePostResults + 1_pInt
           plastic_isotropic_sizePostResult(i,instance) = 1_pInt
           p%outputID = [p%outputID,flowstress_ID]
         case ('strainrate')
           plastic_isotropic_Noutput(instance) = plastic_isotropic_Noutput(instance) + 1_pInt
           plastic_isotropic_output(plastic_isotropic_Noutput(instance),instance) = outputs(i)
           plasticState(phase)%sizePostResults = &
             plasticState(phase)%sizePostResults + 1_pInt
           plastic_isotropic_sizePostResult(i,instance) = 1_pInt
           p%outputID = [p%outputID,strainrate_ID]
       end select
     enddo

!--------------------------------------------------------------------------------------------------
!  sanity checks
     extmsg = ''
     if (p%aTolShear        <= 0.0_pReal) extmsg = trim(extmsg)//"'aTolShear' "
     if (p%tau0              < 0.0_pReal) extmsg = trim(extmsg)//"'tau0' "
     if (p%gdot0            <= 0.0_pReal) extmsg = trim(extmsg)//"'gdot0' "
     if (p%n                <= 0.0_pReal) extmsg = trim(extmsg)//"'n' "
     if (p%tausat           <= p%tau0)    extmsg = trim(extmsg)//"'tausat' "
     if (p%a                <= 0.0_pReal) extmsg = trim(extmsg)//"'a' " 
     if (p%fTaylor          <= 0.0_pReal) extmsg = trim(extmsg)//"'m' "
     if (p%aTolFlowstress   <= 0.0_pReal) extmsg = trim(extmsg)//"'atol_flowstress' "
     if (extmsg /= '') call IO_error(211_pInt,ip=instance,&
                            ext_msg=trim(extmsg)//'('//PLASTICITY_ISOTROPIC_label//')')

!--------------------------------------------------------------------------------------------------
! allocate state arrays
     NipcMyPhase = count(material_phase == phase)                                                   ! number of own material points (including point components ipc)

     sizeDotState   = size(["flowstress       ","accumulated_shear"])
     sizeDeltaState = 0_pInt                                                                         ! no sudden jumps in state
     sizeState      = sizeDotState + sizeDeltaState
     plasticState(phase)%sizeState = sizeState
     plasticState(phase)%sizeDotState = sizeDotState
     plasticState(phase)%sizeDeltaState = sizeDeltaState
     plasticState(phase)%nSlip = 1
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

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState

     state(instance)%flowstress             => plasticState(phase)%state    (1,1:NipcMyPhase)
     dotState(instance)%flowstress          => plasticState(phase)%dotState (1,1:NipcMyPhase)
     plasticState(phase)%state0(1,1:NipcMyPhase) = p%tau0
     plasticState(phase)%aTolState(1)       =  p%aTolFlowstress

     state(instance)%accumulatedShear       => plasticState(phase)%state    (2,1:NipcMyPhase)
     dotState(instance)%accumulatedShear    => plasticState(phase)%dotState (2,1:NipcMyPhase)
     plasticState(phase)%state0 (2,1:NipcMyPhase) = 0.0_pReal
     plasticState(phase)%aTolState(2)       =  p%aTolShear
     ! global alias
     plasticState(phase)%slipRate           => plasticState(phase)%dotState(2:2,1:NipcMyPhase)
     plasticState(phase)%accumulatedSlip    => plasticState(phase)%state   (2:2,1:NipcMyPhase)

endif
 enddo

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
   math_mul33xx33
 use material, only: &
   phasememberAt, &
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

 type(tParameters), pointer :: p
 
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
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 p => param(instance)
 
 Tstar_dev_33 = math_deviatoric33(math_Mandel6to33(Tstar_v))                                        ! deviatoric part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_dev = math_mul33xx33(Tstar_dev_33,Tstar_dev_33)
 norm_Tstar_dev = sqrt(squarenorm_Tstar_dev) 

 if (norm_Tstar_dev <= 0.0_pReal) then                                                              ! Tstar == 0 --> both Lp and dLp_dTstar are zero
   Lp = 0.0_pReal
   dLp_dTstar99 = 0.0_pReal
 else
   gamma_dot = p%gdot0 &
             * ( sqrt(1.5_pReal) * norm_Tstar_dev / p%fTaylor / state(instance)%flowstress(of) ) &
             **p%n

   Lp = Tstar_dev_33/norm_Tstar_dev * gamma_dot/p%fTaylor 

   if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
       .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CONST isotropic >> at el ip g ',el,ip,ipc
     write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Tstar (dev) / MPa', &
                                      transpose(Tstar_dev_33(1:3,1:3))*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> norm Tstar / MPa', norm_Tstar_dev*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> gdot', gamma_dot
   end if
!--------------------------------------------------------------------------------------------------
! Calculation of the tangent of Lp
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,m,n) = (p%n-1.0_pReal) * &
                                      Tstar_dev_33(k,l)*Tstar_dev_33(m,n) / squarenorm_Tstar_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,l,k,l) = dLp_dTstar_3333(k,l,k,l) + 1.0_pReal
   forall (k=1_pInt:3_pInt,m=1_pInt:3_pInt) &
     dLp_dTstar_3333(k,k,m,m) = dLp_dTstar_3333(k,k,m,m) - 1.0_pReal/3.0_pReal
   dLp_dTstar99 = math_Plain3333to99(gamma_dot / p%fTaylor * &
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
   phasememberAt, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Li                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out)  :: &
   dLi_dTstar_3333                                                                                  !< derivative of Li with respect to Tstar as 4th order tensor
 real(pReal), dimension(6),   intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 type(tParameters), pointer :: p
 
 real(pReal), dimension(3,3) :: &
   Tstar_sph_33                                                                                     !< sphiatoric part of the 2nd Piola Kirchhoff stress tensor as 2nd order tensor
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_sph, &                                                                                !< euclidean norm of Tstar_sph
   squarenorm_Tstar_sph                                                                             !< square of the euclidean norm of Tstar_sph
 integer(pInt) :: &
   instance, of, &
   k, l, m, n

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 p => param(instance)
 
 Tstar_sph_33 = math_spherical33(math_Mandel6to33(Tstar_v))                                         ! spherical part of 2nd Piola-Kirchhoff stress
 squarenorm_Tstar_sph = math_mul33xx33(Tstar_sph_33,Tstar_sph_33)
 norm_Tstar_sph = sqrt(squarenorm_Tstar_sph) 

 if (p%dilatation .and. norm_Tstar_sph > 0.0_pReal) then                              ! Tstar == 0 or J2 plascitiy --> both Li and dLi_dTstar are zero
   gamma_dot = p%gdot0 &
               * (sqrt(1.5_pReal) * norm_Tstar_sph / p%fTaylor / state(instance)%flowstress(of) ) &
               **p%n

   Li = Tstar_sph_33/norm_Tstar_sph * gamma_dot/p%fTaylor

   !--------------------------------------------------------------------------------------------------
   ! Calculation of the tangent of Li
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLi_dTstar_3333(k,l,m,n) = (p%n-1.0_pReal) * &
                                      Tstar_sph_33(k,l)*Tstar_sph_33(m,n) / squarenorm_Tstar_sph
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLi_dTstar_3333(k,l,k,l) = dLi_dTstar_3333(k,l,k,l) + 1.0_pReal

   dLi_dTstar_3333 = gamma_dot / p%fTaylor * &
                                      dLi_dTstar_3333 / norm_Tstar_sph
 else
  Li = 0.0_pReal
  dLi_dTstar_3333 = 0.0_pReal
 endif
 end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_dotState(Tstar_v,ipc,ip,el)
 use prec, only: &
   dEq0
 use math, only: &
   math_mul6x6
 use material, only: &
   phasememberAt, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6), intent(in):: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 type(tParameters), pointer :: p
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
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 p => param(instance)
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (p%dilatation) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
!--------------------------------------------------------------------------------------------------
! strain rate 
 gamma_dot = p%gdot0 * ( sqrt(1.5_pReal) * norm_Tstar_v & 
            / &!-----------------------------------------------------------------------------------
           (p%fTaylor*state(instance)%flowstress(of) ))**p%n
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (dEq0(p%tausat_SinhFitA)) then
     saturation = p%tausat
   else
     saturation = p%tausat &
                + asinh( (gamma_dot / p%tausat_SinhFitA&
                         )**(1.0_pReal / p%tausat_SinhFitD)&
                       )**(1.0_pReal / p%tausat_SinhFitC) &
                   / ( p%tausat_SinhFitB &
                       * (gamma_dot / p%gdot0)**(1.0_pReal / p%n) &
                     )
   endif
   hardening = ( p%h0 + p%h0_slopeLnRate * log(gamma_dot) ) &
               * abs( 1.0_pReal - state(instance)%flowstress(of)/saturation )**p%a &
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
   plasticState, &
   material_phase, &
   phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6),  intent(in) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 type(tParameters), pointer :: p
 
 real(pReal), dimension(plasticState(material_phase(ipc,ip,el))%sizePostResults) :: &
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
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 p => param(instance)
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) 2nd Piola-Kirchhoff stress
 if (p%dilatation) then
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_v,Tstar_v))
 else
   Tstar_dev_v(1:3) = Tstar_v(1:3) - sum(Tstar_v(1:3))/3.0_pReal
   Tstar_dev_v(4:6) = Tstar_v(4:6)
   norm_Tstar_v = sqrt(math_mul6x6(Tstar_dev_v,Tstar_dev_v))
 end if
 
 c = 0_pInt
 plastic_isotropic_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,plastic_isotropic_Noutput(instance)
   select case(p%outputID(o))
     case (flowstress_ID)
       plastic_isotropic_postResults(c+1_pInt) = state(instance)%flowstress(of)
       c = c + 1_pInt
     case (strainrate_ID)
       plastic_isotropic_postResults(c+1_pInt) = &
                p%gdot0 * (            sqrt(1.5_pReal) * norm_Tstar_v & 
             / &!----------------------------------------------------------------------------------
              (p%fTaylor * state(instance)%flowstress(of)) ) ** p%n
       c = c + 1_pInt
   end select
 enddo outputsLoop

end function plastic_isotropic_postResults


end module plastic_isotropic
