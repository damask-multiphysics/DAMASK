!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for isotropic plasticity
!> @details Isotropic Plasticity which resembles the phenopowerlaw plasticity without
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
   plastic_isotropic_sizePostResult                                                                 !< size of each post result output
 character(len=64),                   dimension(:,:),   allocatable, target, public :: &
   plastic_isotropic_output                                                                         !< name of each post result output
 
 enum, bind(c) 
   enumerator :: &
     undefined_ID, &
     flowstress_ID, &
     strainrate_ID
 end enum

 type, private :: tParameters
  real(pReal) :: &
     fTaylor, &                                                                                     !< Taylor factor
     tau0, &                                                                                        !< initial critical stress
     gdot0, &                                                                                       !< reference strain rate
     n, &                                                                                           !< stress exponent
     h0, &
     h0_slopeLnRate, &
     tausat, &                                                                                      !< maximum critical stress
     a, &
     tausat_SinhFitA, &
     tausat_SinhFitB, &
     tausat_SinhFitC, &
     tausat_SinhFitD, &
     aTolFlowstress, &
     aTolShear
  integer(kind(undefined_ID)), allocatable, dimension(:) :: & 
     outputID
  logical :: &
     dilatation
 end type

 type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)
 
 type, private :: tIsotropicState
   real(pReal), pointer,     dimension(:) :: &
     flowstress, &
     accumulatedShear
 end type

 type(tIsotropicState), allocatable, dimension(:), private :: &
   dotState, &
   state
   
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
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO, only: &
   IO_error, &
   IO_timeStamp
 use material, only: &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   material_allocatePlasticState, &
   PLASTICITY_ISOTROPIC_label, &
   PLASTICITY_ISOTROPIC_ID, &
   material_phase, &
   plasticState
 use config, only: &
   MATERIAL_partPhase, &
   config_phase
 use lattice  

 implicit none
 
 integer(pInt) :: &
   p, & 
   instance, &
   maxNinstance, &
   sizeDotState, &
   sizeState
 character(len=65536) :: &
   extmsg    = ''
 integer(pInt) :: NipcMyPhase,i
 integer(kind(undefined_ID)) :: &
   outputID                                                                                         !< ID of each post result output

 character(len=65536), dimension(:), allocatable :: outputs

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

 allocate(param(maxNinstance))                                                                      ! one container of parameters per instance
 allocate(state(maxNinstance))                                                                      ! internal state aliases
 allocate(dotState(maxNinstance))

 do p = 1_pInt, size(phase_plasticityInstance)
   if (phase_plasticity(p) /= PLASTICITY_ISOTROPIC_ID) cycle
   instance = phase_plasticityInstance(p)
   associate(prm => param(instance))
   prm%tau0            =  config_phase(p)%getFloat('tau0')
   prm%tausat          =  config_phase(p)%getFloat('tausat')
   prm%gdot0           =  config_phase(p)%getFloat('gdot0')
   prm%n               =  config_phase(p)%getFloat('n')
   prm%h0              =  config_phase(p)%getFloat('h0')
   prm%fTaylor         =  config_phase(p)%getFloat('m')
   prm%h0_slopeLnRate  =  config_phase(p)%getFloat('h0_slopelnrate', defaultVal=0.0_pReal)
   prm%tausat_SinhFitA =  config_phase(p)%getFloat('tausat_sinhfita',defaultVal=0.0_pReal)
   prm%tausat_SinhFitB =  config_phase(p)%getFloat('tausat_sinhfitb',defaultVal=0.0_pReal)
   prm%tausat_SinhFitC =  config_phase(p)%getFloat('tausat_sinhfitc',defaultVal=0.0_pReal)
   prm%tausat_SinhFitD =  config_phase(p)%getFloat('tausat_sinhfitd',defaultVal=0.0_pReal)
   prm%a               =  config_phase(p)%getFloat('a')
   prm%aTolFlowStress  =  config_phase(p)%getFloat('atol_flowstress',defaultVal=1.0_pReal)
   prm%aTolShear       =  config_phase(p)%getFloat('atol_shear',defaultVal=1.0e-6_pReal)
   
   prm%dilatation      = config_phase(p)%keyExists('/dilatation/')

#if defined(__GFORTRAN__)
   outputs = ['GfortranBug86277']
   outputs = config_phase(p)%getStrings('(output)',defaultVal=outputs)
   if (outputs(1) == 'GfortranBug86277') outputs = [character(len=65536)::]
#else
   outputs = config_phase(p)%getStrings('(output)',defaultVal=[character(len=65536)::])
#endif
   allocate(prm%outputID(0))
   do i=1_pInt, size(outputs)
     outputID = undefined_ID
     select case(outputs(i))
       case ('flowstress')
         outputID = flowstress_ID
       case ('strainrate')
         outputID = strainrate_ID
     end select
     
     if (outputID /= undefined_ID) then
      plastic_isotropic_output(i,phase_plasticityInstance(p)) = outputs(i)
      plastic_isotropic_sizePostResult(i,phase_plasticityInstance(p)) = 1_pInt
      prm%outputID = [prm%outputID , outputID]
    endif

   enddo

!--------------------------------------------------------------------------------------------------
!  sanity checks
   extmsg = ''
   if (prm%aTolShear        <= 0.0_pReal) extmsg = trim(extmsg)//"'aTolShear' "
   if (prm%tau0              < 0.0_pReal) extmsg = trim(extmsg)//"'tau0' "
   if (prm%gdot0            <= 0.0_pReal) extmsg = trim(extmsg)//"'gdot0' "
   if (prm%n                <= 0.0_pReal) extmsg = trim(extmsg)//"'n' "
   if (prm%tausat           <= prm%tau0)    extmsg = trim(extmsg)//"'tausat' "
   if (prm%a                <= 0.0_pReal) extmsg = trim(extmsg)//"'a' " 
   if (prm%fTaylor          <= 0.0_pReal) extmsg = trim(extmsg)//"'m' "
   if (prm%aTolFlowstress   <= 0.0_pReal) extmsg = trim(extmsg)//"'atol_flowstress' "
   if (extmsg /= '') call IO_error(211_pInt,ip=instance,&
                          ext_msg=trim(extmsg)//'('//PLASTICITY_ISOTROPIC_label//')')

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   NipcMyPhase = count(material_phase == p)                                                   ! number of own material points (including point components ipc)

   sizeDotState   = size(["flowstress       ","accumulated_shear"])
   sizeState = sizeDotState
   call material_allocatePlasticState(p,NipcMyPhase,sizeState,sizeDotState,0_pInt, &
                                      1_pInt,0_pInt,0_pInt)
   plasticState(p)%sizePostResults = sum(plastic_isotropic_sizePostResult(:,phase_plasticityInstance(p)))
                                      

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState

   state(instance)%flowstress             => plasticState(p)%state    (1,1:NipcMyPhase)
   dotState(instance)%flowstress          => plasticState(p)%dotState (1,1:NipcMyPhase)
   plasticState(p)%state0(1,1:NipcMyPhase) = prm%tau0
   plasticState(p)%aTolState(1)       =  prm%aTolFlowstress

   state(instance)%accumulatedShear       => plasticState(p)%state    (2,1:NipcMyPhase)
   dotState(instance)%accumulatedShear    => plasticState(p)%dotState (2,1:NipcMyPhase)
   plasticState(p)%state0 (2,1:NipcMyPhase) = 0.0_pReal
   plasticState(p)%aTolState(2)       =  prm%aTolShear
   ! global alias
   plasticState(p)%slipRate           => plasticState(p)%dotState(2:2,1:NipcMyPhase)
   plasticState(p)%accumulatedSlip    => plasticState(p)%state   (2:2,1:NipcMyPhase)
  end associate

 enddo

end subroutine plastic_isotropic_init


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LpAndItsTangent(Lp,dLp_dMp,Mp,ipc,ip,el)
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
   math_deviatoric33, &
   math_mul33xx33
 use material, only: &
   phasememberAt, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3),     intent(out) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out) :: &
   dLp_dMp                                                                                          !< derivative of Lp with respect to the Mandel stress

 real(pReal), dimension(3,3), intent(in) :: &
   Mp                
 integer(pInt),               intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element

 
 real(pReal), dimension(3,3) :: &
   Mp_dev                                                                                           !< deviatoric part of the Mandel stress
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Mp_dev, &                                                                                   !< euclidean norm of the Mandel stress
   squarenorm_Mp_dev                                                                                !< square of the euclidean norm of the Mandel stress
 integer(pInt) :: &
   instance, of, &
   k, l, m, n

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 associate(prm => param(instance))
 
 Mp_dev = math_deviatoric33(Mp)
 squarenorm_Mp_dev = math_mul33xx33(Mp_dev,Mp_dev)
 norm_Mp_dev = sqrt(squarenorm_Mp_dev) 

 if (norm_Mp_dev <= 0.0_pReal) then
   Lp = 0.0_pReal
   dLp_dMp = 0.0_pReal
 else
   gamma_dot = prm%gdot0 &
             * ( sqrt(1.5_pReal) * norm_Mp_dev / prm%fTaylor / state(instance)%flowstress(of) ) &
             **prm%n

   Lp = Mp_dev/norm_Mp_dev * gamma_dot/prm%fTaylor 

   if (iand(debug_level(debug_constitutive), debug_levelExtensive) /= 0_pInt &
       .and. ((el == debug_e .and. ip == debug_i .and. ipc == debug_g) &
              .or. .not. iand(debug_level(debug_constitutive),debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CONST isotropic >> at el ip g ',el,ip,ipc
     write(6,'(/,a,/,3(12x,3(f12.4,1x)/))') '<< CONST isotropic >> Tstar (dev) / MPa', &
                                      transpose(Mp_dev)*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> norm Tstar / MPa', norm_Mp_dev*1.0e-6_pReal
     write(6,'(/,a,/,f12.5)') '<< CONST isotropic >> gdot', gamma_dot
   end if

   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLp_dMp(k,l,m,n) = (prm%n-1.0_pReal) * Mp_dev(k,l)*Mp_dev(m,n) / squarenorm_Mp_dev
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLp_dMp(k,l,k,l) = dLp_dMp(k,l,k,l) + 1.0_pReal
   forall (k=1_pInt:3_pInt,m=1_pInt:3_pInt) &
     dLp_dMp(k,k,m,m) = dLp_dMp(k,k,m,m) - 1.0_pReal/3.0_pReal
   dLp_dMp = gamma_dot / prm%fTaylor * dLp_dMp / norm_Mp_dev
 end if
 
 end associate
 
end subroutine plastic_isotropic_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_LiAndItsTangent(Li,dLi_dTstar,Tstar,instance,of)
 use math, only: &
   math_spherical33, &
   math_mul33xx33

 implicit none
 real(pReal), dimension(3,3), intent(out) :: &
   Li                                                                                               !< inleastic velocity gradient
 real(pReal), dimension(3,3,3,3), intent(out)  :: &
   dLi_dTstar                                                                                       !< derivative of Li with respect to the Mandel stress
 
 real(pReal), dimension(3,3),   intent(in) :: &
   Tstar                                                                                            !< Mandel stress ToDo: Mi?
 integer(pInt),                 intent(in) :: &
   instance, &
   of

 real(pReal), dimension(3,3) :: &
   Tstar_sph                                                                                        !< sphiatoric part of the Mandel stress
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   norm_Tstar_sph, &                                                                                !< euclidean norm of Tstar_sph
   squarenorm_Tstar_sph                                                                             !< square of the euclidean norm of Tstar_sph
 integer(pInt) :: &
   k, l, m, n

 associate(prm => param(instance), stt => state(instance))
 
 Tstar_sph = math_spherical33(Tstar)
 squarenorm_Tstar_sph = math_mul33xx33(Tstar_sph,Tstar_sph)
 norm_Tstar_sph = sqrt(squarenorm_Tstar_sph) 

 if (prm%dilatation .and. norm_Tstar_sph > 0.0_pReal) then                                          ! no stress or J2 plastitiy --> Li and its derivative are zero
   gamma_dot = prm%gdot0 * (sqrt(1.5_pReal) * norm_Tstar_sph /(prm%fTaylor*stt%flowstress(of))) **prm%n

   Li = Tstar_sph/norm_Tstar_sph * gamma_dot/prm%fTaylor
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
     dLi_dTstar(k,l,m,n) = (prm%n-1.0_pReal) * Tstar_sph(k,l)*Tstar_sph(m,n) / squarenorm_Tstar_sph
   forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt) &
     dLi_dTstar(k,l,k,l) = dLi_dTstar(k,l,k,l) + 1.0_pReal

   dLi_dTstar = gamma_dot / prm%fTaylor * dLi_dTstar / norm_Tstar_sph
 else
   Li = 0.0_pReal
   dLi_dTstar = 0.0_pReal
 endif
 
 end associate
 
 end subroutine plastic_isotropic_LiAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
subroutine plastic_isotropic_dotState(Mp,ipc,ip,el)
 use prec, only: &
   dEq0
 use math, only: &
   math_mul33xx33, &
   math_deviatoric33
 use material, only: &
   phasememberAt, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress    
 integer(pInt),             intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 real(pReal) :: &
   gamma_dot, &                                                                                     !< strainrate
   hardening, &                                                                                     !< hardening coefficient
   saturation, &                                                                                    !< saturation flowstress
   norm_Mp                                                                                          !< norm of the Mandel stress
 integer(pInt) :: &
   instance, &
   of                                                                                               !< shortcut notation for offset position in state array

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 associate(prm => param(instance), stt => state(instance), dot => dotState(instance))
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) Mandel stress
 if (prm%dilatation) then
   norm_Mp = sqrt(math_mul33xx33(Mp,Mp))
 else
   norm_Mp = sqrt(math_mul33xx33(math_deviatoric33(Mp),math_deviatoric33(Mp)))
 endif
 
 gamma_dot = prm%gdot0 * (sqrt(1.5_pReal) * norm_Mp /(prm%fTaylor*stt%flowstress(of))) **prm%n
 
!--------------------------------------------------------------------------------------------------
! hardening coefficient
 if (abs(gamma_dot) > 1e-12_pReal) then
   if (dEq0(prm%tausat_SinhFitA)) then
     saturation = prm%tausat
   else
     saturation = prm%tausat &
                + asinh( (gamma_dot / prm%tausat_SinhFitA)**(1.0_pReal / prm%tausat_SinhFitD) &
                       )**(1.0_pReal / prm%tausat_SinhFitC) &
                   / prm%tausat_SinhFitB * (gamma_dot / prm%gdot0)**(1.0_pReal / prm%n)
   endif
   hardening = ( prm%h0 + prm%h0_slopeLnRate * log(gamma_dot) ) &
               * abs( 1.0_pReal - stt%flowstress(of)/saturation )**prm%a &
               * sign(1.0_pReal, 1.0_pReal - stt%flowstress(of)/saturation)
 else
   hardening = 0.0_pReal
 endif

 dot%flowstress      (of) = hardening * gamma_dot
 dot%accumulatedShear(of) =             gamma_dot
 
 end associate

end subroutine plastic_isotropic_dotState


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
function plastic_isotropic_postResults(Mp,ipc,ip,el)
 use math, only: &
   math_mul33xx33, &
   math_deviatoric33
 use material, only: &
   plasticState, &
   material_phase, &
   phasememberAt, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(3,3), intent(in) :: &
   Mp                                                                                               !< Mandel stress    
 integer(pInt),              intent(in) :: &
   ipc, &                                                                                           !< component-ID of integration point
   ip, &                                                                                            !< integration point
   el                                                                                               !< element
 
 real(pReal), dimension(sum(plastic_isotropic_sizePostResult(:,phase_plasticityInstance(material_phase(ipc,ip,el))))) :: &
   plastic_isotropic_postResults


 real(pReal) :: &
   norm_Mp                                                                                          !< norm of the Mandel stress
 integer(pInt) :: &
   instance, &                                                                                      !< instance of my instance (unique number of my constitutive model)
   of, &                                                                                            !< shortcut notation for offset position in state array
   c, &
   o

 of = phasememberAt(ipc,ip,el)                                                                      ! phasememberAt should be tackled by material and be renamed to material_phasemember
 instance = phase_plasticityInstance(material_phase(ipc,ip,el))
 associate(prm => param(instance))
 
!--------------------------------------------------------------------------------------------------
! norm of (deviatoric) Mandel stress
 if (prm%dilatation) then
   norm_Mp = sqrt(math_mul33xx33(Mp,Mp))
 else
   norm_Mp = sqrt(math_mul33xx33(math_deviatoric33(Mp),math_deviatoric33(Mp)))
 endif
 
 c = 0_pInt
 plastic_isotropic_postResults = 0.0_pReal

 outputsLoop: do o = 1_pInt,size(prm%outputID)
   select case(prm%outputID(o))
     case (flowstress_ID)
       plastic_isotropic_postResults(c+1_pInt) = state(instance)%flowstress(of)
       c = c + 1_pInt
     case (strainrate_ID)
       plastic_isotropic_postResults(c+1_pInt) = &
                prm%gdot0 * (            sqrt(1.5_pReal) * norm_Mp & 
             / &!----------------------------------------------------------------------------------
              (prm%fTaylor * state(instance)%flowstress(of)) ) ** prm%n
       c = c + 1_pInt
   end select
 enddo outputsLoop
 
 end associate

end function plastic_isotropic_postResults


end module plastic_isotropic
