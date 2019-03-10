!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of slip planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_slipplane_opening
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt), dimension(:), allocatable, private :: kinematics_slipplane_opening_instance
   
 type, private :: tParameters                                                                       !< container type for internal constitutive parameters
   integer(pInt) :: &
     totalNslip
   integer(pInt), dimension(:),   allocatable :: &
     Nslip                                                                                          !< active number of slip systems per family
   real(pReal) :: &
     sdot0, &
     n
   real(pReal),   dimension(:),   allocatable :: &
     critLoad
   real(pReal), dimension(:,:), allocatable     :: &
    slip_direction, &
    slip_normal, &
    slip_transverse
 end type tParameters

     type(tParameters), dimension(:), allocatable, private :: param                                     !< containers of constitutive parameters (len Ninstance)
 public :: &
   kinematics_slipplane_opening_init, &
   kinematics_slipplane_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_init()
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use config, only: &
   config_phase
 use IO, only: &
   IO_error
 use math, only: &
   math_expand
 use material, only: &
   phase_kinematics, &
   KINEMATICS_slipplane_opening_label, &
   KINEMATICS_slipplane_opening_ID
 use lattice

 implicit none

 integer(pInt) :: maxNinstance,p,instance,kinematics

 write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_slipplane_opening_LABEL//' init  -+>>>'

 maxNinstance = count(phase_kinematics == KINEMATICS_slipplane_opening_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_slipplane_opening_instance(size(config_phase)), source=0_pInt)
 do p = 1_pInt, size(config_phase)
   kinematics_slipplane_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_slipplane_opening_ID) ! ToDo: count correct?
 enddo
 
 allocate(param(maxNinstance))
 
 do p = 1_pInt, size(config_phase)
   if (all(phase_kinematics(:,p) /= KINEMATICS_slipplane_opening_ID)) cycle
   associate(prm => param(kinematics_slipplane_opening_instance(p)), &
            config => config_phase(p))
   instance = kinematics_slipplane_opening_instance(p)
   prm%sdot0 = config_phase(p)%getFloat('anisoductile_sdot0')
   prm%n      = config_phase(p)%getFloat('anisoductile_ratesensitivity')
   
   prm%Nslip = config%getInts('nslip') 

   prm%critLoad = config_phase(p)%getFloats('anisoductile_criticalload',requiredSize=size(prm%Nslip ))

   prm%critLoad = math_expand(prm%critLoad, prm%Nslip)
      
prm%slip_direction  = lattice_slip_direction  (prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_normal = lattice_slip_normal (prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))
      prm%slip_transverse     = lattice_slip_transverse(prm%Nslip,config%getString('lattice_structure'),&
                                                    config%getFloat('c/a',defaultVal=0.0_pReal))

   !  if (kinematics_slipplane_opening_sdot_0(instance) <= 0.0_pReal) &
   !    call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//KINEMATICS_slipplane_opening_LABEL//')')
   !  if (any(kinematics_slipplane_opening_critPlasticStrain(:,instance) < 0.0_pReal)) &
   !    call IO_error(211_pInt,el=instance,ext_msg='criticaPlasticStrain ('//KINEMATICS_slipplane_opening_LABEL//')')
   !  if (kinematics_slipplane_opening_N(instance) <= 0.0_pReal) &
   !    call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_slipplane_opening_LABEL//')')
 
   end associate
 enddo
  
end subroutine kinematics_slipplane_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)
 use prec, only: &
   tol_math_check
 use math, only: &
   math_mul33xx33, &
   math_outer
 use material, only: &
   material_phase, &
   material_homogenizationAt, &
   damage, &
   damageMapping

 implicit none
 integer, intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   S
 real(pReal),   intent(out), dimension(3,3) :: &
   Ld                                                                                               !< damage velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLd_dTstar                                                                                       !< derivative of Ld with respect to Tstar (4th-order tensor)
 real(pReal),   dimension(3,3) :: &
   projection_d, projection_t, projection_n                                                         !< projection modes 3x3 tensor
 integer :: &
   instance, phase, &
   homog, damageOffset, &
   i, k, l, m, n
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_slipplane_opening_instance(phase)
 homog = material_homogenizationAt(el)
 damageOffset = damageMapping(homog)%p(ip,el)
 
 associate(prm => param(instance))
 Ld = 0.0_pReal
 dLd_dTstar = 0.0_pReal
 do i = 1, prm%totalNslip

     projection_d = math_outer(prm%slip_direction(1:3,i),prm%slip_normal(1:3,i))
     projection_t = math_outer(prm%slip_transverse(1:3,i),prm%slip_normal(1:3,i))
     projection_n = math_outer(prm%slip_normal(1:3,i),prm%slip_normal(1:3,i))
   
     traction_d    = math_mul33xx33(S,projection_d)
     traction_t    = math_mul33xx33(S,projection_t)
     traction_n    = math_mul33xx33(S,projection_n)
     
     traction_crit = prm%critLoad(i)* damage(homog)%p(damageOffset)                                                        ! degrading critical load carrying capacity by damage 

     udotd = sign(1.0_pReal,traction_d)* &
       prm%sdot0* &
       (abs(traction_d)/traction_crit - &
        abs(traction_d)/prm%critLoad(i))**prm%n
     if (abs(udotd) > tol_math_check) then
       Ld = Ld + udotd*projection_d
       dudotd_dt = udotd*prm%n/traction_d
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudotd_dt*projection_d(k,l)*projection_d(m,n)
     endif                
 
     udott = sign(1.0_pReal,traction_t)* &
       prm%sdot0* &
       (abs(traction_t)/traction_crit - &
        abs(traction_t)/prm%critLoad(i))**prm%n
     if (abs(udott) > tol_math_check) then
       Ld = Ld + udott*projection_t
       dudott_dt = udott*prm%n/traction_t
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudott_dt*projection_t(k,l)*projection_t(m,n)
     endif

     udotn = &
       prm%sdot0* &
       (max(0.0_pReal,traction_n)/traction_crit - &
        max(0.0_pReal,traction_n)/prm%critLoad(i))**prm%n
     if (abs(udotn) > tol_math_check) then
       Ld = Ld + udotn*projection_n
       dudotn_dt = udotn*prm%n/traction_n
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudotn_dt*projection_n(k,l)*projection_n(m,n)
     endif 
 enddo
  
end associate

end subroutine kinematics_slipplane_opening_LiAndItsTangent

end module kinematics_slipplane_opening
