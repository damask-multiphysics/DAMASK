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
     critDip, &
     critPlasticStrain
 end type

! Begin Deprecated
 integer(pInt),                       dimension(:),           allocatable,         private :: &
   kinematics_slipplane_opening_totalNslip                                                                    !< total number of slip systems

 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   kinematics_slipplane_opening_Nslip                                                                         !< number of slip systems per family
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   kinematics_slipplane_opening_sdot_0, &
   kinematics_slipplane_opening_N

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   kinematics_slipplane_opening_critPlasticStrain, &
   kinematics_slipplane_opening_critLoad
! End Deprecated
   
 public :: &
   kinematics_slipplane_opening_init, &
   kinematics_slipplane_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_init()
#if defined(__GFORTRAN__) || __INTEL_COMPILER >= 1800
 use, intrinsic :: iso_fortran_env, only: &
   compiler_version, &
   compiler_options
#endif
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
 use config, only: &
   config_phase
 use IO, only: &
   IO_warning, &
   IO_error, &
   IO_timeStamp
 use material, only: &
   phase_kinematics, &
   KINEMATICS_slipplane_opening_label, &
   KINEMATICS_slipplane_opening_ID
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem

 implicit none
 integer(pInt), allocatable, dimension(:) :: tempInt
 real(pReal), allocatable, dimension(:) :: tempFloat

 integer(pInt) :: maxNinstance,p,instance,kinematics

 write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_slipplane_opening_LABEL//' init  -+>>>'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = int(count(phase_kinematics == KINEMATICS_slipplane_opening_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_slipplane_opening_instance(size(config_phase)), source=0_pInt)
 do p = 1_pInt, size(config_phase)
   kinematics_slipplane_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_slipplane_opening_ID) ! ToDo: count correct?
 enddo
   
 allocate(kinematics_slipplane_opening_critLoad(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_critPlasticStrain(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(kinematics_slipplane_opening_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(kinematics_slipplane_opening_N(maxNinstance),                                   source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_sdot_0(maxNinstance),                              source=0.0_pReal) 

 do p = 1_pInt, size(config_phase)
   if (all(phase_kinematics(:,p) /= KINEMATICS_slipplane_opening_ID)) cycle
   instance = kinematics_slipplane_opening_instance(p)
   kinematics_slipplane_opening_sdot_0(instance) = config_phase(p)%getFloat('anisoductile_sdot0')
   kinematics_slipplane_opening_N(instance)      = config_phase(p)%getFloat('anisoductile_ratesensitivity')
   tempInt = config_phase(p)%getInts('ncleavage')     
   kinematics_slipplane_opening_Nslip(1:size(tempInt),instance) = tempInt

   tempFloat = config_phase(p)%getFloats('anisoductile_criticalplasticstrain',requiredSize=size(tempInt))
   kinematics_slipplane_opening_critPlasticStrain(1:size(tempInt),instance) = tempFloat

   tempFloat = config_phase(p)%getFloats('anisoductile_criticalload',requiredSize=size(tempInt))
   kinematics_slipplane_opening_critLoad(1:size(tempInt),instance) = tempFloat

     kinematics_slipplane_opening_Nslip(1:lattice_maxNslipFamily,instance) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,p),&                                    ! limit active cleavage systems per family to min of available and requested
           kinematics_slipplane_opening_Nslip(1:lattice_maxNslipFamily,instance))
         kinematics_slipplane_opening_totalNslip(instance) = sum(kinematics_slipplane_opening_Nslip(:,instance))
     if (kinematics_slipplane_opening_sdot_0(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//KINEMATICS_slipplane_opening_LABEL//')')
     if (any(kinematics_slipplane_opening_critPlasticStrain(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='criticaPlasticStrain ('//KINEMATICS_slipplane_opening_LABEL//')')
     if (kinematics_slipplane_opening_N(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_slipplane_opening_LABEL//')')
 enddo
  
end subroutine kinematics_slipplane_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)
 use prec, only: &
   tol_math_check
 use math, only: &
   math_mul33xx33
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_sd, &
   lattice_st, &
   lattice_sn
 use material, only: &
   material_phase, &
   material_homog, &
   damage, &
   damageMapping
 use math, only: &
   math_tensorproduct33
 
 implicit none
 integer(pInt), intent(in) :: &
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
 integer(pInt) :: &
   instance, phase, &
   homog, damageOffset, &
   f, i, index_myFamily, k, l, m, n
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_slipplane_opening_instance(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)
 
 Ld = 0.0_pReal
 dLd_dTstar = 0.0_pReal
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                                      ! at which index starts my family
   do i = 1_pInt,kinematics_slipplane_opening_Nslip(f,instance)                                              ! process each (active) slip system in family
     projection_d = math_tensorproduct33(lattice_sd(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_t = math_tensorproduct33(lattice_st(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_n = math_tensorproduct33(lattice_sn(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))

   
     traction_d    = math_mul33xx33(S,projection_d)
     traction_t    = math_mul33xx33(S,projection_t)
     traction_n    = math_mul33xx33(S,projection_n)
     
     traction_crit = kinematics_slipplane_opening_critLoad(f,instance)* &
                     damage(homog)%p(damageOffset)                                                        ! degrading critical load carrying capacity by damage 

     udotd = &
       sign(1.0_pReal,traction_d)* &
       kinematics_slipplane_opening_sdot_0(instance)* &
       (abs(traction_d)/traction_crit - &
        abs(traction_d)/kinematics_slipplane_opening_critLoad(f,instance))**kinematics_slipplane_opening_N(instance)
     if (abs(udotd) > tol_math_check) then
       Ld = Ld + udotd*projection_d
       dudotd_dt = udotd*kinematics_slipplane_opening_N(instance)/traction_d
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudotd_dt*projection_d(k,l)*projection_d(m,n)
     endif                
 
     udott = &
       sign(1.0_pReal,traction_t)* &
       kinematics_slipplane_opening_sdot_0(instance)* &
       (abs(traction_t)/traction_crit - &
        abs(traction_t)/kinematics_slipplane_opening_critLoad(f,instance))**kinematics_slipplane_opening_N(instance)
     if (abs(udott) > tol_math_check) then
       Ld = Ld + udott*projection_t
       dudott_dt = udott*kinematics_slipplane_opening_N(instance)/traction_t
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudott_dt*projection_t(k,l)*projection_t(m,n)
     endif

     udotn = &
       kinematics_slipplane_opening_sdot_0(instance)* &
       (max(0.0_pReal,traction_n)/traction_crit - &
        max(0.0_pReal,traction_n)/kinematics_slipplane_opening_critLoad(f,instance))**kinematics_slipplane_opening_N(instance)
     if (abs(udotn) > tol_math_check) then
       Ld = Ld + udotn*projection_n
       dudotn_dt = udotn*kinematics_slipplane_opening_N(instance)/traction_n
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
           dudotn_dt*projection_n(k,l)*projection_n(m,n)
     endif 
   enddo
 enddo
  
end subroutine kinematics_slipplane_opening_LiAndItsTangent

end module kinematics_slipplane_opening
