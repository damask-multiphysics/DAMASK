!--------------------------------------------------------------------------------------------------
! $Id$
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
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_slipplane_opening_sizePostResults, &                                                            !< cumulative size of post results
   kinematics_slipplane_opening_offset, &                                                                     !< which kinematics is my current damage mechanism?
   kinematics_slipplane_opening_instance                                                                      !< instance of damage kinematics mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   kinematics_slipplane_opening_sizePostResult                                                                !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   kinematics_slipplane_opening_output                                                                        !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   kinematics_slipplane_opening_Noutput                                                                       !< number of outputs per instance of this damage 
   
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
   
 public :: &
   kinematics_slipplane_opening_init, &
   kinematics_slipplane_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level,&
   debug_constitutive,&
   debug_levelBasic
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
   phase_kinematics, &
   phase_Nkinematics, &
   phase_Noutput, &
   KINEMATICS_slipplane_opening_label, &
   KINEMATICS_slipplane_opening_ID, &
   material_Nphase, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,kinematics
 integer(pInt) :: Nchunks_SlipFamilies = 0_pInt, j   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_slipplane_opening_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_kinematics == KINEMATICS_slipplane_opening_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_slipplane_opening_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_slipplane_opening_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   kinematics_slipplane_opening_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_slipplane_opening_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_slipplane_opening_ID) &
       kinematics_slipplane_opening_offset(phase) = kinematics
   enddo    
 enddo
   
 allocate(kinematics_slipplane_opening_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(kinematics_slipplane_opening_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(kinematics_slipplane_opening_output(maxval(phase_Noutput),maxNinstance))
          kinematics_slipplane_opening_output = ''
 allocate(kinematics_slipplane_opening_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(kinematics_slipplane_opening_critLoad(lattice_maxNslipFamily,maxNinstance),     source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_critPlasticStrain(lattice_maxNslipFamily,maxNinstance),source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_Nslip(lattice_maxNslipFamily,maxNinstance),        source=0_pInt)
 allocate(kinematics_slipplane_opening_totalNslip(maxNinstance),                          source=0_pInt)
 allocate(kinematics_slipplane_opening_N(maxNinstance),                                   source=0.0_pReal) 
 allocate(kinematics_slipplane_opening_sdot_0(maxNinstance),                              source=0.0_pReal) 

 rewind(fileUnit)
 phase = 0_pInt
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= MATERIAL_partPhase)         ! wind forward to <phase>
   line = IO_read(fileUnit)
 enddo
 
 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of phase part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif   
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next phase section
     phase = phase + 1_pInt                                                                         ! advance phase section counter
     cycle                                                                                          ! skip to next line
   endif
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == KINEMATICS_slipplane_opening_ID)) then ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_slipplane_opening_instance(phase)                                        ! which instance of my damage is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                              ! extract key
     select case(tag)
       case ('nslip')  !
         Nchunks_SlipFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_SlipFamilies
           kinematics_slipplane_opening_Nslip(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo

       case ('anisoductile_sdot0')
         kinematics_slipplane_opening_sdot_0(instance) = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('anisoductile_criticalplasticstrain')
         do j = 1_pInt, Nchunks_SlipFamilies
           kinematics_slipplane_opening_critPlasticStrain(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
       case ('anisoductile_ratesensitivity')
         kinematics_slipplane_opening_N(instance) = IO_floatValue(line,chunkPos,2_pInt)

       case ('anisoductile_criticalload')
         do j = 1_pInt, Nchunks_SlipFamilies
           kinematics_slipplane_opening_critLoad(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo
         
     end select
   endif; endif
 enddo parsingFile

!--------------------------------------------------------------------------------------------------
!  sanity checks
 sanityChecks: do phase = 1_pInt, material_Nphase   
   myPhase: if (any(phase_kinematics(:,phase) == KINEMATICS_slipplane_opening_ID)) then
     instance = kinematics_slipplane_opening_instance(phase)
     kinematics_slipplane_opening_Nslip(1:lattice_maxNslipFamily,instance) = &
       min(lattice_NslipSystem(1:lattice_maxNslipFamily,phase),&                                    ! limit active cleavage systems per family to min of available and requested
           kinematics_slipplane_opening_Nslip(1:lattice_maxNslipFamily,instance))
         kinematics_slipplane_opening_totalNslip(instance) = sum(kinematics_slipplane_opening_Nslip(:,instance))
     if (kinematics_slipplane_opening_sdot_0(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//KINEMATICS_slipplane_opening_LABEL//')')
     if (any(kinematics_slipplane_opening_critPlasticStrain(:,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='criticaPlasticStrain ('//KINEMATICS_slipplane_opening_LABEL//')')
     if (kinematics_slipplane_opening_N(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_slipplane_opening_LABEL//')')
   endif myPhase
 enddo sanityChecks
  

end subroutine kinematics_slipplane_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar3333, Tstar_v, ipc, ip, el)
 use prec, only: &
   tol_math_check
 use lattice, only: &
   lattice_maxNslipFamily, &
   lattice_NslipSystem, &
   lattice_sd, &
   lattice_st, &
   lattice_sn
 use material, only: &
   mappingConstitutive, &
   material_homog, &
   damage, &
   damageMapping
 use math, only: &
   math_Plain3333to99, &
   math_I3, &
   math_identity4th, &
   math_symmetric33, &
   math_Mandel33to6, &
   math_tensorproduct33, &
   math_det33, &
   math_mul33x33
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(out), dimension(3,3) :: &
   Ld                                                                                               !< damage velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLd_dTstar3333                                                                                   !< derivative of Ld with respect to Tstar (4th-order tensor)
 real(pReal),   dimension(3,3) :: &
   projection_d, projection_t, projection_n                                                         !< projection modes 3x3 tensor
 real(pReal),   dimension(6) :: &
   projection_d_v, projection_t_v, projection_n_v                                                   !< projection modes 3x3 vector
 integer(pInt) :: &
   phase, &
   constituent, &
   instance, &
   homog, damageOffset, &
   f, i, index_myFamily, k, l, m, n
 real(pReal) :: &
   traction_d, traction_t, traction_n, traction_crit, &
   udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt
   
 phase = mappingConstitutive(2,ipc,ip,el)
 constituent = mappingConstitutive(1,ipc,ip,el)
 instance = kinematics_slipplane_opening_instance(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)
 
 Ld = 0.0_pReal
 dLd_dTstar3333 = 0.0_pReal
 do f = 1_pInt,lattice_maxNslipFamily
   index_myFamily = sum(lattice_NslipSystem(1:f-1_pInt,phase))                                      ! at which index starts my family
   do i = 1_pInt,kinematics_slipplane_opening_Nslip(f,instance)                                              ! process each (active) slip system in family
     projection_d = math_tensorproduct33(lattice_sd(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_t = math_tensorproduct33(lattice_st(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))
     projection_n = math_tensorproduct33(lattice_sn(1:3,index_myFamily+i,phase),&
                                       lattice_sn(1:3,index_myFamily+i,phase))

     projection_d_v(1:6) = math_Mandel33to6(math_symmetric33(projection_d(1:3,1:3)))
     projection_t_v(1:6) = math_Mandel33to6(math_symmetric33(projection_t(1:3,1:3)))
     projection_n_v(1:6) = math_Mandel33to6(math_symmetric33(projection_n(1:3,1:3)))
   
     traction_d    = dot_product(Tstar_v,projection_d_v(1:6))
     traction_t    = dot_product(Tstar_v,projection_t_v(1:6))
     traction_n    = dot_product(Tstar_v,projection_n_v(1:6))
     
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
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
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
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
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
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotn_dt*projection_n(k,l)*projection_n(m,n)
     endif 
   enddo
 enddo
  
end subroutine kinematics_slipplane_opening_LiAndItsTangent

end module kinematics_slipplane_opening
