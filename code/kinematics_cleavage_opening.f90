!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut fŸr Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut fŸr Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of cleavage planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_cleavage_opening
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_cleavage_opening_sizePostResults, &                                                                !< cumulative size of post results
   kinematics_cleavage_opening_offset, &                                                                         !< which kinematics is my current damage mechanism?
   kinematics_cleavage_opening_instance                                                                          !< instance of damage kinematics mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   kinematics_cleavage_opening_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   kinematics_cleavage_opening_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   kinematics_cleavage_opening_Noutput                                                                           !< number of outputs per instance of this damage 

 integer(pInt),                       dimension(:),           allocatable,         private :: &
   kinematics_cleavage_opening_totalNcleavage                                                                    !< total number of cleavage systems
   
 integer(pInt),                       dimension(:,:),         allocatable,         private :: &
   kinematics_cleavage_opening_Ncleavage                                                                         !< number of cleavage systems per family
   
 real(pReal),                         dimension(:),           allocatable,         private :: &
   kinematics_cleavage_opening_sdot_0, &
   kinematics_cleavage_opening_N

 real(pReal),                         dimension(:,:),         allocatable,         private :: &
   kinematics_cleavage_opening_critDisp, &
   kinematics_cleavage_opening_critLoad

 public :: &
   kinematics_cleavage_opening_init, &
   kinematics_cleavage_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_init(fileUnit)
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
   KINEMATICS_cleavage_opening_label, &
   KINEMATICS_cleavage_opening_ID, &
   material_Nphase, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank
 use lattice, only: &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,kinematics
 integer(pInt) :: Nchunks_CleavageFamilies = 0_pInt, j   
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_cleavage_opening_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_kinematics == KINEMATICS_cleavage_opening_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_cleavage_opening_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_cleavage_opening_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   kinematics_cleavage_opening_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_cleavage_opening_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_cleavage_opening_ID) &
       kinematics_cleavage_opening_offset(phase) = kinematics
   enddo    
 enddo
   
 allocate(kinematics_cleavage_opening_sizePostResults(maxNinstance),                      source=0_pInt)
 allocate(kinematics_cleavage_opening_sizePostResult(maxval(phase_Noutput),maxNinstance), source=0_pInt)
 allocate(kinematics_cleavage_opening_output(maxval(phase_Noutput),maxNinstance))
          kinematics_cleavage_opening_output = ''
 allocate(kinematics_cleavage_opening_Noutput(maxNinstance),                              source=0_pInt) 
 allocate(kinematics_cleavage_opening_critDisp(lattice_maxNcleavageFamily,maxNinstance),  source=0.0_pReal) 
 allocate(kinematics_cleavage_opening_critLoad(lattice_maxNcleavageFamily,maxNinstance),  source=0.0_pReal) 
 allocate(kinematics_cleavage_opening_Ncleavage(lattice_maxNcleavageFamily,maxNinstance), source=0_pInt)
 allocate(kinematics_cleavage_opening_totalNcleavage(maxNinstance),                       source=0_pInt)
 allocate(kinematics_cleavage_opening_sdot_0(maxNinstance),                               source=0.0_pReal) 
 allocate(kinematics_cleavage_opening_N(maxNinstance),                                    source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == KINEMATICS_cleavage_opening_ID)) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_cleavage_opening_instance(phase)                                                         ! which instance of my damage is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('anisobrittle_sdot0')
         kinematics_cleavage_opening_sdot_0(instance) = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('anisobrittle_ratesensitivity')
         kinematics_cleavage_opening_N(instance) = IO_floatValue(line,chunkPos,2_pInt)
         
       case ('ncleavage')  !
         Nchunks_CleavageFamilies = chunkPos(1) - 1_pInt
         do j = 1_pInt, Nchunks_CleavageFamilies
           kinematics_cleavage_opening_Ncleavage(j,instance) = IO_intValue(line,chunkPos,1_pInt+j)
         enddo

       case ('anisobrittle_criticaldisplacement')
         do j = 1_pInt, Nchunks_CleavageFamilies
           kinematics_cleavage_opening_critDisp(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo

       case ('anisobrittle_criticalload')
         do j = 1_pInt, Nchunks_CleavageFamilies
           kinematics_cleavage_opening_critLoad(j,instance) = IO_floatValue(line,chunkPos,1_pInt+j)
         enddo

     end select
   endif; endif
 enddo parsingFile

!--------------------------------------------------------------------------------------------------
!  sanity checks
 sanityChecks: do phase = 1_pInt, material_Nphase   
   myPhase: if (any(phase_kinematics(:,phase) == KINEMATICS_cleavage_opening_ID)) then
     instance = kinematics_cleavage_opening_instance(phase)
     kinematics_cleavage_opening_Ncleavage(1:lattice_maxNcleavageFamily,instance) = &
       min(lattice_NcleavageSystem(1:lattice_maxNcleavageFamily,phase),&                            ! limit active cleavage systems per family to min of available and requested
           kinematics_cleavage_opening_Ncleavage(1:lattice_maxNcleavageFamily,instance))
     kinematics_cleavage_opening_totalNcleavage(instance)  = sum(kinematics_cleavage_opening_Ncleavage(:,instance)) ! how many cleavage systems altogether
     if (kinematics_cleavage_opening_sdot_0(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='sdot_0 ('//KINEMATICS_cleavage_opening_LABEL//')')
     if (any(kinematics_cleavage_opening_critDisp(1:Nchunks_CleavageFamilies,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='critical_displacement ('//KINEMATICS_cleavage_opening_LABEL//')')
     if (any(kinematics_cleavage_opening_critLoad(1:Nchunks_CleavageFamilies,instance) < 0.0_pReal)) &
       call IO_error(211_pInt,el=instance,ext_msg='critical_load ('//KINEMATICS_cleavage_opening_LABEL//')')
     if (kinematics_cleavage_opening_N(instance) <= 0.0_pReal) &
       call IO_error(211_pInt,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_cleavage_opening_LABEL//')')
   endif myPhase
 enddo sanityChecks
 
end subroutine kinematics_cleavage_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar3333, Tstar_v, ipc, ip, el)
 use prec, only: &
   tol_math_check
 use material, only: &
   mappingConstitutive, &
   material_homog, &
   damage, &
   damageMapping
 use lattice, only: &
   lattice_Scleavage, &
   lattice_Scleavage_v, &
   lattice_maxNcleavageFamily, &
   lattice_NcleavageSystem
 
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
 instance = kinematics_cleavage_opening_instance(phase)
 homog = material_homog(ip,el)
 damageOffset = damageMapping(homog)%p(ip,el)
 
 Ld = 0.0_pReal
 dLd_dTstar3333 = 0.0_pReal
 do f = 1_pInt,lattice_maxNcleavageFamily
   index_myFamily = sum(lattice_NcleavageSystem(1:f-1_pInt,phase))                                   ! at which index starts my family
   do i = 1_pInt,kinematics_cleavage_opening_Ncleavage(f,instance)                                            ! process each (active) cleavage system in family
     traction_d    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,1,index_myFamily+i,phase))
     traction_t    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,2,index_myFamily+i,phase))
     traction_n    = dot_product(Tstar_v,lattice_Scleavage_v(1:6,3,index_myFamily+i,phase))
     traction_crit = kinematics_cleavage_opening_critLoad(f,instance)* &
                     damage(homog)%p(damageOffset)*damage(homog)%p(damageOffset)
     udotd = &
       sign(1.0_pReal,traction_d)* &
       kinematics_cleavage_opening_sdot_0(instance)* &
       (max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
     if (abs(udotd) > tol_math_check) then
       Ld = Ld + udotd*lattice_Scleavage(1:3,1:3,1,index_myFamily+i,phase)
       dudotd_dt = sign(1.0_pReal,traction_d)*udotd*kinematics_cleavage_opening_N(instance)/ &
                   max(0.0_pReal, abs(traction_d) - traction_crit)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotd_dt*lattice_Scleavage(k,l,1,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,1,index_myFamily+i,phase)
     endif                

     udott = &
       sign(1.0_pReal,traction_t)* &
       kinematics_cleavage_opening_sdot_0(instance)* &
       (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
     if (abs(udott) > tol_math_check) then
       Ld = Ld + udott*lattice_Scleavage(1:3,1:3,2,index_myFamily+i,phase)
       dudott_dt = sign(1.0_pReal,traction_t)*udott*kinematics_cleavage_opening_N(instance)/ &
                   max(0.0_pReal, abs(traction_t) - traction_crit)  
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudott_dt*lattice_Scleavage(k,l,2,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,2,index_myFamily+i,phase)
     endif                

     udotn = &
       sign(1.0_pReal,traction_n)* &
       kinematics_cleavage_opening_sdot_0(instance)* &
       (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**kinematics_cleavage_opening_N(instance)
     if (abs(udotn) > tol_math_check) then
       Ld = Ld + udotn*lattice_Scleavage(1:3,1:3,3,index_myFamily+i,phase)
       dudotn_dt = sign(1.0_pReal,traction_n)*udotn*kinematics_cleavage_opening_N(instance)/ &
                   max(0.0_pReal, abs(traction_n) - traction_crit)
       forall (k=1_pInt:3_pInt,l=1_pInt:3_pInt,m=1_pInt:3_pInt,n=1_pInt:3_pInt) &
         dLd_dTstar3333(k,l,m,n) = dLd_dTstar3333(k,l,m,n) + &
           dudotn_dt*lattice_Scleavage(k,l,3,index_myFamily+i,phase)* &
                     lattice_Scleavage(m,n,3,index_myFamily+i,phase)
     endif                

   enddo
 enddo
 
end subroutine kinematics_cleavage_opening_LiAndItsTangent

end module kinematics_cleavage_opening
