!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_thermal_expansion
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_thermal_expansion_sizePostResults, &                                                                !< cumulative size of post results
   kinematics_thermal_expansion_offset, &                                                                         !< which kinematics is my current damage mechanism?
   kinematics_thermal_expansion_instance                                                                          !< instance of damage kinematics mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   kinematics_thermal_expansion_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   kinematics_thermal_expansion_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   kinematics_thermal_expansion_Noutput                                                                           !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,         private :: &
   kinematics_thermal_expansion_coeff

 public :: &
   kinematics_thermal_expansion_init, &
   kinematics_thermal_expansion_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_init(fileUnit)
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
   KINEMATICS_thermal_expansion_label, &
   KINEMATICS_thermal_expansion_ID, &
   material_Nphase, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1+2*MAXNCHUNKS) :: positions
 integer(pInt) :: maxNinstance,phase,instance,kinematics
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_thermal_expansion_LABEL//' init  -+>>>'
   write(6,'(a)')     ' $Id$'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_kinematics == KINEMATICS_thermal_expansion_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_thermal_expansion_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_thermal_expansion_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   kinematics_thermal_expansion_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_thermal_expansion_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_thermal_expansion_ID) &
       kinematics_thermal_expansion_offset(phase) = kinematics
   enddo    
 enddo
   
 allocate(kinematics_thermal_expansion_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(kinematics_thermal_expansion_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(kinematics_thermal_expansion_output(maxval(phase_Noutput),maxNinstance))
          kinematics_thermal_expansion_output = ''
 allocate(kinematics_thermal_expansion_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(kinematics_thermal_expansion_coeff(maxNinstance),                               source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == KINEMATICS_thermal_expansion_ID)) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_thermal_expansion_instance(phase)                                                         ! which instance of my damage is present phase
     positions = IO_stringPos(line,MAXNCHUNKS)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('thermal_expansion_coeff')
         kinematics_thermal_expansion_coeff(instance) = IO_floatValue(line,positions,2_pInt)

     end select
   endif; endif
 enddo parsingFile

end subroutine kinematics_thermal_expansion_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar3333, ipc, ip, el)
 use material, only: &
   material_phase, &
   material_homog, &
   temperatureRate, &
   thermalMapping
 use math, only: &
   math_I3
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),   intent(out), dimension(3,3) :: &
   Li                                                                                               !< thermal velocity gradient
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dLi_dTstar3333                                                                                   !< derivative of Li with respect to Tstar (4th-order tensor)
 integer(pInt) :: &
   phase, &
   instance, &
   homog, offset
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_thermal_expansion_instance(phase)
 homog = material_homog(ip,el)
 offset = thermalMapping(homog)%p(ip,el)
 
 Li = temperatureRate(homog)%p(offset)* &
      kinematics_thermal_expansion_coeff(instance)* &
      math_I3
 dLi_dTstar3333 = 0.0_pReal
  
end subroutine kinematics_thermal_expansion_LiAndItsTangent

end module kinematics_thermal_expansion
