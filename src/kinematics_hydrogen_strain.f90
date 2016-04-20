!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from interstitial hydrogen
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_hydrogen_strain
 use prec, only: &
   pReal, &
   pInt

 implicit none
 private
 integer(pInt),                       dimension(:),           allocatable,         public, protected :: &
   kinematics_hydrogen_strain_sizePostResults, &                                                                !< cumulative size of post results
   kinematics_hydrogen_strain_offset, &                                                                         !< which kinematics is my current damage mechanism?
   kinematics_hydrogen_strain_instance                                                                          !< instance of damage kinematics mechanism

 integer(pInt),                       dimension(:,:),         allocatable, target, public  :: &
   kinematics_hydrogen_strain_sizePostResult                                                                    !< size of each post result output

 character(len=64),                   dimension(:,:),         allocatable, target, public  :: &
   kinematics_hydrogen_strain_output                                                                            !< name of each post result output
   
 integer(pInt),                       dimension(:),           allocatable, target, public  :: &
   kinematics_hydrogen_strain_Noutput                                                                           !< number of outputs per instance of this damage 

 real(pReal),                         dimension(:),           allocatable,         private :: &
   kinematics_hydrogen_strain_coeff

 public :: &
   kinematics_hydrogen_strain_init, &
   kinematics_hydrogen_strain_initialStrain, &
   kinematics_hydrogen_strain_LiAndItsTangent, &
   kinematics_hydrogen_strain_ChemPotAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_hydrogen_strain_init(fileUnit)
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
   KINEMATICS_hydrogen_strain_label, &
   KINEMATICS_hydrogen_strain_ID, &
   material_Nphase, &
   MATERIAL_partPhase
 use numerics,only: &
   worldrank

 implicit none
 integer(pInt), intent(in) :: fileUnit

 integer(pInt), allocatable, dimension(:) :: chunkPos
 integer(pInt) :: maxNinstance,phase,instance,kinematics
 character(len=65536) :: &
   tag  = '', &
   line = ''

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  kinematics_'//KINEMATICS_hydrogen_strain_LABEL//' init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

 maxNinstance = int(count(phase_kinematics == KINEMATICS_hydrogen_strain_ID),pInt)
 if (maxNinstance == 0_pInt) return
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(kinematics_hydrogen_strain_offset(material_Nphase), source=0_pInt)
 allocate(kinematics_hydrogen_strain_instance(material_Nphase), source=0_pInt)
 do phase = 1, material_Nphase
   kinematics_hydrogen_strain_instance(phase) = count(phase_kinematics(:,1:phase) == kinematics_hydrogen_strain_ID)
   do kinematics = 1, phase_Nkinematics(phase)
     if (phase_kinematics(kinematics,phase) == kinematics_hydrogen_strain_ID) &
       kinematics_hydrogen_strain_offset(phase) = kinematics
   enddo    
 enddo
   
 allocate(kinematics_hydrogen_strain_sizePostResults(maxNinstance),                     source=0_pInt)
 allocate(kinematics_hydrogen_strain_sizePostResult(maxval(phase_Noutput),maxNinstance),source=0_pInt)
 allocate(kinematics_hydrogen_strain_output(maxval(phase_Noutput),maxNinstance))
          kinematics_hydrogen_strain_output = ''
 allocate(kinematics_hydrogen_strain_Noutput(maxNinstance),                             source=0_pInt) 
 allocate(kinematics_hydrogen_strain_coeff(maxNinstance),                               source=0.0_pReal) 

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
   if (phase > 0_pInt ) then; if (any(phase_kinematics(:,phase) == KINEMATICS_hydrogen_strain_ID)) then         ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     instance = kinematics_hydrogen_strain_instance(phase)                                                         ! which instance of my damage is present phase
     chunkPos = IO_stringPos(line)
     tag = IO_lc(IO_stringValue(line,chunkPos,1_pInt))                                             ! extract key
     select case(tag)
       case ('hydrogen_strain_coeff')
         kinematics_hydrogen_strain_coeff(instance) = IO_floatValue(line,chunkPos,2_pInt)

     end select
   endif; endif
 enddo parsingFile

end subroutine kinematics_hydrogen_strain_init

!--------------------------------------------------------------------------------------------------
!> @brief  report initial hydrogen strain based on current hydrogen conc deviation from 
!>         equillibrium (0)
!--------------------------------------------------------------------------------------------------
pure function kinematics_hydrogen_strain_initialStrain(ipc, ip, el)
 use math, only: &
   math_I3
 use material, only: &
   material_phase, &
   material_homog, &
   hydrogenConc, &
   hydrogenfluxMapping
 use lattice, only: &
   lattice_equilibriumHydrogenConcentration
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3) :: &
   kinematics_hydrogen_strain_initialStrain                                                       !< initial thermal strain (should be small strain, though)
 integer(pInt) :: &
   phase, &
   homog, offset, instance
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_hydrogen_strain_instance(phase)
 homog = material_homog(ip,el)
 offset = hydrogenfluxMapping(homog)%p(ip,el)
 
 kinematics_hydrogen_strain_initialStrain = &
   (hydrogenConc(homog)%p(offset) - lattice_equilibriumHydrogenConcentration(phase)) * &
   kinematics_hydrogen_strain_coeff(instance)* math_I3
  
end function kinematics_hydrogen_strain_initialStrain

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_hydrogen_strain_LiAndItsTangent(Li, dLi_dTstar3333, ipc, ip, el)
 use material, only: &
   material_phase, &
   material_homog, &
   hydrogenConc, &
   hydrogenConcRate, &
   hydrogenfluxMapping
 use math, only: &
   math_I3
 use lattice, only: &
   lattice_equilibriumHydrogenConcentration
 
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
 real(pReal) :: &
   Ch, ChEq, ChDot  
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_hydrogen_strain_instance(phase)
 homog = material_homog(ip,el)
 offset = hydrogenfluxMapping(homog)%p(ip,el)
 Ch = hydrogenConc(homog)%p(offset)
 ChDot = hydrogenConcRate(homog)%p(offset)
 ChEq = lattice_equilibriumHydrogenConcentration(phase)
 
 Li = ChDot*math_I3* &
      kinematics_hydrogen_strain_coeff(instance)/ &
      (1.0_pReal + kinematics_hydrogen_strain_coeff(instance)*(Ch - ChEq))
 dLi_dTstar3333 = 0.0_pReal
  
end subroutine kinematics_hydrogen_strain_LiAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief  contains the kinematic contribution to hydrogen chemical potential  
!--------------------------------------------------------------------------------------------------
subroutine kinematics_hydrogen_strain_ChemPotAndItsTangent(ChemPot, dChemPot_dCh, Tstar_v, Fi0, Fi, ipc, ip, el)
 use material, only: &
   material_phase
 use math, only: &
   math_inv33, &
   math_mul33x33, &
   math_Mandel6to33, &
   math_transpose33
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v
 real(pReal),  intent(in), dimension(3,3) :: &
   Fi0, Fi
 real(pReal),   intent(out) :: &
   ChemPot, dChemPot_dCh                                                               
 integer(pInt) :: &
   phase, &
   instance
   
 phase = material_phase(ipc,ip,el)
 instance = kinematics_hydrogen_strain_instance(phase)
 
 ChemPot = -kinematics_hydrogen_strain_coeff(instance)* &
            sum(math_mul33x33(Fi,math_Mandel6to33(Tstar_v))* &
                math_mul33x33(math_mul33x33(Fi,math_inv33(Fi0)),Fi))
 dChemPot_dCh = 0.0_pReal
   
end subroutine kinematics_hydrogen_strain_ChemPotAndItsTangent

end module kinematics_hydrogen_strain
