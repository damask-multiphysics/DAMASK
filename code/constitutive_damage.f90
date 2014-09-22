!--------------------------------------------------------------------------------------------------
! $Id: constitutive_damage.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @brief damage internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive_damage
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt), public, protected :: &
   constitutive_damage_maxSizePostResults, &
   constitutive_damage_maxSizeDotState

 public :: & 
   constitutive_damage_init, &
   constitutive_damage_microstructure, &
   constitutive_damage_collectDotState, &
   constitutive_damage_postResults
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_damage_init

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use IO, only: &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_timeStamp
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   FE_Nips, &
   FE_geomtype
 use material, only: &
   material_phase, &
   material_Nphase, &
   material_localFileExt, &    
   material_configFile, &    
   phase_name, &
   phase_damage, &
   phase_damageInstance, &
   phase_Noutput, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   damageState, &
   LOCAL_DAMAGE_NONE_ID, &
   LOCAL_DAMAGE_NONE_label, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   LOCAL_DAMAGE_BRITTLE_label
use damage_none
use damage_local
   
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                              !< grain number
  ph, &
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 logical :: knownDamage
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 
!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_damage == LOCAL_DAMAGE_NONE_ID))       call damage_none_init(FILEUNIT)
 if (any(phase_damage == LOCAL_DAMAGE_BRITTLE_ID))      call damage_local_init(FILEUNIT)
 close(FILEUNIT)
 
 write(6,'(/,a)')   ' <<<+-  constitutive_damage init  -+>>>'
 write(6,'(a)')     ' $Id: constitutive_damage.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputDamage') 
 do ph = 1_pInt,material_Nphase
   instance = phase_damageInstance(ph)                                                           ! which instance of a plasticity is present phase
   knownDamage = .true.
   select case(phase_damage(ph))                                                                 ! split per constititution
     case (LOCAL_DAMAGE_none_ID)
       outputName = LOCAL_DAMAGE_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_DAMAGE_BRITTLE_ID)
       outputName = LOCAL_DAMAGE_BRITTLE_label
       thisOutput => damage_local_output
       thisSize   => damage_local_sizePostResult
     case default
       knownDamage = .false.
   end select   
   write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(ph))//']'
   if (knownDamage) then
     write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
     if (phase_damage(ph) /= LOCAL_DAMAGE_none_ID) then
       do e = 1_pInt,phase_Noutput(ph)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
 enddo
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! allocation of states
 constitutive_damage_maxSizePostResults = 0_pInt
 constitutive_damage_maxSizeDotState = 0_pInt
 PhaseLoop:do ph = 1_pInt,material_Nphase                                                              ! loop over phases
  constitutive_damage_maxSizeDotState = max(constitutive_damage_maxSizeDotState, damageState(ph)%sizeDotState)
  constitutive_damage_maxSizePostResults = max(constitutive_damage_maxSizePostResults, damageState(ph)%sizePostResults)
 enddo PhaseLoop

end subroutine constitutive_damage_init


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_damage_microstructure(Tstar_v, Fe, ipc, ip, el)
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   phase_damage
 use damage_local, only:  &
   damage_local_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Fe
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     call damage_local_microstructure(Tstar_v, Fe, ipc, ip, el)

 end select

end subroutine constitutive_damage_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_damage_collectDotState(Tstar_v, Fe, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   phase_damage
 use damage_local, only:  &
   damage_local_dotState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp, &
   Fe
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     call damage_local_dotState(Tstar_v, Fe, Lp, ipc, ip, el)

 end select

end subroutine constitutive_damage_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_damage_postResults(ipc, ip, el)
 use material, only: &
   damageState, &
   material_phase, &
   LOCAL_DAMAGE_BRITTLE_ID, &
   phase_damage
 use damage_local, only:  &
   damage_local_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(damageState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_damage_postResults

 constitutive_damage_postResults = 0.0_pReal
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (LOCAL_DAMAGE_BRITTLE_ID)
     constitutive_damage_postResults = damage_local_postResults(ipc, ip, el)
 end select
  
end function constitutive_damage_postResults


end module constitutive_damage
