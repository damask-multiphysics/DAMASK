!--------------------------------------------------------------------------------------------------
! $Id: constitutive_thermal.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @brief thermal internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive_thermal
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt), public, protected :: &
   constitutive_thermal_maxSizePostResults, &
   constitutive_thermal_maxSizeDotState
 public :: & 
   constitutive_thermal_init, &
   constitutive_thermal_microstructure, &
   constitutive_thermal_collectDotState, &
   constitutive_thermal_postResults
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_thermal_init

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
   phase_thermal, &
   phase_thermalInstance, &
   phase_Noutput, &
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   thermalState, &
   THERMAL_none_ID, &
   THERMAL_NONE_label, &
   THERMAL_conduction_ID, &
   THERMAL_CONDUCTION_label
 use thermal_none
 use thermal_conduction
   
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e, &                                                                                              !< element number
  cMax, &                                                                                           !< maximum number of grains
  iMax, &                                                                                           !< maximum number of integration points
  eMax, &                                                                                           !< maximum number of elements
  phase, &
  s, &
  p, &
  instance,&
  myNgrains

 integer(pInt), dimension(:,:), pointer :: thisSize
 logical :: knownThermal
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 
!--------------------------------------------------------------------------------------------------
! parse from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_thermal == THERMAL_none_ID))       call thermal_none_init(FILEUNIT)
 if (any(phase_thermal == THERMAL_conduction_ID)) call thermal_conduction_init(FILEUNIT)
 close(FILEUNIT)
 
 write(6,'(/,a)')   ' <<<+-  constitutive_thermal init  -+>>>'
 write(6,'(a)')     ' $Id: constitutive_thermal.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputThermal') 
 do phase = 1_pInt,material_Nphase
   instance = phase_thermalInstance(phase)                                                           ! which instance is present phase
   knownThermal = .true.
   select case(phase_thermal(phase))                                                                 ! split per constititution
     case (THERMAL_none_ID)
       outputName = THERMAL_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (THERMAL_conduction_ID)
       outputName = THERMAL_CONDUCTION_label
       thisOutput => thermal_conduction_output
       thisSize   => thermal_conduction_sizePostResult
     case default
       knownThermal = .false.
   end select   
   write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(phase))//']'
   if (knownThermal) then
     write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
     if (phase_thermal(phase) /= THERMAL_none_ID) then
       do e = 1_pInt,phase_Noutput(phase)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
 enddo
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! allocation of states
 PhaseLoop:do phase = 1_pInt,material_Nphase                                                              ! loop over phases
   instance = phase_thermalInstance(phase)
   select case(phase_thermal(phase))
     case (THERMAL_none_ID) 
       thermalState(phase)%sizePostResults = thermal_none_sizePostResults(instance)

     case (THERMAL_conduction_ID) 
       thermalState(phase)%sizePostResults = thermal_conduction_sizePostResults(instance)
       
   end select
 enddo PhaseLoop
 
 constitutive_thermal_maxSizePostResults = 0_pInt
 constitutive_thermal_maxSizeDotState = 0_pInt
 do p = 1, size(thermalState)
  constitutive_thermal_maxSizeDotState = max(constitutive_thermal_maxSizeDotState, thermalState(p)%sizeDotState)
  constitutive_thermal_maxSizePostResults = max(constitutive_thermal_maxSizePostResults, thermalState(p)%sizePostResults)
 enddo
end subroutine constitutive_thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_thermal_microstructure(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_thermal, &
   THERMAL_conduction_ID
 use thermal_conduction, only: &
   thermal_conduction_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp

 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (THERMAL_conduction_ID)
     call thermal_conduction_microstructure(Tstar_v, Lp, ipc, ip, el)
 end select

end subroutine constitutive_thermal_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_thermal_collectDotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_thermal, &
   THERMAL_adiabatic_ID
! use thermal_conduction, only: &
!   thermal_adiabatic_microstructure

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (THERMAL_adiabatic_ID)
!     call thermal_adiabatic_dotState(Tstar_v, Lp, ipc, ip, el)
 end select

end subroutine constitutive_thermal_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_thermal_postResults(ipc, ip, el)
 use material, only: &
   thermalState, &
   material_phase, &
   phase_thermal, &
   THERMAL_conduction_ID
 use thermal_conduction, only:  &
   thermal_conduction_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(thermalState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_thermal_postResults

 constitutive_thermal_postResults = 0.0_pReal
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
   case (THERMAL_conduction_ID)
     constitutive_thermal_postResults = thermal_conduction_postResults(ipc,ip,el)
 end select
  
end function constitutive_thermal_postResults


end module constitutive_thermal
