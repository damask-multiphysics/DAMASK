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
   constitutive_temperature, &
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
#ifdef NEWSTATE
   LOCAL_THERMAL_none_ID, &
   LOCAL_THERMAL_none_label, &
   LOCAL_THERMAL_heatgen_ID, &
   LOCAL_THERMAL_heatgen_label, &
#else
   THERMAL_none_ID, &
   THERMAL_NONE_label, &
   THERMAL_conduction_ID, &
   THERMAL_CONDUCTION_label, &
#endif
   homogenization_Ngrains, &
   homogenization_maxNgrains, &
   thermalState
 use thermal_none
#ifndef NEWSTATE
 use thermal_conduction
#endif
   
 implicit none
 integer(pInt), parameter :: FILEUNIT = 200_pInt
 integer(pInt) :: &
  e, &                                                                                              !< grain number
  ph, &                                                                                                !< phase
  instance

 integer(pInt), dimension(:,:), pointer :: thisSize
 logical :: knownThermal
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 
!--------------------------------------------------------------------------------------------------
! parse from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
#ifdef NEWSTATE
 if (any(phase_thermal == LOCAL_THERMAL_none_ID))       call thermal_none_init(FILEUNIT)
! if (any(phase_thermal == LOCAL_THERMAL_HEATGEN_ID)) call thermal_heatgen_init(FILEUNIT)
#else
 if (any(phase_thermal == THERMAL_none_ID))       call thermal_none_init(FILEUNIT)
 if (any(phase_thermal == THERMAL_conduction_ID)) call thermal_conduction_init(FILEUNIT)
#endif
 close(FILEUNIT)
 
 write(6,'(/,a)')   ' <<<+-  constitutive_thermal init  -+>>>'
 write(6,'(a)')     ' $Id: constitutive_thermal.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputThermal') 
 do ph = 1_pInt,material_Nphase
   instance = phase_thermalInstance(ph)                                                              ! which instance is present phase
   knownThermal = .true.
   select case(phase_thermal(ph))                                                                 ! split per constititution
#ifdef NEWSTATE
     case (LOCAL_THERMAL_none_ID)
       outputName = LOCAL_THERMAL_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (LOCAL_THERMAL_heatgen_ID)
       outputName = LOCAL_THERMAL_HEATGEN_label
       thisOutput => null()
       thisSize   => null()
#else
     case (THERMAL_none_ID)
       outputName = THERMAL_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (THERMAL_conduction_ID)
       outputName = THERMAL_CONDUCTION_label
       thisOutput => thermal_conduction_output
       thisSize   => thermal_conduction_sizePostResult
#endif
     case default
       knownThermal = .false.
   end select   
   write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(ph))//']'
   if (knownThermal) then
     write(FILEUNIT,'(a)') '(thermal)'//char(9)//trim(outputName)
#ifdef NEWSTATE 
     if (phase_thermal(ph) /= LOCAL_THERMAL_none_ID) then
#else
     if (phase_thermal(ph) /= THERMAL_none_ID) then
#endif
       do e = 1_pInt,phase_Noutput(ph)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
 enddo
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! allocation of states
 constitutive_thermal_maxSizePostResults = 0_pInt
 constitutive_thermal_maxSizeDotState = 0_pInt

 PhaseLoop:do ph = 1_pInt,material_Nphase                                                           ! loop over phases
  constitutive_thermal_maxSizeDotState = max(constitutive_thermal_maxSizeDotState, thermalState(ph)%sizeDotState)
  constitutive_thermal_maxSizePostResults = max(constitutive_thermal_maxSizePostResults, thermalState(ph)%sizePostResults)
 enddo PhaseLoop

end subroutine constitutive_thermal_init


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_thermal_microstructure(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
#ifndef NEWSTATE
   THERMAL_conduction_ID, &
#endif
   phase_thermal
#ifndef NEWSTATE
 use thermal_conduction, only: &
   thermal_conduction_microstructure
#endif

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
#ifndef NEWSTATE
   case (THERMAL_conduction_ID)
     call thermal_conduction_microstructure(Tstar_v, Lp, ipc, ip, el)
#endif   
 end select

end subroutine constitutive_thermal_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_thermal_collectDotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
#ifdef NEWSTATE
   LOCAL_THERMAL_none_ID, &
   LOCAL_THERMAL_HEATGEN_ID, &
#else
   THERMAL_none_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
#endif
   phase_thermal
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
#ifdef NEWSTATE
   case (LOCAL_THERMAL_HEATGEN_ID)
#else
   case (THERMAL_adiabatic_ID)
#endif
!     call thermal_adiabatic_dotState(Tstar_v, Lp, ipc, ip, el)
 end select

end subroutine constitutive_thermal_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief returns temperature based on each thermal model state layout 
!--------------------------------------------------------------------------------------------------
function constitutive_temperature(ipc, ip, el)
 use material, only: &
   material_phase, &
#ifdef NEWSTATE
   LOCAL_THERMAL_none_ID, &
   LOCAL_THERMAL_HEATGEN_ID, &
#else
   THERMAL_none_ID, &
   THERMAL_adiabatic_ID, &
   THERMAL_conduction_ID, &
#endif
   phase_thermal
 use lattice, only: &
   lattice_referenceTemperature
#ifndef NEWSTATE
 use thermal_conduction, only: &
   thermal_conduction_temperature
#endif
! use thermal_adiabatic, only: &
!   thermal_adiabatic_temperature

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal) :: constitutive_temperature
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
#ifdef NEWSTATE
   case (LOCAL_THERMAL_none_ID)
     constitutive_temperature = lattice_referenceTemperature(material_phase(ipc,ip,el))
   
   case (LOCAL_THERMAL_HEATGEN_ID)
     !constitutive_temperature = thermal_heatgen_temperature(ipc, ip, el)
#else
   case (THERMAL_none_ID)
     constitutive_temperature = lattice_referenceTemperature(material_phase(ipc,ip,el))
   
   case (THERMAL_adiabatic_ID)
     !constitutive_temperature = thermal_adiabatic_temperature(ipc, ip, el)

   case (THERMAL_conduction_ID)
     constitutive_temperature = thermal_conduction_temperature(ipc, ip, el)
#endif
 end select

end function constitutive_temperature

!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_thermal_postResults(ipc, ip, el)
 use material, only: &
   thermalState, &
   material_phase, &
#ifndef NEWSTATE
   THERMAL_conduction_ID, &
#endif
   phase_thermal
#ifndef NEWSTATE
 use thermal_conduction, only:  &
   thermal_conduction_postResults
#endif

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(thermalState(material_phase(ipc,ip,el))%sizePostResults) :: &
   constitutive_thermal_postResults

 constitutive_thermal_postResults = 0.0_pReal
 
 select case (phase_thermal(material_phase(ipc,ip,el)))
#ifndef NEWSTATE
   case (THERMAL_conduction_ID)
     constitutive_thermal_postResults = thermal_conduction_postResults(ipc,ip,el)
#endif
 end select
  
end function constitutive_thermal_postResults


end module constitutive_thermal
