!--------------------------------------------------------------------------------------------------
! $Id: constitutive_damage.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @brief damage internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive_damage
 use prec, only: &
   pInt, &
   pReal
 
 implicit none
 private
 integer(pInt), public, dimension(:,:,:), allocatable :: &
   constitutive_damage_sizePostResults                                                                      !< size of postResults array per grain
 integer(pInt), public, protected :: &
   constitutive_damage_maxSizePostResults, &
   constitutive_damage_maxSizeDotState
 public :: & 
   constitutive_damage_init, &
   constitutive_damage_microstructure, &
   constitutive_damage_collectDotState, &
   constitutive_damage_collectDeltaState, &
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
   DAMAGE_none_ID, &
   DAMAGE_NONE_label, &
   DAMAGE_gradient_ID, &
   DAMAGE_GRADIENT_label
use damage_none
use damage_gradient
   
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
 logical :: knownDamage
 character(len=64), dimension(:,:), pointer :: thisOutput
 character(len=32) :: outputName                                                                    !< name of output, intermediate fix until HDF5 output is ready
 
!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(FILEUNIT,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(FILEUNIT,material_configFile)                                                  ! ... open material.config file
 if (any(phase_damage == DAMAGE_none_ID))       call damage_none_init(FILEUNIT)
 if (any(phase_damage == DAMAGE_gradient_ID))   call damage_gradient_init(FILEUNIT)
 close(FILEUNIT)
 
 write(6,'(/,a)')   ' <<<+-  constitutive_damage init  -+>>>'
 write(6,'(a)')     ' $Id: constitutive_damage.f90 3205 2014-06-17 06:54:49Z MPIE\m.diehl $'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(FILEUNIT,'outputDamage') 
 do phase = 1_pInt,material_Nphase
   instance = phase_damageInstance(phase)                                                           ! which instance of a plasticity is present phase
   knownDamage = .true.
   select case(phase_damage(phase))                                                                 ! split per constititution
     case (DAMAGE_none_ID)
       outputName = DAMAGE_NONE_label
       thisOutput => null()
       thisSize   => null()
     case (DAMAGE_gradient_ID)
       outputName = DAMAGE_GRADIENT_label
       thisOutput => damage_gradient_output
       thisSize   => damage_gradient_sizePostResult
     case default
       knownDamage = .false.
   end select   
   write(FILEUNIT,'(/,a,/)') '['//trim(phase_name(phase))//']'
   if (knownDamage) then
     write(FILEUNIT,'(a)') '(damage)'//char(9)//trim(outputName)
     if (phase_damage(phase) /= DAMAGE_none_ID) then
       do e = 1_pInt,phase_Noutput(phase)
         write(FILEUNIT,'(a,i4)') trim(thisOutput(e,instance))//char(9),thisSize(e,instance)
       enddo
     endif
   endif
 enddo
 close(FILEUNIT)
 
!--------------------------------------------------------------------------------------------------
! allocation of states
 cMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 allocate(constitutive_damage_sizePostResults(cMax,iMax,eMax), source=0_pInt)
 
 ElemLoop:do e = 1_pInt,mesh_NcpElems                                                               ! loop over elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e)) 
   IPloop:do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))                                     ! loop over IPs
     GrainLoop:do g = 1_pInt,myNgrains                                                              ! loop over grains
       phase = material_phase(g,i,e)
       instance = phase_damageInstance(phase)
       select case(phase_damage(phase))
         case (DAMAGE_gradient_ID)
           constitutive_damage_sizePostResults(g,i,e) =    damage_gradient_sizePostResults(instance)
           
       end select
     enddo GrainLoop
   enddo IPloop
 enddo ElemLoop
 
 constitutive_damage_maxSizePostResults = maxval(constitutive_damage_sizePostResults)
 constitutive_damage_maxSizeDotState = 0_pInt
 do p = 1, size(damageState)
  constitutive_damage_maxSizeDotState = max(constitutive_damage_maxSizeDotState, damageState(p)%sizeDotState)
 enddo
end subroutine constitutive_damage_init


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_damage_microstructure(Tstar_v, Fe, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_damage, &
   DAMAGE_gradient_ID
 use damage_gradient, only:  &
   damage_gradient_microstructure

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
   case (DAMAGE_gradient_ID)
     call damage_gradient_microstructure(Tstar_v, Fe, ipc, ip, el)

 end select

end subroutine constitutive_damage_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_damage_collectDotState(Tstar_v, Lp, ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_damage, &
   DAMAGE_gradient_ID
 use damage_gradient, only:  &
   damage_gradient_dotState

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 real(pReal),  intent(in), dimension(3,3) :: &
   Lp
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (DAMAGE_gradient_ID)
     call damage_gradient_dotState(Tstar_v, Lp, ipc, ip, el)

 end select

end subroutine constitutive_damage_collectDotState

!--------------------------------------------------------------------------------------------------
!> @brief for constitutive models having an instantaneous change of state (so far, only nonlocal)
!> will return false if delta state is not needed/supported by the constitutive model
!--------------------------------------------------------------------------------------------------
logical function constitutive_damage_collectDeltaState(ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_damage
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number

 select case (phase_damage(material_phase(ipc,ip,el)))

 end select
 constitutive_damage_collectDeltaState = .true.

end function constitutive_damage_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_damage_postResults(ipc, ip, el)
 use material, only: &
   material_phase, &
   phase_damage, &
   DAMAGE_gradient_ID
 use damage_gradient, only:  &
   damage_gradient_postResults

 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(constitutive_damage_sizePostResults(ipc,ip,el)) :: &
   constitutive_damage_postResults

 constitutive_damage_postResults = 0.0_pReal
 
 select case (phase_damage(material_phase(ipc,ip,el)))
   case (DAMAGE_gradient_ID)
     constitutive_damage_postResults = damage_gradient_postResults(ipc,ip,el)
 end select
  
end function constitutive_damage_postResults


end module constitutive_damage
