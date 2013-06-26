! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine for purely elastic material
!--------------------------------------------------------------------------------------------------
module constitutive_none
 use prec, only: &
   pReal, &
   pInt
 
 implicit none
 private
 character (len=*), parameter, public :: constitutive_none_LABEL = 'none'
 
 integer(pInt),   dimension(:), allocatable, public :: &
   constitutive_none_sizeDotState, &
   constitutive_none_sizeState, &
   constitutive_none_sizePostResults
   
 character(len=32), dimension(:), allocatable, private :: &
   constitutive_none_structureName

 integer(pInt), dimension(:,:), allocatable, target, public :: &
   constitutive_none_sizePostResult                                                                 ! size of each post result output

 real(pReal), dimension(:,:,:), allocatable, private :: &
   constitutive_none_Cslip_66

 public :: &
   constitutive_none_init, &
   constitutive_none_stateInit, &
   constitutive_none_aTolState, &
   constitutive_none_homogenizedC, &
   constitutive_none_microstructure, &
   constitutive_none_LpAndItsTangent, &
   constitutive_none_dotState, &
   constitutive_none_deltaState, &
   constitutive_none_dotTemperature, &
   constitutive_none_postResults

 contains


!--------------------------------------------------------------------------------------------------
!> @brief reads in material parameters and allocates arrays
!--------------------------------------------------------------------------------------------------
subroutine constitutive_none_init(myFile)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO, only: &
   IO_read, &
   IO_lc, &
   IO_getTag, &
   IO_isBlank, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_error, &
   IO_timeStamp
 use material
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use lattice, only: lattice_symmetrizeC66

 implicit none
 integer(pInt), intent(in) :: myFile
 
 integer(pInt), parameter :: MAXNCHUNKS = 7_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS) :: positions
 integer(pInt) :: section = 0_pInt, maxNinstance, i
 character(len=65536) :: tag
 character(len=65536) :: line = ''                                                                   ! to start initialized
 
 write(6,'(/,a)')   ' <<<+-  constitutive_'//trim(constitutive_none_LABEL)//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 
 maxNinstance = int(count(phase_plasticity == constitutive_none_label),pInt)
 if (maxNinstance == 0_pInt) return

 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 
 allocate(constitutive_none_sizeDotState(maxNinstance))
          constitutive_none_sizeDotState = 0_pInt
 allocate(constitutive_none_sizeState(maxNinstance))
          constitutive_none_sizeState = 0_pInt
 allocate(constitutive_none_sizePostResults(maxNinstance))
          constitutive_none_sizePostResults = 0_pInt
 allocate(constitutive_none_structureName(maxNinstance))
          constitutive_none_structureName        = ''
 allocate(constitutive_none_Cslip_66(6,6,maxNinstance))
          constitutive_none_Cslip_66 = 0.0_pReal
 
 rewind(myFile)
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= 'phase')                                                  ! wind forward to <phase>
   line = IO_read(myFile)
 enddo
 
 do while (trim(line) /= '#EOF#')                                                                                                                                ! read thru sections of phase part
   line = IO_read(myFile)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt                                                                     ! advance section counter
     cycle
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if statemen). It's not safe in Fortran
     if ( phase_plasticity(section) == constitutive_none_LABEL) then                                ! one of my sections
       i = phase_plasticityInstance(section)                                                        ! which instance of my plasticity is present phase
       positions = IO_stringPos(line,MAXNCHUNKS)
       tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                           ! extract key
       select case(tag)
         case ('plasticity','elasticity')
           cycle
         case ('lattice_structure')
                constitutive_none_structureName(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
         case ('c11')
                constitutive_none_Cslip_66(1,1,i) = IO_floatValue(line,positions,2_pInt)
         case ('c12')
                constitutive_none_Cslip_66(1,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c13')
                constitutive_none_Cslip_66(1,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c22')
                constitutive_none_Cslip_66(2,2,i) = IO_floatValue(line,positions,2_pInt)
         case ('c23')
                constitutive_none_Cslip_66(2,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c33')
                constitutive_none_Cslip_66(3,3,i) = IO_floatValue(line,positions,2_pInt)
         case ('c44')
                constitutive_none_Cslip_66(4,4,i) = IO_floatValue(line,positions,2_pInt)
         case ('c55')
                constitutive_none_Cslip_66(5,5,i) = IO_floatValue(line,positions,2_pInt)
         case ('c66')
                constitutive_none_Cslip_66(6,6,i) = IO_floatValue(line,positions,2_pInt)
         case default
                call IO_error(210_pInt,ext_msg=tag//' ('//constitutive_none_label//')')
       end select
     endif
   endif
 enddo

 do i = 1_pInt,maxNinstance                 
   if (constitutive_none_structureName(i) == '')              call IO_error(205_pInt,e=i)
 enddo

 do i = 1_pInt,maxNinstance
   constitutive_none_sizeDotState(i)    = 1_pInt
   constitutive_none_sizeState(i)       = 1_pInt

   constitutive_none_Cslip_66(:,:,i) = lattice_symmetrizeC66(constitutive_none_structureName(i),&
                                                                      constitutive_none_Cslip_66(:,:,i))
   constitutive_none_Cslip_66(:,:,i) = &
     math_Mandel3333to66(math_Voigt66to3333(constitutive_none_Cslip_66(:,:,i)))

 enddo

end subroutine constitutive_none_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial microstructural state
!--------------------------------------------------------------------------------------------------
pure function constitutive_none_stateInit(myInstance)
  
 implicit none
 integer(pInt), intent(in) :: myInstance

 real(pReal), dimension(1) :: constitutive_none_stateInit
  
 constitutive_none_stateInit = 0.0_pReal

end function constitutive_none_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief relevant microstructural state (ensures convergence as state is always 0.0)
!--------------------------------------------------------------------------------------------------
pure function constitutive_none_aTolState(myInstance)

 implicit none
 integer(pInt), intent(in) ::  myInstance                                                           !< number specifying the current instance of the plasticity

 real(pReal), dimension(constitutive_none_sizeState(myInstance)) :: &
                                                              constitutive_none_aTolState           !< relevant state values for the current instance of this plasticity                                                           
 constitutive_none_aTolState = 1.0_pReal

end function constitutive_none_aTolState


!--------------------------------------------------------------------------------------------------
!> @brief homogenized elacticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_none_homogenizedC(state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
  homogenization_maxNgrains, &
  material_phase, &
  phase_plasticityInstance
 
 implicit none
 real(pReal), dimension(6,6) :: constitutive_none_homogenizedC
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of current integration point
   ip, &                                                                                            !< current integration point
   el                                                                                               !< current element
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state
 integer(pInt)               :: matID

 matID = phase_plasticityInstance(material_phase(ipc,ip,el))
 constitutive_none_homogenizedC = constitutive_none_Cslip_66(1:6,1:6,matID)

end function constitutive_none_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state (not used here)
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_none_microstructure(Temperature,state,ipc,ip,el)
 use prec, only: &
   p_vec
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance
 
 implicit none
 integer(pInt), intent(in) :: &
   ipc, &                                                                                           !< component-ID of current integration point
   ip, &                                                                                            !< current integration point
   el                                                                                               !< current element
 real(pReal), intent(in) :: &
   Temperature                                                                                      !< temperature
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state

end subroutine constitutive_none_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_none_LpAndItsTangent(Lp, dLp_dTstar_99, Tstar_dev_v, Temperature, &
                                                                            state, gr, ip, el)
 use prec, only: &
   p_vec
 use math, only: &
   math_identity2nd
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6), intent(in) :: Tstar_dev_v                                               !< deviatoric part of the 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: Temperature
 integer(pInt),             intent(in) :: &
   gr, &                                                                                            !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state  !< state of the current microstructure

 real(pReal), dimension(3,3), intent(out) ::   Lp                                                   !< plastic velocity gradient
 real(pReal), dimension(9,9), intent(out) ::   dLp_dTstar_99                                        !< derivative of Lp with respect to Tstar (9x9 matrix)

 Lp = 0.0_pReal                                                                                     !< set Lp to zero 
 dLp_dTstar_99 = math_identity2nd(9)                                                                !< set dLp_dTstar to Identity

end subroutine constitutive_none_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
pure function constitutive_none_dotState(Tstar_v, Temperature, state, gr, ip, el)
 use prec, only: &
   p_vec
 use math, only: &
   math_identity2nd
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(1)             :: constitutive_none_dotState
 
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: Temperature
 integer(pInt),             intent(in) :: &
   gr, &                                                                                            !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state  !< state of the current microstructure

 constitutive_none_dotState =  0.0_pReal

end function constitutive_none_dotState


!--------------------------------------------------------------------------------------------------
!> @brief (instantaneous) incremental change of microstructure
!--------------------------------------------------------------------------------------------------
function constitutive_none_deltaState(Tstar_v, Temperature, state, gr, ip, el)
 use prec, only: &
   p_vec
 use math, only: &
   math_identity2nd
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance

 implicit none
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: Temperature
 integer(pInt),             intent(in) :: &
   gr, &                                                                                            !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state  !< state of the current microstructure

 real(pReal), dimension(constitutive_none_sizeDotState(phase_plasticityInstance(&
                                       material_phase(gr,ip,el)))) :: constitutive_none_deltaState 

 constitutive_none_deltaState = 0.0_pReal


end function constitutive_none_deltaState


!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of temperature
!--------------------------------------------------------------------------------------------------
pure real(pReal) function constitutive_none_dotTemperature(Tstar_v, Temperature, state, gr, ip, el)
 use prec, only: &
   p_vec
 use math, only: &
   math_identity2nd
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains
  
 implicit none
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: Temperature
 integer(pInt),             intent(in) :: &
   gr, &                                                                                            !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state  !< state of the current microstructure

 constitutive_none_dotTemperature = 0.0_pReal

end function constitutive_none_dotTemperature


!--------------------------------------------------------------------------------------------------
!> @brief return array of constitutive results
!--------------------------------------------------------------------------------------------------
pure function constitutive_none_postResults(Tstar_v, Temperature, dt, state, gr, ip, el)
 use prec, only: &
   p_vec
 use math, only: &
   math_mul6x6
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase, &
   phase_plasticityInstance, &
   phase_Noutput

 implicit none
 real(pReal), dimension(6), intent(in) :: Tstar_v                                                   !< 2nd Piola Kirchhoff stress tensor in Mandel notation
 real(pReal),               intent(in) :: &
   Temperature, &
   dt                                                                                               !< current time increment
 integer(pInt),             intent(in) :: &
   gr, &                                                                                            !< grain number
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 type(p_vec), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: state  !< state of the current microstructure

 real(pReal), dimension(constitutive_none_sizePostResults(phase_plasticityInstance(&
                                      material_phase(gr,ip,el)))) :: constitutive_none_postResults
 
 constitutive_none_postResults = 0.0_pReal

end function constitutive_none_postResults

end module constitutive_none
