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
!> @brief Isostrain (full constraint Taylor assuption) homogenization scheme
!--------------------------------------------------------------------------------------------------
module homogenization_isostrain
 use prec, only: &
   pInt
 
 implicit none
 private
 character (len=*), parameter,                           public  :: &
   homogenization_isostrain_label = 'isostrain'
 
 integer(pInt), dimension(:),       allocatable,         public  :: &
   homogenization_isostrain_sizeState, &
   homogenization_isostrain_sizePostResults
 integer(pInt), dimension(:,:),     allocatable, target, public  :: &
   homogenization_isostrain_sizePostResult
 character(len=64), dimension(:,:), allocatable, target, public  :: &
  homogenization_isostrain_output                                                                   !< name of each post result output
 
 integer(pInt), dimension(:),       allocatable,         private :: &
   homogenization_isostrain_Ngrains

 public :: &
   homogenization_isostrain_init, &
   homogenization_isostrain_stateInit, &
   homogenization_isostrain_partitionDeformation, &
   homogenization_isostrain_updateState, &
   homogenization_isostrain_averageStressAndItsTangent, &
   homogenization_isostrain_averageTemperature, &
   homogenization_isostrain_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_init(myFile)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: math_Mandel3333to66, math_Voigt66to3333
 use IO
 use material
 integer(pInt), intent(in) :: myFile
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 integer(pInt) section, i, j, output, mySize
 integer :: maxNinstance, k                                                                         ! no pInt (stores a system dependen value from 'count'
 character(len=64)   :: tag
 character(len=1024) :: line = ''                                                                   ! to start initialized
 

 write(6,*)
 write(6,*) '<<<+-  homogenization_',trim(homogenization_isostrain_label),' init  -+>>>'
 write(6,*) '$Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"


 maxNinstance = count(homogenization_type == homogenization_isostrain_label)
 if (maxNinstance == 0) return

 allocate(homogenization_isostrain_sizeState(maxNinstance)) ;      homogenization_isostrain_sizeState = 0_pInt
 allocate(homogenization_isostrain_sizePostResults(maxNinstance)); homogenization_isostrain_sizePostResults = 0_pInt
 allocate(homogenization_isostrain_sizePostResult(maxval(homogenization_Noutput), &
                                                  maxNinstance)); homogenization_isostrain_sizePostResult = 0_pInt
 allocate(homogenization_isostrain_Ngrains(maxNinstance));         homogenization_isostrain_Ngrains = 0_pInt
 allocate(homogenization_isostrain_output(maxval(homogenization_Noutput), &
                                          maxNinstance)) ;         homogenization_isostrain_output = ''
 
 rewind(myFile)
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)                           ! wind forward to <homogenization>
   read(myFile,'(a1024)',END=100) line
 enddo

 do                                                                                                 ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     output = 0_pInt                                                                                ! reset output counter
   endif
   if (section > 0 .and. homogenization_type(section) == homogenization_isostrain_label) then       ! one of my sections
     i = homogenization_typeInstance(section)                                                       ! which instance of my type is present homogenization
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1_pInt
         homogenization_isostrain_output(output,i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('ngrains')
              homogenization_isostrain_Ngrains(i) = IO_intValue(line,positions,2_pInt)
     end select
   endif
 enddo

100 do k = 1,maxNinstance
   homogenization_isostrain_sizeState(i)    = 0_pInt

   do j = 1_pInt,maxval(homogenization_Noutput)
     select case(homogenization_isostrain_output(j,i))
       case('ngrains')
         mySize = 1_pInt
       case default
         mySize = 0_pInt
     end select

     if (mySize > 0_pInt) then                                                                      ! any meaningful output found
       homogenization_isostrain_sizePostResult(j,i) = mySize
       homogenization_isostrain_sizePostResults(i) = &
       homogenization_isostrain_sizePostResults(i) + mySize
     endif
   enddo
 enddo

end subroutine homogenization_isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief sets the initial homogenization stated
!--------------------------------------------------------------------------------------------------
function homogenization_isostrain_stateInit(myInstance)
 use prec, only: &
   pReal
 
 implicit none
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(homogenization_isostrain_sizeState(myInstance)) :: &
              homogenization_isostrain_stateInit

 homogenization_isostrain_stateInit = 0.0_pReal

end function homogenization_isostrain_stateInit


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_partitionDeformation(F,F0,avgF,state,i,e)
 use prec, only: pReal,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: F                           ! partioned def grad per grain
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0                          ! initial partioned def grad per grain
 real(pReal), dimension (3,3), intent(in) :: avgF                                                   ! my average def grad
 type(p_vec), intent(in) :: state                                                                   ! my state
 integer(pInt), intent(in) :: &
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number

 F = spread(avgF,3,homogenization_Ngrains(mesh_element(3,e)))

end subroutine homogenization_isostrain_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and 
! "happy" with result
!--------------------------------------------------------------------------------------------------
function homogenization_isostrain_updateState(state,P,dPdF,i,e)
 use prec, only: &
   pReal,&
   p_vec
 use material, only: &
   homogenization_maxNgrains
 
 implicit none
 type(p_vec), intent(inout) :: state                                                                !< my state 
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P                            !< array of current grain stresses
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF                     !< array of current grain stiffnesses
 integer(pInt), intent(in) :: &
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 logical, dimension(2) :: homogenization_isostrain_updateState

 homogenization_isostrain_updateState = .true.                                                      ! homogenization at material point converged (done and happy)
 
end function homogenization_isostrain_updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,i,e)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: homogenization_maxNgrains, homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (3,3), intent(out) :: avgP                                                  !< average stress at material point
 real(pReal), dimension (3,3,3,3), intent(out) :: dAvgPdAvgF                                        !< average stiffness at material point
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P                            !< array of current grain stresses
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF                     !< array of current grain stiffnesses
 integer(pInt), intent(in) :: &
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 integer(pInt) :: Ngrains

 Ngrains = homogenization_Ngrains(mesh_element(3,e))
 avgP = sum(P,3)/real(Ngrains,pReal)
 dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)

end subroutine homogenization_isostrain_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief derive average temperature from constituent quantities 
!--------------------------------------------------------------------------------------------------
real(pReal) pure function homogenization_isostrain_averageTemperature(Temperature,i,e)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (homogenization_maxNgrains), intent(in) :: Temperature
 integer(pInt), intent(in) :: &
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 integer(pInt) :: Ngrains

 Ngrains = homogenization_Ngrains(mesh_element(3,e))
 homogenization_isostrain_averageTemperature = sum(Temperature(1:Ngrains))/real(Ngrains,pReal)

end function homogenization_isostrain_averageTemperature


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_isostrain_postResults(state,i,e)
 use prec, only: &
   pReal,&
   p_vec
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_typeInstance, &
   homogenization_Noutput
 
 implicit none
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: &
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 integer(pInt) :: homID,o,c
 real(pReal), dimension(homogenization_isostrain_sizePostResults&
        (homogenization_typeInstance(mesh_element(3,e)))) :: homogenization_isostrain_postResults

 c = 0_pInt
 homID = homogenization_typeInstance(mesh_element(3,e))
 homogenization_isostrain_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,e))
   select case(homogenization_isostrain_output(o,homID))
     case ('ngrains')
       homogenization_isostrain_postResults(c+1_pInt) = real(homogenization_isostrain_Ngrains(homID),pReal)
       c = c + 1_pInt
   end select
 enddo

end function homogenization_isostrain_postResults

end module homogenization_isostrain
