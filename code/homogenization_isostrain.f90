! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
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
!##############################################################
!* $Id$
!*****************************************************
!*      Module: HOMOGENIZATION_ISOSTRAIN             *
!*****************************************************
!* contains:                                         *
!*****************************************************

! [isostrain]
! type            isostrain
! Ngrains         6
!   (output)        Ngrains

module homogenization_isostrain

 use prec, only: pInt
 
 implicit none
 character (len=*), parameter :: &
   homogenization_isostrain_label = 'isostrain'
 
 integer(pInt),dimension(:), allocatable :: &
   homogenization_isostrain_sizeState, &
   homogenization_isostrain_Ngrains, &
   homogenization_isostrain_sizePostResults
 
 integer(pInt), dimension(:,:), allocatable, target :: &
   homogenization_isostrain_sizePostResult
 
 character(len=64), dimension(:,:), allocatable, target :: &
  homogenization_isostrain_output             ! name of each post result output


contains
!****************************************
!* - homogenization_isostrain_init
!* - homogenization_isostrain_stateInit
!* - homogenization_isostrain_deformationPartititon
!* - homogenization_isostrain_stateUpdate
!* - homogenization_isostrain_averageStressAndItsTangent
!* - homogenization_isostrain_postResults
!****************************************


!**************************************
!*      Module initialization         *
!**************************************
subroutine homogenization_isostrain_init(myFile)                           ! file pointer to material configuration
 use, intrinsic :: iso_fortran_env                                ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: pInt
 use math, only: math_Mandel3333to66, math_Voigt66to3333
 use IO
 use material
 integer(pInt), intent(in) :: myFile
 integer(pInt), parameter :: maxNchunks = 2_pInt
 integer(pInt), dimension(1_pInt+2_pInt*maxNchunks) :: positions
 integer(pInt) section, i, j, output, mySize
 integer :: maxNinstance, k                                                                         !no pInt (stores a system dependen value from 'count'
 character(len=64)   :: tag
 character(len=1024) :: line
 
 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  homogenization_',trim(homogenization_isostrain_label),' init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
 !$OMP END CRITICAL (write2out)

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
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)     ! wind forward to <homogenization>
   read(myFile,'(a1024)',END=100) line
 enddo

 do                                                       ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                ! next section
     section = section + 1_pInt
     output = 0_pInt                                           ! reset output counter
   endif
   if (section > 0 .and. homogenization_type(section) == homogenization_isostrain_label) then  ! one of my sections
     i = homogenization_typeInstance(section)             ! which instance of my type is present homogenization
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))        ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1_pInt
         homogenization_isostrain_output(output,i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       case ('ngrains')
              homogenization_isostrain_Ngrains(i) = IO_intValue(line,positions,2_pInt)
     end select
   endif
 enddo

100 do k = 1,maxNinstance                                        ! sanity checks
 enddo

 do k = 1,maxNinstance
   homogenization_isostrain_sizeState(i)    = 0_pInt

   do j = 1_pInt,maxval(homogenization_Noutput)
     select case(homogenization_isostrain_output(j,i))
       case('ngrains')
         mySize = 1_pInt
       case default
         mySize = 0_pInt
     end select

     if (mySize > 0_pInt) then                               ! any meaningful output found
       homogenization_isostrain_sizePostResult(j,i) = mySize
       homogenization_isostrain_sizePostResults(i) = &
       homogenization_isostrain_sizePostResults(i) + mySize
     endif
   enddo
 enddo

end subroutine homogenization_isostrain_init


!*********************************************************************
!* initial homogenization state                                      *
!*********************************************************************
function homogenization_isostrain_stateInit(myInstance)
 use prec, only: pReal,pInt
 
 implicit none
 integer(pInt), intent(in) :: myInstance
 real(pReal), dimension(homogenization_isostrain_sizeState(myInstance)) :: &
              homogenization_isostrain_stateInit

 homogenization_isostrain_stateInit = 0.0_pReal

endfunction homogenization_isostrain_stateInit


!********************************************************************
! partition material point def grad onto constituents
!********************************************************************
subroutine homogenization_isostrain_partitionDeformation(&
   F, &             ! partioned def grad per grain
!
   F0, &            ! initial partioned def grad per grain
   avgF, &          ! my average def grad
   state, &         ! my state
   ip, &            ! my integration point
   el  &            ! my element
  )
 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains,homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(out) :: F
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in)  :: F0
 real(pReal), dimension (3,3), intent(in) :: avgF
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: ip,el
 integer(pInt) i
 
! homID = homogenization_typeInstance(mesh_element(3,el))
 forall (i = 1_pInt:homogenization_Ngrains(mesh_element(3,el))) &
   F(1:3,1:3,i)  = avgF

end subroutine homogenization_isostrain_partitionDeformation


!********************************************************************
! update the internal state of the homogenization scheme
! and tell whether "done" and "happy" with result
!********************************************************************
function homogenization_isostrain_updateState(&
   state, &         ! my state
!
   P, &             ! array of current grain stresses
   dPdF, &          ! array of current grain stiffnesses
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use material, only: homogenization_maxNgrains
 implicit none

!* Definition of variables
 type(p_vec), intent(inout) :: state
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 integer(pInt), intent(in) :: ip,el
! integer(pInt) homID
 logical, dimension(2) :: homogenization_isostrain_updateState

! homID = homogenization_typeInstance(mesh_element(3,el))
 homogenization_isostrain_updateState = .true.                      ! homogenization at material point converged (done and happy)
 
end function homogenization_isostrain_updateState


!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
subroutine homogenization_isostrain_averageStressAndItsTangent(&
   avgP, &          ! average stress at material point
   dAvgPdAvgF, &    ! average stiffness at material point
!
   P, &             ! array of current grain stresses
   dPdF, &          ! array of current grain stiffnesses
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains, homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (3,3), intent(out) :: avgP
 real(pReal), dimension (3,3,3,3), intent(out) :: dAvgPdAvgF
 real(pReal), dimension (3,3,homogenization_maxNgrains), intent(in) :: P
 real(pReal), dimension (3,3,3,3,homogenization_maxNgrains), intent(in) :: dPdF
 integer(pInt), intent(in) :: ip,el
 integer(pInt) Ngrains

! homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 avgP = sum(P,3)/real(Ngrains,pReal)
 dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)

end subroutine homogenization_isostrain_averageStressAndItsTangent


!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
function homogenization_isostrain_averageTemperature(&
   Temperature, &   ! temperature
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_maxNgrains, homogenization_Ngrains
 
 implicit none
 real(pReal), dimension (homogenization_maxNgrains), intent(in) :: Temperature
 integer(pInt), intent(in) :: ip,el
 real(pReal) homogenization_isostrain_averageTemperature
 integer(pInt) Ngrains

! homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 homogenization_isostrain_averageTemperature = sum(Temperature(1:Ngrains))/real(Ngrains,pReal)

end function homogenization_isostrain_averageTemperature


!********************************************************************
! return array of homogenization results for post file inclusion
!********************************************************************
pure function homogenization_isostrain_postResults(&
   state, &         ! my state
   ip, &            ! my integration point
   el  &            ! my element
  )

 use prec, only: pReal,pInt,p_vec
 use mesh, only: mesh_element
 use material, only: homogenization_typeInstance,homogenization_Noutput
 
 implicit none
 type(p_vec), intent(in) :: state
 integer(pInt), intent(in) :: ip,el
 integer(pInt) :: homID,o,c
 real(pReal), dimension(homogenization_isostrain_sizePostResults&
        (homogenization_typeInstance(mesh_element(3,el)))) :: homogenization_isostrain_postResults
 c = 0_pInt
 homID = homogenization_typeInstance(mesh_element(3,el))
 homogenization_isostrain_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_isostrain_output(o,homID))
     case ('ngrains')
       homogenization_isostrain_postResults(c+1_pInt) = real(homogenization_isostrain_Ngrains(homID),pReal)
       c = c + 1_pInt
   end select
 enddo
 
 return

end function homogenization_isostrain_postResults

end module homogenization_isostrain
