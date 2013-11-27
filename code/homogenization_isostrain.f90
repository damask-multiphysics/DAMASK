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
 integer(pInt),        dimension(:),    allocatable,         public, protected :: &
   homogenization_isostrain_sizeState, &
   homogenization_isostrain_sizePostResults
 integer(pInt),        dimension(:,:),  allocatable, target, public :: &
   homogenization_isostrain_sizePostResult
 
 character(len=64), dimension(:,:), allocatable, target, public :: &
  homogenization_isostrain_output                                                                   !< name of each post result output
 character(len=64), dimension(:),   allocatable,         private :: &
   homogenization_isostrain_mapping
 integer(pInt),        dimension(:),    allocatable,         private :: &
   homogenization_isostrain_Ngrains

 public :: &
   homogenization_isostrain_init, &
   homogenization_isostrain_partitionDeformation, &
   homogenization_isostrain_averageStressAndItsTangent, &
   homogenization_isostrain_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_init(myUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_Mandel3333to66, &
   math_Voigt66to3333
 use IO
 use material
 
 implicit none
 integer(pInt),                                      intent(in) :: myUnit
 integer(pInt),                                      parameter  :: MAXNCHUNKS = 2_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS)             :: positions
 integer(pInt) :: &
   section, i, j, output, mySize
 integer :: &
   maxNinstance, k                                                                                  ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''                                                                                        ! to start initialized
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_ISOSTRAIN_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)
 if (maxNinstance == 0) return

 allocate(homogenization_isostrain_sizeState(maxNinstance))
          homogenization_isostrain_sizeState = 0_pInt
 allocate(homogenization_isostrain_sizePostResults(maxNinstance))
          homogenization_isostrain_sizePostResults = 0_pInt
 allocate(homogenization_isostrain_sizePostResult(maxval(homogenization_Noutput),maxNinstance))
          homogenization_isostrain_sizePostResult = 0_pInt
 allocate(homogenization_isostrain_Ngrains(maxNinstance))
          homogenization_isostrain_Ngrains = 0_pInt
 allocate(homogenization_isostrain_mapping(maxNinstance))
          homogenization_isostrain_mapping = 'avg'
 allocate(homogenization_isostrain_output(maxval(homogenization_Noutput),maxNinstance))
          homogenization_isostrain_output = ''
 
 rewind(myUnit)
 section = 0_pInt
 
 do while (trim(line) /= '#EOF#' .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization) ! wind forward to <homogenization>
   line = IO_read(myUnit)
 enddo

 do while (trim(line) /= '#EOF#')
   line = IO_read(myUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     output = 0_pInt                                                                                ! reset output counter
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (homogenization_type(section) == HOMOGENIZATION_ISOSTRAIN_ID) then                          ! one of my sections
       i = homogenization_typeInstance(section)                                                     ! which instance of my type is present homogenization
       positions = IO_stringPos(line,MAXNCHUNKS)
       tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                           ! extract key
       select case(tag)
         case ('(output)')
           output = output + 1_pInt
           homogenization_isostrain_output(output,i) = IO_lc(IO_stringValue(line,positions,2_pInt))
         case ('ngrains','ncomponents')
                homogenization_isostrain_Ngrains(i) = IO_intValue(line,positions,2_pInt)
         case ('mapping')
                homogenization_isostrain_mapping(i) = IO_lc(IO_stringValue(line,positions,2_pInt))
       end select
     endif
   endif
 enddo

 do k = 1,maxNinstance
   homogenization_isostrain_sizeState(i)    = 0_pInt

   do j = 1_pInt,maxval(homogenization_Noutput)
     select case(homogenization_isostrain_output(j,i))
       case('ngrains','ncomponents','temperature')
         mySize = 1_pInt
       case('ipcoords')
         mySize = 3_pInt
       case('avgdefgrad','avgf','avgp','avgfirstpiola','avg1stpiola')
         mySize = 9_pInt
       case default
         mySize = 0_pInt
     end select

     outputFound: if (mySize > 0_pInt) then
       homogenization_isostrain_sizePostResult(j,i) = mySize
       homogenization_isostrain_sizePostResults(i) = &
         homogenization_isostrain_sizePostResults(i) + mySize
     endif outputFound
   enddo
 enddo

end subroutine homogenization_isostrain_init


!--------------------------------------------------------------------------------------------------
!> @brief partitions the deformation gradient onto the constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_partitionDeformation(F,avgF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains
 
 implicit none
 real(pReal),   dimension (3,3,homogenization_maxNgrains), intent(out) :: F                         !< partioned def grad per grain
 real(pReal),   dimension (3,3),                           intent(in)  :: avgF                      !< my average def grad
 integer(pInt),                                            intent(in)  :: &
   el                                                                                               !< element number

 F = spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))

end subroutine homogenization_isostrain_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities 
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_averageStressAndItsTangent(avgP,dAvgPdAvgF,P,dPdF,el)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_maxNgrains, &
   homogenization_Ngrains, &
   homogenization_typeInstance
 
 implicit none
 real(pReal),   dimension (3,3),                               intent(out) :: avgP                  !< average stress at material point
 real(pReal),   dimension (3,3,3,3),                           intent(out) :: dAvgPdAvgF            !< average stiffness at material point
 real(pReal),   dimension (3,3,homogenization_maxNgrains),     intent(in)  :: P                     !< array of current grain stresses
 real(pReal),   dimension (3,3,3,3,homogenization_maxNgrains), intent(in)  :: dPdF                  !< array of current grain stiffnesses
 integer(pInt),                                                intent(in)  :: el                    !< element number
 integer(pInt) :: &
   homID, & 
   Ngrains

 homID = homogenization_typeInstance(mesh_element(3,el))
 Ngrains = homogenization_Ngrains(mesh_element(3,el))

 select case (homogenization_isostrain_mapping(homID))
   case ('parallel','sum')
     avgP       = sum(P,3)
     dAvgPdAvgF = sum(dPdF,5)
   case ('average','mean','avg')
     avgP       = sum(P,3)   /real(Ngrains,pReal)
     dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)
   case default
     avgP       = sum(P,3)   /real(Ngrains,pReal)
     dAvgPdAvgF = sum(dPdF,5)/real(Ngrains,pReal)
 end select

end subroutine homogenization_isostrain_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion 
!--------------------------------------------------------------------------------------------------
pure function homogenization_isostrain_postResults(ip,el,avgP,avgF)
 use prec, only: &
   pReal
 use mesh, only: &
   mesh_element, &
   mesh_ipCoordinates
 use material, only: &
   homogenization_typeInstance, &
   homogenization_Noutput
 use crystallite, only: &
   crystallite_temperature
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 real(pReal), dimension(3,3), intent(in) :: &
   avgP, &                                                                                          !< average stress at material point
   avgF                                                                                             !< average deformation gradient at material point
 real(pReal),  dimension(homogenization_isostrain_sizePostResults &
                         (homogenization_typeInstance(mesh_element(3,el)))) :: &
   homogenization_isostrain_postResults
 
 integer(pInt) :: &
   homID, &
   o, c
   
 c = 0_pInt
 homID = homogenization_typeInstance(mesh_element(3,el))
 homogenization_isostrain_postResults = 0.0_pReal
 
 do o = 1_pInt,homogenization_Noutput(mesh_element(3,el))
   select case(homogenization_isostrain_output(o,homID))
     case ('ngrains','ncomponents')
       homogenization_isostrain_postResults(c+1_pInt) = real(homogenization_isostrain_Ngrains(homID),pReal)
       c = c + 1_pInt
     case ('temperature')
       homogenization_isostrain_postResults(c+1_pInt) = crystallite_temperature(ip,el)
       c = c + 1_pInt
     case ('avgdefgrad','avgf')
       homogenization_isostrain_postResults(c+1_pInt:c+9_pInt) = reshape(avgF,[9])
       c = c + 9_pInt
     case ('avgp','avgfirstpiola','avg1stpiola')
       homogenization_isostrain_postResults(c+1_pInt:c+9_pInt) = reshape(avgP,[9])
       c = c + 9_pInt
     case ('ipcoords')
       homogenization_isostrain_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_isostrain_postResults

end module homogenization_isostrain
