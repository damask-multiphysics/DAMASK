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
 integer(pInt),               dimension(:),   allocatable,         public, protected :: &
   homogenization_isostrain_sizePostResults
 integer(pInt),               dimension(:,:), allocatable, target, public :: &
   homogenization_isostrain_sizePostResult
 
 character(len=64),           dimension(:,:), allocatable, target, public :: &
  homogenization_isostrain_output                                                                   !< name of each post result output
 integer(pInt),               dimension(:),   allocatable, target, public :: &
   homogenization_isostrain_Noutput                                                                 !< number of outputs per homog instance
 integer(pInt),               dimension(:),   allocatable,         private :: &
   homogenization_isostrain_Ngrains
 enum, bind(c) 
   enumerator :: undefined_ID, &
                 nconstituents_ID, &
                 ipcoords_ID, &
                 avgdefgrad_ID, &
                 avgfirstpiola_ID
 end enum
 enum, bind(c) 
   enumerator :: parallel_ID, &
                 average_ID
 end enum
 integer(kind(undefined_ID)), dimension(:,:), allocatable,         private :: &
  homogenization_isostrain_outputID                                                                 !< ID of each post result output
 integer(kind(average_ID)),   dimension(:),   allocatable,         private :: &
  homogenization_isostrain_mapping                                                                  !< mapping type


 public :: &
   homogenization_isostrain_init, &
   homogenization_isostrain_partitionDeformation, &
   homogenization_isostrain_averageStressAndItsTangent, &
   homogenization_isostrain_postResults

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields, reads information from material configuration file
!--------------------------------------------------------------------------------------------------
subroutine homogenization_isostrain_init(fileUnit)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pReal
 use debug, only: &
   debug_HOMOGENIZATION, &
   debug_level, &
   debug_levelBasic
 use IO
 use material
 
 implicit none
 integer(pInt),                                      intent(in) :: fileUnit
 integer(pInt),                                      parameter  :: MAXNCHUNKS = 2_pInt
 integer(pInt), dimension(1_pInt+2_pInt*MAXNCHUNKS)             :: positions
 integer(pInt) :: &
   section = 0_pInt, i, mySize, o
 integer :: &
   maxNinstance, &
   homog, &
   instance
 integer :: &
   NofMyHomog                                                                                       ! no pInt (stores a system dependen value from 'count'
 character(len=65536) :: &
   tag  = '', &
   line = ''
 
 write(6,'(/,a)')   ' <<<+-  homogenization_'//HOMOGENIZATION_ISOSTRAIN_label//' init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 maxNinstance = count(homogenization_type == HOMOGENIZATION_ISOSTRAIN_ID)
 if (maxNinstance == 0) return
 
 if (iand(debug_level(debug_HOMOGENIZATION),debug_levelBasic) /= 0_pInt) &
   write(6,'(a16,1x,i5,/)') '# instances:',maxNinstance
 allocate(homogenization_isostrain_sizePostResults(maxNinstance),          source=0_pInt)
 allocate(homogenization_isostrain_sizePostResult(maxval(homogenization_Noutput),maxNinstance), &
                                                                           source=0_pInt)
 allocate(homogenization_isostrain_Noutput(maxNinstance),                  source=0_pInt)
 allocate(homogenization_isostrain_Ngrains(maxNinstance),                  source=0_pInt)
 allocate(homogenization_isostrain_mapping(maxNinstance),                  source=average_ID)
 allocate(homogenization_isostrain_output(maxval(homogenization_Noutput),maxNinstance))
          homogenization_isostrain_output = ''
 allocate(homogenization_isostrain_outputID(maxval(homogenization_Noutput),maxNinstance), &
                                                                           source=undefined_ID)

 rewind(fileUnit)
 do while (trim(line) /= IO_EOF .and. IO_lc(IO_getTag(line,'<','>')) /= material_partHomogenization)! wind forward to <homogenization>
   line = IO_read(fileUnit)
 enddo

 parsingFile: do while (trim(line) /= IO_EOF)                                                       ! read through sections of homogenization part
   line = IO_read(fileUnit)
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') then                                                          ! stop at next part
     line = IO_read(fileUnit, .true.)                                                               ! reset IO_read
     exit                                                                                           
   endif
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     cycle
   endif
   if (section > 0_pInt ) then                                                                      ! do not short-circuit here (.and. with next if-statement). It's not safe in Fortran
     if (homogenization_type(section) == HOMOGENIZATION_ISOSTRAIN_ID) then                          ! one of my sections
       i = homogenization_typeInstance(section)                                                     ! which instance of my type is present homogenization
       positions = IO_stringPos(line,MAXNCHUNKS)
       tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                           ! extract key
       select case(tag)
         case ('(output)')
           select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
             case('nconstituents','ngrains')
               homogenization_isostrain_Noutput(i) = homogenization_isostrain_Noutput(i) + 1_pInt
               homogenization_isostrain_outputID(homogenization_isostrain_Noutput(i),i) = nconstituents_ID
               homogenization_isostrain_output(homogenization_isostrain_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,positions,2_pInt))
             case('ipcoords')
               homogenization_isostrain_Noutput(i) = homogenization_isostrain_Noutput(i) + 1_pInt
               homogenization_isostrain_outputID(homogenization_isostrain_Noutput(i),i) = ipcoords_ID
               homogenization_isostrain_output(homogenization_isostrain_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,positions,2_pInt))
             case('avgdefgrad','avgf')
               homogenization_isostrain_Noutput(i) = homogenization_isostrain_Noutput(i) + 1_pInt
               homogenization_isostrain_outputID(homogenization_isostrain_Noutput(i),i) = avgdefgrad_ID
               homogenization_isostrain_output(homogenization_isostrain_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,positions,2_pInt))
             case('avgp','avgfirstpiola','avg1stpiola')
               homogenization_isostrain_Noutput(i) = homogenization_isostrain_Noutput(i) + 1_pInt
               homogenization_isostrain_outputID(homogenization_isostrain_Noutput(i),i) = avgfirstpiola_ID
               homogenization_isostrain_output(homogenization_isostrain_Noutput(i),i) = &
                 IO_lc(IO_stringValue(line,positions,2_pInt))

           end select
         case ('nconstituents','ngrains')
           homogenization_isostrain_Ngrains(i) = IO_intValue(line,positions,2_pInt)
         case ('mapping')
           select case(IO_lc(IO_stringValue(line,positions,2_pInt)))
             case ('parallel','sum')
               homogenization_isostrain_mapping(i) = parallel_ID
             case ('average','mean','avg')
               homogenization_isostrain_mapping(i) = average_ID
             case default
               call IO_error(211_pInt,ext_msg=trim(tag)//' ('//HOMOGENIZATION_isostrain_label//')')
           end select

       end select
     endif
   endif
 enddo parsingFile

 initializeInstances: do homog = 1_pInt, material_Nhomogenization
   myHomog: if (homogenization_type(homog) == HOMOGENIZATION_ISOSTRAIN_ID) then
     NofMyHomog = count(material_homog == homog)
     instance = homogenization_typeInstance(homog)

! *  Determine size of postResults array
     outputsLoop: do o = 1_pInt, homogenization_isostrain_Noutput(instance)
       select case(homogenization_isostrain_outputID(o,instance))
        case(nconstituents_ID)
          mySize = 1_pInt
        case(ipcoords_ID)
          mySize = 3_pInt
        case(avgdefgrad_ID, avgfirstpiola_ID)
          mySize = 9_pInt
        case default
          mySize = 0_pInt
       end select

       outputFound: if (mySize > 0_pInt) then
        homogenization_isostrain_sizePostResult(o,instance) = mySize
        homogenization_isostrain_sizePostResults(instance) = &
          homogenization_isostrain_sizePostResults(instance) + mySize
       endif outputFound
     enddo outputsLoop

! allocate state arrays
     homogState(homog)%sizeState = 0_pInt
     homogState(homog)%sizePostResults = homogenization_isostrain_sizePostResults(instance)
     allocate(homogState(homog)%state0   (0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%subState0(0_pInt,NofMyHomog), source=0.0_pReal)
     allocate(homogState(homog)%state    (0_pInt,NofMyHomog), source=0.0_pReal)

   endif myHomog
 enddo initializeInstances

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
 F=0.0_pReal
 F(1:3,1:3,1:homogenization_Ngrains(mesh_element(3,el)))= &
   spread(avgF,3,homogenization_Ngrains(mesh_element(3,el)))

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
   case (parallel_ID)
     avgP       = sum(P,3)
     dAvgPdAvgF = sum(dPdF,5)
   case (average_ID)
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
   select case(homogenization_isostrain_outputID(o,homID))
     case (nconstituents_ID)
       homogenization_isostrain_postResults(c+1_pInt) = real(homogenization_isostrain_Ngrains(homID),pReal)
       c = c + 1_pInt
     case (avgdefgrad_ID)
       homogenization_isostrain_postResults(c+1_pInt:c+9_pInt) = reshape(avgF,[9])
       c = c + 9_pInt
     case (avgfirstpiola_ID)
       homogenization_isostrain_postResults(c+1_pInt:c+9_pInt) = reshape(avgP,[9])
       c = c + 9_pInt
     case (ipcoords_ID)
       homogenization_isostrain_postResults(c+1_pInt:c+3_pInt) = mesh_ipCoordinates(1:3,ip,el)                       ! current ip coordinates
       c = c + 3_pInt
    end select
 enddo

end function homogenization_isostrain_postResults

end module homogenization_isostrain
