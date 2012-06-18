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
!##############################################################
module IO   
!##############################################################
 use prec, only: pInt, pReal
 
 implicit none
 private
 public ::  IO_init, &
            IO_checkAndRewind, &
            IO_open_file_stat, &
            IO_open_jobFile_stat, &
            IO_open_file, &
            IO_open_jobFile, &
            IO_write_jobFile, &
            IO_write_jobBinaryFile, &
            IO_read_jobBinaryFile, &
            IO_hybridIA, &
            IO_isBlank, &
            IO_getTag, &
            IO_countSections, &
            IO_countTagInPart, &
            IO_spotTagInPart, &
            IO_stringPos, &
            IO_stringValue, &
            IO_fixedStringValue ,&
            IO_floatValue, &
            IO_fixedNoEFloatValue, &
            IO_intValue, &
            IO_fixedIntValue, &
            IO_lc, &
            IO_skipChunks, &
            IO_extractValue, &
            IO_countDataLines, &
            IO_countContinuousIntValues, &
            IO_continuousIntValues, &
            IO_error, &
            IO_warning
#ifndef Spectral
 public ::  IO_open_inputFile, &
            IO_open_logFile
#endif
#ifdef Abaqus  
 public ::  IO_abaqus_hasNoPart
#endif

 private :: IO_fixedFloatValue, &
            IO_lcInplace ,&
            hybridIA_reps
#ifdef Abaqus 
 private :: abaqus_assembleInputFile
#endif

contains


!--------------------------------------------------------------------------------------------------
!> @brief only output of revision number
!--------------------------------------------------------------------------------------------------
subroutine IO_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)

  write(6,*)
  write(6,*) '<<<+-  IO init  -+>>>'
  write(6,*) '$Id$'
#include "compilation_info.f90"
  flush(6)

end subroutine IO_init


!--------------------------------------------------------------------------------------------------
!> @brief Checks if unit is opened for reading, if true rewinds. Otherwise stops with
!! error message 102
!--------------------------------------------------------------------------------------------------
subroutine IO_checkAndRewind(myUnit)
 
implicit none
 integer(pInt), intent(in) :: myUnit
 logical :: fileOpened
 character(len=15) :: fileRead
 inquire(unit=myUnit, opened=fileOpened, read = fileRead) 
 if (fileOpened .neqv. .true. .or. trim(fileRead)/='YES') call IO_error(102_pInt)
 rewind(myUnit)

end subroutine IO_checkAndRewind


!--------------------------------------------------------------------------------------------------
!> @brief Open existing file to given unit path to file is relative to working directory
!--------------------------------------------------------------------------------------------------
logical function IO_open_file_stat(myUnit,relPath)
 use DAMASK_interface, &
   only: getSolverWorkingDirectoryName
 
 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: relPath

 integer(pInt)                  :: myStat
 character(len=1024)            :: path
 
 path = trim(getSolverWorkingDirectoryName())//relPath
 open(myUnit,status='old',iostat=myStat,file=path)
 IO_open_file_stat = (myStat == 0_pInt)
 
end function IO_open_file_stat


!--------------------------------------------------------------------------------------------------
!> @brief Open (write) file related to current job but with different extension to given unit
!--------------------------------------------------------------------------------------------------
logical function IO_open_jobFile_stat(myUnit,newExt)
 use DAMASK_interface, only: &
   getSolverWorkingDirectoryName, &
   getSolverJobName

 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: newExt

 integer(pInt)                  :: myStat
 character(len=1024)            :: path

 path = trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.'//newExt
 open(myUnit,status='old',iostat=myStat,file=path)
 IO_open_jobFile_stat = (myStat == 0_pInt)

end function IO_open_JobFile_stat


!--------------------------------------------------------------------------------------------------
!> @brief Open existing file to given unit path to file is relative to working directory
!--------------------------------------------------------------------------------------------------
subroutine IO_open_file(myUnit,relPath)

 use DAMASK_interface,       only: getSolverWorkingDirectoryName
 
 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: relPath

 integer(pInt)                  :: myStat
 character(len=1024)            :: path
 
 path = trim(getSolverWorkingDirectoryName())//relPath
 open(myUnit,status='old',iostat=myStat,file=path)
 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
 
end subroutine IO_open_file


!--------------------------------------------------------------------------------------------------
!> @brief Open (write) file related to current job but with different extension to given unit
!--------------------------------------------------------------------------------------------------
subroutine IO_open_jobFile(myUnit,newExt)

 use DAMASK_interface,       only: getSolverWorkingDirectoryName, &
                                   getSolverJobName

 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: newExt

 integer(pInt)                  :: myStat
 character(len=1024)            :: path


 path = trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.'//newExt
 open(myUnit,status='old',iostat=myStat,file=path)
 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
 
end subroutine IO_open_jobFile

#ifndef Spectral
!********************************************************************
! open FEM inputfile to given myUnit
! AP: 12.07.10 
!   : changed the function to open *.inp_assembly, which is basically 
!     the input file without comment lines and possibly assembled includes
!********************************************************************
subroutine IO_open_inputFile(myUnit,model)

 use DAMASK_interface, only: &
   getSolverWorkingDirectoryName,&
   getSolverJobName, &
   inputFileExtension

 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: model

 integer(pInt)                  :: myStat
 character(len=1024)            :: path
 
#ifdef Abaqus
 path = trim(getSolverWorkingDirectoryName())//trim(model)//InputFileExtension
 open(myUnit+1,status='old',iostat=myStat,file=path)
 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
   
 path = trim(getSolverWorkingDirectoryName())//trim(model)//InputFileExtension//'_assembly'
 open(myUnit,iostat=myStat,file=path)
 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
    if (.not.abaqus_assembleInputFile(myUnit,myUnit+1_pInt)) call IO_error(103_pInt)     ! strip comments and concatenate any "include"s
 close(myUnit+1_pInt) 
#endif
#ifdef Marc
   path = trim(getSolverWorkingDirectoryName())//trim(model)//InputFileExtension
   open(myUnit,status='old',iostat=myStat,file=path)
   if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
#endif

end subroutine IO_open_inputFile


!********************************************************************
! open FEM logfile to given myUnit
!********************************************************************
subroutine IO_open_logFile(myUnit)

 use DAMASK_interface, only: &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   LogFileExtension

 implicit none
 integer(pInt),      intent(in) :: myUnit

 integer(pInt)                  :: myStat
 character(len=1024)            :: path

 path = trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//LogFileExtension
 open(myUnit,status='old',iostat=myStat,file=path)
 if (myStat /= 0) call IO_error(100_pInt,ext_msg=path)

end subroutine IO_open_logFile
#endif

!********************************************************************
! open (write) file related to current job
! but with different extension to given myUnit
!********************************************************************
subroutine IO_write_jobFile(myUnit,newExt)

 use DAMASK_interface,       only: getSolverWorkingDirectoryName,&
                                   getSolverJobName
 
 implicit none
 integer(pInt),      intent(in) :: myUnit
 character(len=*),   intent(in) :: newExt

 integer(pInt)                  :: myStat
 character(len=1024)            :: path

 path = trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.'//newExt
 open(myUnit,status='replace',iostat=myStat,file=path)
 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
 
end subroutine IO_write_jobFile


!********************************************************************
! open (write) binary file related to current job
! but with different extension to given myUnit
!********************************************************************
subroutine IO_write_jobBinaryFile(myUnit,newExt,recMultiplier)

 use DAMASK_interface,                 only: getSolverWorkingDirectoryName, &
                                             getSolverJobName
 
 implicit none
 integer(pInt),      intent(in)           :: myUnit
 integer(pInt),      intent(in), optional :: recMultiplier
 character(len=*),   intent(in)           :: newExt

 integer(pInt)                            :: myStat
 character(len=1024)                      :: path

 path = trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.'//newExt
 if (present(recMultiplier)) then
   open(myUnit,status='replace',form='unformatted',access='direct', &
                                                   recl=pReal*recMultiplier,iostat=myStat,file=path)
 else
   open(myUnit,status='replace',form='unformatted',access='direct', &
                                                   recl=pReal,iostat=myStat,file=path)
 endif

 if (myStat /= 0_pInt) call IO_error(100_pInt,ext_msg=path)
 
end subroutine IO_write_jobBinaryFile


!********************************************************************
! open (read) binary file related to restored job
! and with different extension to given myUnit
!********************************************************************
subroutine IO_read_jobBinaryFile(myUnit,newExt,jobName,recMultiplier)

 use DAMASK_interface,                 only: getSolverWorkingDirectoryName
 
 implicit none
 integer(pInt),      intent(in)           :: myUnit
 integer(pInt),      intent(in), optional :: recMultiplier
 character(len=*),   intent(in)           :: newExt, jobName

 integer(pInt)                            :: myStat
 character(len=1024)                      :: path

 path = trim(getSolverWorkingDirectoryName())//trim(jobName)//'.'//newExt
 if (present(recMultiplier)) then
   open(myUnit,status='old',form='unformatted',access='direct', & 
                                               recl=pReal*recMultiplier,iostat=myStat,file=path)
 else
   open(myUnit,status='old',form='unformatted',access='direct', &
                                               recl=pReal,iostat=myStat,file=path)
 endif
 if (myStat /= 0) call IO_error(100_pInt,ext_msg=path)
 
end subroutine IO_read_jobBinaryFile

#ifdef Abaqus
!***********************************************************
! check if the input file for Abaqus contains part info
!***********************************************************
logical function IO_abaqus_hasNoPart(myUnit)

 implicit none
 integer(pInt),    intent(in)                :: myUnit

 integer(pInt),    parameter                 :: maxNchunks = 1_pInt
 integer(pInt),    dimension(1+2*maxNchunks) :: myPos
 character(len=300)                          :: line
 
 IO_abaqus_hasNoPart = .true.
 
610 FORMAT(A300)
 rewind(myUnit)
 do
   read(myUnit,610,END=620) line
   myPos = IO_stringPos(line,maxNchunks)
   if (IO_lc(IO_stringValue(line,myPos,1_pInt)) == '*part' ) then
     IO_abaqus_hasNoPart = .false.
     exit
   endif
 enddo
 
620 end function IO_abaqus_hasNoPart
#endif

!********************************************************************
! hybrid IA sampling of ODFfile
!********************************************************************
function IO_hybridIA(Nast,ODFfileName)

 implicit none
 integer(pInt),    intent(in)   :: Nast
 real(pReal), dimension(3,Nast) :: IO_hybridIA

 character(len=*), intent(in)   :: ODFfileName

 real(pReal),      parameter  :: PI = 3.14159265358979323846264338327950288419716939937510_pReal
 real(pReal),      parameter  :: INRAD = PI/180.0_pReal
 character(len=*), parameter  :: fileFormat = '(A80)'

 integer(pInt) :: i,j,bin,NnonZero,Nset,Nreps,reps,phi1,Phi,phi2
 integer(pInt), dimension(7)                :: myPos
 integer(pInt), dimension(3)                :: steps
 integer(pInt), dimension(:), allocatable   :: binSet
 real(pReal) :: center,sum_dV_V,prob,dg_0,C,lowerC,upperC,rnd
 real(pReal), dimension(3)                  :: limits, &
                                               deltas

 real(pReal), dimension(:,:,:), allocatable :: dV_V
 character(len=80) :: line
 
 call IO_open_file(999_pInt,ODFfileName)
 
!--- parse header of ODF file ---
!--- limits in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 myPos = IO_stringPos(line,3_pInt)
 if (myPos(1).ne.3) goto 100
 do i=1_pInt,3_pInt
   limits(i) = IO_floatValue(line,myPos,i)*INRAD
 enddo

!--- deltas in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 myPos = IO_stringPos(line,3_pInt)
 if (myPos(1).ne.3) goto 100
 do i=1_pInt,3_pInt
   deltas(i) = IO_floatValue(line,myPos,i)*INRAD
 enddo
 steps = nint(limits/deltas,pInt)
 allocate(dV_V(steps(3),steps(2),steps(1)))

!--- box boundary/center at origin? ---
 read(999,fmt=fileFormat,end=100) line
 if (index(IO_lc(line),'bound')>0) then
   center = 0.5_pReal
 else
   center = 0.0_pReal
 endif
 
!--- skip blank line ---
 read(999,fmt=fileFormat,end=100) line

 sum_dV_V = 0.0_pReal
 dV_V = 0.0_pReal
 dg_0 = deltas(1)*deltas(3)*2.0_pReal*sin(deltas(2)/2.0_pReal)
 NnonZero = 0_pInt
 
 do phi1=1_pInt,steps(1)
   do Phi=1_pInt,steps(2)
     do phi2=1_pInt,steps(3)
       read(999,fmt=*,end=100) prob
       if (prob > 0.0_pReal) then
         NnonZero = NnonZero+1_pInt
         sum_dV_V = sum_dV_V+prob
       else
         prob = 0.0_pReal
       endif
       dV_V(phi2,Phi,phi1) = prob*dg_0*sin((Phi-1.0_pReal+center)*deltas(2))
     enddo
   enddo
 enddo  

 dV_V = dV_V/sum_dV_V  ! normalize to 1
 
!--- now fix bounds ---
 Nset = max(Nast,NnonZero)                             ! if less than non-zero voxel count requested, sample at least that much
 lowerC = 0.0_pReal
 upperC = real(Nset, pReal)
 
 do while (hybridIA_reps(dV_V,steps,upperC) < Nset)
   lowerC = upperC
   upperC = upperC*2.0_pReal
 enddo
!--- binary search for best C ---
 do
   C = (upperC+lowerC)/2.0_pReal
   Nreps = hybridIA_reps(dV_V,steps,C)
   if (abs(upperC-lowerC) < upperC*1.0e-14_pReal) then
     C = upperC
     Nreps = hybridIA_reps(dV_V,steps,C)
     exit
   elseif (Nreps < Nset) then
     lowerC = C
   elseif (Nreps > Nset) then
     upperC = C
   else
     exit
   endif
 enddo

 allocate(binSet(Nreps))
 bin = 0_pInt ! bin counter
 i = 1_pInt ! set counter
 do phi1=1_pInt,steps(1)
   do Phi=1_pInt,steps(2)
     do phi2=1_pInt,steps(3)
       reps = nint(C*dV_V(phi2,Phi,phi1), pInt)
       binSet(i:i+reps-1) = bin
       bin = bin+1_pInt ! advance bin
       i = i+reps ! advance set
     enddo
   enddo
 enddo

 do i=1_pInt,Nast
   if (i < Nast) then
     call random_number(rnd)
     j = nint(rnd*(Nreps-i)+i+0.5_pReal,pInt)
   else
     j = i
   endif
   bin = binSet(j)
   IO_hybridIA(1,i) = deltas(1)*(real(mod(bin/(steps(3)*steps(2)),steps(1)),pReal)+center)  ! phi1
   IO_hybridIA(2,i) = deltas(2)*(real(mod(bin/ steps(3)          ,steps(2)),pReal)+center)  ! Phi
   IO_hybridIA(3,i) = deltas(3)*(real(mod(bin                    ,steps(3)),pReal)+center)  ! phi2
   binSet(j) = binSet(i)
 enddo
 close(999)
 return

! on error
100 IO_hybridIA = -1.0_pReal
 close(999)
 
end function IO_hybridIA


!********************************************************************
! identifies lines without content
!********************************************************************
logical pure function IO_isBlank(line)

 implicit none
 character(len=*), intent(in) :: line

 character(len=*),  parameter :: blankChar = achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 character(len=*),  parameter :: comment = achar(35)                               ! comment id '#'

 integer :: posNonBlank, posComment                                                ! no pInt
 
 posNonBlank = verify(line,blankChar)
 posComment  = scan(line,comment)
 IO_isBlank = posNonBlank == 0 .or. posNonBlank == posComment
 
end function IO_isBlank


!********************************************************************
! get tagged content of line
!********************************************************************
pure function IO_getTag(line,openChar,closeChar)

 implicit none
 character(len=*), intent(in)  :: line
 character(len=len_trim(line)) :: IO_getTag
 
 character(len=*), intent(in)  :: openChar, & 
                                  closeChar

 character(len=*), parameter   :: sep=achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces

 integer :: left,right                                                          ! no pInt

 IO_getTag = ''
 left = scan(line,openChar)
 right = scan(line,closeChar)
 
 if (left == verify(line,sep) .and. right > left) & ! openChar is first and closeChar occurs
   IO_getTag = line(left+1:right-1)

end function IO_getTag

!*********************************************************************
!
!*********************************************************************
integer(pInt) function IO_countSections(myFile,part)

 implicit none
 integer(pInt),      intent(in) :: myFile
 character(len=*),   intent(in) :: part

 character(len=1024)            :: line

 line = ''
 IO_countSections = 0_pInt
 rewind(myFile)

 do while (IO_getTag(line,'<','>') /= part)      ! search for part
   read(myFile,'(a1024)',END=100) line
 enddo

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     IO_countSections = IO_countSections + 1_pInt
 enddo

100 end function IO_countSections
 

!*********************************************************************
! return array of myTag counts within <part> for at most N[sections]
!*********************************************************************
function IO_countTagInPart(myFile,part,myTag,Nsections)

 implicit none
 integer(pInt),   intent(in)                :: Nsections
 integer(pInt),   dimension(Nsections)      :: IO_countTagInPart
 
 integer(pInt),   intent(in)                :: myFile
 character(len=*),intent(in)                :: part, &
                                               myTag

 integer(pInt),   parameter                 :: maxNchunks = 1_pInt

 integer(pInt),   dimension(Nsections)      :: counter
 integer(pInt),   dimension(1+2*maxNchunks) :: positions
 integer(pInt)                              :: section
 character(len=1024)                        :: line, &
                                               tag
 line = ''
 counter = 0_pInt
 section = 0_pInt

 rewind(myFile) 
 do while (IO_getTag(line,'<','>') /= part)               ! search for part
   read(myFile,'(a1024)',END=100) line
 enddo

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     section = section + 1_pInt
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))        ! extract key
     if (tag == myTag) &                                  ! match
       counter(section) = counter(section) + 1_pInt
   endif   
 enddo

100 IO_countTagInPart = counter

end function IO_countTagInPart


!*********************************************************************
! return array of myTag presence within <part> for at most N[sections]
!*********************************************************************
function IO_spotTagInPart(myFile,part,myTag,Nsections)

 implicit none
 integer(pInt),    intent(in)  :: Nsections
 logical, dimension(Nsections) :: IO_spotTagInPart
 
 integer(pInt),    intent(in)  :: myFile
 character(len=*), intent(in)  :: part, &
                                  myTag

 integer(pInt), parameter     :: maxNchunks = 1_pInt

 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt)                            :: section
 character(len=1024)                      :: line, &
                                             tag

 IO_spotTagInPart = .false.                               ! assume to nowhere spot tag
 section = 0_pInt
 line =''

 rewind(myFile)
 do while (IO_getTag(line,'<','>') /= part)               ! search for part
   read(myFile,'(a1024)',END=100) line
 enddo

 do
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     section = section + 1_pInt
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))        ! extract key
     if (tag == myTag) &                                  ! match
       IO_spotTagInPart(section) = .true.
   endif   
 enddo

100 end function IO_spotTagInPart


!********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts in line and
! the left/right positions of at most N to be used by IO_xxxVal
!********************************************************************
pure function IO_stringPos(line,N)

 implicit none
 integer(pInt),    intent(in) :: N
 integer(pInt)                ::  IO_stringPos(1_pInt+N*2_pInt)
 
 character(len=*), intent(in) :: line
 
 character(len=*), parameter  :: sep=achar(44)//achar(32)//achar(9)//achar(10)//achar(13) ! comma and whitespaces

 integer                      :: left, right                      !no pInt (verify and scan return default integer)


 IO_stringPos = -1_pInt
 IO_stringPos(1) = 0_pInt
 right = 0
 
 do while (verify(line(right+1:),sep)>0)
   left  = right + verify(line(right+1:),sep)
   right = left + scan(line(left:),sep) - 2
   if ( line(left:left) == '#' ) then
     exit
   endif
   if ( IO_stringPos(1)<N ) then
     IO_stringPos(1_pInt+IO_stringPos(1)*2_pInt+1_pInt) = int(left, pInt)
     IO_stringPos(1_pInt+IO_stringPos(1)*2_pInt+2_pInt) = int(right, pInt)
   endif
   IO_stringPos(1) = IO_stringPos(1)+1_pInt
 enddo

end function IO_stringPos


!********************************************************************
! read string value at myPos from line
!********************************************************************
 pure function IO_stringValue(line,positions,myPos)
 
 implicit none

 integer(pInt),                                intent(in) :: positions(*), &
                                                             myPos

 character(len=1+positions(myPos*2+1)-positions(myPos*2)) :: IO_stringValue

 character(len=*),                             intent(in) :: line

 if (positions(1) < myPos) then
   IO_stringValue = ''
 else
   IO_stringValue = line(positions(myPos*2):positions(myPos*2+1))
 endif

end function IO_stringValue


!********************************************************************
! read string value at myPos from fixed format line
!********************************************************************
pure function IO_fixedStringValue (line,ends,myPos)
 
 implicit none

 integer(pInt),                intent(in) :: ends(*), &
                                             myPos
                                             
 character(len=ends(myPos+1)-ends(myPos)) :: IO_fixedStringValue
 
 character(len=*),             intent(in) :: line

 IO_fixedStringValue = line(ends(myPos)+1:ends(myPos+1))

end function IO_fixedStringValue


!********************************************************************
! read float value at myPos from line
!********************************************************************
real(pReal) pure function IO_floatValue (line,positions,myPos)
 
 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: positions(*), &
                                 myPos

 if (positions(1) < myPos) then
   IO_floatValue = 0.0_pReal
 else
   read(UNIT=line(positions(myPos*2):positions(myPos*2+1)),ERR=100,FMT=*) IO_floatValue
 endif
 return
100 IO_floatValue = huge(1.0_pReal)

end function IO_floatValue


!********************************************************************
! read float value at myPos from fixed format line
!********************************************************************
real(pReal) pure function IO_fixedFloatValue (line,ends,myPos)
 
 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: ends(*), &
                                 myPos

 read(UNIT=line(ends(myPos-1)+1:ends(myPos)),ERR=100,FMT=*) IO_fixedFloatValue
 return
100 IO_fixedFloatValue = huge(1.0_pReal)

end function IO_fixedFloatValue


!********************************************************************
! read float x.y+z value at myPos from format line line
!********************************************************************
real(pReal) pure function IO_fixedNoEFloatValue (line,ends,myPos)

 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: ends(*), &
                                 myPos

 integer(pInt)                :: expon
 integer                      :: pos_exp
 real(pReal)                  :: base
 
 pos_exp = scan(line(ends(myPos)+1:ends(myPos+1)),'+-',back=.true.)
 if (pos_exp > 1) then
   read(UNIT=line(ends(myPos)+1:ends(myPos)+pos_exp-1),ERR=100,FMT=*) base
   read(UNIT=line(ends(myPos)+pos_exp:ends(myPos+1)),ERR=100,FMT=*) expon
 else
   read(UNIT=line(ends(myPos)+1:ends(myPos+1)),ERR=100,FMT=*) base
   expon = 0_pInt
 endif
 IO_fixedNoEFloatValue = base*10.0_pReal**expon
 return
100 IO_fixedNoEFloatValue = huge(1.0_pReal)

end function IO_fixedNoEFloatValue


!********************************************************************
! read int value at myPos from line
!********************************************************************
integer(pInt) pure function IO_intValue(line,positions,myPos)

 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: positions(*), &
                                 myPos

 if (positions(1) < myPos) then
   IO_intValue = 0_pInt
 else
   read(UNIT=line(positions(myPos*2):positions(myPos*2+1)),ERR=100,FMT=*) IO_intValue
 endif
 return
100 IO_intValue = huge(1_pInt)

end function IO_intValue


!********************************************************************
! read int value at myPos from fixed format line
!********************************************************************
integer(pInt) pure function IO_fixedIntValue(line,ends,myPos)
 
 implicit none
 character(len=*), intent(in) :: line
 integer(pInt),    intent(in) :: ends(*), &
                                 myPos

 read(UNIT=line(ends(myPos)+1:ends(myPos+1)),ERR=100,FMT=*) IO_fixedIntValue
 return
100 IO_fixedIntValue = huge(1_pInt)

end function IO_fixedIntValue


!********************************************************************
! change character in line to lower case
!********************************************************************
pure function IO_lc(line)

 implicit none
 character(26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
 character(26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
 character(len=*), intent(in) :: line
 character(len=len(line))     :: IO_lc
 
 integer                      :: i,n                      ! no pInt (len returns default integer)

 IO_lc = line
 do i=1,len(line)
   n = index(upper,IO_lc(i:i))
   if (n/=0) IO_lc(i:i) = lower(n:n)
 enddo

end function IO_lc


!********************************************************************
! in place change of character in line to lower case
!********************************************************************
subroutine IO_lcInplace(line)

 implicit none
 character(26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
 character(26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
 character(len=*), intent(inout) :: line
 character(len=len(line))        :: IO_lc

 integer                         :: i,n                   ! no pInt (len returns default integer)

 do i=1,len(line)
   n = index(upper,line(i:i))
   if (n/=0) then 
     IO_lc(i:i) = lower(n:n)
   else
     IO_lc(i:i) = line(i:i)
   endif 
 enddo

 end subroutine IO_lcInplace


!********************************************************************
! read on in file to skip (at least) N chunks (may be over multiple lines)
!********************************************************************
subroutine IO_skipChunks(myUnit,N)

 implicit none
 integer(pInt), intent(in)                :: myUnit, &
                                             N

 integer(pInt), parameter                 :: maxNchunks = 64_pInt
 
 integer(pInt)                            :: remainingChunks
 integer(pInt), dimension(1+2*maxNchunks) :: myPos
 character(len=300)                       :: line

 remainingChunks = N
 do while (remainingChunks > 0)
   read(myUnit,'(A300)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   remainingChunks = remainingChunks - myPos(1)
 enddo
100  end subroutine IO_skipChunks


!********************************************************************
! extract value from key=value pair and check whether key matches
!********************************************************************
character(len=300) pure function IO_extractValue(line,key)
 
 implicit none
 character(len=*), intent(in) :: line, &
                                 key

 character(len=*), parameter  :: sep = achar(61)         ! '='

 integer                      :: myPos                                          ! no pInt (scan returns default integer)

 IO_extractValue = ''

 myPos = scan(line,sep)
 if (myPos > 0 .and. line(:myPos-1) == key(:myPos-1)) &       ! key matches expected key
   IO_extractValue = line(myPos+1:)                       ! extract value

end function IO_extractValue


!********************************************************************
! count lines containig data up to next *keyword
! AP: changed the function to neglect comment lines between keyword definitions.
!   : is not changed back to the original version since *.inp_assembly does not
!   : contain any comment lines (12.07.2010)
!********************************************************************
integer(pInt) function IO_countDataLines(myUnit)

 implicit none
 integer(pInt), intent(in)                :: myUnit
 
 integer(pInt), parameter                 :: maxNchunks = 1_pInt

 integer(pInt), dimension(1+2*maxNchunks) :: myPos
 character(len=300)                       :: line, &
                                             tmp

 IO_countDataLines = 0_pInt

 do
   read(myUnit,'(A300)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   tmp = IO_lc(IO_stringValue(line,myPos,1_pInt))
   if (tmp(1:1) == '*' .and. tmp(2:2) /= '*') then  ! found keyword
     exit
   else
     if (tmp(2:2) /= '*') IO_countDataLines = IO_countDataLines + 1_pInt
   endif
 enddo
100 backspace(myUnit)

end function IO_countDataLines

 
!********************************************************************
! count items in consecutive lines
! Marc:      ints concatenated by "c" as last char or range of values a "to" b
! Abaqus:    triplet of start,stop,inc
! Spectral:  ints concatenated range of a "to" b, multiple entries with a "copies of" b
!********************************************************************
integer(pInt) function IO_countContinuousIntValues(myUnit)

 implicit none
 integer(pInt), intent(in) :: myUnit
 
 integer(pInt), parameter :: maxNchunks = 8192_pInt
#ifdef Abaqus 
 integer(pInt)                            :: l,c
#endif
 integer(pInt), dimension(1+2*maxNchunks) :: myPos
 character(len=65536)                     :: line

 IO_countContinuousIntValues = 0_pInt

#ifndef Abaqus
 do
   read(myUnit,'(A300)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   if (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'to' ) then                  ! found range indicator
     IO_countContinuousIntValues = 1_pInt + IO_intValue(line,myPos,3_pInt) - IO_intValue(line,myPos,1_pInt)
     exit                                                                       ! only one single range indicator allowed
   else if (IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'copies' .and. &
            IO_lc(IO_stringValue(line,myPos,3_pInt)) == 'of'           ) then   ! found multiple entries indicator
     IO_countContinuousIntValues = IO_intValue(line,myPos,1_pInt)
     exit                                                                       ! only one single multiplier allowed
   else
     IO_countContinuousIntValues = IO_countContinuousIntValues+myPos(1)-1_pInt  ! add line's count when assuming 'c'
     if ( IO_lc(IO_stringValue(line,myPos,myPos(1))) /= 'c' ) then              ! line finished, read last value
       IO_countContinuousIntValues = IO_countContinuousIntValues+1_pInt
       exit                                                                     ! data ended
     endif
   endif
 enddo
#else
 c = IO_countDataLines(myUnit)
 do l = 1_pInt,c
   backspace(myUnit)
 enddo
     
 do l = 1_pInt,c
   read(myUnit,'(A300)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   IO_countContinuousIntValues = IO_countContinuousIntValues + 1_pInt + &    ! assuming range generation
                                (IO_intValue(line,myPos,2_pInt)-IO_intValue(line,myPos,1_pInt))/&
                                                     max(1_pInt,IO_intValue(line,myPos,3_pInt))
 enddo
#endif

100 end function IO_countContinuousIntValues


!********************************************************************
! return integer list corrsponding to items in consecutive lines.
! First integer in array is counter
! Marc:      ints concatenated by "c" as last char, range of a "to" b, or named set
! Abaqus:    triplet of start,stop,inc or named set
! Spectral:  ints concatenated range of a "to" b, multiple entries with a "copies of" b
!********************************************************************
function IO_continuousIntValues(myUnit,maxN,lookupName,lookupMap,lookupMaxN)

 implicit none
 integer(pInt),                     intent(in) :: maxN
 integer(pInt),     dimension(1+maxN)          :: IO_continuousIntValues
 
 integer(pInt),                     intent(in) :: myUnit, &
                                                  lookupMaxN
 integer(pInt),     dimension(:,:), intent(in) :: lookupMap
 character(len=64), dimension(:),   intent(in) :: lookupName
 integer(pInt), parameter :: maxNchunks = 8192_pInt
 integer(pInt) :: i
#ifdef Abaqus
 integer(pInt) :: j,l,c,first,last
#endif

 integer(pInt), dimension(1+2*maxNchunks) :: myPos
 character(len=65536) line
 logical rangeGeneration

 IO_continuousIntValues = 0_pInt
 rangeGeneration = .false.

#ifndef Abaqus
 do
   read(myUnit,'(A65536)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   if (verify(IO_stringValue(line,myPos,1_pInt),'0123456789') > 0) then     ! a non-int, i.e. set name
     do i = 1_pInt, lookupMaxN                                             ! loop over known set names
       if (IO_stringValue(line,myPos,1_pInt) == lookupName(i)) then         ! found matching name
         IO_continuousIntValues = lookupMap(:,i)                      ! return resp. entity list
         exit
       endif
     enddo
     exit
   else if (myPos(1) > 2_pInt .and. IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'to' ) then         ! found range indicator
     do i = IO_intValue(line,myPos,1_pInt),IO_intValue(line,myPos,3_pInt)
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1_pInt
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = i
     enddo
     exit
   else if (myPos(1) > 3_pInt .and. IO_lc(IO_stringValue(line,myPos,2_pInt)) == 'copies' &
                              .and. IO_lc(IO_stringValue(line,myPos,3_pInt)) == 'of' ) then         ! found multiple entries indicator
     IO_continuousIntValues(1) = IO_intValue(line,myPos,1_pInt)
     IO_continuousIntValues(2:IO_continuousIntValues(1)+1) = IO_intValue(line,myPos,4_pInt)
     exit
   else
     do i = 1_pInt,myPos(1)-1_pInt  ! interpret up to second to last value
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1_pInt
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,myPos,i)
     enddo
     if ( IO_lc(IO_stringValue(line,myPos,myPos(1))) /= 'c' ) then       ! line finished, read last value
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1_pInt
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,myPos,myPos(1))
       exit
     endif
   endif
 enddo
#else
 c = IO_countDataLines(myUnit)
 do l = 1_pInt,c
   backspace(myUnit)
 enddo
 
   !heck if the element values in the elset are auto generated
 backspace(myUnit)
 read(myUnit,'(A65536)',end=100) line
 myPos = IO_stringPos(line,maxNchunks)
 do i = 1_pInt,myPos(1)
   if (IO_lc(IO_stringValue(line,myPos,i)) == 'generate') rangeGeneration = .true.
 enddo
 
 do l = 1_pInt,c
   read(myUnit,'(A65536)',end=100) line
   myPos = IO_stringPos(line,maxNchunks)
   if (verify(IO_stringValue(line,myPos,1_pInt),'0123456789') > 0) then     ! a non-int, i.e. set names follow on this line
     do i = 1_pInt,myPos(1)                                                 ! loop over set names in line
       do j = 1_pInt,lookupMaxN                                      ! look thru known set names
         if (IO_stringValue(line,myPos,i) == lookupName(j)) then       ! found matching name
           first = 2_pInt + IO_continuousIntValues(1)                      ! where to start appending data
           last  = first + lookupMap(1,j) - 1_pInt                        ! up to where to append data
           IO_continuousIntValues(first:last) = lookupMap(2:1+lookupMap(1,j),j)    ! add resp. entity list
           IO_continuousIntValues(1) = IO_continuousIntValues(1) + lookupMap(1,j)   ! count them
         endif
       enddo
     enddo
   else if (rangeGeneration) then                                    ! range generation
     do i = IO_intValue(line,myPos,1_pInt),IO_intValue(line,myPos,2_pInt),max(1_pInt,IO_intValue(line,myPos,3_pInt))
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1_pInt
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = i
     enddo
   else                                                              ! read individual elem nums
     do i = 1_pInt,myPos(1)
      !  write(*,*)'IO_CIV-int',IO_intValue(line,myPos,i)
       IO_continuousIntValues(1) = IO_continuousIntValues(1) + 1_pInt
       IO_continuousIntValues(1+IO_continuousIntValues(1)) = IO_intValue(line,myPos,i)
     enddo
   endif
 enddo
#endif

100 end function IO_continuousIntValues


!********************************************************************
! write error statements to standard out
! and terminate the Marc run with exit #9xxx
! in ABAQUS either time step is reduced or execution terminated
!********************************************************************
subroutine IO_error(error_ID,e,i,g,ext_msg)
 implicit none
 integer(pInt),              intent(in) :: error_ID
 integer(pInt),    optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 
 character(len=1024)                    :: msg

 select case (error_ID)

 !* file handling errors
 
 case (100_pInt)
   msg = 'could not open file:'
 case (101_pInt)
   msg = 'write error for file:'
 case (102_pInt)
   msg = 'could not read file:'
 case (103_pInt)
   msg = 'could not assemble input files'
 

 !* material error messages and related messages in mesh
 
 case (150_pInt)
   msg = 'crystallite index out of bounds'
 case (151_pInt)
   msg = 'phase index out of bounds'
 case (152_pInt)
   msg = 'texture index out of bounds'
 case (153_pInt)
   msg = 'sum of phase fractions differs from 1'
 case (154_pInt)
   msg = 'homogenization index out of bounds'
 case (155_pInt)
   msg = 'microstructure index out of bounds'
 case (156_pInt)
   msg = 'reading from ODF file'
 case (160_pInt)
   msg = 'no entries in config part'
 case (170_pInt)
   msg = 'no homogenization specified via State Variable 2'
 case (180_pInt)
   msg = 'no microstructure specified via State Variable 3'


 !* plasticity error messages

 case (200_pInt)
   msg = 'unknown elasticity specified:' 
 case (201_pInt)
   msg = 'unknown plasticity specified:' 
 case (205_pInt)
   msg = 'unknown lattice structure encountered'

 case (210_pInt)
   msg = 'unknown material parameter for j2 plasticity phase:'
 case (211_pInt)
   msg = 'material parameter for j2 plasticity phase out of bounds:'
 case (212_pInt)
   msg = 'unknown plasticity output for j2 plasticity:'

 case (220_pInt)
   msg = 'unknown material parameter for phenopowerlaw plasticity phase:'
 case (221_pInt)
   msg = 'material parameter for phenopowerlaw plasticity phase out of bounds:'
 case (222_pInt)
   msg = 'unknown plasticity output for phenopowerlaw plasticity:'

 case (230_pInt)
   msg = 'unknown material parameter for titanmod plasticity phase:'
 case (231_pInt)
   msg = 'material parameter for titanmod plasticity phase out of bounds:'
 case (232_pInt)
   msg = 'unknown plasticity output for titanmod plasticity:'

 case (240_pInt)
   msg = 'unknown material parameter for dislotwin plasticity phase:'
 case (241_pInt)
   msg = 'material parameter for dislotwin plasticity phase out of bounds:'
 case (242_pInt)
   msg = 'unknown plasticity output for dislotwin plasticity:'
 case (243_pInt)
   msg = 'zero stacking fault energy'

 case (250_pInt)
   msg = 'unknown material parameter for nonlocal plasticity phase:'
 case (251_pInt)
   msg = 'material parameter for nonlocal plasticity phase out of bounds:'
 case (252_pInt)
   msg = 'unknown plasticity output for nonlocal plasticity:'
 case (253_pInt)
   msg = 'element type not supported for nonlocal plasticity'

 
 !* numerics error messages 

 case (300_pInt)
   msg = 'unknown numerics parameter:'
 case (301_pInt)
   msg = 'numerics parameter out of bounds:'
 
 
 !* math errors
 
 case (400_pInt)
   msg = 'matrix inversion error'
 case (401_pInt)
   msg = 'math_check: quat -> axisAngle -> quat failed'
 case (402_pInt)
   msg = 'math_check: quat -> R -> quat failed'
 case (403_pInt)
   msg = 'math_check: quat -> euler -> quat failed'
 case (404_pInt)
   msg = 'math_check: R -> euler -> R failed'
 case (405_pInt)
   msg = 'I_TO_HALTON-error: An input base BASE is <= 1'
 case (406_pInt)
   msg = 'Prime-error: N must be between 0 and PRIME_MAX'
 case (407_pInt)
   msg = 'Dimension in nearest neigbor search wrong'
 case (408_pInt)
   msg = 'Polar decomposition error'
 case (450_pInt)
   msg = 'unknown symmetry type specified'
 case (460_pInt)
   msg = 'kdtree2 error'

 !* homogenization errors

 case (500_pInt)
   msg = 'unknown homogenization specified'


 !* DAMASK_marc errors
 
 case (700_pInt)
   msg = 'invalid materialpoint result requested'


 !* errors related to spectral solver

 case (808_pInt)
   msg = 'precision not suitable for FFTW'
 case (809_pInt)
   msg = 'initializing FFTW'
 case (831_pInt)
   msg = 'mask consistency violated in spectral loadcase'
 case (832_pInt)
   msg = 'ill-defined L (each line should be either fully or not at all defined) in spectral loadcase'
 case (834_pInt)
   msg = 'negative time increment in spectral loadcase'
 case (835_pInt)
   msg = 'non-positive increments in spectral loadcase'
 case (836_pInt)
   msg = 'non-positive result frequency in spectral loadcase'
 case (837_pInt)
   msg = 'incomplete loadcase'
 case (838_pInt)
   msg = 'mixed boundary conditions allow rotation'
 case (841_pInt)
   msg = 'missing header length info in spectral mesh'
 case (842_pInt)
   msg = 'homogenization in spectral mesh'
 case (843_pInt)
   msg = 'resolution in spectral mesh'
 case (844_pInt)
   msg = 'dimension in spectral mesh'
 case (845_pInt)
   msg = 'incomplete information in spectral mesh header'
 case (846_pInt)
   msg = 'not a rotation defined for loadcase rotation'
 case (847_pInt)
   msg = 'updating of gamma operator not possible if it is pre calculated'
 case (880_pInt)
   msg = 'mismatch of microstructure count and a*b*c in geom file'


 !* Error messages related to parsing of Abaqus input file

 case (900_pInt)
   msg = 'PARSE ERROR: Improper definition of nodes in input file (Nnodes < 2)'
 case (901_pInt)
   msg = 'PARSE ERROR: No Elements defined in input file (Nelems = 0)'
 case (902_pInt)
   msg = 'PARSE ERROR: No Element sets defined in input file (Atleast one *Elset must exist)'
 case (903_pInt)
   msg = 'PARSE ERROR: No Materials defined in input file (Look into section assigments)'
 case (904_pInt)
   msg = 'PARSE ERROR: No elements could be assigned for Elset: '
 case (905_pInt)
   msg = 'PARSE ERROR: Error in mesh_abaqus_map_materials'
 case (906_pInt)
   msg = 'PARSE ERROR: Error in mesh_abaqus_count_cpElements'
 case (907_pInt)
   msg = 'PARSE ERROR: Incorrect size of mesh_mapFEtoCPelem in mesh_abaqus_map_elements; Size cannot be zero'
 case (908_pInt)
   msg = 'PARSE ERROR: Incorrect size of mesh_mapFEtoCPnode in mesh_abaqus_map_nodes; Size cannot be zero'
 case (909_pInt)
   msg = 'PARSE ERROR: Incorrect size of mesh_node in mesh_abaqus_build_nodes; must be equal to mesh_Nnodes'
 case (910_pInt)
   msg = 'PARSE ERROR: Incorrect element type mapping in '
 
 
 !* general error messages
 
 case (666_pInt)
   msg = 'memory leak detected'
 case default
   msg = 'Unknown error number...'

 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,'(a38)')        '+------------------------------------+'
 write(6,'(a38)')        '+               error                +'
 write(6,'(a17,i3,a18)') '+                ',error_ID,'                 +'
 write(6,'(a38)')        '+                                    +'
 write(6,'(a2,a)')       '+ ', trim(msg)
 if (present(ext_msg))  write(6,'(a2,a)') '+ ', trim(ext_msg)
 if (present(e)) then
   if (present(i) .and. present(g)) then
     write(6,'(a13,i6,a4,i2,a7,i4,a2)') '+ at element ',e,' IP ',i,' grain ',g,' +'
   else
     write(6,'(a18,i6,a14)') '+              at ',e,'             +'
   endif
 endif
 write(6,'(a38)') '+------------------------------------+'
 flush(6)
 call quit(9000_pInt+error_ID)
 !$OMP END CRITICAL (write2out)

! ABAQUS returns in some cases

end subroutine IO_error


!********************************************************************
! write warning statements to standard out
!********************************************************************
subroutine IO_warning(warning_ID,e,i,g,ext_msg)

 implicit none
 integer(pInt),              intent(in) :: warning_ID
 integer(pInt),    optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 
 character(len=1024)                    :: msg

 select case (warning_ID)
 case (34_pInt)
   msg = 'invalid restart increment given'
 case (35_pInt)
   msg = 'could not get $DAMASK_NUM_THREADS'
 case (40_pInt)
   msg = 'Found Spectral solver parameter '
 case (47_pInt)
   msg = 'No valid parameter for FFTW given, using FFTW_PATIENT'
 case (101_pInt)
   msg = '+    crystallite debugging off...    +'
 case (600_pInt)
   msg = '+  crystallite responds elastically  +'
 case (601_pInt)
   msg = '+      stiffness close to zero       +'
 case (650_pInt)
   msg = '+     polar decomposition failed     +'
 case (700_pInt)
   msg = '+      unknown crystal symmetry      +'
 case default
   msg = '+     unknown warning number...      +'
 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,'(a38)')        '+------------------------------------+'
 write(6,'(a38)')        '+              warning               +'
 write(6,'(a38)')        '+                                    +'
 write(6,'(a17,i3,a18)') '+                ',warning_ID,'                 +'
 write(6,'(a2,a)') '+ ', trim(msg)
 if (present(ext_msg))  write(6,'(a2,a)') '+ ', trim(ext_msg)
 if (present(e)) then
   if (present(i)) then
     if (present(g)) then
       write(6,'(a12,1x,i6,1x,a2,1x,i2,1x,a5,1x,i4,a2)') '+ at element',e,'IP',i,'grain',g,' +'
     else
       write(6,'(a12,1x,i6,1x,a2,1x,i2,a13)') '+ at element',e,'IP',i,'            +'
     endif
   else
     write(6,'(a12,1x,i6,a19)') '+ at element',e,'             +'
   endif
 endif
 write(6,'(a38)') '+------------------------------------+'
 flush(6)
 !$OMP END CRITICAL (write2out)

end subroutine IO_warning


! INTERNAL (HELPER) FUNCTIONS:

#ifdef Abaqus 
!********************************************************************
! AP: 12.07.10
!    create a new input file for abaqus simulations
!    by removing all comment lines and including "include"s
!********************************************************************
recursive function abaqus_assembleInputFile(unit1,unit2) result(createSuccess)

 use DAMASK_interface, only: getSolverWorkingDirectoryName

 implicit none
 integer(pInt), intent(in)                :: unit1, &
                                             unit2
 
 integer(pInt), parameter                 :: maxNchunks = 6_pInt

 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=300)                       :: line,fname
 logical                                  :: createSuccess,fexist


 
 do
   read(unit2,'(A300)',END=220) line
   positions = IO_stringPos(line,maxNchunks)

!   call IO_lcInPlace(line)
   if (IO_lc(IO_StringValue(line,positions,1_pInt))=='*include') then
     fname = trim(getSolverWorkingDirectoryName())//trim(line(9+scan(line(9:),'='):))
     inquire(file=fname, exist=fexist)
     if (.not.(fexist)) then
       !$OMP CRITICAL (write2out)
         write(6,*)'ERROR: file does not exist error in abaqus_assembleInputFile'
         write(6,*)'filename: ', trim(fname)
       !$OMP END CRITICAL (write2out)
       createSuccess = .false.
       return
     endif
     open(unit2+1,err=200,status='old',file=fname)
     if (abaqus_assembleInputFile(unit1,unit2+1_pInt)) then
       createSuccess=.true.
       close(unit2+1)
     else
       createSuccess=.false.
       return
     endif
   else if (line(1:2) /= '**' .OR. line(1:8)=='**damask') then
     write(unit1,'(A)') trim(line)
   endif
 enddo
 
220 createSuccess = .true.
 return
 
200 createSuccess =.false.

end function abaqus_assembleInputFile
#endif

!********************************************************************
! hybrid IA repetition counter
!********************************************************************
integer(pInt) function hybridIA_reps(dV_V,steps,C)

 implicit none
  integer(pInt), intent(in), dimension(3) :: &
   steps
 real(pReal),   intent(in), dimension(steps(3),steps(2),steps(1)) :: &
   dV_V
 real(pReal),   intent(in) :: &
   C
 
  integer(pInt) :: phi1,Phi,phi2
 
 hybridIA_reps = 0_pInt
 do phi1=1_pInt,steps(1)
   do Phi =1_pInt,steps(2)
     do phi2=1_pInt,steps(3)
       hybridIA_reps = hybridIA_reps+nint(C*dV_V(phi2,Phi,phi1), pInt)
     enddo
   enddo
 enddo
 
end function hybridIA_reps
 
end module IO
