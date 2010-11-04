!* $Id$
!##############################################################
 MODULE IO   
!##############################################################

 CONTAINS
!---------------------------
! function IO_abaqus_assembleInputFile
! function IO_open_file(unit,relPath)
! function IO_open_inputFile(unit)
! function IO_hybridIA(Nast,ODFfileName)
! private function hybridIA_reps(dV_V,steps,C)
! function IO_stringPos(line,maxN)
! function IO_stringValue(line,positions,pos)
! function IO_floatValue(line,positions,pos)
! function IO_intValue(line,positions,pos)
! function IO_fixedStringValue(line,ends,pos)
! function IO_fixedFloatValue(line,ends,pos)
! function IO_fixedFloatNoEValue(line,ends,pos)
! function IO_fixedIntValue(line,ends,pos)
! function IO_continousIntValues(unit,maxN)
! function IO_lc(line)
! subroutine IO_lcInplace(line)
! subroutine IO_error(ID)
! subroutine IO_warning(ID)
!---------------------------

!********************************************************************
! output version number
!********************************************************************
subroutine IO_init ()

  !$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  IO init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
  call flush(6)
  !$OMP END CRITICAL (write2out)
 
 return
endsubroutine



!********************************************************************
! AP: 12.07.10
!    create a new input file for abaqus simulations
!    by removing all comment lines and including "include"s
!********************************************************************
recursive function IO_abaqus_assembleInputFile(unit1,unit2) result(createSuccess)
 use prec
 use mpie_interface
 implicit none
 
 character(len=300) line,fname
 integer(pInt), intent(in) :: unit1, unit2
 logical createSuccess,fexist
 integer(pInt) i
 
 do
   read(unit2,'(A300)',END=220) line
   line = IO_lc(trim(line))
!   call IO_lcInPlace(line)
   if (line(1:8)=='*include') then
     fname = trim(getSolverWorkingDirectoryName())//trim(line(9+scan(line(9:),'='):))
     inquire(file=fname, exist=fexist)
     if (.not.(fexist)) then
       write(6,*)'ERROR: file does not exist error in IO_abaqus_assembleInputFile'
       write(6,*)'filename: ', trim(fname)
       createSuccess = .false.
       return
     endif
     open(unit2+1,err=200,status='old',file=fname)
     if (IO_abaqus_assembleInputFile(unit1,unit2+1)) then
       createSuccess=.true.
       close(unit2+1)
     else
       createSuccess=.false.
       return
     endif
   else if (line(1:2) /= '**') then 
     write(unit1,'(A)') trim(line)
   endif
 enddo
 
220 createSuccess = .true.
 return
 
200 createSuccess =.false.
 return
 
end function

!***********************************************************
! check if the input file for Abaqus contains part info
!***********************************************************
 function IO_abaqus_hasNoPart(unit)
 
 use prec, only: pInt
 implicit none
 
 integer(pInt) unit
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension(1+2*maxNchunks) :: pos
 logical IO_abaqus_hasNoPart
 character(len=300) line
 
 IO_abaqus_hasNoPart = .true.
 
610 FORMAT(A300)
 rewind(unit)
 do
   read(unit,610,END=620) line
   pos = IO_stringPos(line,maxNchunks)
   if (IO_lc(IO_stringValue(line,pos,1)) == '*part' ) then
     IO_abaqus_hasNoPart = .false.
     exit
   endif
 enddo
 
620 return 
 
 endfunction



!********************************************************************
! open existing file to given unit
! path to file is relative to working directory
!********************************************************************
 logical function IO_open_file(unit,relPath)

 use prec, only: pInt
 use mpie_interface
 implicit none

 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forward and backward slash
 character(len=*) relPath
 integer(pInt) unit

 IO_open_file = .false.
 
 open(unit,status='old',err=100,file=trim(getSolverWorkingDirectoryName())//relPath)
 IO_open_file = .true.
 
100 return
 
 endfunction


!********************************************************************
! open FEM inputfile to given unit
! AP: 12.07.10 
!   : changed the function to open *.inp_assembly, which is basically 
!     the input file without comment lines and possibly assembled includes
!********************************************************************
 logical function IO_open_inputFile(unit)

 use prec, only: pReal, pInt
 use mpie_interface
 implicit none

 integer(pInt), intent(in) :: unit

 IO_open_inputFile = .false.
 
 if (FEsolver == 'Abaqus') then
   open(unit+1,status='old',err=100,&
               file=trim(getSolverWorkingDirectoryName())//&
                    trim(getSolverJobName())//InputFileExtension)
   open(unit,err=100,file=trim(getSolverWorkingDirectoryName())//&
                          trim(getSolverJobName())//InputFileExtension//'_assembly')
   IO_open_inputFile = IO_abaqus_assembleInputFile(unit,unit+1)          ! strip comments and concatenate any "include"s
   close(unit+1) 
 else
   open(unit,status='old',err=100,file=trim(getSolverWorkingDirectoryName())//&
                                       trim(getSolverJobName())//InputFileExtension)
   IO_open_inputFile = .true.
 endif

100 return

 endfunction


!********************************************************************
! open FEM logfile to given unit
!********************************************************************
 logical function IO_open_logFile(unit)

 use prec, only: pReal, pInt
 use mpie_interface
 implicit none

 integer(pInt), intent(in) :: unit

 IO_open_logFile = .false.
 
 open(unit,status='old',err=100,file=trim(getSolverWorkingDirectoryName())//&
                                     trim(getSolverJobName())//LogFileExtension)
 IO_open_logFile = .true.

100 return

 endfunction


!********************************************************************
! open (write) file related to current job
! but with different extension to given unit
!********************************************************************
 logical function IO_open_jobFile(unit,newExt)

 use prec, only: pReal, pInt
 use mpie_interface
 implicit none

 integer(pInt), intent(in) :: unit
 character(*), intent(in) :: newExt

 IO_open_jobFile = .false.
 
 open(unit,status='replace',err=100,file=trim(getSolverWorkingDirectoryName())//&
                                     trim(getSolverJobName())//'.'//newExt)
 IO_open_jobFile = .true.
 
100 return

 endfunction


!********************************************************************
! open (write) binary file related to current job
! but with different extension to given unit
!********************************************************************
 logical function IO_write_jobBinaryFile(unit,newExt,recMultiplier)

 use prec, only: pReal, pInt
 use mpie_interface
 implicit none

 integer(pInt), intent(in) :: unit
 integer(pInt), intent(in), optional :: recMultiplier
 character(*), intent(in) :: newExt

 IO_write_jobBinaryFile = .false.
 if (present(recMultiplier)) then
   open(unit,status='replace',form='unformatted',access='direct',recl=pReal*recMultiplier, &
             err=100,file=trim(getSolverWorkingDirectoryName())//&
                                       trim(getSolverJobName())//'.'//newExt)
  else
   open(unit,status='replace',form='unformatted',access='direct',recl=pReal, &
             err=100,file=trim(getSolverWorkingDirectoryName())//&
                                       trim(getSolverJobName())//'.'//newExt)
 endif
 IO_write_jobBinaryFile = .true.
 
100 return

 endfunction


!********************************************************************
! open (read) binary file related to restored job
! and with different extension to given unit
!********************************************************************
 logical function IO_read_jobBinaryFile(unit,newExt,jobName,recMultiplier)

 use prec, only: pReal, pInt
 use mpie_interface
 implicit none

 integer(pInt), intent(in) :: unit
 integer(pInt), intent(in), optional :: recMultiplier
 character(*), intent(in) :: newExt, jobName

 IO_read_jobBinaryFile = .false.
 if (present(recMultiplier)) then
   open(unit,status='old',form='unformatted',access='direct',recl=pReal*recMultiplier, &
             err=100,file=trim(getSolverWorkingDirectoryName())//&
                                                  trim(jobName)//'.'//newExt)
  else
   open(unit,status='old',form='unformatted',access='direct',recl=pReal, &
             err=100,file=trim(getSolverWorkingDirectoryName())//&
                                                  trim(jobName)//'.'//newExt)
 endif
 IO_read_jobBinaryFile = .true.
 
100 return

 endfunction


!********************************************************************
! hybrid IA repetition counter
!********************************************************************
 function hybridIA_reps(dV_V,steps,C)

 use prec, only: pReal, pInt
 implicit none
 
 integer(pInt), intent(in), dimension(3) :: steps
 integer(pInt) hybridIA_reps, phi1,Phi,phi2
 real(pReal), intent(in), dimension(steps(3),steps(2),steps(1)) :: dV_V
 real(pReal), intent(in) :: C
 
 hybridIA_reps = 0_pInt
 do phi1=1,steps(1)
   do Phi =1,steps(2)
     do phi2=1,steps(3)
       hybridIA_reps = hybridIA_reps+nint(C*dV_V(phi2,Phi,phi1), pInt)
     enddo
   enddo
 enddo
 return
 
 endfunction


!********************************************************************
! hybrid IA sampling of ODFfile
!********************************************************************
 function IO_hybridIA(Nast,ODFfileName)

 use prec, only: pReal, pInt
 implicit none
 
 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal
 real(pReal), parameter :: inRad = pi/180.0_pReal

 character(len=*) ODFfileName
 character(len=80) line
 character(len=*), parameter :: fileFormat = '(A80)'
 integer(pInt) i,j,bin,Nast,NnonZero,Nset,Nreps,reps,phi1,Phi,phi2
 integer(pInt), dimension(7) :: pos
 integer(pInt), dimension(3) :: steps
 integer(pInt), dimension(:), allocatable :: binSet
 real(pReal) center,sum_dV_V,prob,dg_0,C,lowerC,upperC,rnd
 real(pReal), dimension(3) :: limits,deltas
 real(pReal), dimension(:,:,:), allocatable :: dV_V
 real(pReal), dimension(3,Nast) :: IO_hybridIA

 if (.not. IO_open_file(999,ODFfileName)) goto 100
 
!--- parse header of ODF file ---
!--- limits in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 pos = IO_stringPos(line,3)
 if (pos(1).ne.3) goto 100
 do i=1,3
   limits(i) = IO_intValue(line,pos,i)*inRad
 enddo

!--- deltas in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 pos = IO_stringPos(line,3)
 if (pos(1).ne.3) goto 100
 do i=1,3
   deltas(i) = IO_intValue(line,pos,i)*inRad
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
 
 do phi1=1,steps(1)
   do Phi=1,steps(2)
     do phi2=1,steps(3)
       read(999,fmt=*,end=100) prob
       if (prob > 0.0_pReal) then
         NnonZero = NnonZero+1
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
 bin = 0 ! bin counter
 i = 1 ! set counter
 do phi1=1,steps(1)
   do Phi=1,steps(2)
     do phi2=1,steps(3)
       reps = nint(C*dV_V(phi2,Phi,phi1), pInt)
       binSet(i:i+reps-1) = bin
       bin = bin+1 ! advance bin
       i = i+reps ! advance set
     enddo
   enddo
 enddo

 do i=1,Nast
   if (i < Nast) then
     call random_number(rnd)
     j = nint(rnd*(Nreps-i)+i+0.5_pReal,pInt)
   else
     j = i
   endif
   bin = binSet(j)
   IO_hybridIA(1,i) = deltas(1)*(mod(bin/(steps(3)*steps(2)),steps(1))+center)  ! phi1
   IO_hybridIA(2,i) = deltas(2)*(mod(bin/ steps(3)          ,steps(2))+center)  ! Phi
   IO_hybridIA(3,i) = deltas(3)*(mod(bin                    ,steps(3))+center)  ! phi2
   binSet(j) = binSet(i)
 enddo
 close(999)
 return

! on error
100 IO_hybridIA = -1
 close(999)
 return
 
 endfunction 


!********************************************************************
! identifies lines without content
!********************************************************************
 pure function IO_isBlank (line)

 use prec, only: pInt
 implicit none

 character(len=*), intent(in) :: line
 character(len=*), parameter :: blank = achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 character(len=*), parameter :: comment = achar(35)                               ! comment id '#'
 integer(pInt) posNonBlank, posComment
 logical IO_isBlank
 
 posNonBlank = verify(line,blank)
 posComment  = scan(line,comment)
 IO_isBlank = posNonBlank == 0 .or. posNonBlank == posComment

 return
 
 endfunction

!********************************************************************
! get tagged content of line
!********************************************************************
 pure function IO_getTag (line,openChar,closechar)

 use prec, only: pInt
 implicit none

 character(len=*), intent(in) :: line,openChar,closeChar
 character(len=*), parameter :: sep=achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 character(len=len_trim(line)) IO_getTag
 integer(pInt)  left,right

 IO_getTag = ''
 left = scan(line,openChar)
 right = scan(line,closeChar)
 
 if (left == verify(line,sep) .and. right > left) & ! openChar is first and closeChar occurs
   IO_getTag = line(left+1:right-1)

 return
 
 endfunction


!*********************************************************************
 function IO_countSections(file,part)
!*********************************************************************
 use prec, only: pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: file
 character(len=*), intent(in) :: part
 integer(pInt) IO_countSections
 character(len=1024) line

 IO_countSections = 0
 line = ''
 rewind(file)

 do while (IO_getTag(line,'<','>') /= part)      ! search for part
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     IO_countSections = IO_countSections + 1
 enddo

100 return

 endfunction
 

!*********************************************************************
! return array of myTag counts within <part> for at most N[sections]
!*********************************************************************
 function IO_countTagInPart(file,part,myTag,Nsections)

 use prec, only: pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: file, Nsections
 character(len=*), intent(in) :: part, myTag
 integer(pInt), dimension(Nsections) :: IO_countTagInPart, counter
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section
 character(len=1024) line,tag

 counter = 0_pInt
 section = 0_pInt
 line = ''
 rewind(file)

 do while (IO_getTag(line,'<','>') /= part)               ! search for part
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     section = section + 1
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     if (tag == myTag) &                                  ! match
       counter(section) = counter(section) + 1
   endif   
 enddo

100 IO_countTagInPart = counter
 return

endfunction


!*********************************************************************
! return array of myTag presence within <part> for at most N[sections]
!*********************************************************************
 function IO_spotTagInPart(file,part,myTag,Nsections)

 use prec, only: pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: file, Nsections
 character(len=*), intent(in) :: part, myTag
 logical, dimension(Nsections) :: IO_spotTagInPart
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section
 character(len=1024) line,tag

 IO_spotTagInPart = .false.                               ! assume to nowhere spot tag
 section = 0_pInt
 line = ''
 rewind(file)

 do while (IO_getTag(line,'<','>') /= part)               ! search for part
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     section = section + 1
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     if (tag == myTag) &                                  ! match
       IO_spotTagInPart(section) = .true.
   endif   
 enddo

100 return

endfunction


!********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts in line and
! the left/right positions of at most N to be used by IO_xxxVal
!********************************************************************
! pure function IO_stringPos (line,N)
 function IO_stringPos (line,N)

 use prec, only: pReal,pInt
 implicit none

 character(len=*), intent(in) :: line
 character(len=*), parameter :: sep=achar(44)//achar(32)//achar(9)//achar(10)//achar(13) ! comma and whitespaces
 integer(pInt), intent(in) :: N
 integer(pInt)  left,right
 integer(pInt) IO_stringPos(1+N*2)

 IO_stringPos = -1
 IO_stringPos(1) = 0
 right = 0

 do while (verify(line(right+1:),sep)>0)
   left  = right + verify(line(right+1:),sep)
   right = left + scan(line(left:),sep) - 2
   if ( IO_stringPos(1)<N ) then
     IO_stringPos(1+IO_stringPos(1)*2+1) = left
     IO_stringPos(1+IO_stringPos(1)*2+2) = right
   endif
   IO_stringPos(1) = IO_stringPos(1)+1
 enddo

 return

 endfunction


!********************************************************************
! read string value at pos from line
!********************************************************************
 pure function IO_stringValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 character(len=1+positions(pos*2+1)-positions(pos*2)) IO_stringValue

 if (positions(1) < pos) then
   IO_stringValue = ''
 else
   IO_stringValue = line(positions(pos*2):positions(pos*2+1))
 endif
 return

 endfunction


!********************************************************************
! read string value at pos from fixed format line
!********************************************************************
 pure function IO_fixedStringValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 character(len=ends(pos+1)-ends(pos)) IO_fixedStringValue

 IO_fixedStringValue = line(ends(pos)+1:ends(pos+1))
 return

 endfunction


!********************************************************************
! read float value at pos from line
!********************************************************************
 pure function IO_floatValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 real(pReal) IO_floatValue

 if (positions(1) < pos) then
   IO_floatValue = 0.0_pReal
 else
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_floatValue
 endif
 return
100 IO_floatValue = huge(1.0_pReal)
 return

 endfunction


!********************************************************************
! read float value at pos from fixed format line
!********************************************************************
 pure function IO_fixedFloatValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 real(pReal) IO_fixedFloatValue

 read(UNIT=line(ends(pos-1)+1:ends(pos)),ERR=100,FMT=*) IO_fixedFloatValue
 return
100 IO_fixedFloatValue = huge(1.0_pReal)
 return

 endfunction


!********************************************************************
! read float x.y+z value at pos from format line line
!********************************************************************
 pure function IO_fixedNoEFloatValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 integer(pInt) pos_exp,expon
 real(pReal) IO_fixedNoEFloatValue,base
 
 pos_exp = scan(line(ends(pos)+1:ends(pos+1)),'+-',back=.true.)
 if (pos_exp > 1) then
   read(UNIT=line(ends(pos)+1:ends(pos)+pos_exp-1),ERR=100,FMT=*) base
   read(UNIT=line(ends(pos)+pos_exp:ends(pos+1)),ERR=100,FMT=*) expon
 else
   read(UNIT=line(ends(pos)+1:ends(pos+1)),ERR=100,FMT=*) base
   expon = 0_pInt
 endif
 IO_fixedNoEFloatValue = base*10.0_pReal**expon
 return
100 IO_fixedNoEFloatValue = huge(1.0_pReal)
 return

 endfunction


!********************************************************************
! read int value at pos from line
!********************************************************************
 pure function IO_intValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 integer(pInt) IO_intValue

 if (positions(1) < pos) then
   IO_intValue = 0_pInt
 else
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_intValue
 endif
 return
100 IO_intValue = huge(1_pInt)
 return

 endfunction


!********************************************************************
! read int value at pos from fixed format line
!********************************************************************
 pure function IO_fixedIntValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 integer(pInt) IO_fixedIntValue

 read(UNIT=line(ends(pos)+1:ends(pos+1)),ERR=100,FMT=*) IO_fixedIntValue
 return
100 IO_fixedIntValue = huge(1_pInt)
 return

 endfunction


!********************************************************************
! change character in line to lower case
!********************************************************************
 pure function IO_lc (line)

 use prec, only: pInt
 implicit none

 character (len=*), intent(in) :: line
 character (len=len(line)) IO_lc
 integer(pInt) i

 IO_lc = line
 do i=1,len(line)
    if(64<iachar(line(i:i)) .and. iachar(line(i:i))<91) IO_lc(i:i)=achar(iachar(line(i:i))+32)
 enddo
 return 

 endfunction


!********************************************************************
! in place change of character in line to lower case
!********************************************************************
 subroutine IO_lcInplace (line)

 use prec, only: pInt
 implicit none

 character (len=*) line
 character (len=len(line)) IO_lc
 integer(pInt) i

 IO_lc = line
 do i=1,len(line)
    if(64<iachar(line(i:i)) .and. iachar(line(i:i))<91) IO_lc(i:i)=achar(iachar(line(i:i))+32)
 enddo
 line = IO_lc
 return 

 endsubroutine


!********************************************************************
! read on in file to skip (at least) N chunks (may be over multiple lines)
!********************************************************************
 subroutine IO_skipChunks (unit,N)

 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  remainingChunks,unit,N
 integer(pInt), parameter :: maxNchunks = 64
 integer(pInt), dimension(1+2*maxNchunks) :: pos
 character(len=300) line

 remainingChunks = N
 do while (remainingChunks > 0)
   read(unit,'(A300)',end=100) line
   pos = IO_stringPos(line,maxNchunks)
   remainingChunks = remainingChunks - pos(1)
 enddo
100 return

 endsubroutine


!********************************************************************
! extract value from key=value pair and check whether key matches
!********************************************************************
 pure function IO_extractValue (line,key)
 
 use prec, only: pReal,pInt
 implicit none

 character(len=*), intent(in) :: line,key
 character(len=*), parameter :: sep = achar(61)         ! '='
 integer(pInt) pos
 character(len=300) IO_extractValue

 IO_extractValue = ''

 pos = scan(line,sep)
 if (pos > 0 .and. line(:pos-1) == key(:pos-1)) &       ! key matches expected key
   IO_extractValue = line(pos+1:)                       ! extract value

 return

 endfunction


!********************************************************************
! count lines containig data up to next *keyword
! AP: changed the function to neglect comment lines between keyword definitions.
!   : is not changed back to the original version since *.inp_assembly does not
!   : contain any comment lines (12.07.2010)
!********************************************************************
 function IO_countDataLines (unit)

 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  IO_countDataLines,unit
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension(1+2*maxNchunks) :: pos
 character(len=300) line,tmp

 IO_countDataLines = 0
 do
   read(unit,'(A300)',end=100) line
   pos = IO_stringPos(line,maxNchunks)
   tmp = IO_lc(IO_stringValue(line,pos,1))
   if (tmp(1:1) == '*' .and. tmp(2:2) /= '*') then  ! found keyword
     exit
   else
     if (tmp(2:2) /= '*') IO_countDataLines = IO_countDataLines + 1_pInt
   endif
 enddo
100 backspace(unit)
 return

 endfunction

 
!********************************************************************
! count items in consecutive lines
! Marc:   ints concatenated by "c" as last char or range of values a "to" b
! Abaqus: triplet of start,stop,inc
!********************************************************************
 function IO_countContinousIntValues (unit)

 use mpie_interface
 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  unit,l,count
 integer(pInt)  IO_countContinousIntValues
 integer(pInt), parameter :: maxNchunks = 64
 integer(pInt), dimension(1+2*maxNchunks) :: pos
 character(len=300) line

 IO_countContinousIntValues = 0_pInt

 select case (FEsolver)
   case ('Marc')
   
     do
       read(unit,'(A300)',end=100) line
       pos = IO_stringPos(line,maxNchunks)
       if (IO_lc(IO_stringValue(line,pos,2)) == 'to' ) then               ! found range indicator
         IO_countContinousIntValues = 1 + IO_intValue(line,pos,3) - IO_intValue(line,pos,1)
         exit                                                             ! only one single range indicator allowed
       else
         IO_countContinousIntValues = IO_countContinousIntValues+pos(1)-1 ! add line's count when assuming 'c'
         if ( IO_lc(IO_stringValue(line,pos,pos(1))) /= 'c' ) then        ! line finished, read last value
           IO_countContinousIntValues = IO_countContinousIntValues+1
           exit                                                           ! data ended
         endif
       endif
     enddo
   
   case('Abaqus')

     count = IO_countDataLines(unit)
     do l = 1,count
       backspace(unit)
     enddo
     
     do l = 1,count
       read(unit,'(A300)',end=100) line
       pos = IO_stringPos(line,maxNchunks)
       IO_countContinousIntValues = IO_countContinousIntValues + 1 + &    ! assuming range generation
                                    (IO_intValue(line,pos,2)-IO_intValue(line,pos,1))/max(1,IO_intValue(line,pos,3))
     enddo
 
 endselect

100 return

 endfunction


!********************************************************************
! return integer list corrsponding to items in consecutive lines
! Marc:   ints concatenated by "c" as last char, range of a "to" b, or named set
! Abaqus: triplet of start,stop,inc or named set
!********************************************************************
 function IO_continousIntValues (unit,maxN,lookupName,lookupMap,lookupMaxN)

 use mpie_interface
 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  unit,maxN,i,j,l,count,first,last
 integer(pInt), dimension(1+maxN) :: IO_continousIntValues
 integer(pInt), parameter :: maxNchunks = 64
 integer(pInt), dimension(1+2*maxNchunks) :: pos
 character(len=64), dimension(:) :: lookupName
 integer(pInt) :: lookupMaxN
 integer(pInt), dimension(:,:) :: lookupMap
 character(len=300) line
 logical rangeGeneration

 IO_continousIntValues = 0
 rangeGeneration = .false.

 select case (FEsolver)
   case ('Marc')
   
     do
       read(unit,'(A300)',end=100) line
       pos = IO_stringPos(line,maxNchunks)
       if (verify(IO_stringValue(line,pos,1),"0123456789") > 0) then     ! a non-int, i.e. set name
         do i = 1,lookupMaxN                                             ! loop over known set names
           if (IO_stringValue(line,pos,1) == lookupName(i)) then         ! found matching name
             IO_continousIntValues = lookupMap(:,i)                      ! return resp. entity list
             exit
           endif
         enddo
         exit
       else if (IO_lc(IO_stringValue(line,pos,2)) == 'to' ) then         ! found range indicator
         do i = IO_intValue(line,pos,1),IO_intValue(line,pos,3)
           IO_continousIntValues(1) = IO_continousIntValues(1) + 1
           IO_continousIntValues(1+IO_continousIntValues(1)) = i
         enddo
         exit
       else
         do i = 1,pos(1)-1  ! interpret up to second to last value
           IO_continousIntValues(1) = IO_continousIntValues(1) + 1
           IO_continousIntValues(1+IO_continousIntValues(1)) = IO_intValue(line,pos,i)
         enddo
         if ( IO_lc(IO_stringValue(line,pos,pos(1))) /= 'c' ) then       ! line finished, read last value
           IO_continousIntValues(1) = IO_continousIntValues(1)+1
           IO_continousIntValues(1+IO_continousIntValues(1)) = IO_intValue(line,pos,pos(1))
           exit
         endif
       endif
     enddo
   
   case('Abaqus')

     count = IO_countDataLines(unit)
     do l = 1,count
       backspace(unit)
     enddo
     
!      check if the element values in the elset are auto generated
     backspace(unit)
     read(unit,'(A300)',end=100) line
     pos = IO_stringPos(line,maxNchunks)
     do i = 1,pos(1)
       if (IO_lc(IO_stringValue(line,pos,i)) == 'generate') rangeGeneration = .true.
     enddo
     
     do l = 1,count
       read(unit,'(A300)',end=100) line
       pos = IO_stringPos(line,maxNchunks)
       if (verify(IO_stringValue(line,pos,1),"0123456789") > 0) then     ! a non-int, i.e. set names follow on this line
         do i = 1,pos(1)                                                 ! loop over set names in line
           do j = 1,lookupMaxN                                           ! look thru known set names
             if (IO_stringValue(line,pos,i) == lookupName(j)) then       ! found matching name
               first = 2 + IO_continousIntValues(1)                      ! where to start appending data
               last  = first + lookupMap(1,j) - 1                        ! up to where to append data
               IO_continousIntValues(first:last) = lookupMap(2:1+lookupMap(1,j),j)    ! add resp. entity list
               IO_continousIntValues(1) = IO_continousIntValues(1) + lookupMap(1,j)   ! count them
             endif
           enddo
         enddo
       else if (rangeGeneration) then                                    ! range generation
         do i = IO_intValue(line,pos,1),IO_intValue(line,pos,2),max(1,IO_intValue(line,pos,3))
           IO_continousIntValues(1) = IO_continousIntValues(1) + 1
           IO_continousIntValues(1+IO_continousIntValues(1)) = i
         enddo
       else                                                              ! read individual elem nums
         do i = 1,pos(1)
!            write(*,*)'IO_CIV-int',IO_intValue(line,pos,i)
           IO_continousIntValues(1) = IO_continousIntValues(1) + 1
           IO_continousIntValues(1+IO_continousIntValues(1)) = IO_intValue(line,pos,i)
         enddo
       endif
     enddo
 
 endselect

100 return

 endfunction



!********************************************************************
! write error statements to standard out
! and terminate the Marc run with exit #9xxx
! in ABAQUS either time step is reduced or execution terminated
!********************************************************************
 subroutine IO_error(ID,e,i,g,ext_msg)

 use prec, only: pInt
 implicit none

 integer(pInt), intent(in) :: ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=120) msg

 select case (ID)
 case (40)
   msg = 'path rectification error'
 case (41)
   msg = 'path too long'
 case (42)
   msg = 'missing descriptive information in spectral mesh'
 case (43)
   msg = 'resolution error in spectral mesh'
 case (44)
   msg = 'dimension error in spectral mesh'
 case (45)
   msg = 'error opening spectral loadcase'
 case (46)
   msg = 'missing parameter in spectral loadcase'
 case (47)
   msg = 'mask consistency violated in spectral loadcase'
 case (48)
   msg = 'Non-positive relative tolerance for defGrad correction'
 case (50)
   msg = 'Error writing constitutive output description'
 case (100)
   msg = 'Cannot open config file'
 case (101)
   msg = 'Cannot open input file'
 case (102)
   msg = 'argument count error (mesh and loadcase) for mpie_spectral'
 case (105)
   msg = 'Error reading from ODF file'
 case (110)
   msg = 'No homogenization specified via State Variable 2'
 case (120)
   msg = 'No microstructure specified via State Variable 3'
 case (125)
   msg = 'No entries in config part'
 case (130)
   msg = 'Homogenization index out of bounds'
 case (140)
   msg = 'Microstructure index out of bounds'
 case (150)
   msg = 'Crystallite index out of bounds'
 case (155)
   msg = 'Phase index out of bounds'
 case (160)
   msg = 'Texture index out of bounds'
 case (170)
   msg = 'Sum of phase fractions differs from 1'
 case (200)
   msg = 'Unknown constitution specified'
 case (201)
   msg = 'Unknown homogenization specified'
 case (205)
   msg = 'Unknown lattice structure encountered'
 case (210)
   msg = 'Negative initial resistance'
 case (211)
   msg = 'Non-positive reference shear rate'
 case (212)
   msg = 'Non-positive stress exponent'
 case (213)
   msg = 'Non-positive saturation stress'
 case (214)
   msg = 'Zero hardening exponent'
 case (220)
   msg = 'Negative initial dislocation density'
 case (221)
   msg = 'Negative Bugers vector'
 case (222)
   msg = 'Negative activation energy for edge dislocation glide'
 case (223)
   msg = 'Negative self diffusion energy'
 case (224)
   msg = 'Non-positive diffusion prefactor'
 case (225)
   msg = 'No slip systems specified'
 case (226)
   msg = 'Non-positive prefactor for dislocation velocity'
 case (227)
   msg = 'Non-positive prefactor for mean free path'
 case (228)
   msg = 'Non-positive minimum stable dipole distance'
 case (229)
   msg = 'Non-positive hardening interaction coefficients'
 case (230)
   msg = 'Non-positive atomic volume'
 case (231)
   msg = 'Non-positive prefactor for self-diffusion coefficient'
 case (232)
   msg = 'Non-positive activation energy for dislocation climb'   
 case (233)
   msg = 'Non-positive relevant dislocation density'   
 case (240)
   msg = 'Non-positive Taylor factor'
 case (241)
   msg = 'Non-positive hardening exponent'
 case (242)
   msg = 'Non-positive relevant slip resistance'   
 case (260)
   msg = 'Non-positive relevant strain'
 case (261)
   msg = 'Frequency for Stiffness update smaller than zero'
 case (262)
   msg = 'Frequency for Jacobian update of Lp residuum smaller than zero'
 case (263)
   msg = 'Non-positive perturbation value'
 case (264)
   msg = 'Limit for homogenization loop too small'
 case (265)
   msg = 'Limit for crystallite loop too small'
 case (266)
   msg = 'Limit for state loop too small'
 case (267)
   msg = 'Limit for stress loop too small'
 case (268)
   msg = 'Non-positive minimum substep size'
 case (269)
   msg = 'Non-positive relative tolerance for state'
 case (270)
   msg = 'Non-positive relative tolerance for stress'
 case (271)
   msg = 'Non-positive absolute tolerance for stress'

!* Error messages related to RGC numerical parameters <<<updated 31.07.2009>>>
 case (272)
   msg = 'Non-positive relative tolerance of residual in RGC'
 case (273)
   msg = 'Non-positive absolute tolerance of residual in RGC'
 case (274)
   msg = 'Non-positive relative maximum of residual in RGC'
 case (275)
   msg = 'Non-positive absolute maximum of residual in RGC'
 case (276)
   msg = 'Non-positive penalty perturbation in RGC'
 case (277)
   msg = 'Non-positive relevant mismatch in RGC'
 case (278)
   msg = 'Non-positive definite viscosity model in RGC'
 case (288)
   msg = 'Non-positive maximum threshold of relaxation change in RGC'
 case (289)
   msg = 'Non-positive definite volume discrepancy penalty in RGC'

 case (294)
   msg = 'Non-positive tolerance for deformation gradient'

 case (298)
   msg = 'Chosen integration method does not exist'
 case (299)
   msg = 'Chosen perturbation method does not exist'

 case (300)
   msg = 'This material can only be used with elements with three direct stress components'
 case (500)
   msg = 'Unknown lattice type specified'
 case (550)
   msg = 'Unknown symmetry type specified'
 case (600)
   msg = 'Convergence not reached'
 case (610)
   msg = 'Stress loop not converged'

 case (666)
   msg = 'Memory leak detected'

 case (670)
   msg = 'math_check: quat -> axisAngle -> quat failed'
 case (671)
   msg = 'math_check: quat -> R -> quat failed'
 case (672)
   msg = 'math_check: quat -> euler -> quat failed'
 case (673)
   msg = 'math_check: R -> euler -> R failed'

 case (700)
   msg = 'Singular matrix in stress iteration'

!    Error messages related to parsing of Abaqus input file
 case (900)
   msg = 'PARSE ERROR: Improper definition of nodes in input file (Nnodes < 2)'
 case (901)
   msg = 'PARSE ERROR: No Elements defined in input file (Nelems = 0)'
 case (902)
   msg = 'PARSE ERROR: No Element sets defined in input file (Atleast one *Elset must exist)'
 case (903)
   msg = 'PARSE ERROR: No Materials defined in input file (Look into section assigments)'
 case (904)
   msg = 'PARSE ERROR: No elements could be assigned for Elset: '//ext_msg
 case (905)
   msg = 'PARSE ERROR: Error in mesh_abaqus_map_materials'
 case (906)
   msg = 'PARSE ERROR: Error in mesh_abaqus_count_cpElements'
 case (907)
   msg = 'PARSE ERROR: Incorrect size of mesh_mapFEtoCPelem in mesh_abaqus_map_elements; Size cannot be zero'
 case (908)
   msg = 'PARSE ERROR: Incorrect size of mesh_mapFEtoCPnode in mesh_abaqus_map_nodes; Size cannot be zero'
 case (909)
   msg = 'PARSE ERROR: Incorrect size of mesh_node in mesh_abaqus_build_nodes; must be equal to mesh_Nnodes'
 case(910)
   msg = 'PARSE ERROR: Incorrect element type mapping in '//ext_msg
 
 
 case default
   msg = 'Unknown error number...'
 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,'(a38)') '+------------------------------------+'
 write(6,'(a38)') '+               error                +'
 write(6,'(a38)') '+                                    +'
 write(6,'(a2,a)')       '+ ',msg
 if (present(ext_msg))  write(6,*) '+ ',ext_msg
 if (present(e)) then
   if (present(i) .and. present(g)) then
     write(6,'(a12,x,i6,x,a2,x,i2,x,a5,x,i4,a2)') '+ at element',e,'IP',i,'grain',g,' +'
   else
     write(6,'(a17,i6,a14)') '+              at ',e,'             +'
   endif
 endif
 write(6,'(a38)') '+------------------------------------+'
 
 call flush(6)
 call quit(9000+ID)
 !$OMP END CRITICAL (write2out)

! ABAQUS returns in some cases
 return

 endsubroutine


!********************************************************************
! write warning statements to standard out
!********************************************************************
 subroutine IO_warning(ID,e,i,g,ext_msg)

 use prec, only: pInt
 implicit none

 integer(pInt), intent(in) :: ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=80) msg

 select case (ID)
 case (101)
   msg = '+    crystallite debugging off...    +'
 case (600)
   msg = '+  crystallite responds elastically  +'
 case (601)
   msg = '+      stiffness close to zero       +'
 case (650)
   msg = '+     polar decomposition failed     +'
 case (700)
   msg = '+      unknown crystal symmetry      +'
 case default
   msg = '+     unknown warning number...      +'
 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,'(a38)') '+------------------------------------+'
 write(6,'(a38)') '+              warning               +'
 write(6,'(a38)') '+                                    +'
 write(6,'(a38)') msg
 if (present(ext_msg))  write(6,*) '+ ',ext_msg
 if (present(e)) then
   if (present(i)) then
     if (present(g)) then
       write(6,'(a12,x,i6,x,a2,x,i2,x,a5,x,i4,a2)') '+ at element',e,'IP',i,'grain',g,' +'
     else
       write(6,'(a12,x,i6,x,a2,x,i2,a13)') '+ at element',e,'IP',i,'            +'
     endif
   else
     write(6,'(a12,x,i6,a19)') '+ at element',e,'             +'
   endif
 endif
 write(6,'(a38)') '+------------------------------------+'
 call flush(6)
 !$OMP END CRITICAL (write2out)

 endsubroutine
 
 END MODULE IO
