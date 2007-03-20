
!##############################################################
 MODULE IO   
!##############################################################

 CONTAINS
!---------------------------
! function IO_open_file(unit,relPath)
! function IO_open_inputFile(unit)
! function IO_stringPos(line,N)
! function IO_stringValue(line,positions,pos)
! function IO_floatValue(line,positions,pos)
! function IO_intValue(line,positions,pos)
! function IO_lowercase(line)
! subroutine IO_error(ID)
!---------------------------



!********************************************************************
! open existing file to given unit
! path to file is relative to working directory
!********************************************************************
 logical FUNCTION IO_open_file(unit,relPath)

 use prec, only: pInt
 implicit none

 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 character(len=*) relPath
 integer(pInt) unit
 character(256) path

 inquire(6, name=path) ! determine outputfile
 open(unit,status='old',err=100,file=path(1:scan(path,pathSep,back=.true.))//relPath)
 IO_open_file = .true.
 return
100 IO_open_file = .false.
 return
 END FUNTION


!********************************************************************
! open FEM inputfile to given unit
!********************************************************************
 logical FUNCTION IO_open_inputFile(unit)

 use prec, only: pReal, pInt
 implicit none

 character(256) outName
 integer(pInt) unit, extPos
 character(3) ext

 inquire(6, name=outName) ! determine outputfileName
 extPos = len_trim(outName)-2
 if(outName(extPos:extPos+2)=='out') then
     ext='dat' ! MARC
 else
     ext='inp' ! ABAQUS
 end if
 open(unit,status='old',err=100,file=outName(1:extPos-1)//ext)
 IO_open_inputFile = .true.
 return
100 IO_open_inputFile = .false.
 return

 END FUNCTION


!********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts found and
! their left/right positions to be used by IO_xxxVal
!********************************************************************
 FUNCTION IO_stringPos (line,N)

 use prec, only: pReal,pInt
 implicit none

 character(len=*) line
 character(len=*), parameter :: sep=achar(32)//achar(9) ! whitespaces
 integer(pInt) N, part
 integer(pInt) IO_stringPos(1+N*2)

 IO_stringPos = -1
 IO_stringPos(1) = 0
 part = 1
 do while ((N<1 .or. part<=N) .and. verify(line(IO_stringPos(part*2-1)+1:),sep)>0)
   IO_stringPos(part*2) = IO_stringPos(part*2-1)+verify(line(IO_stringPos(part*2-1)+1:),sep)
   IO_stringPos(part*2+1) = IO_stringPos(part*2)+scan(line(IO_stringPos(part*2):),sep)-2
   part = part+1
 end do
 IO_stringPos(1) = part-1
 return

 END FUNCTION


!********************************************************************
! read string value at pos from line
!********************************************************************
 FUNCTION IO_stringValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*) line
 integer(pInt) positions(*),pos
 character(len=1+positions(pos*2+1)-positions(pos*2)) IO_stringValue

 IO_stringValue = line(positions(pos*2):positions(pos*2+1))
 return

 END FUNCTION


!********************************************************************
! read float value at pos from line
!********************************************************************
 FUNCTION IO_floatValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*) line
 real(pReal) IO_floatValue
 integer(pInt) positions(*),pos

 READ(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT='(F)') IO_floatValue
 return
100 IO_floatValue = -1.0_pReal
 return

 END FUNCTION


!********************************************************************
! read int value at pos from line
!********************************************************************
 FUNCTION IO_intValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*) line
 integer(pInt) IO_intValue
 integer(pInt) positions(*),pos

 READ(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT='(I)') IO_intValue
 return
100 IO_intValue = -1_pInt
 return

 END FUNCTION


!********************************************************************
! change character in line to lower case
!********************************************************************
 FUNCTION IO_lowercase (line)

 use prec, only: pInt
 implicit none

 character (len=*) line
 character (len=len(line)) IO_lowercase
 integer(pInt) i

 IO_lowercase = line
 forall (i=1:len(line),64<ichar(line(i:i)).and.ichar(line(i:i))<91) IO_lowercase(i:i)=achar(ichar(line(i:i))+32)
 return 

 END FUNCTION


!********************************************************************
! in place change character in line to lower case
!********************************************************************
 SUBROUTINE IO_lowercaseInplace (line)

 use prec, only: pInt
 implicit none

 character (len=*) line
 integer(pInt) i

 forall (i=1:len(line),64<ichar(line(i:i)).and.ichar(line(i:i))<91) line(i:i)=achar(ichar(line(i:i))+32)
 return 

 END SUBROUTINE


!********************************************************************
! write error statements to standard out
! and terminate the Marc run with exit #9xxx
! in ABAQUS either time step is reduced or execution terminated
!********************************************************************
 SUBROUTINE IO_error(ID)

 use prec, only: pInt
 implicit none

 integer(pInt) ID
 character(len=80) msg

 select case (ID)
 case (100)
   msg='File material.mpie can not be opened'
 case (110)
   msg='File material99.mpie can not be opened'
 case (120)
   msg='File with c-coefficience can not be opened'
 case (130)
   msg='File with single orientations can not be opened'
 case (200)
   msg='Error reading from file material.mpie'
 case (210)
   msg='Error reading from file material99.mpie'
 case (220)
   msg='Error reading from file containing c-coefficiences'
 case (230)
   msg='Error reading from file containing single orientations'
 case (300)
   msg='This material can only be used with &
  &elements with three direct stress components'
 case (400)
   msg='Unknown alloy number specified'
 case (500)
   msg='Unknown lattice number specified'
 case (510)
   msg='Unknown component type specified'
 case (520)
   msg='Unknown component symmetry specified'
 case (600)
   msg='Stress iteration did not converge'
 case (700)
   msg='Singular matrix in stress iteration'
 case default
   msg='Unknown error number'
 end select
 
 write(6,*) 'MPIE Material Routine Ver. 0.7 by Dr. F. Roters'
 write(6,*)
 write(6,*) msg
 call flush(6)
!  call quit(9000+ID)
! ABAQUS returns in some cases
 return

 END SUBROUTINE


 END MODULE IO
