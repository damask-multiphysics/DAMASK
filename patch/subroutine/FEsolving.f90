
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) cycleCounter
 integer(pInt) theInc,theCycle,theLovl
 real(pReal)   theTime
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false.,outdatedFFN1 = .false.
 logical :: symmetricSolver = .false. 
 logical :: parallelExecution = .true. 


 CONTAINS

!***********************************************************
! determine wether a symmetric solver is used 
!***********************************************************
 subroutine FE_init()
 
 use prec, only: pInt
 use IO
 implicit none
 
 integer(pInt), parameter :: fileunit = 222
 integer(pInt), dimension (1+2*2) :: pos
 character(len=1024) line

 if (IO_open_inputFile(fileunit)) then
 
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=100) line
     pos = IO_stringPos(line,1)
     if( IO_lc(IO_stringValue(line,pos,1)) == 'solver' ) then
       read (fileunit,'(a1024)',END=100) line  ! Garbage line
       pos = IO_stringPos(line,2)
       symmetricSolver = (IO_intValue(line,pos,2) /= 1_pInt)
       exit
     endif
   enddo
 else
   call IO_error(100) ! cannot open input file
 endif

100 close(fileunit)

 return

 end subroutine

 END MODULE FEsolving
