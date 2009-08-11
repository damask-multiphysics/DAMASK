
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) cycleCounter
 integer(pInt) theInc,theCycle,theLovl
 real(pReal)   theTime
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false.,outdatedFFN1 = .false.,terminallyIll = .false.
 logical :: symmetricSolver = .false. 
 logical :: parallelExecution = .true. 
 integer(pInt), dimension(:,:), allocatable :: FEsolving_execIP
 integer(pInt), dimension(2) :: FEsolving_execElem


 CONTAINS

!***********************************************************
! determine wether a symmetric solver is used 
!***********************************************************
 subroutine FE_init()
 
 use prec, only: pInt
 use IO
 implicit none
 
 integer(pInt), parameter :: fileunit = 222
 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=1024) line

 write(6,*)
 write(6,*) '<<<+-  FEsolving init  -+>>>'
 write(6,*)
 
 if (IO_open_inputFile(fileunit)) then
 
   rewind(fileunit)
   do
     read (fileunit,'(a1024)',END=100) line
     positions = IO_stringPos(line,1)
     if( IO_lc(IO_stringValue(line,positions,1)) == 'solver' ) then
       read (fileunit,'(a1024)',END=100) line  ! Garbage line
       positions = IO_stringPos(line,2)
       symmetricSolver = (IO_intValue(line,positions,2) /= 1_pInt)
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
