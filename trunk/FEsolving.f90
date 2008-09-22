
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) cycleCounter
 integer(pInt) theInc,theCycle,theLovl
 real(pReal)   theTime
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false., outdatedFFN1 = .false.
 logical :: symmetricSolver = .false. 

 CONTAINS

!***********************************************************
! determine wether a symmetric solver is used 
!***********************************************************
 subroutine FE_get_solverSymmetry(unit)
 
 use prec, only: pInt
 use IO
 implicit none
 
 integer(pInt) unit
 integer(pInt), dimension (133) :: pos
 character*300 line
 
610 FORMAT(A300)
 
 rewind(unit)
 do
   read (unit,610,END=630) line
   pos = IO_stringPos(line,1)
   if( IO_lc(IO_stringValue(line,pos,1)) == 'solver' ) then
     read (unit,610,END=630) line  ! Garbage line
     pos = IO_stringPos(line,2)  ! limit to 64 nodes max (plus ID, type)
     if(IO_intValue(line,pos,2) /= 1_pInt) symmetricSolver = .true.
!$OMP CRITICAL (write2out)
     write (6,*)
     write (6,*) 'Symmetric solver detected. d-Matrix will be symmetrized!'
!$OMP END CRITICAL (write2out)
   endif
 enddo

630 return

 end subroutine

 END MODULE FEsolving
