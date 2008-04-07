
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) cycleCounter
 integer(pInt) theInc,theCycle,theLovl
 real(pReal)   theTime
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false.

 END MODULE FEsolving
