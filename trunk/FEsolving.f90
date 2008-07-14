
!##############################################################
 MODULE FEsolving
!##############################################################

 use prec, only: pInt,pReal
 implicit none

 integer(pInt) cycleCounter
 integer(pInt) theInc,theCycle,theLovl
 real(pReal)   theTime
 logical :: lastIncConverged = .false.,outdatedByNewInc = .false., outdatedFFN1 = .false.

 END MODULE FEsolving
