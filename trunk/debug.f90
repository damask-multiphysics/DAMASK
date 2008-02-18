
!##############################################################
 MODULE debug
!##############################################################
 use prec

 implicit none
 integer(pInt), dimension(nCutback+1) :: debug_cutbackDistribution
 integer(pInt), dimension(nInner) :: debug_innerLoopDistribution
 integer(pInt), dimension(nOuter) :: debug_outerLoopDistribution

 CONTAINS


!********************************************************************
! write debug statements to standard out
!********************************************************************
 SUBROUTINE debug_info()

 use prec
 implicit none

 integer(pInt) i

 write(6,*) 'DEBUG Info'
 write(6,*)	'distribution_cutback :'
 do i=0,nCutback
   if (debug_cutbackDistribution(i+1) > 0) write(6,*) i,debug_cutbackDistribution(i+1)
 enddo
 write(6,*) 'total',sum(debug_cutbackDistribution)
 write(6,*)
 
 write(6,*)	'distribution_innerLoop :'
 do i=1,nInner
   if (debug_innerLoopDistribution(i) > 0) write(6,*) i,debug_innerLoopDistribution(i)
 enddo
 write(6,*) 'total',sum(debug_innerLoopDistribution)
 write(6,*)
 
 write(6,*)	'distribution_outerLoop :'
 do i=1,nOuter
   if (debug_outerLoopDistribution(i) > 0) write(6,*) i,debug_outerLoopDistribution(i)
 enddo
 write(6,*) 'total',sum(debug_outerLoopDistribution)
 write(6,*)

 END SUBROUTINE
 
 END MODULE debug
