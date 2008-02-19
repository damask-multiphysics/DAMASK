
!##############################################################
 MODULE debug
!##############################################################
 use prec


 implicit none
 integer(pInt), dimension(nCutback+1) :: debug_cutbackDistribution
 integer(pInt), dimension(nInner) :: debug_InnerLoopDistribution
 integer(pInt), dimension(nOuter) :: debug_OuterLoopDistribution
 logical :: debugger = .false.

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
   if (debug_cutbackDistribution(i+1) /= 0) write(6,*) i,debug_cutbackDistribution(i+1)
 enddo
 write(6,*) 'total',sum(debug_cutbackDistribution)
 write(6,*)
 
 write(6,*)	'distribution_InnerLoop :'
 do i=1,nInner
   if (debug_InnerLoopDistribution(i) /= 0) write(6,*) i,debug_InnerLoopDistribution(i)
 enddo
 write(6,*) 'total',sum(debug_InnerLoopDistribution)
 write(6,*)
 
 write(6,*)	'distribution_OuterLoop :'
 do i=1,nOuter
   if (debug_OuterLoopDistribution(i) /= 0) write(6,*) i,debug_OuterLoopDistribution(i)
 enddo
 write(6,*) 'total',sum(debug_OuterLoopDistribution)
 write(6,*)

 END SUBROUTINE
 
 END MODULE debug
