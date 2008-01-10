
!##############################################################
 MODULE debug
!##############################################################
 use prec

 implicit none
 integer(pInt), dimension(nCutback) :: debug_cutbackDistribution
 integer(pInt), dimension(nStress) :: debug_stressLoopDistribution
 integer(pInt), dimension(nState) :: debug_stateLoopDistribution

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
 do i=1,nCutback
   if (debug_cutbackDistribution(i) > 0) write(6,*) i,debug_cutbackDistribution(i)
 enddo
 write(6,*)
 
 write(6,*)	'distribution_stressLoop :'
 do i=1,nStress
   if (debug_stressLoopDistribution(i) > 0) write(6,*) i,debug_stressLoopDistribution(i)
 enddo
 write(6,*)
 
 write(6,*)	'distribution_stateLoop :'
 do i=1,nState
   if (debug_stateLoopDistribution(i) > 0) write(6,*) i,debug_stateLoopDistribution(i)
 enddo
 write(6,*)

 END SUBROUTINE
 
 END MODULE debug
