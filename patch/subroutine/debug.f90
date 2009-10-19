
!##############################################################
 MODULE debug
!##############################################################
 use prec


 implicit none
 integer(pInt), dimension(nCutback+1) :: debug_cutbackDistribution = 0_pInt
 integer(pInt), dimension(nInner) :: debug_InnerLoopDistribution = 0_pInt
 integer(pInt), dimension(nOuter) :: debug_OuterLoopDistribution = 0_pInt
 integer(pLongInt) :: debug_cumLpTicks = 0_pInt
 integer(pLongInt) :: debug_cumDotStateTicks = 0_pInt
 integer(pInt) :: debug_cumLpCalls = 0_pInt
 integer(pInt) :: debug_cumDotStateCalls = 0_pInt
 logical :: debugger = .false.
 logical :: distribution_init = .false.

 CONTAINS


!********************************************************************
! write debug statements to standard out
!********************************************************************
 SUBROUTINE debug_info()

 use prec
 implicit none

 integer(pInt) i,integral
 integer(pLongInt) tickrate

 write(6,*)
 write(6,*) 'DEBUG Info'
 write(6,*)
 write(6,'(a33,x,i12)')	'total calls to LpAndItsTangent  :',debug_cumLpCalls
 if (debug_cumLpCalls > 0_pInt) then
   call system_clock(count_rate=tickrate)
   write(6,'(a33,x,f12.6)') 'avg CPU time/microsecs per call :',dble(debug_cumLpTicks)/tickrate/1.0e-6_pReal/debug_cumLpCalls
   write(6,'(a33,x,i12)')   'total CPU ticks                 :',debug_cumLpTicks
 endif
 write(6,*)
 write(6,'(a33,x,i12)')	'total calls to dotState             :',debug_cumDotStateCalls
 if (debug_cumdotStateCalls > 0_pInt) then
   call system_clock(count_rate=tickrate)
   write(6,'(a33,x,f12.6)') 'avg CPU time/microsecs per call :',&
     dble(debug_cumDotStateTicks)/tickrate/1.0e-6_pReal/debug_cumDotStateCalls
   write(6,'(a33,x,i12)')   'total CPU ticks                 :',debug_cumDotStateTicks
 endif
 write(6,*)
 write(6,*)	'distribution_cutback :'
 do i=0,nCutback
   if (debug_cutbackDistribution(i+1) /= 0) write(6,*) i,debug_cutbackDistribution(i+1)
 enddo
 write(6,*) 'total',sum(debug_cutbackDistribution)
 write(6,*)
 
 integral = 0_pInt
 write(6,*)	'distribution_InnerLoop :'
 do i=1,nInner
   if (debug_InnerLoopDistribution(i) /= 0) then
     integral = integral + i*debug_InnerLoopDistribution(i)
     write(6,*) i,debug_InnerLoopDistribution(i)
   endif
 enddo
 write(6,*) 'total',sum(debug_InnerLoopDistribution),integral
 write(6,*)
 
 integral = 0_pInt
 write(6,*)	'distribution_OuterLoop :'
 do i=1,nOuter
   if (debug_OuterLoopDistribution(i) /= 0) then
     integral = integral + i*debug_OuterLoopDistribution(i)
     write(6,*) i,debug_OuterLoopDistribution(i)
   endif
 enddo
 write(6,*) 'total',sum(debug_OuterLoopDistribution),integral
 write(6,*)

 END SUBROUTINE
 
 END MODULE debug
