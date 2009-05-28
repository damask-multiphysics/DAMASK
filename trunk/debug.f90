
!##############################################################
 MODULE debug
!##############################################################
 use prec


 implicit none
 integer(pInt), dimension(nStress) :: debug_StressLoopDistribution = 0_pInt
 integer(pInt), dimension(nState) :: debug_StateLoopDistribution = 0_pInt
 integer(pInt), dimension(nState) :: debug_StiffnessStateLoopDistribution = 0_pInt
 integer(pInt), dimension(nCryst) :: debug_CrystalliteLoopDistribution = 0_pInt
 integer(pLongInt) :: debug_cumLpTicks = 0_pInt
 integer(pLongInt) :: debug_cumDotStateTicks = 0_pInt
 integer(pInt) :: debug_cumLpCalls = 0_pInt
 integer(pInt) :: debug_cumDotStateCalls = 0_pInt
 logical :: debugger = .false.
 logical :: distribution_init = .false.

 CONTAINS

!********************************************************************
! reset debug distributions
!********************************************************************
SUBROUTINE debug_reset()

 use prec
 implicit none

 debug_StressLoopDistribution         = 0_pInt ! initialize debugging data
 debug_StateLoopDistribution          = 0_pInt
 debug_StiffnessStateLoopDistribution = 0_pInt
 debug_CrystalliteLoopDistribution    = 0_pInt
 debug_cumLpTicks       = 0_pInt
 debug_cumDotStateTicks = 0_pInt
 debug_cumLpCalls       = 0_pInt
 debug_cumDotStateCalls = 0_pInt

END SUBROUTINE

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

 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_StressLoop :'
 do i=1,nStress
   if (debug_StressLoopDistribution(i) /= 0) then
     integral = integral + i*debug_StressLoopDistribution(i)
     write(6,'(i25,i10)') i,debug_StressLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,i10)') '          total',sum(debug_StressLoopDistribution),integral
 
 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_StateLoop :'
 do i=1,nState
   if (debug_StateLoopDistribution(i) /= 0) then
     integral = integral + i*debug_StateLoopDistribution(i)
     write(6,'(i25,i10)') i,debug_StateLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,i10)') '          total',sum(debug_StateLoopDistribution),integral

 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_StiffnessStateLoop :'
 do i=1,nState
   if (debug_StiffnessStateLoopDistribution(i) /= 0) then
     integral = integral + i*debug_StiffnessStateLoopDistribution(i)
     write(6,'(i25,i10)') i,debug_StiffnessStateLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,i10)') '          total',sum(debug_StiffnessStateLoopDistribution),integral
 
 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_CrystalliteLoop :'
 do i=1,nCryst
   if (debug_CrystalliteLoopDistribution(i) /= 0) then
     integral = integral + i*debug_CrystalliteLoopDistribution(i)
     write(6,'(i25,i10)') i,debug_CrystalliteLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,i10)') '          total',sum(debug_CrystalliteLoopDistribution),integral
 write(6,*)

 END SUBROUTINE
 
 END MODULE debug
