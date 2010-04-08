!* $Id$
!##############################################################
 MODULE debug
!##############################################################
 use prec

 implicit none
 integer(pInt), dimension(:), allocatable :: debug_StressLoopDistribution
 integer(pInt), dimension(:), allocatable :: debug_CrystalliteStateLoopDistribution
 integer(pInt), dimension(:), allocatable :: debug_StiffnessStateLoopDistribution
 integer(pInt), dimension(:), allocatable :: debug_CrystalliteLoopDistribution
 integer(pInt), dimension(:), allocatable :: debug_MaterialpointStateLoopDistribution
 integer(pInt), dimension(:), allocatable :: debug_MaterialpointLoopDistribution
 integer(pLongInt) :: debug_cumLpTicks = 0_pInt
 integer(pLongInt) :: debug_cumDotStateTicks = 0_pInt
 integer(pLongInt) :: debug_cumDotTemperatureTicks = 0_pInt
 integer(pInt) :: debug_cumLpCalls = 0_pInt
 integer(pInt) :: debug_cumDotStateCalls = 0_pInt
 integer(pInt) :: debug_cumDotTemperatureCalls = 0_pInt
 integer(pInt) :: debug_e = 1_pInt
 integer(pInt) :: debug_i = 1_pInt
 integer(pInt) :: debug_g = 1_pInt
 logical :: selectiveDebugger = .false.
 logical :: verboseDebugger = .false.
 logical :: debugger = .false.
 logical :: distribution_init = .false.

 CONTAINS


!********************************************************************
! initialize the debugging capabilities
!********************************************************************
subroutine debug_init()
  
  use prec,     only: pInt  
  use numerics, only: nStress, &
                      nState, &
                      nCryst, &
                      nMPstate, &
                      nHomog
  implicit none
  
  write(6,*)
  write(6,*) '<<<+-  debug init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
 
  allocate(debug_StressLoopDistribution(nStress)) ;         debug_StressLoopDistribution = 0_pInt
  allocate(debug_CrystalliteStateLoopDistribution(nState)); debug_CrystalliteStateLoopDistribution = 0_pInt
  allocate(debug_StiffnessStateLoopDistribution(nState)) ;  debug_StiffnessStateLoopDistribution = 0_pInt
  allocate(debug_CrystalliteLoopDistribution(nCryst+1)) ;   debug_CrystalliteLoopDistribution = 0_pInt
  allocate(debug_MaterialpointStateLoopDistribution(nMPstate)) ; debug_MaterialpointStateLoopDistribution = 0_pInt
  allocate(debug_MaterialpointLoopDistribution(nHomog+1)) ;      debug_MaterialpointLoopDistribution = 0_pInt
endsubroutine
 
!********************************************************************
! reset debug distributions
!********************************************************************
subroutine debug_reset()

  use prec
  implicit none

  debug_StressLoopDistribution              = 0_pInt ! initialize debugging data
  debug_CrystalliteStateLoopDistribution    = 0_pInt
  debug_StiffnessStateLoopDistribution      = 0_pInt
  debug_CrystalliteLoopDistribution         = 0_pInt
  debug_MaterialpointStateLoopDistribution  = 0_pInt
  debug_MaterialpointLoopDistribution       = 0_pInt
  debug_cumLpTicks             = 0_pInt
  debug_cumDotStateTicks       = 0_pInt
  debug_cumDotTemperatureTicks = 0_pInt
  debug_cumLpCalls             = 0_pInt
  debug_cumDotStateCalls       = 0_pInt
  debug_cumDotTemperatureCalls = 0_pInt

endsubroutine

!********************************************************************
! write debug statements to standard out
!********************************************************************
 subroutine debug_info()

 use prec
 use numerics, only: nStress, &
                     nState, &
                     nCryst, &
                     nMPstate, &
                     nHomog
 implicit none

 integer(pInt)       i,integral
 integer(pLongInt)   tickrate
 
 call system_clock(count_rate=tickrate)

 write(6,*)
 write(6,*) 'DEBUG Info'
 write(6,*)
 write(6,'(a33,x,i12)')	     'total calls to LpAndItsTangent  :',debug_cumLpCalls
 if (debug_cumLpCalls > 0_pInt) then
   write(6,'(a33,x,f12.3)')  'total CPU time/s                :',dble(debug_cumLpTicks)/tickrate
   write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
     dble(debug_cumLpTicks)*1.0e6_pReal/tickrate/debug_cumLpCalls
 endif
 write(6,*)
 write(6,'(a33,x,i12)')	     'total calls to collectDotState  :',debug_cumDotStateCalls
 if (debug_cumdotStateCalls > 0_pInt) then
   write(6,'(a33,x,f12.3)')  'total CPU time/s                :',dble(debug_cumDotStateTicks)/tickrate
   write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
     dble(debug_cumDotStateTicks)*1.0e6_pReal/tickrate/debug_cumDotStateCalls
 endif
 write(6,*)
 write(6,'(a33,x,i12)')	     'total calls to dotTemperature   :',debug_cumDotTemperatureCalls
 if (debug_cumdotTemperatureCalls > 0_pInt) then
   write(6,'(a33,x,f12.3)')  'total CPU time/s                :', dble(debug_cumDotTemperatureTicks)/tickrate
   write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
     dble(debug_cumDotTemperatureTicks)*1.0e6_pReal/tickrate/debug_cumDotTemperatureCalls
 endif

 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_StressLoop :'
 do i=1,nStress
   if (debug_StressLoopDistribution(i) /= 0) then
     integral = integral + i*debug_StressLoopDistribution(i)
     write(6,'(i25,x,i10)') i,debug_StressLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_StressLoopDistribution)
 
 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_CrystalliteStateLoop :'
 do i=1,nState
   if (debug_CrystalliteStateLoopDistribution(i) /= 0) then
     integral = integral + i*debug_CrystalliteStateLoopDistribution(i)
     write(6,'(i25,x,i10)') i,debug_CrystalliteStateLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_CrystalliteStateLoopDistribution)

 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_StiffnessStateLoop :'
 do i=1,nState
   if (debug_StiffnessStateLoopDistribution(i) /= 0) then
     integral = integral + i*debug_StiffnessStateLoopDistribution(i)
     write(6,'(i25,x,i10)') i,debug_StiffnessStateLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_StiffnessStateLoopDistribution)
 
 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_CrystalliteLoop :'
 do i=1,nCryst+1
   if (debug_CrystalliteLoopDistribution(i) /= 0) then
     integral = integral + i*debug_CrystalliteLoopDistribution(i)
     if (i <= nCryst) then
       write(6,'(i25,x,i10)') i,debug_CrystalliteLoopDistribution(i)
     else
       write(6,'(i25,a1,i10)') i-1,'+',debug_CrystalliteLoopDistribution(i)
     endif
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_CrystalliteLoopDistribution)
 write(6,*)

!* Material point loop counter <<<updated 31.07.2009>>>
 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_MaterialpointStateLoop :'
 do i=1,nMPstate
   if (debug_MaterialpointStateLoopDistribution(i) /= 0) then
     integral = integral + i*debug_MaterialpointStateLoopDistribution(i)
     write(6,'(i25,x,i10)') i,debug_MaterialpointStateLoopDistribution(i)
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_MaterialpointStateLoopDistribution)
 write(6,*)

 integral = 0_pInt
 write(6,*)
 write(6,*)	'distribution_MaterialpointLoop :'
 do i=1,nHomog+1
   if (debug_MaterialpointLoopDistribution(i) /= 0) then
     integral = integral + i*debug_MaterialpointLoopDistribution(i)
     if (i <= nHomog) then
       write(6,'(i25,x,i10)') i,debug_MaterialpointLoopDistribution(i)
     else
       write(6,'(i25,a1,i10)') i-1,'+',debug_MaterialpointLoopDistribution(i)
     endif
   endif
 enddo
 write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_MaterialpointLoopDistribution)
 write(6,*)


 endsubroutine
 
 END MODULE debug
