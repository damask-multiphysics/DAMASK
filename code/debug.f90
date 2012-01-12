! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
!* $Id$
!##############################################################
MODULE debug
!##############################################################
use prec

implicit none
character(len=64), parameter :: debug_configFile = 'debug.config' ! name of configuration file
integer(pInt), parameter :: debug_spectralGeneral    = 1_pInt, &
                            debug_spectralDivergence = 2_pInt, &
                            debug_spectralRestart    = 4_pInt

integer(pInt), dimension(:,:), allocatable :: debug_StressLoopDistribution
integer(pInt), dimension(:,:), allocatable :: debug_LeapfrogBreakDistribution
integer(pInt), dimension(:,:), allocatable :: debug_StateLoopDistribution
integer(pInt), dimension(:), allocatable ::   debug_CrystalliteLoopDistribution
integer(pInt), dimension(:), allocatable ::   debug_MaterialpointStateLoopDistribution
integer(pInt), dimension(:), allocatable ::   debug_MaterialpointLoopDistribution
integer(pLongInt) :: debug_cumLpTicks             = 0_pInt
integer(pLongInt) :: debug_cumDotStateTicks       = 0_pInt
integer(pLongInt) :: debug_cumDotTemperatureTicks = 0_pInt
integer(pInt) :: debug_cumLpCalls             = 0_pInt
integer(pInt) :: debug_cumDotStateCalls       = 0_pInt
integer(pInt) :: debug_cumDotTemperatureCalls = 0_pInt
integer(pInt) :: debug_e = 1_pInt
integer(pInt) :: debug_i = 1_pInt
integer(pInt) :: debug_g = 1_pInt
integer(pInt), dimension(2) :: debug_stressMaxLocation = 0_pInt
integer(pInt), dimension(2) :: debug_stressMinLocation = 0_pInt
integer(pInt), dimension(2) :: debug_jacobianMaxLocation = 0_pInt
integer(pInt), dimension(2) :: debug_jacobianMinLocation = 0_pInt
real(pReal) :: debug_stressMax
real(pReal) :: debug_stressMin
real(pReal) :: debug_jacobianMax
real(pReal) :: debug_jacobianMin
logical :: debug_selectiveDebugger = .true.
integer(pInt) :: debug_verbosity = 1_pInt
integer(pInt) :: debug_spectral = 0_pInt

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
  use IO,       only: IO_error, &
                      IO_open_file, &
                      IO_isBlank, &
                      IO_stringPos, &
                      IO_stringValue, &
                      IO_lc, &
                      IO_floatValue, &
                      IO_intValue
  implicit none
  
  !*** input variables ***!
  
  !*** output variables ***!
  
  !*** local variables ***!
  integer(pInt), parameter ::                 fileunit = 300  
  integer(pInt), parameter ::                 maxNchunks = 2
  integer(pInt), dimension(1+2*maxNchunks) :: positions
  character(len=64)                           tag
  character(len=1024)                         line
  
  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '<<<+-  debug init  -+>>>'
    write(6,*) '$Id$'
    write(6,*)
  !$OMP END CRITICAL (write2out)
  
  allocate(debug_StressLoopDistribution(nStress,2)) ;            debug_StressLoopDistribution             = 0_pInt
  allocate(debug_LeapfrogBreakDistribution(nStress,2)) ;         debug_LeapfrogBreakDistribution          = 0_pInt
  allocate(debug_StateLoopDistribution(nState,2)) ;              debug_StateLoopDistribution              = 0_pInt
  allocate(debug_CrystalliteLoopDistribution(nCryst+1)) ;        debug_CrystalliteLoopDistribution        = 0_pInt
  allocate(debug_MaterialpointStateLoopDistribution(nMPstate)) ; debug_MaterialpointStateLoopDistribution = 0_pInt
  allocate(debug_MaterialpointLoopDistribution(nHomog+1)) ;      debug_MaterialpointLoopDistribution      = 0_pInt
  
  ! try to open the config file
  if(IO_open_file(fileunit,debug_configFile)) then 
  
    line = ''
    ! read variables from config file and overwrite parameters
    do
      read(fileunit,'(a1024)',END=100) line
      if (IO_isBlank(line)) cycle                           ! skip empty lines
      positions = IO_stringPos(line,maxNchunks)
      tag = IO_lc(IO_stringValue(line,positions,1))         ! extract key
      select case(tag)
        case ('element','e','el')
              debug_e = IO_intValue(line,positions,2)
        case ('integrationpoint','i','ip')
              debug_i = IO_intValue(line,positions,2)
        case ('grain','g','gr')
              debug_g = IO_intValue(line,positions,2)
        case ('selective')
              debug_selectiveDebugger = IO_intValue(line,positions,2) > 0_pInt
        case ('verbosity')
              debug_verbosity = IO_intValue(line,positions,2)
        case ('(spectral)')
            select case(IO_lc(IO_stringValue(line,positions,2)))
              case('general')
                   debug_spectral = ior(debug_spectral, debug_spectralGeneral)
              case('divergence')
                   debug_spectral = ior(debug_spectral, debug_spectralDivergence)
              case('restart')
                   debug_spectral = ior(debug_spectral, debug_spectralRestart)
            endselect
      endselect
    enddo
    100 close(fileunit)
    
    if (debug_verbosity > 0) then
      !$OMP CRITICAL (write2out)
        write(6,*) '   ... using values from config file'
        write(6,*)
      !$OMP END CRITICAL (write2out)
    endif
    
  ! no config file, so we use standard values
  else 

    if (debug_verbosity > 0) then
      !$OMP CRITICAL (write2out)
        write(6,*) '   ... using standard values'
        write(6,*)
      !$OMP END CRITICAL (write2out)
    endif
    
  endif  

  if (debug_verbosity > 0) then
    !$OMP CRITICAL (write2out)
      write(6,'(a24,x,i1)')   'verbose:                ',debug_verbosity
      write(6,'(a24,x,l)')    'selective:              ',debug_selectiveDebugger
    !$OMP END CRITICAL (write2out)
  endif
  if (debug_selectiveDebugger) then
    if (debug_verbosity > 0) then
      !$OMP CRITICAL (write2out)
        write(6,'(a24,x,i8)') 'element:              ',debug_e
        write(6,'(a24,x,i8)') 'ip:                   ',debug_i
        write(6,'(a24,x,i8)') 'grain:                ',debug_g
      !$OMP END CRITICAL (write2out)
    endif
  else
    debug_e = 0_pInt                                                            ! switch off selective debugging
    debug_i = 0_pInt
    debug_g = 0_pInt
  endif
  !$OMP CRITICAL (write2out)                                                   ! bitwise coded
  if (iand(debug_spectral,debug_spectralGeneral)    > 0_pInt)  write(6,'(a)') ' spectral general    debugging'
  if (iand(debug_spectral,debug_spectralDivergence) > 0_pInt)  write(6,'(a)') ' spectral divergence debugging'
  if (iand(debug_spectral,debug_spectalRestart)     > 0_pInt)  write(6,'(a)') ' spectral restart    debugging'
  !$OMP END CRITICAL (write2out)

endsubroutine
 
!********************************************************************
! reset debug distributions
!********************************************************************
subroutine debug_reset()

  use prec
  implicit none

  debug_StressLoopDistribution              = 0_pInt ! initialize debugging data
  debug_LeapfrogBreakDistribution           = 0_pInt
  debug_StateLoopDistribution               = 0_pInt
  debug_CrystalliteLoopDistribution         = 0_pInt
  debug_MaterialpointStateLoopDistribution  = 0_pInt
  debug_MaterialpointLoopDistribution       = 0_pInt
  debug_cumLpTicks             = 0_pInt
  debug_cumDotStateTicks       = 0_pInt
  debug_cumDotTemperatureTicks = 0_pInt
  debug_cumLpCalls             = 0_pInt
  debug_cumDotStateCalls       = 0_pInt
  debug_cumDotTemperatureCalls = 0_pInt
  debug_stressMaxLocation = 0_pInt
  debug_stressMinLocation = 0_pInt
  debug_jacobianMaxLocation = 0_pInt
  debug_jacobianMinLocation = 0_pInt
  debug_stressMax = -huge(1.0_pReal)
  debug_stressMin = huge(1.0_pReal)
  debug_jacobianMax = -huge(1.0_pReal)
  debug_jacobianMin = huge(1.0_pReal)


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

  if (debug_verbosity > 4) then
    !$OMP CRITICAL (write2out)
  
      write(6,*)
      write(6,*) 'DEBUG Info (from previous cycle)'
      write(6,*)
      write(6,'(a33,x,i12)')      'total calls to LpAndItsTangent  :',debug_cumLpCalls
      if (debug_cumLpCalls > 0_pInt) then
        write(6,'(a33,x,f12.3)')  'total CPU time/s                :',dble(debug_cumLpTicks)/tickrate
        write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
          dble(debug_cumLpTicks)*1.0e6_pReal/tickrate/debug_cumLpCalls
      endif
      write(6,*)
      write(6,'(a33,x,i12)')      'total calls to collectDotState  :',debug_cumDotStateCalls
      if (debug_cumdotStateCalls > 0_pInt) then
        write(6,'(a33,x,f12.3)')  'total CPU time/s                :',dble(debug_cumDotStateTicks)/tickrate
        write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
          dble(debug_cumDotStateTicks)*1.0e6_pReal/tickrate/debug_cumDotStateCalls
      endif
      write(6,*)
      write(6,'(a33,x,i12)')      'total calls to dotTemperature   :',debug_cumDotTemperatureCalls
      if (debug_cumdotTemperatureCalls > 0_pInt) then
        write(6,'(a33,x,f12.3)')  'total CPU time/s                :', dble(debug_cumDotTemperatureTicks)/tickrate
        write(6,'(a33,x,f12.6)')  'avg CPU time/microsecs per call :',&
          dble(debug_cumDotTemperatureTicks)*1.0e6_pReal/tickrate/debug_cumDotTemperatureCalls
      endif
    
      integral = 0_pInt
      write(6,*)
      write(6,*)
      write(6,*) 'distribution_StressLoop :    stress  frogbreak  stiffness  frogbreak'
      do i=1,nStress
        if (any(debug_StressLoopDistribution(i,:)     /= 0_pInt ) .or. &
            any(debug_LeapfrogBreakDistribution(i,:)  /= 0_pInt ) ) then
          integral = integral + i*debug_StressLoopDistribution(i,1) + i*debug_StressLoopDistribution(i,2)
          write(6,'(i25,x,i10,x,i10,x,i10,x,i10)')   i,debug_StressLoopDistribution(i,1),debug_LeapfrogBreakDistribution(i,1), &
                                                 debug_StressLoopDistribution(i,2),debug_LeapfrogBreakDistribution(i,2)
        endif
      enddo
      write(6,'(a15,i10,x,i10,12x,i10)') '          total',integral,&
                                                           sum(debug_StressLoopDistribution(:,1)), &
                                                           sum(debug_StressLoopDistribution(:,2))
      
      integral = 0_pInt
      write(6,*)
      write(6,*) 'distribution_CrystalliteStateLoop :'
      do i=1,nState
        if (any(debug_StateLoopDistribution(i,:) /= 0)) then
          integral = integral + i*debug_StateLoopDistribution(i,1) + i*debug_StateLoopDistribution(i,2)
          write(6,'(i25,x,i10,12x,i10)') i,debug_StateLoopDistribution(i,1),debug_StateLoopDistribution(i,2)
        endif
      enddo
      write(6,'(a15,i10,x,i10,12x,i10)') '          total',integral,&
                                                         sum(debug_StateLoopDistribution(:,1)), &
                                                         sum(debug_StateLoopDistribution(:,2))
     
      integral = 0_pInt
      write(6,*)
      write(6,*) 'distribution_CrystalliteCutbackLoop :'
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

    !$OMP END CRITICAL (write2out)
  endif
      
  if (debug_verbosity > 2) then
    !$OMP CRITICAL (write2out)
  
      integral = 0_pInt
      write(6,*)
      write(6,*) 'distribution_MaterialpointStateLoop :'
      do i=1,nMPstate
        if (debug_MaterialpointStateLoopDistribution(i) /= 0) then
          integral = integral + i*debug_MaterialpointStateLoopDistribution(i)
          write(6,'(i25,x,i10)') i,debug_MaterialpointStateLoopDistribution(i)
        endif
      enddo
      write(6,'(a15,i10,x,i10)') '          total',integral,sum(debug_MaterialpointStateLoopDistribution) 
     
      integral = 0_pInt
      write(6,*)
      write(6,*) 'distribution_MaterialpointCutbackLoop :'
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
      write(6,*)
      write(6,*) 'Extreme values of returned stress and jacobian'
      write(6,*)
      write(6,'(a39)') '                      value     el   ip'
      write(6,'(a14,x,e12.3,x,i6,x,i4)') 'stress   min :', debug_stressMin, debug_stressMinLocation
      write(6,'(a14,x,e12.3,x,i6,x,i4)') '         max :', debug_stressMax, debug_stressMaxLocation
      write(6,'(a14,x,e12.3,x,i6,x,i4)') 'jacobian min :', debug_jacobianMin, debug_jacobianMinLocation
      write(6,'(a14,x,e12.3,x,i6,x,i4)') '         max :', debug_jacobianMax, debug_jacobianMaxLocation  
      write(6,*)

    !$OMP END CRITICAL (write2out)
  endif

endsubroutine
 
END MODULE debug
