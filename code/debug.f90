! Copyright 2011,2012 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
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
!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Reading in and interpretating the debugging settings for the various modules
!--------------------------------------------------------------------------------------------------
module debug
 use prec, only: &
   pInt, &
   pReal, &
   pLongInt

 implicit none
 private
 
 integer(pInt), parameter, public :: &
   debug_levelSelective     = 2_pInt**0_pInt, &
   debug_levelBasic         = 2_pInt**1_pInt, &
   debug_levelExtensive     = 2_pInt**2_pInt
 integer(pInt), parameter, private :: &
   debug_maxGeneral         = debug_levelExtensive                                                 ! must be set to the last bitcode used by (potentially) all debug types
 integer(pInt), parameter, public :: &
   debug_spectralRestart    = debug_maxGeneral*2_pInt**1_pInt, &
   debug_spectralFFTW       = debug_maxGeneral*2_pInt**2_pInt, &
   debug_spectralDivergence = debug_maxGeneral*2_pInt**3_pInt, &
   debug_spectralRotation   = debug_maxGeneral*2_pInt**4_pInt, &
   debug_spectralPETSc      = debug_maxGeneral*2_pInt**5_pInt

 integer(pInt), parameter, public :: &
   debug_debug                   =  1_pInt, &
   debug_math                    =  2_pInt, &
   debug_FEsolving               =  3_pInt, &
   debug_mesh                    =  4_pInt, &                                                       !< stores debug level for mesh part of DAMASK bitwise coded
   debug_material                =  5_pInt, &                                                       !< stores debug level for material part of DAMASK bitwise coded
   debug_lattice                 =  6_pInt, &                                                       !< stores debug level for lattice part of DAMASK bitwise coded
   debug_constitutive            =  7_pInt, &                                                       !< stores debug level for constitutive part of DAMASK bitwise coded
   debug_crystallite             =  8_pInt, &
   debug_homogenization          =  9_pInt, &
   debug_CPFEM                   = 10_pInt, &
   debug_spectral                = 11_pInt, &
   debug_abaqus                  = 12_pInt
 integer(pInt), parameter, private :: &
   debug_maxNtype                = debug_abaqus                                                     ! must be set to the maximum defined debug type
   
 integer(pInt),protected, dimension(debug_maxNtype+2_pInt),  public :: &                            ! specific ones, and 2 for "all" and "other"
   debug_level                    = 0_pInt

 integer(pInt), public :: &
   debug_cumLpCalls              = 0_pInt, &                                                        ! total number of calls to LpAndItsTangent
   debug_cumDeltaStateCalls      = 0_pInt, &                                                        ! total number of calls to deltaState
   debug_cumDotStateCalls        = 0_pInt, &                                                        ! total number of calls to dotState
   debug_cumDotTemperatureCalls  = 0_pInt, &                                                        ! total number of calls to dotTemprature
   debug_e                       = 1_pInt, &
   debug_i                       = 1_pInt, &
   debug_g                       = 1_pInt

 integer(pLongInt), public :: &
   debug_cumLpTicks              = 0_pLongInt, &                                                    ! total cpu ticks spent in LpAndItsTangent
   debug_cumDeltaStateTicks      = 0_pLongInt, &                                                    ! total cpu ticks spent in deltaState
   debug_cumDotStateTicks        = 0_pLongInt, &                                                    ! total cpu ticks spent in dotState
   debug_cumDotTemperatureTicks  = 0_pLongInt                                                       ! total cpu ticks spent in dotTemperature
 
 integer(pInt), dimension(2), public :: &
   debug_stressMaxLocation       = 0_pInt, &
   debug_stressMinLocation       = 0_pInt, &
   debug_jacobianMaxLocation     = 0_pInt, &
   debug_jacobianMinLocation     = 0_pInt

 integer(pInt), dimension(:), allocatable, public :: &
   debug_CrystalliteLoopDistribution, &                                                             ! distribution of crystallite cutbacks
   debug_MaterialpointStateLoopDistribution, &
   debug_MaterialpointLoopDistribution

 integer(pInt), dimension(:,:), allocatable, public :: &
   debug_StressLoopDistribution, &                                                                  ! distribution of stress iterations until convergence
   debug_StateLoopDistribution                                                                      ! distribution of state iterations until convergence
 
 real(pReal), public :: &
   debug_stressMax               = -huge(1.0_pReal), &
   debug_stressMin               =  huge(1.0_pReal), &
   debug_jacobianMax             = -huge(1.0_pReal), &
   debug_jacobianMin             =  huge(1.0_pReal)
 
 character(len=64), parameter, private ::  &
   debug_configFile         = 'debug.config'                                                        ! name of configuration file

#ifdef PETSc 
 character(len=1024), parameter, public :: &
   PETScDebug = ' -snes_view -snes_monitor ' 
#endif
 public :: debug_init, &
           debug_reset, &
           debug_info

contains


!********************************************************************
! initialize the debugging capabilities
!********************************************************************
subroutine debug_init

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use numerics, only: nStress, &
                     nState, &
                     nCryst, &
                     nMPstate, &
                     nHomog
 use IO,       only: IO_error, &
                     IO_open_file_stat, &
                     IO_isBlank, &
                     IO_stringPos, &
                     IO_stringValue, &
                     IO_lc, &
                     IO_floatValue, &
                     IO_intValue

 implicit none
 integer(pInt), parameter                 :: fileunit    = 300_pInt  
 integer(pInt), parameter                 :: maxNchunks  = 7_pInt  
 
 integer(pInt)                            :: i, what
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 character(len=64)                        :: tag
 character(len=1024)                      :: line
 !$OMP CRITICAL (write2out)
   write(6,*)
   write(6,*) '<<<+-  debug init  -+>>>'
   write(6,*) '$Id$'
#include "compilation_info.f90"
 !$OMP END CRITICAL (write2out)
 
 if (allocated(debug_StressLoopDistribution)) &
    deallocate(debug_StressLoopDistribution)
      allocate(debug_StressLoopDistribution(nStress+1,2))
               debug_StressLoopDistribution = 0_pInt
 if (allocated(debug_StateLoopDistribution)) &
    deallocate(debug_StateLoopDistribution)
      allocate(debug_StateLoopDistribution(nState+1,2))
               debug_StateLoopDistribution = 0_pInt
 if (allocated(debug_CrystalliteLoopDistribution)) &
    deallocate(debug_CrystalliteLoopDistribution)
      allocate(debug_CrystalliteLoopDistribution(nCryst+1))
               debug_CrystalliteLoopDistribution = 0_pInt
 if (allocated(debug_MaterialpointStateLoopDistribution)) &
    deallocate(debug_MaterialpointStateLoopDistribution)
      allocate(debug_MaterialpointStateLoopDistribution(nMPstate))
               debug_MaterialpointStateLoopDistribution = 0_pInt
 if (allocated(debug_MaterialpointLoopDistribution)) &
    deallocate(debug_MaterialpointLoopDistribution)
      allocate(debug_MaterialpointLoopDistribution(nHomog+1))
               debug_MaterialpointLoopDistribution = 0_pInt
 
 
 ! try to open the config file
 if(IO_open_file_stat(fileunit,debug_configFile)) then
 
   ! read variables from config file and overwrite parameters
   do
     read(fileunit,'(a1024)',END=100) line
     if (IO_isBlank(line)) cycle                                                                    ! skip empty lines
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('element','e','el')
         debug_e = IO_intValue(line,positions,2_pInt)
       case ('integrationpoint','i','ip')
         debug_i = IO_intValue(line,positions,2_pInt)
       case ('grain','g','gr')
         debug_g = IO_intValue(line,positions,2_pInt)
     end select
     
     what = 0_pInt
     select case(tag)
       case ('debug')
         what = debug_debug
       case ('math')
         what = debug_math
       case ('fesolving', 'fe')
         what = debug_FEsolving
       case ('mesh')
         what = debug_mesh
       case ('material')
         what = debug_material
       case ('lattice')
         what = debug_lattice
       case ('constitutive')
         what = debug_constitutive
       case ('crystallite')
         what = debug_crystallite
       case ('homogenization')
         what = debug_homogenization
       case ('cpfem')
         what = debug_CPFEM
       case ('spectral')
         what = debug_spectral
       case ('abaqus')
         what = debug_abaqus
       case ('all')
         what = debug_maxNtype + 1_pInt
       case ('other')
         what = debug_maxNtype + 2_pInt
     end select
     if(what /= 0) then
       do i = 2_pInt, maxNchunks
         select case(IO_lc(IO_stringValue(line,positions,i)))
           case('basic')
             debug_level(what) = ior(debug_level(what), debug_levelBasic)
           case('extensive')
             debug_level(what) = ior(debug_level(what), debug_levelExtensive)
           case('selective')
             debug_level(what) = ior(debug_level(what), debug_levelSelective)
           case('restart')
             debug_level(what) = ior(debug_level(what), debug_spectralRestart)
           case('fft','fftw')
             debug_level(what) = ior(debug_level(what), debug_spectralFFTW)
           case('divergence')
             debug_level(what) = ior(debug_level(what), debug_spectralDivergence)
           case('rotation')
             debug_level(what) = ior(debug_level(what), debug_spectralRotation)
           case('petsc')
             debug_level(what) = ior(debug_level(what), debug_spectralPETSc)
         end select
       enddo
      endif
   enddo
   100 close(fileunit)
 
   do i = 1_pInt, debug_maxNtype
     if (debug_level(i) == 0) &
       debug_level(i) = ior(debug_level(i), debug_level(debug_maxNtype + 2_pInt))                         ! fill undefined debug types with levels specified by "other" 

       debug_level(i) = ior(debug_level(i), debug_level(debug_maxNtype + 1_pInt))                         ! fill all debug types with levels specified by "all" 
   enddo
  
   if (iand(debug_level(debug_debug),debug_levelBasic) /= 0) then
     !$OMP CRITICAL (write2out)
       write(6,*) 'using values from config file'
       write(6,*)
     !$OMP END CRITICAL (write2out)
   endif

 ! no config file, so we use standard values
 else 
   if (iand(debug_level(debug_debug),debug_levelBasic) /= 0) then
     !$OMP CRITICAL (write2out)
       write(6,*) 'using standard values'
       write(6,*)
     !$OMP END CRITICAL (write2out)
   endif
 endif

 !output switched on (debug level for debug must be extensive)
 if (iand(debug_level(debug_debug),debug_levelExtensive) /= 0) then
     do i = 1_pInt, debug_maxNtype
       select case(i)
         case (debug_debug)
           tag = 'Debug'
         case (debug_math)
           tag = 'Math'
         case (debug_FEsolving)
           tag = 'FEsolving'
         case (debug_mesh)
           tag = 'Mesh'
         case (debug_material)
           tag = 'Material'
         case (debug_lattice)
           tag = 'Lattice'
         case (debug_constitutive)
           tag = 'Constitutive'
         case (debug_crystallite)
           tag = 'Crystallite'
         case (debug_homogenization)
           tag = 'Homogenizaiton'
         case (debug_CPFEM)
           tag = 'CPFEM'
         case (debug_spectral)
           tag = 'Spectral solver'
         case (debug_abaqus)
           tag = 'ABAQUS FEM solver'
       end select
           
       if(debug_level(i) /= 0) then
   !$OMP CRITICAL (write2out)
         write(6,'(a,a)') tag,' debugging:'
         if(iand(debug_level(i),debug_levelBasic)        /= 0) write(6,'(a)') ' basic'
         if(iand(debug_level(i),debug_levelExtensive)    /= 0) write(6,'(a)') ' extensive'
         if(iand(debug_level(i),debug_levelSelective)    /= 0) then
           write(6,'(a)') 'selective on:'
           write(6,'(a24,1x,i8)') 'element:              ',debug_e
           write(6,'(a24,1x,i8)') 'ip:                   ',debug_i
           write(6,'(a24,1x,i8)') 'grain:                ',debug_g
         endif
         if(iand(debug_level(i),debug_spectralRestart)   /= 0) write(6,'(a)') ' restart'
         if(iand(debug_level(i),debug_spectralFFTW)      /= 0) write(6,'(a)') ' FFTW'
         if(iand(debug_level(i),debug_spectralDivergence)/= 0) write(6,'(a)') ' divergence'
         if(iand(debug_level(i),debug_spectralRotation)  /= 0) write(6,'(a)') ' rotation'
         if(iand(debug_level(i),debug_spectralPETSc)     /= 0) write(6,'(a)') ' PETSc'
   !$OMP END CRITICAL (write2out)
       endif
     enddo
 endif

end subroutine debug_init
 
!********************************************************************
! reset debug distributions
!********************************************************************
subroutine debug_reset

 implicit none

 debug_StressLoopDistribution              = 0_pInt ! initialize debugging data
 debug_StateLoopDistribution               = 0_pInt
 debug_CrystalliteLoopDistribution         = 0_pInt
 debug_MaterialpointStateLoopDistribution  = 0_pInt
 debug_MaterialpointLoopDistribution       = 0_pInt
 debug_cumLpTicks                          = 0_pLongInt
 debug_cumDeltaStateTicks                  = 0_pLongInt
 debug_cumDotStateTicks                    = 0_pLongInt
 debug_cumDotTemperatureTicks              = 0_pLongInt
 debug_cumLpCalls                          = 0_pInt
 debug_cumDeltaStateCalls                  = 0_pInt
 debug_cumDotStateCalls                    = 0_pInt
 debug_cumDotTemperatureCalls              = 0_pInt
 debug_stressMaxLocation                   = 0_pInt
 debug_stressMinLocation                   = 0_pInt
 debug_jacobianMaxLocation                 = 0_pInt
 debug_jacobianMinLocation                 = 0_pInt
 debug_stressMax                           = -huge(1.0_pReal)
 debug_stressMin                           =  huge(1.0_pReal)
 debug_jacobianMax                         = -huge(1.0_pReal)
 debug_jacobianMin                         =  huge(1.0_pReal)

end subroutine debug_reset

!********************************************************************
! write debug statements to standard out
!********************************************************************
subroutine debug_info

 use numerics, only: nStress, &
                     nState, &
                     nCryst, &
                     nMPstate, &
                     nHomog

 implicit none
 integer(pInt)     :: i,integral
 integer(pLongInt) :: tickrate
 character(len=1)  :: exceed

 call system_clock(count_rate=tickrate)

 !$OMP CRITICAL (write2out)
   if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0) then
     write(6,*)
     write(6,*) 'DEBUG Info (from previous cycle)'
     write(6,*)
     write(6,'(a33,1x,i12)')      'total calls to LpAndItsTangent  :',debug_cumLpCalls
     if (debug_cumLpCalls > 0_pInt) then
       write(6,'(a33,1x,f12.3)')  'total CPU time/s                :',real(debug_cumLpTicks,pReal)&
                                                                            /real(tickrate,pReal)
       write(6,'(a33,1x,f12.6)')  'avg CPU time/microsecs per call :',&
         real(debug_cumLpTicks,pReal)*1.0e6_pReal/real(tickrate,pReal)/real(debug_cumLpCalls,pReal)
     endif
     write(6,*)
     write(6,'(a33,1x,i12)')      'total calls to collectDotState  :',debug_cumDotStateCalls
     if (debug_cumdotStateCalls > 0_pInt) then
       write(6,'(a33,1x,f12.3)')  'total CPU time/s                :',real(debug_cumDotStateTicks,pReal)&
                                                                            /real(tickrate,pReal)
       write(6,'(a33,1x,f12.6)')  'avg CPU time/microsecs per call :',&
         real(debug_cumDotStateTicks,pReal)*1.0e6_pReal/real(tickrate,pReal)&
                                                                     /real(debug_cumDotStateCalls,pReal)
     endif
     write(6,*)
     write(6,'(a33,1x,i12)')      'total calls to collectDeltaState:',debug_cumDeltaStateCalls
     if (debug_cumDeltaStateCalls > 0_pInt) then
       write(6,'(a33,1x,f12.3)')  'total CPU time/s                :',real(debug_cumDeltaStateTicks,pReal)&
                                                                            /real(tickrate,pReal)
       write(6,'(a33,1x,f12.6)')  'avg CPU time/microsecs per call :',&
         real(debug_cumDeltaStateTicks,pReal)*1.0e6_pReal/real(tickrate,pReal)&
                                                                     /real(debug_cumDeltaStateCalls,pReal)
     endif
     write(6,*)
     write(6,'(a33,1x,i12)')      'total calls to dotTemperature   :',debug_cumDotTemperatureCalls
     if (debug_cumdotTemperatureCalls > 0_pInt) then
       write(6,'(a33,1x,f12.3)')  'total CPU time/s                :',real(debug_cumDotTemperatureTicks,pReal)&
                                                                            /real(tickrate,pReal)
       write(6,'(a33,1x,f12.6)')  'avg CPU time/microsecs per call :',&
         real(debug_cumDotTemperatureTicks,pReal)*1.0e6_pReal/real(tickrate,pReal)&
                                                                 /real(debug_cumDotTemperatureCalls,pReal)
     endif
   
     integral = 0_pInt
     write(6,*)
     write(6,*)
     write(6,*) 'distribution_StressLoop :    stress  stiffness'
     do i=1_pInt,nStress+1_pInt
       if (any(debug_StressLoopDistribution(i,:)     /= 0_pInt )) then
         integral = integral + i*(debug_StressLoopDistribution(i,1) + debug_StressLoopDistribution(i,2))
         exceed = ' '
         if (i > nStress) exceed = '+'                                                                    ! last entry gets "+"
         write(6,'(i25,a1,i10,1x,i10)') min(nStress,i),exceed,debug_StressLoopDistribution(i,1),&
                                                              debug_StressLoopDistribution(i,2)
       endif
     enddo
     write(6,'(a15,i10,2(1x,i10))') '          total',integral,sum(debug_StressLoopDistribution(:,1)), &
                                                               sum(debug_StressLoopDistribution(:,2))
     
     integral = 0_pInt
     write(6,*)
     write(6,*) 'distribution_CrystalliteStateLoop :'
     do i=1_pInt,nState+1_pInt
       if (any(debug_StateLoopDistribution(i,:) /= 0)) then
         integral = integral + i*(debug_StateLoopDistribution(i,1) + debug_StateLoopDistribution(i,2))
         exceed = ' '
         if (i > nState) exceed = '+'                                                                    ! last entry gets "+"
         write(6,'(i25,a1,i10,1x,i10)') min(nState,i),exceed,debug_StateLoopDistribution(i,1),&
                                                             debug_StateLoopDistribution(i,2)
       endif
     enddo
     write(6,'(a15,i10,2(1x,i10))') '          total',integral,sum(debug_StateLoopDistribution(:,1)), &
                                                               sum(debug_StateLoopDistribution(:,2))
    
     integral = 0_pInt
     write(6,*)
     write(6,*) 'distribution_CrystalliteCutbackLoop :'
     do i=1_pInt,nCryst+1_pInt
       if (debug_CrystalliteLoopDistribution(i) /= 0) then
         integral = integral + i*debug_CrystalliteLoopDistribution(i)
         exceed = ' '
         if (i > nCryst) exceed = '+'
         write(6,'(i25,a1,i10)') min(nCryst,i),exceed,debug_CrystalliteLoopDistribution(i)
       endif
     enddo
     write(6,'(a15,i10,1x,i10)') '          total',integral,sum(debug_CrystalliteLoopDistribution)
   endif
     
   if (iand(debug_level(debug_homogenization),debug_levelBasic) /= 0) then
     integral = 0_pInt
     write(6,*)
     write(6,*) 'distribution_MaterialpointStateLoop :'
     do i=1_pInt,nMPstate
       if (debug_MaterialpointStateLoopDistribution(i) /= 0) then
         integral = integral + i*debug_MaterialpointStateLoopDistribution(i)
         write(6,'(i25,1x,i10)') i,debug_MaterialpointStateLoopDistribution(i)
       endif
     enddo
     write(6,'(a15,i10,1x,i10)') '          total',integral,sum(debug_MaterialpointStateLoopDistribution) 
    
     integral = 0_pInt
     write(6,*)
     write(6,*) 'distribution_MaterialpointCutbackLoop :'
     do i=1_pInt,nHomog+1_pInt
       if (debug_MaterialpointLoopDistribution(i) /= 0) then
         integral = integral + i*debug_MaterialpointLoopDistribution(i)
         exceed = ' '
         if (i > nHomog) exceed = '+'
         write(6,'(i25,a1,i10)') min(nHomog,i),exceed,debug_MaterialpointLoopDistribution(i)
       endif
     enddo
     write(6,'(a15,i10,1x,i10)') '          total',integral,sum(debug_MaterialpointLoopDistribution)    
   endif
     
   if (iand(debug_level(debug_CPFEM),debug_levelBasic) /= 0) then
     write(6,*)
     write(6,*)
     write(6,*) 'Extreme values of returned stress and jacobian'
     write(6,*)
     write(6,'(a39)')                      '                      value     el   ip'
     write(6,'(a14,1x,e12.3,1x,i6,1x,i4)') 'stress   min :', debug_stressMin, debug_stressMinLocation
     write(6,'(a14,1x,e12.3,1x,i6,1x,i4)') '         max :', debug_stressMax, debug_stressMaxLocation
     write(6,'(a14,1x,e12.3,1x,i6,1x,i4)') 'jacobian min :', debug_jacobianMin, debug_jacobianMinLocation
     write(6,'(a14,1x,e12.3,1x,i6,1x,i4)') '         max :', debug_jacobianMax, debug_jacobianMaxLocation  
     write(6,*)
   endif
 !$OMP END CRITICAL (write2out)
 
end subroutine debug_info
 
end module debug
