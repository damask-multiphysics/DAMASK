! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the various spectral solvers
!> @details doing cutbacking, forwarding in case of restart, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_spectral_Driver
 use, intrinsic :: &
   iso_fortran_env                                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use prec, only: &
   pInt, &
   pReal
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   geometryFile, &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   appendToOutFile
 use IO, only: &
   IO_read, &
   IO_isBlank, &
   IO_open_file, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_error, &
   IO_lc, &
   IO_read_jobBinaryFile, &
   IO_write_jobBinaryFile, &
   IO_intOut, &
   IO_warning, &
   IO_timeStamp
 use debug, only: &
   debug_level, &
   debug_spectral, &
   debug_levelBasic
 use math                                                                                           ! need to include the whole module for FFTW
 use CPFEM, only: &
   CPFEM_initAll
 use FEsolving, only: &
   restartWrite, &
   restartInc
 use numerics, only: &
   maxCutBack, &
   rotation_tol, &
   mySpectralSolver, &
   regridMode
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results
 use DAMASK_spectral_Utilities, only: &
   grid, &
   geomSize, &
   tBoundaryCondition, &
   tSolutionState, &
   cutBack
 use DAMASK_spectral_SolverBasic
#ifdef PETSc
 use DAMASK_spectral_SolverBasicPETSC
 use DAMASK_spectral_SolverAL
 use DAMASK_spectral_SolverPolarisation
#endif

 implicit none
 type tLoadCase
   real(pReal), dimension (3,3) :: rotation               = math_I3                                 !< rotation of BC
   type(tBoundaryCondition) ::     P, &                                                             !< stress BC
                                   deformation                                                      !< deformation BC (Fdot or L)
   real(pReal) ::                  time                   = 0.0_pReal, &                            !< length of increment
                                   temperature            = 300.0_pReal                             !< isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               !< number of increments
                                   outputfrequency        = 1_pInt, &                               !< frequency of result writes
                                   restartfrequency       = 0_pInt, &                               !< frequency of restart writes
                                   logscale               = 0_pInt                                  !< linear/logarithmic time inc flag
   logical ::                      followFormerTrajectory = .true.                                  !< follow trajectory of former loadcase 
 end type tLoadCase

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: temp_valueVector = 0.0_pReal                                          !< temporarily from loadcase file when reading in tensors (initialize to 0.0)
 logical,     dimension(9) :: temp_maskVector  = .false.                                            !< temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: maxNchunks       = (1_pInt + 9_pInt)*3_pInt + &                       ! deformation, rotation, and stress
                                                 (1_pInt + 1_pInt)*5_pInt + &                       ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                  1_pInt, &                                         ! dropguessing
                              myUnit           = 234_pInt                                           !< file unit, DAMASK IO does not support newunit feature
 integer(pInt), dimension(1_pInt + maxNchunks*2_pInt) :: positions                                  ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_t    = 0_pInt, &                                                                               !< # of time indicators found in load case file 
   N_n    = 0_pInt, &                                                                               !< # of increment specifiers found in load case file
   N_def = 0_pInt                                                                                   !< # of rate of deformation specifiers found in load case file
 character(len=65536) :: &
   line

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal), dimension(3,3), parameter :: &
   ones  = 1.0_pReal, &
   zeros = 0.0_pReal 
 integer(pInt), parameter :: &
   subStepFactor = 2_pInt                                                                           !< for each substep, divide the last time increment by 2.0
 real(pReal) :: &
   time = 0.0_pReal, &                                                                              !< elapsed time
   time0 = 0.0_pReal, &                                                                             !< begin of interval
   timeinc = 1.0_pReal, &                                                                           !< current time interval
   timeIncOld = 0.0_pReal, &                                                                        !< previous time interval
   remainingLoadCaseTime = 0.0_pReal                                                                !< remaining time of current load case
 logical :: &
   guess                                                                                            !< guess along former trajectory
 integer(pInt) :: &
   i, j, k, l, &
   errorID, &
   cutBackLevel = 0_pInt, &                                                                         !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
   stepFraction = 0_pInt                                                                            !< fraction of current time interval
 integer(pInt) :: &
   currentLoadcase = 0_pInt, &                                                                      !< current load case
   inc, &                                                                                           !< current increment in current load case
   totalIncsCounter = 0_pInt, &                                                                     !< total No. of increments
   convergedCounter = 0_pInt, &                                                                     !< No. of converged increments
   notConvergedCounter = 0_pInt, &                                                                  !< No. of non-converged increments
   resUnit = 0_pInt, &                                                                              !< file unit for results writing
   statUnit = 0_pInt, &                                                                             !< file unit for statistics output
   lastRestartWritten = 0_pInt                                                                      !< total increment No. at which last restart information was written
 character(len=6)  :: loadcase_string
 character(len=1024)  :: incInfo                                                                    !< string parsed to solution with information about current load case
 type(tLoadCase), allocatable, dimension(:) :: loadCases                                            !< array of all load cases
 type(tSolutionState) solres
 external :: quit
!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
 call CPFEM_initAll(temperature = 300.0_pReal, element = 1_pInt, IP= 1_pInt)
 write(6,'(/,a)')   ' <<<+-  DAMASK_spectral_driver init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

!--------------------------------------------------------------------------------------------------
! reading basic information from load case file and allocate data structure containing load cases
 call IO_open_file(myUnit,trim(loadCaseFile))
 rewind(myUnit)
 do
   line = IO_read(myUnit)
   if (trim(line) == '#EOF#') exit
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   positions = IO_stringPos(line,maxNchunks)
   do i = 1_pInt, positions(1)                                                                      ! reading compulsory parameters for loadcase
     select case (IO_lc(IO_stringValue(line,positions,i)))
       case('l','velocitygrad','velgrad','velocitygradient','fdot','dotf','f')
         N_def = N_def + 1_pInt
       case('t','time','delta')
         N_t = N_t + 1_pInt
       case('n','incs','increments','steps','logincs','logincrements','logsteps')
         N_n = N_n + 1_pInt
     end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

 if ((N_def /= N_n) .or. (N_n /= N_t)) &                                                            ! sanity check
   call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                                    ! error message for incomplete loadcase
 allocate (loadCases(N_n))                                                                          ! array of load cases
 loadCases%P%myType='p'

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)
 do
   line = IO_read(myUnit)
   if (trim(line) == '#EOF#') exit
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   currentLoadCase = currentLoadCase + 1_pInt
   positions = IO_stringPos(line,maxNchunks)
   do i = 1_pInt, positions(1)
     select case (IO_lc(IO_stringValue(line,positions,i)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient','f')                      ! assign values for the deformation BC matrix
         temp_valueVector = 0.0_pReal
         if (IO_lc(IO_stringValue(line,positions,i)) == 'fdot'.or. &                                ! in case of Fdot, set type to fdot
             IO_lc(IO_stringValue(line,positions,i)) == 'dotf') then
           loadCases(currentLoadCase)%deformation%myType = 'fdot'
         else if (IO_lc(IO_stringValue(line,positions,i)) == 'f') then
           loadCases(currentLoadCase)%deformation%myType = 'f'
         else
           loadCases(currentLoadCase)%deformation%myType = 'l'
         endif
         do j = 1_pInt, 9_pInt
           temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'  ! true if not a *
         enddo
         do j = 1_pInt,9_pInt 
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)          ! read value where applicable
         enddo
         loadCases(currentLoadCase)%deformation%maskLogical = &                                     ! logical mask in 3x3 notation
               transpose(reshape(temp_maskVector,[ 3,3]))  
         loadCases(currentLoadCase)%deformation%maskFloat   = &                                     ! float (1.0/0.0) mask in 3x3 notation
               merge(ones,zeros,loadCases(currentLoadCase)%deformation%maskLogical)
         loadCases(currentLoadCase)%deformation%values = math_plain9to33(temp_valueVector)          ! values in 3x3 notation
       case('p','pk1','piolakirchhoff','stress', 's')
         temp_valueVector = 0.0_pReal
         do j = 1_pInt, 9_pInt
           temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'                           ! true if not an asterisk
         enddo
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)          ! read value where applicable
         enddo
         loadCases(currentLoadCase)%P%maskLogical = transpose(reshape(temp_maskVector,[ 3,3]))
         loadCases(currentLoadCase)%P%maskFloat   = merge(ones,zeros,&
                                                        loadCases(currentLoadCase)%P%maskLogical)
         loadCases(currentLoadCase)%P%values      = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         loadCases(currentLoadCase)%time = IO_floatValue(line,positions,i+1_pInt)
       case('temp','temperature')                                                                   ! starting temperature
         loadCases(currentLoadCase)%temperature = IO_floatValue(line,positions,i+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         loadCases(currentLoadCase)%incs = IO_intValue(line,positions,i+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         loadCases(currentLoadCase)%incs = IO_intValue(line,positions,i+1_pInt)
         loadCases(currentLoadCase)%logscale = 1_pInt
       case('freq','frequency','outputfreq')                                                        ! frequency of result writings
         loadCases(currentLoadCase)%outputfrequency = IO_intValue(line,positions,i+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         loadCases(currentLoadCase)%restartfrequency = &
               max(0_pInt,IO_intValue(line,positions,i+1_pInt))                
       case('guessreset','dropguessing')
         loadCases(currentLoadCase)%followFormerTrajectory = .false.                                ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of currentLoadCase given in euler angles
         temp_valueVector = 0.0_pReal
         l = 1_pInt                                                                                 ! assuming values given in degrees
         k = 1_pInt                                                                                 ! assuming keyword indicating degree/radians present
         select case (IO_lc(IO_stringValue(line,positions,i+1_pInt)))
           case('deg','degree')
           case('rad','radian')                                                                     ! don't convert from degree to radian           
             l = 0_pInt
           case default   
             k = 0_pInt           
         end select
         do j = 1_pInt, 3_pInt
           temp_valueVector(j) = IO_floatValue(line,positions,i+k+j)
         enddo
         if (l == 1_pInt) temp_valueVector(1:3) = temp_valueVector(1:3) * inRad                     ! convert to rad
         loadCases(currentLoadCase)%rotation = math_EulerToR(temp_valueVector(1:3))                 ! convert rad Eulers to rotation matrix
       case('rotation','rot')                                                                       ! assign values for the rotation of currentLoadCase matrix
         temp_valueVector = 0.0_pReal
         do j = 1_pInt, 9_pInt
           temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         loadCases(currentLoadCase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
 close(myUnit)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 loadCases(1)%followFormerTrajectory = .false.                                                      ! cannot guess along trajectory for first inc of first currentLoadCase
 errorID = 0_pInt
 checkLoadcases: do currentLoadCase = 1_pInt, size(loadCases)
   write (loadcase_string, '(i6)' ) currentLoadCase
   write(6,'(1x,a,i6)') 'load case: ', currentLoadCase
   if (.not. loadCases(currentLoadCase)%followFormerTrajectory) &
     write(6,'(2x,a)') 'drop guessing along trajectory'
   if (loadCases(currentLoadCase)%deformation%myType=='l') then
     do j = 1_pInt, 3_pInt
       if (any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .true.) .and. &
           any(loadCases(currentLoadCase)%deformation%maskLogical(j,1:3) .eqv. .false.)) &
                                                                  errorID = 832_pInt                ! each row should be either fully or not at all defined
     enddo
     write(6,'(2x,a)') 'velocity gradient:'
   else if (loadCases(currentLoadCase)%deformation%myType=='f') then
     write(6,'(2x,a)') 'deformation gradient at end of load case:'
   else
     write(6,'(2x,a)') 'deformation gradient rate:'
   endif
   write(6,'(3(3(3x,f12.7,1x)/))',advance='no') &
              merge(math_transpose33(loadCases(currentLoadCase)%deformation%values), &
              reshape(spread(huge(1.0_pReal),1,9),[ 3,3]), &                                        ! print *** (huge) for undefined
              transpose(loadCases(currentLoadCase)%deformation%maskLogical))
   if (any(loadCases(currentLoadCase)%P%maskLogical .eqv. &
            loadCases(currentLoadCase)%deformation%maskLogical)) errorID = 831_pInt                 ! exclusive or masking only
   if (any(loadCases(currentLoadCase)%P%maskLogical .and. &                                   
           transpose(loadCases(currentLoadCase)%P%maskLogical) .and. &
           reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
           errorID = 838_pInt                                                                       ! no rotation is allowed by stress BC
   write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'stress / GPa:',&
              1e-9_pReal*merge(math_transpose33(loadCases(currentLoadCase)%P%values),&
              reshape(spread(huge(1.0_pReal),1,9),[ 3,3]),&
              transpose(loadCases(currentLoadCase)%P%maskLogical))
  if (any(abs(math_mul33x33(loadCases(currentLoadCase)%rotation, &
              math_transpose33(loadCases(currentLoadCase)%rotation))-math_I3) >&
              reshape(spread(rotation_tol,1,9),[ 3,3]))&
              .or. abs(math_det33(loadCases(currentLoadCase)%rotation)) > &
              1.0_pReal + rotation_tol) errorID = 846_pInt                                          ! given rotation matrix contains strain
   if (any(loadCases(currentLoadCase)%rotation /= math_I3)) &
     write(6,'(2x,a,/,3(3(3x,f12.7,1x)/))',advance='no') 'rotation of loadframe:',&
              math_transpose33(loadCases(currentLoadCase)%rotation)
   write(6,'(2x,a,f12.6)') 'temperature:', loadCases(currentLoadCase)%temperature
   if (loadCases(currentLoadCase)%time < 0.0_pReal)          errorID = 834_pInt                     ! negative time increment
   write(6,'(2x,a,f12.6)') 'time:       ', loadCases(currentLoadCase)%time
   if (loadCases(currentLoadCase)%incs < 1_pInt)             errorID = 835_pInt                     ! non-positive incs count
   write(6,'(2x,a,i5)')    'increments: ', loadCases(currentLoadCase)%incs
   if (loadCases(currentLoadCase)%outputfrequency < 1_pInt)  errorID = 836_pInt                     ! non-positive result frequency
   write(6,'(2x,a,i5)')    'output  frequency:  ', &
              loadCases(currentLoadCase)%outputfrequency
   write(6,'(2x,a,i5,/)')    'restart frequency:  ', &
              loadCases(currentLoadCase)%restartfrequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)               ! exit with error message
 enddo checkLoadcases

!--------------------------------------------------------------------------------------------------
! doing initialization depending on selected solver 
 select case (myspectralsolver)
   case (DAMASK_spectral_SolverBasic_label)
     call basic_init(loadCases(1)%temperature)
#ifdef PETSc
   case (DAMASK_spectral_SolverBasicPETSc_label)
     call basicPETSc_init(loadCases(1)%temperature)
   case (DAMASK_spectral_SolverAL_label)
     if(iand(debug_level(debug_spectral),debug_levelBasic)/= 0) &
       call IO_warning(42_pInt, ext_msg='debug Divergence')
     call AL_init(loadCases(1)%temperature)
   case (DAMASK_spectral_SolverPolarisation_label)
     if(iand(debug_level(debug_spectral),debug_levelBasic)/= 0) &
       call IO_warning(42_pInt, ext_msg='debug Divergence')
     call Polarisation_init(loadCases(1)%temperature)
#endif
   case default
      call IO_error(error_ID = 891, ext_msg = trim(myspectralsolver))
 end select 
 
!--------------------------------------------------------------------------------------------------
! write header of output file
 if (appendToOutFile) then                                                                          ! after restart, append to existing results file
   open(newunit=resUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.spectralOut',form='UNFORMATTED', position='APPEND', status='OLD')
                             !   '.spectralOut',access='STREAM', position='APPEND', status='OLD')
   open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.sta',form='FORMATTED', position='APPEND', status='OLD')
 else                                                                                               ! open new files ...
   open(newunit=resUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.spectralOut',form='UNFORMATTED',status='REPLACE')
                             !   '.spectralOut',access='STREAM',status='REPLACE')
   write(resUnit) 'load',       trim(loadCaseFile)                                                  ! ... and write header
   write(resUnit) 'workingdir', trim(getSolverWorkingDirectoryName())
   write(resUnit) 'geometry',   trim(geometryFile)
   !write(resUnit) 'grid', grid
   !write(resUnit) 'size',  geomSize
   write(resUnit) 'resolution', grid
   write(resUnit) 'dimension',  geomSize
   write(resUnit) 'materialpoint_sizeResults', materialpoint_sizeResults
   write(resUnit) 'loadcases',  size(loadCases)
   write(resUnit) 'frequencies', loadCases%outputfrequency                                          ! one entry per currentLoadCase
   write(resUnit) 'times', loadCases%time                                                           ! one entry per currentLoadCase
   write(resUnit) 'logscales',  loadCases%logscale         
   write(resUnit) 'increments', loadCases%incs                                                      ! one entry per currentLoadCase
   write(resUnit) 'startingIncrement', restartInc - 1_pInt                                          ! start with writing out the previous inc
   write(resUnit) 'eoh'                                                                             ! end of header
   write(resUnit) materialpoint_results                                                             ! initial (non-deformed or read-in) results
   open(newunit=statUnit,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//&
                               '.sta',form='FORMATTED',status='REPLACE')
   write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                   ! statistics file
   if (iand(debug_level(debug_spectral),debug_levelBasic) /= 0) &
     write(6,'(/,a)') ' header of result file written out'
   flush(6)
 endif

!--------------------------------------------------------------------------------------------------
! loopping over loadcases
 loadCaseLooping: do currentLoadCase = 1_pInt, size(loadCases)
   time0 = time                                                                                     ! currentLoadCase start time                
   if (loadCases(currentLoadCase)%followFormerTrajectory) then
     guess = .true.
   else
     guess = .false.                                                                                ! change of load case, homogeneous guess for the first inc
   endif

!--------------------------------------------------------------------------------------------------
! loop oper incs defined in input file for current currentLoadCase
   incLooping: do inc = 1_pInt, loadCases(currentLoadCase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeIncOld = timeinc
     if (loadCases(currentLoadCase)%logscale == 0_pInt) then                                        ! linear scale
       timeinc = loadCases(currentLoadCase)%time/loadCases(currentLoadCase)%incs                    ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (currentLoadCase == 1_pInt) then                                                          ! 1st currentLoadCase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st currentLoadCase of logarithmic scale
           timeinc = loadCases(1)%time*(2.0_pReal**real(    1_pInt-loadCases(1)%incs ,pReal))       ! assume 1st inc is equal to 2nd 
         else                                                                                       ! not-1st inc of 1st currentLoadCase of logarithmic scale
           timeinc = loadCases(1)%time*(2.0_pReal**real(inc-1_pInt-loadCases(1)%incs ,pReal))
         endif
       else                                                                                         ! not-1st currentLoadCase of logarithmic scale
         timeinc = time0 * &
              ( (1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real(          inc,pReal)/&
                                                    real(loadCases(currentLoadCase)%incs ,pReal))&
               -(1.0_pReal + loadCases(currentLoadCase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                     real(loadCases(currentLoadCase)%incs ,pReal)))
       endif
     endif
     timeinc = timeinc / 2.0_pReal**real(cutBackLevel,pReal)                                        ! depending on cut back level, decrease time step

     forwarding: if(totalIncsCounter >= restartInc) then
       stepFraction = 0_pInt

!--------------------------------------------------------------------------------------------------
! loop over sub incs 
       subIncLooping: do while (stepFraction/subStepFactor**cutBackLevel <1_pInt)
         time = time + timeinc                                                                      ! forward time
         stepFraction = stepFraction + 1_pInt 
         remainingLoadCaseTime = time0 - time + loadCases(currentLoadCase)%time + timeInc
!--------------------------------------------------------------------------------------------------
! report begin of new increment
         write(6,'(/,a)') ' ###########################################################################'
         write(6,'(1x,a,es12.5'//&
                 ',a,'//IO_intOut(inc)//',a,'//IO_intOut(loadCases(currentLoadCase)%incs)//&
                 ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(subStepFactor**cutBackLevel)//&
                 ',a,'//IO_intOut(currentLoadCase)//',a,'//IO_intOut(size(loadCases))//')') &
                 'Time', time, &
                 's: Increment ', inc, '/', loadCases(currentLoadCase)%incs,&
                 '-', stepFraction, '/', subStepFactor**cutBackLevel,&
                 ' of load case ', currentLoadCase,'/',size(loadCases)
         flush(6)
         write(incInfo,'(a,'//IO_intOut(totalIncsCounter)//',a,'//IO_intOut(sum(loadCases%incs))//&
               ',a,'//IO_intOut(stepFraction)//',a,'//IO_intOut(subStepFactor**cutBackLevel)//')') &
               'Increment ',totalIncsCounter,'/',sum(loadCases%incs),&
               '-',stepFraction, '/', subStepFactor**cutBackLevel
         select case(myspectralsolver)
         
!--------------------------------------------------------------------------------------------------
! calculate solution 
           case (DAMASK_spectral_SolverBasic_label)
             solres = basic_solution (&
                 incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                 P_BC               = loadCases(currentLoadCase)%P, &
                 F_BC               = loadCases(currentLoadCase)%deformation, &
                 temperature_bc     = loadCases(currentLoadCase)%temperature, &
                 rotation_BC        = loadCases(currentLoadCase)%rotation)
#ifdef PETSc
           case (DAMASK_spectral_SolverBasicPETSC_label)
             solres = BasicPETSC_solution (&
                 incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                 P_BC               = loadCases(currentLoadCase)%P, &
                 F_BC               = loadCases(currentLoadCase)%deformation, &
                 temperature_bc     = loadCases(currentLoadCase)%temperature, &
                 rotation_BC        = loadCases(currentLoadCase)%rotation)
           case (DAMASK_spectral_SolverAL_label)
             solres = AL_solution (&
                 incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                 P_BC               = loadCases(currentLoadCase)%P, &
                 F_BC               = loadCases(currentLoadCase)%deformation, &
                 temperature_bc     = loadCases(currentLoadCase)%temperature, &
                 rotation_BC        = loadCases(currentLoadCase)%rotation)
           case (DAMASK_spectral_SolverPolarisation_label)
             solres = Polarisation_solution (&
                 incInfo,guess,timeinc,timeIncOld,remainingLoadCaseTime, &
                 P_BC               = loadCases(currentLoadCase)%P, &
                 F_BC               = loadCases(currentLoadCase)%deformation, &
                 temperature_bc     = loadCases(currentLoadCase)%temperature, &
                 rotation_BC        = loadCases(currentLoadCase)%rotation)
#endif
         end select 

!--------------------------------------------------------------------------------------------------
! check solution 
         cutBack = .False.                                                                   
         if(solres%termIll .or. .not. solres%converged) then                                        ! no solution found
           if (cutBackLevel < maxCutBack) then                                                      ! do cut back
             write(6,'(/,a)') ' cut back detected'
             cutBack = .True.
             stepFraction = (stepFraction - 1_pInt) * subStepFactor                                 ! adjust to new denominator
             cutBackLevel = cutBackLevel + 1_pInt
             time    = time - timeinc                                                               ! rewind time
             timeIncOld = timeinc
             timeinc = timeinc/2.0_pReal
           elseif (solres%termIll) then                                                             ! material point model cannot find a solution
             if(regridMode > 0_pInt) call quit(-1_pInt*(lastRestartWritten+1_pInt))                 ! regrid requested (mode 1 or 2)
             call IO_error(850_pInt)                                                                ! no regrid (give up)
           else
             if(regridMode == 2_pInt) call quit(-1_pInt*(lastRestartWritten+1_pInt))                ! regrid also if BVP solver do not converge
             guess = .true.                                                                         ! continue from non-converged solution and start guessing after accepted (sub)inc
           endif
         else
           guess = .true.                                                                           ! start guessing after first converged (sub)inc
         endif
         if (.not. cutBack) &
           write(statUnit,*) totalIncsCounter, time, cutBackLevel, &
                             solres%converged, solres%iterationsNeeded                              ! write statistics about accepted solution
       enddo subIncLooping
       cutBackLevel = max(0_pInt, cutBackLevel - 1_pInt)                                            ! try half number of subincs next inc
       if(solres%converged) then                                                                    ! report converged inc
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &
                                     ' increment ', totalIncsCounter, ' converged'
       else
         write(6,'(/,a,'//IO_intOut(totalIncsCounter)//',a)') &                                     ! report non-converged inc
                                     ' increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif; flush(6)

       if (mod(inc,loadCases(currentLoadCase)%outputFrequency) == 0_pInt) then                      ! at output frequency
         write(6,'(1/,a)') ' ... writing results to file ......................................'
         write(resUnit)  materialpoint_results                                                      ! write result to file
       endif
       if( loadCases(currentLoadCase)%restartFrequency > 0_pInt .and. &                             ! at frequency of writing restart information set restart parameter for FEsolving 
                      mod(inc,loadCases(currentLoadCase)%restartFrequency) == 0_pInt) then          ! ToDo first call to CPFEM_general will write? 
         restartWrite = .true.
         lastRestartWritten = inc
       endif 
     else forwarding
       time = time + timeinc
       guess = .true.
     endif forwarding

    enddo incLooping
 enddo loadCaseLooping
 
 select case (myspectralsolver)
   case (DAMASK_spectral_SolverBasic_label)
     call basic_destroy()
#ifdef PETSc
   case (DAMASK_spectral_SolverBasicPETSC_label)
     call BasicPETSC_destroy()
   case (DAMASK_spectral_SolverAL_label)
     call AL_destroy()
   case (DAMASK_spectral_SolverPolarisation_label)
     call Polarisation_destroy()
#endif 
 end select
 
!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
 write(6,'(/,a)') ' ##################################################################'
 write(6,'(1x,i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                   notConvergedCounter + convergedCounter, ' (', &
                                   real(convergedCounter, pReal)/&
                                   real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                   ' %) increments converged!'
 close(resUnit)
 if (notConvergedCounter > 0_pInt) call quit(3_pInt)                                                ! error if some are not converged
 call quit(0_pInt)                                                                                  ! no complains ;)

end program DAMASK_spectral_Driver


!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief quit subroutine to mimic behavior of FEM solvers
!> @details exits the Spectral solver and reports time and duration. Exit code 0 signals
!> everything went fine. Exit code 1 signals an error, message according to IO_error. Exit code 
!> 2 signals request for regridding, increment of last saved restart information is written to
!> stderr. Exit code 3 signals no severe problems, but some increments did not converge
!--------------------------------------------------------------------------------------------------
subroutine quit(stop_id)
 use prec, only: &
   pInt
   
 implicit none
 integer(pInt), intent(in) :: stop_id
 integer, dimension(8) :: dateAndTime                                                               ! type default integer

 call date_and_time(values = dateAndTime)
 write(6,'(/,a)') 'DAMASK terminated on:'
 write(6,'(a,2(i2.2,a),i4.4)') 'Date:               ',dateAndTime(3),'/',&
                                                      dateAndTime(2),'/',&
                                                      dateAndTime(1) 
 write(6,'(a,2(i2.2,a),i2.2)') 'Time:               ',dateAndTime(5),':',&
                                                      dateAndTime(6),':',&
                                                      dateAndTime(7)  
 if (stop_id == 0_pInt) stop 0                                                                      ! normal termination
 if (stop_id <  0_pInt) then                                                                        ! trigger regridding
   write(0,'(a,i6)') 'restart at ', stop_id*(-1_pInt)
   stop 2
 endif
 if (stop_id == 3_pInt) stop 3                                                                      ! not all incs converged
 stop 1                                                                                             ! error (message from IO_error)

end subroutine quit
