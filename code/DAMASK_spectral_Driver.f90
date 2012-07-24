! Copyright 2012 Max-Planck-Institut für Eisenforschung GmbH
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
!##################################################################################################
!* $Id$
!##################################################################################################
! Material subroutine for BVP solution using spectral method
!
! Run 'DAMASK_spectral.exe --help' to get usage hints
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts,
!            D.D. Tjahjanto,
!            C. Kords,
!            M. Diehl,
!            R. Lebensohn
!
! MPI fuer Eisenforschung, Duesseldorf


program DAMASK_spectral_Driver
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 
 use DAMASK_interface, only: &
   DAMASK_interface_init, &
   loadCaseFile, &
   geometryFile, &
   getSolverWorkingDirectoryName, &
   getSolverJobName, &
   appendToOutFile
   
 use prec, only: &
   pInt, &
   pReal, &
   DAMASK_NaN
   
 use IO, only: &
   IO_isBlank, &
   IO_open_file, &
   IO_stringPos, &
   IO_stringValue, &
   IO_floatValue, &
   IO_intValue, &
   IO_error, &
   IO_lc, &
   IO_read_jobBinaryFile, &
   IO_write_jobBinaryFile
      
 use math
 
 use mesh,  only : &
   mesh_spectral_getResolution, &
   mesh_spectral_getDimension, &
   mesh_spectral_getHomogenization
 
 use CPFEM, only: &
   CPFEM_initAll
   
 use FEsolving, only: &
   restartWrite, &
   restartInc
   
 use numerics, only: &
   rotation_tol
   
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results
   
 use DAMASK_spectral_SolverAL
 use DAMASK_spectral_SolverBasic
 use DAMASK_spectral_Utilities
 
 implicit none
 
 type loadcase
   real(pReal), dimension (3,3) :: deformation            = 0.0_pReal, &                            ! applied velocity gradient or time derivative of deformation gradient
                                   stress                 = 0.0_pReal, &                            ! stress BC (if applicable)
                                   rotation               = math_I3                                 ! rotation of BC (if applicable)
   real(pReal) ::                  time                   = 0.0_pReal, &                            ! length of increment
                                   temperature            = 300.0_pReal                             ! isothermal starting conditions
   integer(pInt) ::                incs                   = 0_pInt, &                               ! number of increments
                                   outputfrequency        = 1_pInt, &                               ! frequency of result writes
                                   restartfrequency       = 0_pInt, &                               ! frequency of restart writes
                                   logscale               = 0_pInt                                  ! linear/logaritmic time inc flag
   logical ::                      followFormerTrajectory = .true., &                               ! follow trajectory of former loadcase
                                   velGradApplied         = .false.                                 ! decide wether velocity gradient or fdot is given 
   logical, dimension(3,3) ::      maskDeformation        = .false., &                              ! mask of deformation boundary conditions
                                   maskStress             = .false.                                 ! mask of stress boundary conditions
   logical, dimension(9) ::        maskStressVector       = .false.                                 ! linear mask of boundary conditions    
 end type

!--------------------------------------------------------------------------------------------------
! variables related to information from load case and geom file
 real(pReal), dimension(9) :: & 
   temp_valueVector                                                                                 !> temporarily from loadcase file when reading in tensors
 logical,     dimension(9) :: &
   temp_maskVector                                                                                  !> temporarily from loadcase file when reading in tensors
 integer(pInt), parameter  :: maxNchunksLoadcase = (1_pInt + 9_pInt)*3_pInt +&                      ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt +&                      ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                    1_pInt, &                                       ! dropguessing
                              maxNchunksGeom     = 7_pInt, &                                        ! 4 identifiers, 3 values
                              myUnit             = 234_pInt
 integer(pInt), dimension(1_pInt + maxNchunksLoadcase*2_pInt) :: positions                          ! this is longer than needed for geometry parsing
 
 integer(pInt) :: &
   N_l    = 0_pInt, &
   N_t    = 0_pInt, &
   N_n    = 0_pInt, &
   N_Fdot = 0_pInt                                                                                 ! number of Fourier points

 character(len=1024) :: &
   line

 type(loadcase), allocatable, dimension(:) ::  bc
 type(solutionState) solres

           
!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval 
 real(pReal) :: guessmode             
 real(pReal),    dimension(3,3) ::          temp33_Real
 integer(pInt) :: i, j, k, l, errorID
 integer(pInt) :: currentLoadcase = 0_pInt, inc, &
                  totalIncsCounter = 0_pInt,&
                  notConvergedCounter = 0_pInt, convergedCounter = 0_pInt
 character(len=6)  :: loadcase_string

 call DAMASK_interface_init
 
 write(6,'(a)') ''
 write(6,'(a)') ' <<<+-  DAMASK_spectral init  -+>>>'
 write(6,'(a)') ' $Id$'
#include "compilation_info.f90"
 write(6,'(a)') ' Working Directory:    ',trim(getSolverWorkingDirectoryName())
 write(6,'(a)') ' Solver Job Name:      ',trim(getSolverJobName())
 write(6,'(a)') ''

!--------------------------------------------------------------------------------------------------
! reading the load case file and allocate data structure containing load cases
 call IO_open_file(myUnit,trim(loadCaseFile))
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt, maxNchunksLoadcase, 1_pInt                                                        ! reading compulsory parameters for loadcase
       select case (IO_lc(IO_stringValue(line,positions,i)))
            case('l','velocitygrad','velgrad','velocitygradient')
                 N_l = N_l + 1_pInt
            case('fdot','dotf')
                 N_Fdot = N_Fdot + 1_pInt
            case('t','time','delta')
                 N_t = N_t + 1_pInt
            case('n','incs','increments','steps','logincs','logincrements','logsteps')
                 N_n = N_n + 1_pInt
        end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

100 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                     ! sanity check
       call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                               ! error message for incomplete loadcase
 allocate (bc(N_n))

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   currentLoadcase = currentLoadcase + 1_pInt
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,positions,i)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient')                          ! assign values for the deformation BC matrix
         bc(currentLoadcase)%velGradApplied = &
                     (IO_lc(IO_stringValue(line,positions,i)) == 'l'.or. &                          ! in case of given L, set flag to true
                      IO_lc(IO_stringValue(line,positions,i)) == 'velocitygrad'.or.&
                      IO_lc(IO_stringValue(line,positions,i)) == 'velgrad'.or.&
                      IO_lc(IO_stringValue(line,positions,i)) == 'velocitygradient')
         temp_valueVector = 0.0_pReal
         temp_maskVector = .false.
         forall (j = 1_pInt:9_pInt) temp_maskVector(j) = IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (temp_maskVector(j)) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         enddo
         bc(currentLoadcase)%maskDeformation = transpose(reshape(temp_maskVector,[ 3,3]))
         bc(currentLoadcase)%deformation = math_plain9to33(temp_valueVector)
       case('p','pk1','piolakirchhoff','stress')
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) bc(currentLoadcase)%maskStressVector(j) =&
                                                          IO_stringValue(line,positions,i+j) /= '*'
         do j = 1_pInt,9_pInt
           if (bc(currentLoadcase)%maskStressVector(j)) temp_valueVector(j) =&
                                                          IO_floatValue(line,positions,i+j)         ! assign values for the bc(currentLoadcase)%stress matrix
         enddo
         bc(currentLoadcase)%maskStress = transpose(reshape(bc(currentLoadcase)%maskStressVector,[ 3,3]))
         bc(currentLoadcase)%stress = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         bc(currentLoadcase)%time = IO_floatValue(line,positions,i+1_pInt)
       case('temp','temperature')                                                                   ! starting temperature
         bc(currentLoadcase)%temperature = IO_floatValue(line,positions,i+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         bc(currentLoadcase)%incs = IO_intValue(line,positions,i+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         bc(currentLoadcase)%incs = IO_intValue(line,positions,i+1_pInt)
         bc(currentLoadcase)%logscale = 1_pInt
       case('f','freq','frequency','outputfreq')                                                    ! frequency of result writings
         bc(currentLoadcase)%outputfrequency = IO_intValue(line,positions,i+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         bc(currentLoadcase)%restartfrequency = max(0_pInt,IO_intValue(line,positions,i+1_pInt))                
       case('guessreset','dropguessing')
         bc(currentLoadcase)%followFormerTrajectory = .false.                                              ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of currentLoadcase given in euler angles
         l = 0_pInt                                                                                 ! assuming values given in radians
         k = 1_pInt                                                                                 ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,positions,i+1_pInt)))
           case('deg','degree')
             l = 1_pInt                                                                             ! for conversion from degree to radian           
           case('rad','radian') 
           case default               
             k = 0_pInt                                                                             ! immediately reading in angles, assuming radians
         end select
         forall(j = 1_pInt:3_pInt)  temp33_Real(j,1) = &
                                        IO_floatValue(line,positions,i+k+j) * real(l,pReal) * inRad
         bc(currentLoadcase)%rotation = math_EulerToR(temp33_Real(:,1))
       case('rotation','rot')                                                                       ! assign values for the rotation of currentLoadcase matrix
         temp_valueVector = 0.0_pReal
         forall (j = 1_pInt:9_pInt) temp_valueVector(j) = IO_floatValue(line,positions,i+j)
         bc(currentLoadcase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
101 close(myUnit)
print*, 'my Unit closed'
!-------------------------------------------------------------------------------------------------- ToDo: if temperature at CPFEM is treated properly, move this up immediately after interface init
! initialization of all related DAMASK modules (e.g. mesh.f90 reads in geometry)
 call CPFEM_initAll(bc(1)%temperature,1_pInt,1_pInt)
 
!--------------------------------------------------------------------------------------------------
! output of geometry information
 write(6,'(a)')          ''
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'DAMASK spectral:'
 write(6,'(a)')          'The spectral method boundary value problem solver for'
 write(6,'(a)')          'the Duesseldorf Advanced Material Simulation Kit'
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'geometry file:        ',trim(geometryFile)
 write(6,'(a)')          '============================================================='
 write(6,'(a,3(i12  ))') 'resolution a b c:',      mesh_spectral_getResolution()
 write(6,'(a,3(f12.5))') 'dimension  x y z:',      mesh_spectral_getDimension()
 write(6,'(a,i5)')       'homogenization:       ', mesh_spectral_getHomogenization()
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'currentLoadcase file:        ',trim(loadCaseFile)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 bc(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first currentLoadcase
 errorID = 0_pInt
 checkLoadcases: do currentLoadcase = 1_pInt, size(bc)
   write (loadcase_string, '(i6)' ) currentLoadcase

   write(6,'(a)') '============================================================='
   write(6,'(a,i6)') 'currentLoadcase:            ', currentLoadcase

   if (.not. bc(currentLoadcase)%followFormerTrajectory) write(6,'(a)') 'drop guessing along trajectory'
   if (bc(currentLoadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(currentLoadcase)%maskDeformation(j,1:3) .eqv. .true.) .and. &
           any(bc(currentLoadcase)%maskDeformation(j,1:3) .eqv. .false.)) errorID = 832_pInt               ! each row should be either fully or not at all defined
     enddo
     write(6,'(a)')'velocity gradient:'
   else
     write(6,'(a)')'deformation gradient rate:'
   endif
   write (6,'(3(3(f12.7,1x)/))',advance='no') merge(math_transpose33(bc(currentLoadcase)%deformation),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(currentLoadcase)%maskDeformation))
   write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' stress / GPa:',&
        1e-9_pReal*merge(math_transpose33(bc(currentLoadcase)%stress),&
                         reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(currentLoadcase)%maskStress))
   if (any(bc(currentLoadcase)%rotation /= math_I3)) &
     write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' rotation of loadframe:',&
                                                          math_transpose33(bc(currentLoadcase)%rotation)
   write(6,'(a,f12.6)') 'temperature:', bc(currentLoadcase)%temperature
   write(6,'(a,f12.6)') 'time:       ', bc(currentLoadcase)%time
   write(6,'(a,i5)')    'increments: ', bc(currentLoadcase)%incs
   write(6,'(a,i5)')    'output  frequency:  ', bc(currentLoadcase)%outputfrequency
   write(6,'(a,i5)')    'restart frequency:  ', bc(currentLoadcase)%restartfrequency

   if (any(bc(currentLoadcase)%maskStress .eqv. bc(currentLoadcase)%maskDeformation)) errorID = 831_pInt          ! exclusive or masking only
   if (any(bc(currentLoadcase)%maskStress .and. transpose(bc(currentLoadcase)%maskStress) .and. &
     reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
                                               errorID = 838_pInt                                   ! no rotation is allowed by stress BC
   if (any(abs(math_mul33x33(bc(currentLoadcase)%rotation,math_transpose33(bc(currentLoadcase)%rotation))&
                                      -math_I3) > reshape(spread(rotation_tol,1,9),[ 3,3]))&
                    .or. abs(math_det33(bc(currentLoadcase)%rotation)) > 1.0_pReal + rotation_tol)&
                                               errorID = 846_pInt                                   ! given rotation matrix contains strain
   if (bc(currentLoadcase)%time < 0.0_pReal)          errorID = 834_pInt                                   ! negative time increment
   if (bc(currentLoadcase)%incs < 1_pInt)             errorID = 835_pInt                                   ! non-positive incs count
   if (bc(currentLoadcase)%outputfrequency < 1_pInt)  errorID = 836_pInt                                   ! non-positive result frequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo checkLoadcases

!--------------------------------------------------------------------------------------------------
! write header of output file
 if (appendToOutFile) then
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED', position='APPEND', status='OLD')
 else
   open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())//'.spectralOut',&
                                   form='UNFORMATTED',status='REPLACE')
   write(538) 'load',       trim(loadCaseFile)
   write(538) 'workingdir', trim(getSolverWorkingDirectoryName())
   write(538) 'geometry',   trim(geometryFile)
   write(538) 'resolution', mesh_spectral_getResolution()
   write(538) 'dimension',  mesh_spectral_getDimension()
   write(538) 'materialpoint_sizeResults', materialpoint_sizeResults
   write(538) 'loadcases',        size(bc)
   write(538) 'frequencies', bc%outputfrequency                                      ! one entry per currentLoadcase
   write(538) 'times', bc%time                                                       ! one entry per currentLoadcase
   write(538) 'logscales',  bc%logscale         
   write(538) 'increments', bc%incs                                                  ! one entry per currentLoadcase
   write(538) 'startingIncrement', restartInc - 1_pInt                                              ! start with writing out the previous inc
   write(538) 'eoh'                                                                                 ! end of header
   write(538) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)                ! initial (non-deformed or read-in) results
   if (debugGeneral) write(6,'(a)') 'Header of result file written out'
 endif

 call Basic_init()

!##################################################################################################
! Loop over loadcases defined in the currentLoadcase file
!##################################################################################################
 loadCaseLooping: do currentLoadcase = 1_pInt, size(bc)
   time0 = time                                                                                     ! currentLoadcase start time                
   if (bc(currentLoadcase)%followFormerTrajectory) then
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                                          ! change of load case, homogeneous guess for the first inc
   endif

!##################################################################################################
! loop oper incs defined in input file for current currentLoadcase
!##################################################################################################
   incLooping: do inc = 1_pInt,  bc(currentLoadcase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeinc_old = timeinc
     if (bc(currentLoadcase)%logscale == 0_pInt) then                                                      ! linear scale
       timeinc = bc(currentLoadcase)%time/bc(currentLoadcase)%incs                                                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (currentLoadcase == 1_pInt) then                                                                 ! 1st currentLoadcase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st currentLoadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(    1_pInt-bc(1)%incs ,pReal))                     ! assume 1st inc is equal to 2nd 
         else                                                                                       ! not-1st inc of 1st currentLoadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(inc-1_pInt-bc(1)%incs ,pReal))
         endif
       else                                                                                         ! not-1st currentLoadcase of logarithmic scale
           timeinc = time0 *( (1.0_pReal + bc(currentLoadcase)%time/time0 )**(real(          inc,pReal)/&
                                                                  real(bc(currentLoadcase)%incs ,pReal))&
                             -(1.0_pReal + bc(currentLoadcase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                                   real(bc(currentLoadcase)%incs ,pReal)) )
       endif
     endif
     time = time + timeinc

     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       write(6,'(a)') '##################################################################'
       write(6,'(A,I5.5,A,es12.5)') 'Increment ', totalIncsCounter, ' Time ',time
       
       solres =basic_solution (&
               guessmode,timeinc,timeinc_old, &
                P_BC              = bc(currentLoadcase)%stress, &
                F_BC              = bc(currentLoadcase)%deformation, &
               ! temperature_bc       = bc(currentLoadcase)%temperature, &
                mask_stressVector = bc(currentLoadcase)%maskStressVector, &
                velgrad           = bc(currentLoadcase)%velGradApplied, &
                rotation_BC       = bc(currentLoadcase)%rotation)
 
       write(6,'(a)') ''
       write(6,'(a)') '=================================================================='
       if(solres%converged) then
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' converged'
       else
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       endif

       if (mod(inc,bc(currentLoadcase)%outputFrequency) == 0_pInt) then                                    ! at output frequency
         write(6,'(a)') ''
         write(6,'(a)') '... writing results to file ......................................'
         write(538)  materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)       ! write result to file
       endif
       
     endif ! end calculation/forwarding
     guessmode = 1.0_pReal                                                                          ! keep guessing along former trajectory during same currentLoadcase

    enddo incLooping
 enddo loadCaseLooping
 write(6,'(a)') ''
 write(6,'(a)') '##################################################################'
 write(6,'(i6.6,a,i6.6,a,f5.1,a)') convergedCounter, ' out of ', &
                                   notConvergedCounter + convergedCounter, ' (', &
                                   real(convergedCounter, pReal)/&
                                   real(notConvergedCounter + convergedCounter,pReal)*100.0_pReal, &
                                   ' %) increments converged!'
 close(538)
 if (notConvergedCounter > 0_pInt) call quit(3_pInt)
 call quit(0_pInt)

end program DAMASK_spectral_Driver




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
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
 if (stop_id == 3_pInt) stop 3                                                                      ! not all steps converged
 stop 1                                                                                             ! error (message from IO_error)
end subroutine
