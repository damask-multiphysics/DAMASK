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

#include "spectral_quit.f90"

program DAMASK_spectral
 use prec, only: &
   pInt, &
   pReal
   
 use IO, only: &
   IO_error,&
   IO_write_jobBinaryFile
      
 use math
 
 use FEsolving, only: &
   restartWrite, &
   restartInc
    
 use homogenization, only: &
   materialpoint_sizeResults, &
   materialpoint_results
   
 use DAMASK_spectralSovler
 
 implicit none

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 integer(pInt) :: i, j, k, l, m, n, p, errorID


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

100 N_Loadcases = N_n
 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                     ! sanity check
   call IO_error(error_ID=837_pInt,ext_msg = trim(loadCaseFile))                               ! error message for incomplete loadcase
 allocate (bc(N_Loadcases))

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)

 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   loadcase = loadcase + 1_pInt
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do j = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,positions,j)))
       case('fdot','dotf','l','velocitygrad','velgrad','velocitygradient')                          ! assign values for the deformation BC matrix
         bc(loadcase)%velGradApplied = &
                     (IO_lc(IO_stringValue(line,positions,j)) == 'l'.or. &                          ! in case of given L, set flag to true
                      IO_lc(IO_stringValue(line,positions,j)) == 'velocitygrad'.or.&
                      IO_lc(IO_stringValue(line,positions,j)) == 'velgrad'.or.&
                      IO_lc(IO_stringValue(line,positions,j)) == 'velocitygradient')
         temp_valueVector = 0.0_pReal
         temp_maskVector = .false.
         forall (k = 1_pInt:9_pInt) temp_maskVector(k) = IO_stringValue(line,positions,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (temp_maskVector(k)) temp_valueVector(k) = IO_floatValue(line,positions,j+k)
         enddo
         bc(loadcase)%maskDeformation = transpose(reshape(temp_maskVector,[ 3,3]))
         bc(loadcase)%deformation = math_plain9to33(temp_valueVector)
       case('p','pk1','piolakirchhoff','stress')
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) bc(loadcase)%maskStressVector(k) =&
                                                          IO_stringValue(line,positions,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (bc(loadcase)%maskStressVector(k)) temp_valueVector(k) =&
                                                          IO_floatValue(line,positions,j+k)         ! assign values for the bc(loadcase)%stress matrix
         enddo
         bc(loadcase)%maskStress = transpose(reshape(bc(loadcase)%maskStressVector,[ 3,3]))
         bc(loadcase)%stress = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                                     ! increment time
         bc(loadcase)%time = IO_floatValue(line,positions,j+1_pInt)
       case('temp','temperature')                                                                   ! starting temperature
         bc(loadcase)%temperature = IO_floatValue(line,positions,j+1_pInt)
       case('n','incs','increments','steps')                                                        ! number of increments
         bc(loadcase)%incs = IO_intValue(line,positions,j+1_pInt)
       case('logincs','logincrements','logsteps')                                                   ! number of increments (switch to log time scaling)
         bc(loadcase)%incs = IO_intValue(line,positions,j+1_pInt)
         bc(loadcase)%logscale = 1_pInt
       case('f','freq','frequency','outputfreq')                                                    ! frequency of result writings
         bc(loadcase)%outputfrequency = IO_intValue(line,positions,j+1_pInt)                
       case('r','restart','restartwrite')                                                           ! frequency of writing restart information
         bc(loadcase)%restartfrequency = max(0_pInt,IO_intValue(line,positions,j+1_pInt))                
       case('guessreset','dropguessing')
         bc(loadcase)%followFormerTrajectory = .false.                                              ! do not continue to predict deformation along former trajectory
       case('euler')                                                                                ! rotation of loadcase given in euler angles
         p = 0_pInt                                                                                 ! assuming values given in radians
         l = 1_pInt                                                                                 ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,positions,j+1_pInt)))
           case('deg','degree')
             p = 1_pInt                                                                             ! for conversion from degree to radian           
           case('rad','radian') 
           case default               
             l = 0_pInt                                                                             ! immediately reading in angles, assuming radians
         end select
         forall(k = 1_pInt:3_pInt)  temp33_Real(k,1) = &
                                        IO_floatValue(line,positions,j+l+k) * real(p,pReal) * inRad
         bc(loadcase)%rotation = math_EulerToR(temp33_Real(:,1))
       case('rotation','rot')                                                                       ! assign values for the rotation of loadcase matrix
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) temp_valueVector(k) = IO_floatValue(line,positions,j+k)
         bc(loadcase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo
101 close(myUnit)

!-------------------------------------------------------------------------------------------------- ToDo: if temperature at CPFEM is treated properly, move this up immediately after interface init
! initialization of all related DAMASK modules (e.g. mesh.f90 reads in geometry)
 call CPFEM_initAll(bc(1)%temperature,1_pInt,1_pInt)
 
!--------------------------------------------------------------------------------------------------
! get resolution, dimension, homogenization and variables derived from resolution
 res     = mesh_spectral_getResolution()
 geomdim = mesh_spectral_getDimension()
 homog   = mesh_spectral_getHomogenization()
 res1_red = res(1)/2_pInt + 1_pInt                                                                  ! size of complex array in first dimension (c2r, r2c)
 Npoints = res(1)*res(2)*res(3)
 wgt = 1.0_pReal/real(Npoints, pReal)

!--------------------------------------------------------------------------------------------------
! output of geometry
 write(6,'(a)')          ''
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'DAMASK spectral:'
 write(6,'(a)')          'The spectral method boundary value problem solver for'
 write(6,'(a)')          'the Duesseldorf Advanced Material Simulation Kit'
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'geometry file:        ',trim(geometryFile)
 write(6,'(a)')          '============================================================='
 write(6,'(a,3(i12  ))') 'resolution a b c:', res
 write(6,'(a,3(f12.5))') 'dimension  x y z:', geomdim
 write(6,'(a,i5)')       'homogenization:       ',homog
 write(6,'(a)')          '#############################################################'
 write(6,'(a)')          'loadcase file:        ',trim(loadCaseFile)

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 bc(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first loadcase
 errorID = 0_pInt
 do loadcase = 1_pInt, N_Loadcases
   write (loadcase_string, '(i6)' ) loadcase

   write(6,'(a)') '============================================================='
   write(6,'(a,i6)') 'loadcase:            ', loadcase

   if (.not. bc(loadcase)%followFormerTrajectory) write(6,'(a)') 'drop guessing along trajectory'
   if (bc(loadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .true.) .and. &
           any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .false.)) errorID = 832_pInt               ! each row should be either fully or not at all defined
     enddo
     write(6,'(a)')'velocity gradient:'
   else
     write(6,'(a)')'deformation gradient rate:'
   endif
   write (6,'(3(3(f12.7,1x)/))',advance='no') merge(math_transpose33(bc(loadcase)%deformation),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskDeformation))
   write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' stress / GPa:',&
        1e-9_pReal*merge(math_transpose33(bc(loadcase)%stress),&
                         reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskStress))
   if (any(bc(loadcase)%rotation /= math_I3)) &
     write (6,'(a,/,3(3(f12.7,1x)/))',advance='no') ' rotation of loadframe:',&
                                                          math_transpose33(bc(loadcase)%rotation)
   write(6,'(a,f12.6)') 'temperature:', bc(loadcase)%temperature
   write(6,'(a,f12.6)') 'time:       ', bc(loadcase)%time
   write(6,'(a,i5)')    'increments: ', bc(loadcase)%incs
   write(6,'(a,i5)')    'output  frequency:  ', bc(loadcase)%outputfrequency
   write(6,'(a,i5)')    'restart frequency:  ', bc(loadcase)%restartfrequency

   if (any(bc(loadcase)%maskStress .eqv. bc(loadcase)%maskDeformation)) errorID = 831_pInt          ! exclusive or masking only
   if (any(bc(loadcase)%maskStress .and. transpose(bc(loadcase)%maskStress) .and. &
     reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
                                               errorID = 838_pInt                                   ! no rotation is allowed by stress BC
   if (any(abs(math_mul33x33(bc(loadcase)%rotation,math_transpose33(bc(loadcase)%rotation))&
                                      -math_I3) > reshape(spread(rotation_tol,1,9),[ 3,3]))&
                    .or. abs(math_det33(bc(loadcase)%rotation)) > 1.0_pReal + rotation_tol)&
                                               errorID = 846_pInt                                   ! given rotation matrix contains strain
   if (bc(loadcase)%time < 0.0_pReal)          errorID = 834_pInt                                   ! negative time increment
   if (bc(loadcase)%incs < 1_pInt)             errorID = 835_pInt                                   ! non-positive incs count
   if (bc(loadcase)%outputfrequency < 1_pInt)  errorID = 836_pInt                                   ! non-positive result frequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo
!##################################################################################################
! Loop over loadcases defined in the loadcase file
!##################################################################################################
 do loadcase = 1_pInt,  N_Loadcases
   time0 = time                                                                                     ! loadcase start time                
   if (bc(loadcase)%followFormerTrajectory .and. &
       (restartInc < totalIncsCounter .or. &
        restartInc > totalIncsCounter+bc(loadcase)%incs) ) then                                     ! continue to guess along former trajectory where applicable
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                                          ! change of load case, homogeneous guess for the first inc
   endif

!--------------------------------------------------------------------------------------------------
! arrays for mixed boundary conditions
   mask_defgrad = merge(ones,zeroes,bc(loadcase)%maskDeformation)                                   
   mask_stress  = merge(ones,zeroes,bc(loadcase)%maskStress)
   size_reduced = int(count(bc(loadcase)%maskStressVector), pInt)
   allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)

!##################################################################################################
! loop oper incs defined in input file for current loadcase
!##################################################################################################
   do inc = 1_pInt,  bc(loadcase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt                                                 

!--------------------------------------------------------------------------------------------------
! forwarding time
     timeinc_old = timeinc
     if (bc(loadcase)%logscale == 0_pInt) then                                                      ! linear scale
       timeinc = bc(loadcase)%time/bc(loadcase)%incs                                                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
     else
       if (loadcase == 1_pInt) then                                                                 ! 1st loadcase of logarithmic scale            
         if (inc == 1_pInt) then                                                                    ! 1st inc of 1st loadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(    1_pInt-bc(1)%incs ,pReal))                     ! assume 1st inc is equal to 2nd 
         else                                                                                       ! not-1st inc of 1st loadcase of logarithmic scale
           timeinc = bc(1)%time*(2.0_pReal**real(inc-1_pInt-bc(1)%incs ,pReal))
         endif
       else                                                                                         ! not-1st loadcase of logarithmic scale
           timeinc = time0 *( (1.0_pReal + bc(loadcase)%time/time0 )**(real(          inc,pReal)/&
                                                                  real(bc(loadcase)%incs ,pReal))&
                             -(1.0_pReal + bc(loadcase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                                   real(bc(loadcase)%incs ,pReal)) )
       endif
     endif
     time = time + timeinc

     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 
       if (bc(loadcase)%velGradApplied) then                                                        ! calculate deltaF_aim from given L and current F
         deltaF_aim = timeinc * mask_defgrad * math_mul33x33(bc(loadcase)%deformation, F_aim)
       else                                                                                         ! deltaF_aim = fDot *timeinc where applicable
         deltaF_aim = timeinc * mask_defgrad * bc(loadcase)%deformation
       endif

!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
       temp33_Real = F_aim                                            
       F_aim = F_aim &                                                                         
                  + guessmode * mask_stress * (F_aim - F_aim_lastInc)*timeinc/timeinc_old &      
                  + deltaF_aim
       F_aim_lastInc = temp33_Real

!--------------------------------------------------------------------------------------------------
! update local deformation gradient and coordinates
       deltaF_aim = math_rotate_backward33(deltaF_aim,bc(loadcase)%rotation)
       call 

       call deformed_fft(res,geomdim,math_rotate_backward33(F_aim,bc(loadcase)%rotation),&          ! calculate current coordinates
                                                          1.0_pReal,F_lastInc,coordinates)

!--------------------------------------------------------------------------------------------------
! calculate reduced compliance
       if(size_reduced > 0_pInt) then                                                               ! calculate compliance in case stress BC is applied
         C_lastInc = math_rotate_forward3333(C,bc(loadcase)%rotation)                               ! calculate stiffness from former inc
         temp99_Real = math_Plain3333to99(C_lastInc)
         k = 0_pInt                                                                                 ! build reduced stiffness
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
               if(bc(loadcase)%maskStressVector(m)) then
                 j = j + 1_pInt
                 c_reduced(k,j) = temp99_Real(n,m)
         endif; enddo; endif; enddo
         call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)                         ! invert reduced stiffness
         if(errmatinv) call IO_error(error_ID=400_pInt)
         temp99_Real = 0.0_pReal                                                                    ! build full compliance
         k = 0_pInt
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
             if(bc(loadcase)%maskStressVector(m)) then
                   j = j + 1_pInt
                   temp99_Real(n,m) = s_reduced(k,j)
         endif; enddo; endif; enddo
         S_lastInc = (math_Plain99to3333(temp99_Real))
       endif

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       write(6,'(a)') '##################################################################'
       write(6,'(A,I5.5,A,es12.5)') 'Increment ', totalIncsCounter, ' Time ',time
       
       guessmode = 1.0_pReal                                                                        ! keep guessing along former trajectory during same loadcase
       iter = 0_pInt
       err_div = huge(err_div_tol)                                                                  ! go into loop 

converged =  solution(mySolver,ForwardFields(solver,deltaF_aim,timeinc/timeinc_old,guessmode))
           
       CPFEM_mode = 1_pInt                                                                          ! winding forward
       C = C * wgt
       write(6,'(a)') ''
       write(6,'(a)') '=================================================================='
       if(err_div > err_div_tol .or. err_stress > err_stress_tol) then
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       else
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' converged'
       endif

       if (mod(inc,bc(loadcase)%outputFrequency) == 0_pInt) then                                    ! at output frequency
         write(6,'(a)') ''
         write(6,'(a)') '... writing results to file ......................................'
         write(538)  materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)       ! write result to file
         flush(538)
       endif
       
       if( bc(loadcase)%restartFrequency > 0_pInt .and. &
                      mod(inc,bc(loadcase)%restartFrequency) == 0_pInt) then                        ! at frequency of writing restart information set restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
         restartInc=totalIncsCounter
         restartWrite = .true.
         write(6,'(a)') 'writing converged results for restart'
         call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F))                        ! writing deformation gradient field to file
         write (777,rec=1) F
         close (777)
         call IO_write_jobBinaryFile(777,'C',size(C))
         write (777,rec=1) C
         close(777)
       endif 
       
     endif ! end calculation/forwarding
   enddo  ! end looping over incs in current loadcase
   deallocate(c_reduced)
   deallocate(s_reduced)
   enddo    ! end looping over loadcases
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
end program DAMASK_spectral
