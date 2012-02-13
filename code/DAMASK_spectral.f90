! Copyright 2012 Max-Planck-Institut fuer Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Duesseldorf Advanced Material Simulation Kit.
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
!##################################################################################################
! used modules
!##################################################################################################
program DAMASK_spectral

 use DAMASK_interface
 use prec,             only: pInt, pReal, DAMASK_NaN
 use IO 
 use debug,            only: debug_spectral, &
                             debug_spectralGeneral, &
                             debug_spectralDivergence, &
                             debug_spectralRestart, &
                             debug_spectralFFTW
 use math
 use kdtree2_module
 use CPFEM,            only: CPFEM_general, CPFEM_initAll
 use FEsolving,        only: restartWrite, restartInc
 use numerics,         only: err_div_tol, err_stress_tolrel, rotation_tol, itmax, &
                             memory_efficient, update_gamma, &
                             simplified_algorithm, divergence_correction, &
                             cut_off_value, &
                             DAMASK_NumThreadsInt, &
                             fftw_planner_flag, fftw_timelimit
 use homogenization,   only: materialpoint_sizeResults, materialpoint_results
 !$ use OMP_LIB                                                                                     ! the openMP function library
!##################################################################################################
! variable declaration
!##################################################################################################
 implicit none

!--------------------------------------------------------------------------------------------------
! variables to read from load case and geom file
 real(pReal), dimension(9) :: temp_valueVector                                                      ! stores information temporarily from loadcase file
 logical,     dimension(9) :: temp_maskVector
 integer(pInt), parameter  :: maxNchunksLoadcase = (1_pInt + 9_pInt)*3_pInt +&                      ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt +&                      ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                    1_pInt, &                                       ! dropguessing
                              maxNchunksGeom     = 7_pInt, &                                        ! 4 identifiers, 3 values
                              myUnit             = 234_pInt
 integer(pInt), dimension(1_pInt + maxNchunksLoadcase*2_pInt) :: positions                          ! this is longer than needed for geometry parsing
 integer(pInt) :: headerLength,&
                  N_l    = 0_pInt,&
                  N_t    = 0_pInt,&
                  N_n    = 0_pInt,&
                  N_Fdot = 0_pInt
 character(len=1024) :: path, line, keyword
 logical ::  gotResolution     = .false.,&
             gotDimension      = .false.,&
             gotHomogenization = .false.
 type bc_type
   real(pReal), dimension (3,3) :: deformation            = 0.0_pReal, &                            ! applied velocity gradient or time derivative of deformation gradient
                                   stress                 = 0.0_pReal, &                            ! stress BC (if applicable)
                                   rotation               = math_I3                                 ! rotation of BC (if applicable)
   real(pReal) ::                  time                   = 0.0_pReal, &                            ! length of increment
                                   temperature            = 300_pReal                               ! isothermal starting conditions
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
 
 type(bc_type), allocatable, dimension(:) ::  bc
 character(len=6) ::                          loadcase_string

!--------------------------------------------------------------------------------------------------
! variables storing information from geom file
 real(pReal) ::                               wgt
 real(pReal), dimension(3) ::                 geomdim = 0.0_pReal                                   ! physical dimension of volume element per direction
 integer(pInt) ::                             Npoints,&                                             ! number of Fourier points
                                              homog                                                 ! homogenization scheme used
 integer(pInt), dimension(3) ::               res = 1_pInt                                          ! resolution (number of Fourier points) in each direction
 integer(pInt)               ::               res1_red                                              ! to store res(1)/2 +1

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) ::           pstress, pstress_av, &
                                          defgradAim = math_I3, defgradAimOld = math_I3,&
                                          mask_stress, mask_defgrad, deltaF, &
                                          pstress_av_lab, defgradAim_lab, defgrad_av_lab            ! quantities rotated to other coordinate system
 real(pReal), dimension(3,3,3,3) ::       dPdF, c0_reference, c_current = 0.0_pReal, s_prev, c_prev,& ! stiffness and compliance
                                          s0_reference
 real(pReal), dimension(6) ::             cstress                                                   ! cauchy stress
 real(pReal), dimension(6,6) ::           dsde, c0_66, s0_66                                        ! small strain stiffness
 real(pReal), dimension(9,9) ::           s_prev99, c_prev99, c0_99, s0_99                          ! compliance and stiffness in matrix notation 
 real(pReal), dimension(:,:), allocatable ::  s_reduced, c_reduced                                  ! reduced compliance and stiffness (only for stress BC)
 real(pReal), dimension(6,6)   :: mask_inversion = reshape([&
                             1.0_pReal, 1.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal,&
                             1.0_pReal, 1.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal,&
                             1.0_pReal, 1.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal,&
                             0.0_pReal, 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal,&
                             0.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal,&
                             0.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal, 1.0_pReal],&
                                                                            [ 6_pInt, 6_pInt])
 real(pReal), dimension(3,3,3,3) :: temp_3333 = 0.0_pReal
 integer(pInt) ::                         size_reduced = 0.0_pReal                                  ! number of stress BCs

!--------------------------------------------------------------------------------------------------
! pointwise data 
 type(C_PTR) :: tensorField, tau                                                                    ! fields in real an fourier space
 real(pReal),    dimension(:,:,:,:,:), pointer :: tensorField_real                                  ! fields in real space (pointer)
 real(pReal),    dimension(:,:,:,:,:), pointer :: tau_real
 complex(pReal), dimension(:,:,:,:,:), pointer :: tensorField_fourier                               ! fields in fourier space (pointer)
 complex(pReal), dimension(:,:,:,:,:), pointer :: tau_fourier
 real(pReal),    dimension(:,:,:,:,:), allocatable ::  defgrad, defgradold
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 type(C_PTR) ::  plan_stress, plan_correction, plan_tau                                             ! plans for fftw
 real(pReal), dimension(3,3) ::                         xiDyad                                      ! product of wave vectors
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                   ! gamma operator (field) for spectral method
 real(pReal), dimension(:,:,:,:), allocatable ::        xi                                          ! wave vector field for divergence and for gamma operator
 integer(pInt), dimension(3) ::                         k_s, cutting_freq

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval 
 real(pReal) :: guessmode, err_div, err_stress, err_stress_tol              
 real(pReal), dimension(3,3), parameter ::  ones = 1.0_pReal, zeroes = 0.0_pReal
 complex(pReal), dimension(3)   ::          temp3_Complex
 complex(pReal), dimension(3,3) ::          temp33_Complex
 real(pReal),    dimension(3,3) ::          temp33_Real
 integer(pInt) :: i, j, k, l, m, n, p, errorID
 integer(pInt) :: N_Loadcases, loadcase, inc, iter, ielem, CPFEM_mode, &
                  ierr, totalIncsCounter = 0_pInt,&
                  notConvergedCounter = 0_pInt, convergedCounter = 0_pInt
 logical :: errmatinv
 real(pReal) :: defgradDet, correctionFactor

!--------------------------------------------------------------------------------------------------
!variables controlling debugging
 logical :: debugGeneral, debugDivergence, debugRestart, debugFFTW

!--------------------------------------------------------------------------------------------------
!variables for additional output due to general debugging
 real(pReal) :: defgradDetMax, defgradDetMin, maxCorrectionSym, maxCorrectionSkew, max_diag, max_offdiag

!--------------------------------------------------------------------------------------------------
! variables for additional output of divergence calculations
 type(C_PTR) :: divergence, plan_divergence
 real(pReal),    dimension(:,:,:,:), pointer :: divergence_real
 complex(pReal), dimension(:,:,:,:), pointer :: divergence_fourier
 real(pReal), dimension(:,:,:,:), allocatable  :: divergence_postProc
 real(pReal) :: pstress_av_L2, err_div_RMS, err_real_div_RMS,&
                err_div_max, err_real_div_max,&
                max_div_error

!--------------------------------------------------------------------------------------------------
! variables for debugging fft using a scalar field
 type(C_PTR) :: scalarField_realC, scalarField_fourierC,&
                plan_scalarField_forth, plan_scalarField_back
 complex(pReal), dimension(:,:,:), pointer :: scalarField_real
 complex(pReal), dimension(:,:,:), pointer :: scalarField_fourier
 integer(pInt) :: row, column

!##################################################################################################
! reading of information from load case file and geometry file
!##################################################################################################
 !$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                  ! set number of threads for parallel execution set by DAMASK_NUM_THREADS
 call DAMASK_interface_init()

 print '(a)', ''
 print '(a)', ' <<<+-  DAMASK_spectral init  -+>>>'
 print '(a)', ' $Id$'
#include "compilation_info.f90"
 print '(a,a)', ' Working Directory:    ',trim(getSolverWorkingDirectoryName())
 print '(a,a)', ' Solver Job Name:      ',trim(getSolverJobName())
 print '(a)', ''
!--------------------------------------------------------------------------------------------------
! reading the load case file and allocate data structure containing load cases
 path = getLoadcaseName()
 if (.not. IO_open_file(myUnit,path)) call IO_error(error_ID = 30_pInt,ext_msg = trim(path))
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
            case('n','incs','increments','steps','logincs','logsteps')
                 N_n = N_n + 1_pInt
        end select
   enddo                                                                                            ! count all identifiers to allocate memory and do sanity check
 enddo

100 N_Loadcases = N_n
 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                                                     ! sanity check
   call IO_error(error_ID=37_pInt,ext_msg = trim(path))                                             ! error message for incomplete loadcase

 allocate (bc(N_Loadcases))

!--------------------------------------------------------------------------------------------------
! reading the load case and assign values to the allocated data structure
 rewind(myUnit)
 loadcase = 0_pInt
 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   loadcase = loadcase + 1_pInt
   positions = IO_stringPos(line,maxNchunksLoadcase)
   do j = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,positions,j)))
       case('fdot','l','velocitygrad','velgrad','velocitygradient')                                 ! assign values for the deformation BC matrix
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
 if (update_gamma .and. .not. memory_efficient) call IO_error(error_ID = 47_pInt)

!--------------------------------------------------------------------------------------------------
! read header of geom file to get size information. complete geom file is intepretated by mesh.f90
 path = getModelName()

 if (.not. IO_open_file(myUnit,trim(path)//InputFileExtension))&
        call IO_error(error_ID=101_pInt,ext_msg = trim(path)//InputFileExtension)
 rewind(myUnit)
 read(myUnit,'(a1024)') line
 positions = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,positions,2_pInt))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,positions,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=42_pInt)
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   positions = IO_stringPos(line,maxNchunksGeom)             
   select case ( IO_lc(IO_StringValue(line,positions,1)) )
     case ('dimension')
       gotDimension = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,positions,j)))
           case('x')
              geomdim(1) = IO_floatValue(line,positions,j+1_pInt)
           case('y')
              geomdim(2) = IO_floatValue(line,positions,j+1_pInt)
           case('z')
              geomdim(3) = IO_floatValue(line,positions,j+1_pInt)
         end select
       enddo
     case ('homogenization')
       gotHomogenization = .true.
       homog = IO_intValue(line,positions,2_pInt)
     case ('resolution')
       gotResolution = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,positions,j)))
           case('a')
             res(1) = IO_intValue(line,positions,j+1_pInt)
           case('b')
             res(2) = IO_intValue(line,positions,j+1_pInt)
           case('c')
             res(3) = IO_intValue(line,positions,j+1_pInt)
         end select
       enddo
   end select
 enddo
 close(myUnit)

!--------------------------------------------------------------------------------------------------
! sanity checks of geometry parameters
 if (.not.(gotDimension .and. gotHomogenization .and. gotResolution))&
                                                                call IO_error(error_ID = 45_pInt)
 if (any(geomdim<=0.0_pReal)) call IO_error(error_ID = 102_pInt)
 if(mod(res(1),2_pInt)/=0_pInt .or.&
    mod(res(2),2_pInt)/=0_pInt .or.&
   (mod(res(3),2_pInt)/=0_pInt .and. res(3)/= 1_pInt))&
                                                                call IO_error(error_ID = 103_pInt)

!--------------------------------------------------------------------------------------------------
! variables derived from resolution
 res1_red = res(1)/2_pInt + 1_pInt                                                                  ! size of complex array in first dimension (c2r, r2c)
 Npoints = res(1)*res(2)*res(3)
 wgt = 1.0_pReal/real(Npoints, pReal)
 if (cut_off_value <0.0_pReal .or. cut_off_value >0.9_pReal) stop
 cutting_freq = nint(real(res,pReal)*cut_off_value,pInt)                                            ! for cut_off_value=0.0 just the highest freq. is removed

!--------------------------------------------------------------------------------------------------
! output of geometry
 print '(a)',  ''
 print '(a)',   '#############################################################'
 print '(a)',   'DAMASK spectral:'
 print '(a)',   'The spectral method boundary value problem solver for'
 print '(a)',   'the Duesseldorf Advanced Material Simulation Kit'
 print '(a)',   '#############################################################'
 print '(a,a)', 'geometry file:        ',trim(path)//'.geom'
 print '(a)',   '============================================================='
 print '(a,3(i12  ))','resolution a b c:', res
 print '(a,3(f12.5))','dimension  x y z:', geomdim
 print '(a,i5)','homogenization:       ',homog
 if(cut_off_value/=0.0_pReal) print '(a,3(i12),a)', 'cutting away    ', cutting_freq, ' frequencies'
 print '(a)',   '#############################################################'
 print '(a,a)', 'loadcase file:        ',trim(getLoadcaseName())

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 bc(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first loadcase
 errorID = 0_pInt
 do loadcase = 1_pInt, N_Loadcases
   write (loadcase_string, '(i6)' ) loadcase

   print '(a)', '============================================================='
   print '(a,i6)', 'loadcase:            ', loadcase

   if (.not. bc(loadcase)%followFormerTrajectory) print '(a)', 'drop guessing along trajectory'
   if (bc(loadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .true.) .and. &
           any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .false.)) errorID = 32_pInt                ! each row should be either fully or not at all defined
     enddo
     print '(a)','velocity gradient:'
   else
     print '(a)','deformation gradient rate:'
   endif
   write (*,'(3(3(f12.7,1x)/))',advance='no') merge(math_transpose33(bc(loadcase)%deformation),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskDeformation))
   write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') ' stress / GPa:',&
        1e-9*merge(math_transpose33(bc(loadcase)%stress),reshape(spread(DAMASK_NaN,1,9),[ 3,3])&
                                                        ,transpose(bc(loadcase)%maskStress))
   if (any(bc(loadcase)%rotation /= math_I3)) &
     write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') ' rotation of loadframe:',&
                                                          math_transpose33(bc(loadcase)%rotation)
   print '(a,f12.6)','temperature:',bc(loadcase)%temperature
   print '(a,f12.6)','time:       ',bc(loadcase)%time
   print '(a,i5)'   ,'increments: ',bc(loadcase)%incs
   print '(a,i5)','output  frequency:  ',bc(loadcase)%outputfrequency
   print '(a,i5)','restart frequency:  ',bc(loadcase)%restartfrequency

   if (any(bc(loadcase)%maskStress .eqv. bc(loadcase)%maskDeformation)) errorID = 31                ! exclusive or masking only
   if (any(bc(loadcase)%maskStress .and. transpose(bc(loadcase)%maskStress) .and. &
     reshape([ .false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false.],[ 3,3]))) &
                                               errorID = 38_pInt                                    ! no rotation is allowed by stress BC
   if (any(abs(math_mul33x33(bc(loadcase)%rotation,math_transpose33(bc(loadcase)%rotation))&
                                      -math_I3) > reshape(spread(rotation_tol,1,9),[ 3,3]))&
                    .or. abs(math_det33(bc(loadcase)%rotation)) > 1.0_pReal + rotation_tol)&
                                               errorID = 46_pInt                                    ! given rotation matrix contains strain
   if (bc(loadcase)%time < 0.0_pReal)          errorID = 34_pInt                                    ! negative time increment
   if (bc(loadcase)%incs < 1_pInt)             errorID = 35_pInt                                    ! non-positive incs count
   if (bc(loadcase)%outputfrequency < 1_pInt)  errorID = 36_pInt                                    ! non-positive result frequency
   if (errorID > 0_pInt) call IO_error(error_ID = errorID, ext_msg = loadcase_string)
 enddo

!--------------------------------------------------------------------------------------------------
! debugging parameters
 debugGeneral    = iand(debug_spectral,debug_spectralGeneral)    > 0_pInt
 debugDivergence = iand(debug_spectral,debug_spectralDivergence) > 0_pInt
 debugRestart    = iand(debug_spectral,debug_spectralRestart)    > 0_pInt
 debugFFTW       = iand(debug_spectral,debug_spectralFFTW)       > 0_pInt
 
!##################################################################################################
! initialization 
!##################################################################################################

 allocate (defgrad    (  res(1),  res(2),res(3),3,3));  defgrad     = 0.0_pReal
 allocate (defgradold (  res(1),  res(2),res(3),3,3));  defgradold  = 0.0_pReal
 allocate (coordinates(  res(1),  res(2),res(3),3));    coordinates = 0.0_pReal
 allocate (temperature(  res(1),  res(2),res(3)));      temperature = bc(1)%temperature       ! start out isothermally
 allocate (xi         (3,res1_red,res(2),res(3)));  xi          = 0.0_pReal
 tensorField = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))                ! allocate continous data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, tensorField_real,    [ res(1)+2_pInt,res(2),res(3),3,3])       ! place a pointer for the real representation
 call c_f_pointer(tensorField, tensorField_fourier, [ res1_red,     res(2),res(3),3,3])       ! place a pointer for the complex representation
 
!--------------------------------------------------------------------------------------------------
! init fields to no deformation
 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt 
   defgrad(i,j,k,1:3,1:3) = math_I3
   defgradold(i,j,k,1:3,1:3) = math_I3
   coordinates(i,j,k,1:3) = geomdim/real(res, pReal)*[i,j,k] - geomdim/real(2_pInt*res,pReal)
   call CPFEM_general(2_pInt,coordinates(i,j,k,1:3),math_I3,math_I3,temperature(i,j,k),&
                                       0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
   c_current = c_current + dPdF
 enddo; enddo; enddo

 c0_reference = c_current * wgt                                                              ! linear reference material stiffness
 c0_66 = math_Mandel3333to66(c0_reference)
 call math_invert(6_pInt, c0_66, s0_66, i, errmatinv)                                        ! invert in mandel notation
 if(errmatinv) call IO_error(error_ID=800_pInt)
 s0_reference = math_Mandel66to3333(s0_66)

!--------------------------------------------------------------------------------------------------
! possible restore deformation gradient from saved state
 if (restartInc > 1_pInt) then                                                             ! using old values from file                                                      
   if (debugRestart) print '(a,i6,a)' , 'Reading values of increment ',&
                                             restartInc - 1_pInt,' from file' 
   if (IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                trim(getSolverJobName()),size(defgrad))) then
     read (777,rec=1) defgrad
     close (777)
   endif
   defgradold = defgrad
   defgradAim = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     defgradAim = defgradAim + defgrad(i,j,k,1:3,1:3)                                         ! calculating old average deformation
   enddo; enddo; enddo
   defgradAim = defgradAim * wgt
   defgradAimOld = defgradAim
 endif
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around) and remove the given highest frequencies
 do k = 1_pInt, res(3)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1, res1_red
           k_s(1) = i - 1_pInt
           xi(1:3,i,j,k) = real(k_s, pReal)/geomdim
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! calculate the gamma operator
 if(memory_efficient) then                                                                           ! allocate just single fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3)); gamma_hat = 0.0_pReal
 else                                                                                                ! precalculation of gamma_hat field
   allocate (gamma_hat(res1_red ,res(2),res(3),3,3,3,3)); gamma_hat = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
   ! if(k==res(3)/2 .or. k==res(3)/2+2 .or.&
   !    j==res(2)/2 .or. j==res(2)/2+2 .or.&
   !    i==res(1)/2 .or. i==res(1)/2+2) then 
   !      gamma_hat(i,j,k,1:3,1:3,1:3,1:3) = s0_reference
   ! else
     forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
       xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
     forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
       temp33_Real(l,m) = sum(c0_reference(l,1:3,m,1:3)*xiDyad)
     temp33_Real = math_inv33(temp33_Real)
     forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, p=1_pInt:3_pInt)&
       gamma_hat(i,j,k, l,m,n,p) =  temp33_Real(l,n)*xiDyad(m,p)
   ! endif  
   enddo; enddo; enddo
   gamma_hat(1,1,1, 1:3,1:3,1:3,1:3) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif

!--------------------------------------------------------------------------------------------------
! general initialization of fftw (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=108_pInt)                     ! check for correct precision in C
#ifdef _OPENMP
    if(DAMASK_NumThreadsInt > 0_pInt) then
      ierr = fftw_init_threads()
      if (ierr == 0_pInt) call IO_error(error_ID = 109_pInt)
      call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
    endif
#endif
 call fftw_set_timelimit(fftw_timelimit)                                                      ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans
 plan_stress =    fftw_plan_many_dft_r2c(3,[ res(3),res(2) ,res(1)],9,&                     ! dimensions , length in each dimension in reversed order
                          tensorField_real,[ res(3),res(2) ,res(1)+2_pInt],&                ! input data , physical length in each dimension in reversed order
                                         1,  res(3)*res(2)*(res(1)+2_pInt),&                 ! striding   , product of physical lenght in the 3 dimensions
                       tensorField_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 plan_correction =fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],9,&
                       tensorField_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,&
                          tensorField_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                         1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)

!--------------------------------------------------------------------------------------------------
! depending on (debug) options, allocate more memory and create additional plans 
 if (.not. simplified_algorithm) then
   print*, 'using polarization field based algorithm'
   tau = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))
   call c_f_pointer(tau, tau_real,    [ res(1)+2_pInt,res(2),res(3),3,3])
   call c_f_pointer(tau, tau_fourier, [ res1_red,     res(2),res(3),3,3])
   plan_tau         = fftw_plan_many_dft_r2c(3,[ res(3),res(2) ,res(1)],9,&
                                      tau_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                             1,  res(3)*res(2)*(res(1)+2_pInt),&
                                   tau_fourier,[ res(3),res(2) ,res1_red],&
                                             1,  res(3)*res(2)* res1_red,fftw_planner_flag)   
 endif

 if (debugDivergence) then
   divergence = fftw_alloc_complex(int(res1_red*res(2)*res(3)*3_pInt,C_SIZE_T))
   call c_f_pointer(divergence, divergence_real,    [ res(1)+2_pInt,res(2),res(3),3])
   call c_f_pointer(divergence, divergence_fourier, [ res1_red,     res(2),res(3),3])
   allocate (divergence_postProc(res(1),res(2),res(3),3));  divergence_postProc= 0.0_pReal
   plan_divergence = fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],3,&
                           divergence_fourier,[ res(3),res(2) ,res1_red],&
                                            1,  res(3)*res(2)* res1_red,&
                              divergence_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                            1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 endif

 if (debugFFTW) then
   scalarField_realC    = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))                    ! do not do an inplace transform  
   scalarField_fourierC = fftw_alloc_complex(int(res(1)*res(2)*res(3),C_SIZE_T))
   call c_f_pointer(scalarField_realC,    scalarField_real,    [res(1),res(2),res(3)])
   call c_f_pointer(scalarField_fourierC, scalarField_fourier, [res(1),res(2),res(3)])
   plan_scalarField_forth = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_real,scalarField_fourier,-1,fftw_planner_flag)
   plan_scalarField_back  = fftw_plan_dft_3d(res(3),res(2),res(1),&                                 !reversed order
                                      scalarField_fourier,scalarField_real,+1,fftw_planner_flag)
 endif 

 if (debugGeneral) print '(a)' , 'FFTW initialized'

!--------------------------------------------------------------------------------------------------
! do not correct divergence criterion (usefull to kill dimension and resolution dependenc)
 correctionFactor = 1.0_pReal

!--------------------------------------------------------------------------------------------------
! write header of output file
 open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())&
                                        //'.spectralOut',form='UNFORMATTED',status='REPLACE')
 write(538) 'load',       trim(getLoadcaseName())
 write(538) 'workingdir', trim(getSolverWorkingDirectoryName())
 write(538) 'geometry',   trim(getSolverJobName())//InputFileExtension
 write(538) 'resolution', res
 write(538) 'dimension',  geomdim
 write(538) 'materialpoint_sizeResults', materialpoint_sizeResults
 write(538) 'loadcases',        N_Loadcases
 write(538) 'frequencies', bc(1:N_Loadcases)%outputfrequency                                 ! one entry per loadcase
 write(538) 'times', bc(1:N_Loadcases)%time                                                  ! one entry per loadcase
 write(538) 'logscales',  bc(1:N_Loadcases)%logscale         
 write(538) 'increments', bc(1:N_Loadcases)%incs                                             ! one entry per loadcase
 write(538) 'startingIncrement', restartInc - 1_pInt                                         ! start with writing out the previous inc
 write(538) 'eoh'                                                                            ! end of header
 write(538) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)         ! initial (non-deformed or read-in) results
 if (debugGeneral) print '(a)' , 'Header of result file written out'

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
   size_reduced = count(bc(loadcase)%maskStressVector)
   allocate (c_reduced(size_reduced,size_reduced));          c_reduced = 0.0_pReal
   allocate (s_reduced(size_reduced,size_reduced));          s_reduced = 0.0_pReal


!##################################################################################################
! loop oper incs defined in input file for current loadcase
!##################################################################################################
   do inc = 1_pInt,  bc(loadcase)%incs
     totalIncsCounter = totalIncsCounter + 1_pInt
     if(totalIncsCounter >= restartInc)  then                                                       ! do calculations (otherwise just forwarding) 
     
!--------------------------------------------------------------------------------------------------
! forwarding time
       timeinc_old = timeinc
       if (bc(loadcase)%logscale == 0_pInt) then                                                    ! linear scale
         timeinc = bc(loadcase)%time/bc(loadcase)%incs                                              ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
       else
         if (loadcase == 1_pInt) then                                                               ! 1st loadcase of logarithmic scale            
           if (inc == 1_pInt) then                                                                  ! 1st inc of 1st loadcase of logarithmic scale
             timeinc = bc(1)%time*(2.0_pReal**real(    1_pInt-bc(1)%incs ,pReal))                   ! assume 1st inc is equal to 2nd 
           else                                                                                     ! not-1st inc of 1st loadcase of logarithmic scale
             timeinc = bc(1)%time*(2.0_pReal**real(inc-1_pInt-bc(1)%incs ,pReal))
           endif
         else                                                                                       ! not-1st loadcase of logarithmic scale
             timeinc = time0 *( (1.0_pReal + bc(loadcase)%time/time0 )**(real(          inc,pReal)/&
                                                                    real(bc(loadcase)%incs ,pReal))&
                               -(1.0_pReal + bc(loadcase)%time/time0 )**(real( (inc-1_pInt),pReal)/&
                                                                    real(bc(loadcase)%incs ,pReal)) )
         endif
       endif
       time = time + timeinc

       if (bc(loadcase)%velGradApplied) then                                                        ! calculate deltaF from given L and current F
         deltaF = timeinc * mask_defgrad * math_mul33x33(bc(loadcase)%deformation, defgradAim)
       else                                                                                         ! deltaF = fDot *timeinc where applicable
         deltaF = timeinc * mask_defgrad * bc(loadcase)%deformation
       endif

!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
       temp33_Real = defgradAim                                            
       defgradAim = defgradAim &                                                                         
                  + guessmode * mask_stress * (defgradAim - defgradAimOld)*timeinc/timeinc_old &      
                  + deltaF
       defgradAimOld = temp33_Real

!--------------------------------------------------------------------------------------------------
! update local deformation gradient
       deltaF = math_rotate_backward33(deltaF,bc(loadcase)%rotation)
       do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
         temp33_Real = defgrad(i,j,k,1:3,1:3)
         defgrad(i,j,k,1:3,1:3) = defgrad(i,j,k,1:3,1:3) &                                          ! decide if guessing along former trajectory or apply homogeneous addon
                                + guessmode * (defgrad(i,j,k,1:3,1:3) - defgradold(i,j,k,1:3,1:3))& ! guessing... 
                                            *timeinc/timeinc_old &
                                + (1.0_pReal-guessmode) * deltaF                                    ! if not guessing, use prescribed average deformation where applicable
         defgradold(i,j,k,1:3,1:3) = temp33_Real 
       enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! calculate reduced compliance
       c_prev = math_rotate_forward3333(c_current*wgt,bc(loadcase)%rotation)                        ! calculate stiffness from former inc
       if(size_reduced > 0_pInt) then                                                               ! calculate compliance in case stress BC is applied
         c_prev99 = math_Plain3333to99(c_prev)
         k = 0_pInt                                                                                 ! build reduced stiffness
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
               if(bc(loadcase)%maskStressVector(m)) then
                 j = j + 1_pInt
                 c_reduced(k,j) = c_prev99(n,m)
         endif; enddo; endif; enddo
         call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)                         ! invert reduced stiffness
         if(errmatinv) call IO_error(error_ID=800_pInt)
         s_prev99 = 0.0_pReal                                                                       ! build full compliance
         k = 0_pInt
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
             if(bc(loadcase)%maskStressVector(m)) then
                   j = j + 1_pInt
                   s_prev99(n,m) = s_reduced(k,j)
         endif; enddo; endif; enddo
         s_prev = (math_Plain99to3333(s_prev99))
       endif

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       print '(a)', '##################################################################'
       print '(A,I5.5,A,es12.6)', 'Increment ', totalIncsCounter, ' Time ',time
       
       guessmode = 1.0_pReal                                                                        ! keep guessing along former trajectory during same loadcase
       CPFEM_mode = 1_pInt                                                                          ! winding forward
       iter = 0_pInt
       err_div = 2.0_pReal * err_div_tol                                                            ! go into loop 

!##################################################################################################
! convergence loop (looping over iterations)
!##################################################################################################
       do while(iter < itmax .and. &
               (err_div     > err_div_tol    .or. &
                err_stress  > err_stress_tol)) 
         iter = iter + 1_pInt

!--------------------------------------------------------------------------------------------------
! report begin of new iteration
         print '(a)', ''
         print '(a)', '=================================================================='
         print '(5(a,i6.6))', 'Loadcase ',loadcase,' Increment ',inc,'/',bc(loadcase)%incs,&
                                                                  ' @ Iteration ',iter,'/',itmax
         do n = 1_pInt,3_pInt; do m = 1_pInt,3_pInt
           defgrad_av_lab(m,n) = sum(defgrad(1:res(1),1:res(2),1:res(3),m,n)) * wgt
         enddo; enddo
         write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'deformation gradient:',&
                    math_transpose33(math_rotate_forward33(defgrad_av_lab,bc(loadcase)%rotation))
         print '(a)', ''
         print '(a)', '... update stress field P(F) .....................................'

!--------------------------------------------------------------------------------------------------
! evaluate constitutive response
         call deformed_fft(res,geomdim,defgrad_av_lab,1.0_pReal,defgrad,coordinates)                ! calculate current coordinates
         ielem = 0_pInt
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(3_pInt,&                                                              ! collect cycle
                              coordinates(i,j,k,1:3), defgradold(i,j,k,1:3,1:3), defgrad(i,j,k,1:3,1:3),&
                              temperature(i,j,k),timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo; enddo

         tensorField_real = 0.0_pReal                                                               ! needed because of the padding for FFTW
         c_current = 0.0_pReal
         ielem = 0_pInt 
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(CPFEM_mode,&                                                          ! first element in first iteration retains CPFEM_mode 1, 
                              coordinates(i,j,k,1:3),&
                              defgradold(i,j,k,1:3,1:3), defgrad(i,j,k,1:3,1:3),&                   ! others get 2 (saves winding forward effort)
                              temperature(i,j,k),timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress,dPdF)
           CPFEM_mode = 2_pInt
           tensorField_real(i,j,k,1:3,1:3) = pstress
           c_current = c_current + dPdF
         enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! copy one component of the stress field to to a single FT and check for mismatch
         if (debugFFTW) then
           row =    (mod(totalIncsCounter+iter-2_pInt,9_pInt))/3_pInt + 1_pInt                      ! go through the elements of the tensors, controlled by totalIncsCounter and iter, starting at 1
           column = (mod(totalIncsCounter+iter-2_pInt,3_pInt))        + 1_pInt
           scalarField_real(1:res(1),1:res(2),1:res(3)) =&                                          ! store the selected component
                  tensorField_real(1:res(1),1:res(2),1:res(3),row,column)
         endif

!--------------------------------------------------------------------------------------------------
! build polarization field
         if (.not. simplified_algorithm) then
           tau_real = 0.0_pReal                                                                     ! padding
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
            tau_real(i,j,k,1:3,1:3)&
                            = tensorField_real(i,j,k,1:3,1:3) &
                            - math_mul3333xx33(c0_reference,defgrad(i,j,k,1:3,1:3))
           enddo; enddo; enddo
           call fftw_execute_dft_r2c(plan_tau,tau_real,tau_fourier)
         endif

!--------------------------------------------------------------------------------------------------
! call function to calculate divergence from math (for post processing) to check results
         if (debugDivergence) &
              call divergence_fft(res,geomdim,3_pInt,&
              tensorField_real(1:res(1),1:res(2),1:res(3),1:3,1:3),divergence_postProc)             !padding
              
!--------------------------------------------------------------------------------------------------
! doing the FT because it simplifies calculation of average stress in real space also
         call fftw_execute_dft_r2c(plan_stress,tensorField_real,tensorField_fourier)

         pstress_av_lab = real(tensorField_fourier(1,1,1,1:3,1:3),pReal)*wgt
         pstress_av = math_rotate_forward33(pstress_av_lab,bc(loadcase)%rotation)
         write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'Piola-Kirchhoff stress / MPa:',&
                                                          math_transpose33(pstress_av)/1.e6

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 FT results
         if (debugFFTW) then
           call fftw_execute_dft(plan_scalarField_forth,scalarField_real,scalarField_fourier)
           print '(a,i1,1x,i1)', 'checking FT results of compontent ', row, column
           print '(a,2(es10.4,1x))',  'max FT relative error ',&
             maxval( real((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-& 
                           tensorField_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                           scalarField_fourier(1:res1_red,1:res(2),1:res(3)))), &
             maxval(aimag((scalarField_fourier(1:res1_red,1:res(2),1:res(3))-&
                           tensorField_fourier(1:res1_red,1:res(2),1:res(3),row,column))/&
                           scalarField_fourier(1:res1_red,1:res(2),1:res(3))))
         endif

!--------------------------------------------------------------------------------------------------
! removing highest frequencies
         tensorField_fourier  (  res1_red,1:res(2) ,             1:res(3)              ,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         tensorField_fourier  (1:res1_red,  res(2)/2_pInt+1_pInt,1:res(3)              ,1:3,1:3)& 
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         if(res(3)>1_pInt) &
          tensorField_fourier (1:res1_red,1:res(2),                res(3)/2_pInt+1_pInt,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! stress BC handling
         if(size_reduced > 0_pInt) then                                                              ! calculate stress BC if applied
           err_stress = maxval(abs(mask_stress * (pstress_av - bc(loadcase)%stress)))                ! maximum deviaton (tensor norm not applicable)
           err_stress_tol = maxval(abs(pstress_av)) * err_stress_tolrel                              ! don't use any tensor norm because the comparison should be coherent
           print '(a)', '' 
           print '(a)', '... correcting deformation gradient to fulfill BCs ...............'
           print '(a,f6.2,a,es10.4,a)', 'error stress = ', err_stress/err_stress_tol, &
                                                                            ' (',err_stress,' Pa)'  
           defgradAim = defgradAim - math_mul3333xx33(s_prev, ((pstress_av - bc(loadcase)%stress))) ! residual on given stress components
           if(debugGeneral)  write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'new deformation aim:',&
                                                            math_transpose33(defgradAim)
           print '(a,1x,es10.4)'         , 'determinant of new deformation: ', math_det33(defgradAim)
         else
           err_stress_tol = +huge(1.0_pReal)
         endif
                                            
         defgradAim_lab = math_rotate_backward33(defgradAim,bc(loadcase)%rotation)                  ! boundary conditions from load frame into lab (Fourier) frame

!--------------------------------------------------------------------------------------------------
! actual spectral method         
         print '(a)', ''
         print '(a)', '... calculating equilibrium with spectral method .................'

!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
         pstress_av_L2 = sqrt(maxval (math_eigenvalues33(math_mul33x33(pstress_av_lab,&                 ! L_2 norm of average stress
                                                  math_transpose33(pstress_av_lab)))))
         err_div_RMS = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2)
           do i = 2_pInt, res1_red -1_pInt                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
             err_div_RMS = err_div_RMS &
                   + 2.0_pReal*(sum (real(math_mul33x3_complex(tensorField_fourier(i,j,k,1:3,1:3),& ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                   xi(1:3,i,j,k))*two_pi_img)**2.0_pReal)&          ! --> sum squared L_2 norm of vector 
                               +sum(aimag(math_mul33x3_complex(tensorField_fourier(i,j,k,1:3,1:3),& 
                                                   xi(1:3,i,j,k))*two_pi_img)**2.0_pReal))           
           enddo
           err_div_RMS = err_div_RMS &                                                              ! Those two layers (DC and Nyquist) do not have a conjugate complex counterpart
                      + sum(real(math_mul33x3_complex(tensorField_fourier(1       ,j,k,1:3,1:3),&
                                       xi(1:3,1       ,j,k))*two_pi_img)**2.0_pReal)&
                      + sum(aimag(math_mul33x3_complex(tensorField_fourier(1       ,j,k,1:3,1:3),&
                                       xi(1:3,1       ,j,k))*two_pi_img)**2.0_pReal)&
                      + sum(real(math_mul33x3_complex(tensorField_fourier(res1_red,j,k,1:3,1:3),&
                                       xi(1:3,res1_red,j,k))*two_pi_img)**2.0_pReal)&
                      + sum(aimag(math_mul33x3_complex(tensorField_fourier(res1_red,j,k,1:3,1:3),&
                                       xi(1:3,res1_red,j,k))*two_pi_img)**2.0_pReal)
         enddo; enddo
         err_div_RMS = sqrt(err_div_RMS)*wgt                                                        ! RMS in real space calculated with Parsevals theorem from Fourier space
         if(err_div_RMS/pstress_av_L2*sqrt(wgt) * correctionFactor>err_div&
                       .and.iter >2_pInt&
                       .and.err_stress  < err_stress_tol) then
           print*, 'Increasing divergence, stopping iterations'
           iter = itmax
         endif
         err_div = err_div_RMS/pstress_av_L2*sqrt(wgt) * correctionFactor                               ! criterion to stop iterations

!--------------------------------------------------------------------------------------------------
! calculate additional divergence criteria and report
         if(debugDivergence) then                                                                   ! calculate divergence again
           err_div_max = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
             temp3_Complex = math_mul33x3_complex(tensorField_fourier(i,j,k,1:3,1:3),&
                                                xi(1:3,i,j,k))*two_pi_img
             err_div_max = max(err_div_max,sqrt(sum(abs(temp3_Complex)**2.0_pReal)))
             divergence_fourier(i,j,k,1:3) = temp3_Complex                                          ! need divergence NOT squared
           enddo; enddo; enddo

           call fftw_execute_dft_c2r(plan_divergence,divergence_fourier,divergence_real)
           divergence_real = divergence_real*wgt
           err_real_div_RMS = 0.0_pReal
           err_real_div_max = 0.0_pReal
           max_div_error = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             max_div_error= max(max_div_error,maxval((divergence_real(i,j,k,1:3)&
                        -divergence_postProc(i,j,k,1:3))/divergence_real(i,j,k,1:3)))
             err_real_div_RMS = err_real_div_RMS    +      sum(divergence_real(i,j,k,1:3)**2.0_pReal)        ! avg of L_2 norm of div(stress) in real space
             err_real_div_max = max(err_real_div_max, sqrt(sum(divergence_real(i,j,k,1:3)**2.0_pReal)))      ! maximum of L two norm of div(stress) in real space
           enddo; enddo; enddo
           err_real_div_RMS = sqrt(wgt*err_real_div_RMS)                                   ! RMS in real space
           err_div_max      = err_div_max*sqrt(wgt) 
           
           print '(a,es10.4)',        'error divergence  FT  RMS = ',err_div_RMS
           print '(a,es10.4)',        'error divergence  FT  max = ',err_div_max
           print '(a,es10.4)',        'error divergence Real RMS = ',err_real_div_RMS
           print '(a,es10.4)',        'error divergence Real max = ',err_real_div_max
           print '(a,es10.4)',        'divergence RMS FT/real    = ',err_div_RMS/err_real_div_RMS
           print '(a,es10.4)',        'divergence max FT/real    = ',err_div_max/err_real_div_max
           print '(a,es10.4)',        'max deviat. from postProc = ',max_div_error
         endif
         print '(a,f6.2,a,es10.4,a)', 'error divergence = ', err_div/err_div_tol, ' (',err_div,' 1/m)'    
         
!--------------------------------------------------------------------------------------------------
! divergence is calculated from FT(stress), depending on algorithm use field for spectral method
         if (.not. simplified_algorithm) tensorField_fourier = tau_fourier
 max_diag = tiny(1.0_pReal)
 max_offdiag = tiny(1.0_pReal)
!--------------------------------------------------------------------------------------------------
! to the actual spectral method calculation (mechanical equilibrium)
         if(memory_efficient) then                                                                           ! memory saving version, on-the-fly calculation of gamma_hat
           
           do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res1_red
             ! if(k==res(3)/2 .or. k==res(3)/2+2 .or.&
                ! j==res(2)/2 .or. j==res(2)/2+2 .or.&
                ! i==res(1)/2 .or. i==res(1)/2+2) then 
               ! forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt)&
                 ! temp33_Complex(m,n) = sum(s0_reference(m,n, 1:3,1:3)* tensorField_fourier(i,j,k,1:3,1:3))
             ! else
               forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                 xiDyad(l,m) = xi(l, i,j,k)*xi(m, i,j,k)
               forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                 temp33_Real(l,m) = sum(c0_reference(l,1:3,m,1:3)*xiDyad)
               temp33_Real = math_inv33(temp33_Real)
               forall(l=1_pInt:3_pInt, m=1_pInt:3_pInt, n=1_pInt:3_pInt, p=1_pInt:3_pInt)&
                 gamma_hat(1,1,1, l,m,n,p) =  temp33_Real(l,n)*xiDyad(m,p)
               forall(l = 1_pInt:3_pInt, m = 1_pInt:3_pInt) &
                 temp33_Complex(l,m) = sum(gamma_hat(1,1,1, l,m, 1:3,1:3) * tensorField_fourier(i,j,k,1:3,1:3))
               tensorField_fourier(i,j,k,1:3,1:3) = temp33_Complex 
            ! endif             
           enddo; enddo; enddo
   
         else                                                                                           ! use precalculated gamma-operator
           
           do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
             forall( m = 1_pInt:3_pInt, n = 1_pInt:3_pInt) &
               temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n, 1:3,1:3) * tensorField_fourier(i,j,k,1:3,1:3))
             tensorField_fourier(i,j,k, 1:3,1:3) = temp33_Complex
           enddo; enddo; enddo

         endif

         if (simplified_algorithm) then                                                               ! do not use the polarization field based algorithm
           tensorField_fourier(1,1,1,1:3,1:3) = (defgrad_av_lab - defgradAim_lab) &                   ! assign (negative) average deformation gradient change to zero frequency (real part)
                                                * real(Npoints,pReal)                                 ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1
         else
           tensorField_fourier(1,1,1,1:3,1:3) = defgradAim_lab  * real(Npoints,pReal)                 ! assign deformation aim to zero frequency (real part)
         endif

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
         if (debugFFTW) then
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
              scalarField_fourier(i,j,k) = tensorField_fourier(i,j,k,row,column)
           enddo; enddo; enddo
           do i = 0_pInt, res(1)/2_pInt-2_pInt !unpack fft data for conj complex symmetric part. can be directly used in calculation of cstress_field
            m = 1_pInt
            do k = 1_pInt, res(3)
              n = 1_pInt
              do j = 1_pInt, res(2)
                scalarField_fourier(res(1)-i,j,k) = conjg(scalarField_fourier(2+i,n,m))
                if(n == 1_pInt) n = res(2) + 1_pInt
                n = n-1_pInt
             enddo
             if(m == 1_pInt) m = res(3) + 1_pInt
             m = m -1_pInt
           enddo; enddo
         endif
!--------------------------------------------------------------------------------------------------
! doing the inverse FT
         call fftw_execute_dft_c2r(plan_correction,tensorField_fourier,tensorField_real)            ! back transform of fluct deformation gradient

!--------------------------------------------------------------------------------------------------
! comparing 1 and 3x3 inverse FT results
         if (debugFFTW) then
           print '(a,i1,1x,i1)', 'checking iFT results of compontent ', row, column
           call fftw_execute_dft(plan_scalarField_back,scalarField_fourier,scalarField_real)
           print '(a,es10.4)', 'max iFT relative error ',&
               maxval((real(scalarField_real(1:res(1),1:res(2),1:res(3)))-&
                       tensorField_real(1:res(1),1:res(2),1:res(3),row,column))/&
                       real(scalarField_real(1:res(1),1:res(2),1:res(3))))
         endif

!--------------------------------------------------------------------------------------------------
! calculate some additional output
         if(debugGeneral) then
           maxCorrectionSkew = 0.0_pReal
           maxCorrectionSym  = 0.0_pReal
           temp33_Real = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             maxCorrectionSym  = max(maxCorrectionSym,&
                                     maxval(math_symmetric33(tensorField_real(i,j,k,1:3,1:3))))
             maxCorrectionSkew = max(maxCorrectionSkew,&
                                     maxval(math_skew33(tensorField_real(i,j,k,1:3,1:3))))
             temp33_Real = temp33_Real + tensorField_real(i,j,k,1:3,1:3)
           enddo; enddo; enddo
           print '(a,1x,es10.4)'       , 'max symmetrix correction of deformation:',&
                                         maxCorrectionSym*wgt
           print '(a,1x,es10.4)'       , 'max skew      correction of deformation:',&
                                         maxCorrectionSkew*wgt
           print '(a,1x,es10.4)'       , 'max sym/skew of avg correction:         ',&
                                         maxval(math_symmetric33(temp33_real))/&
                                         maxval(math_skew33(temp33_real))
         endif

!--------------------------------------------------------------------------------------------------
! updated deformation gradient
         defgrad = defgrad - tensorField_real(1:res(1),1:res(2),1:res(3),1:3,1:3)*wgt               ! F(x)^(n+1) = F(x)^(n) + correction;  *wgt: correcting for missing normalization

!--------------------------------------------------------------------------------------------------
! updated deformation gradient in case of fluctuation field algorithm
         if (.not.simplified_algorithm) then
           defgrad = tensorField_real(1:res(1),1:res(2),1:res(3),1:3,1:3) * wgt
         endif
                                                    
!--------------------------------------------------------------------------------------------------
! calculate bounds of det(F) and report
         if(debugGeneral) then
           defgradDetMax = -huge(1.0_pReal)
           defgradDetMin = +huge(1.0_pReal)
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             defgradDet = math_det33(defgrad(i,j,k,1:3,1:3))
             defgradDetMax = max(defgradDetMax,defgradDet)
             defgradDetMin = min(defgradDetMin,defgradDet) 
           enddo; enddo; enddo

           print '(a,1x,es10.4)'       , 'max determinant of deformation:', defgradDetMax
           print '(a,1x,es10.4)'       , 'min determinant of deformation:', defgradDetMin
         endif
         
       enddo    ! end looping when convergency is achieved 
           
       print '(a)', ''
       print '(a)', '=================================================================='
       if(err_div > err_div_tol .or. err_stress > err_stress_tol) then
         print '(A,I5.5,A)', 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       else
         convergedCounter = convergedCounter + 1_pInt
         print '(A,I5.5,A)', 'increment ', totalIncsCounter, ' converged'
       endif

       if (mod(totalIncsCounter -1_pInt,bc(loadcase)%outputfrequency) == 0_pInt) then                ! at output frequency
         print '(a)', ''
         print '(a)', '... writing results to file ......................................'
         write(538)  materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)        ! write result to file
       endif
       
       if( bc(loadcase)%restartFrequency > 0_pInt .and. &
                      mod(inc - 1_pInt,bc(loadcase)%restartFrequency) == 0_pInt) then                ! at frequency of writing restart information set restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
         restartWrite = .true.
         print '(A)', 'writing converged results for restart'
         if(IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(defgrad))) then               ! writing deformation gradient field to file
           write (777,rec=1) defgrad
           close (777)
         endif
         restartInc=totalIncsCounter
       endif 
       
       if (update_gamma) then
         print*, 'update c0_reference '
         c0_reference = c_current*wgt
         ! s0_reference = math_Plain99to3333(s0_99)
 
 !c0_99 = math_Plain3333to99(c0_reference)
 ! call math_invert(9_pInt, s0_99, c0_99, i, errmatinv)                                        ! invert reduced stiffness
 ! if(errmatinv) call IO_error(error_ID=800_pInt)
 ! print*, (c0_reference - math_Plain99to3333(c0_99))/c0_reference
! pause
       endif
       
     endif ! end calculation/forwarding
   enddo  ! end looping over incs in current loadcase
   deallocate(c_reduced)
   deallocate(s_reduced)
   enddo    ! end looping over loadcases
   print '(a)', ''
   print '(a)', '##################################################################'
   print '(i6.6,a,i6.6,a)', notConvergedCounter, ' out of ', &
                            notConvergedCounter + convergedCounter, ' increments did not converge!'
 close(538)
 call fftw_destroy_plan(plan_stress); call fftw_destroy_plan(plan_correction)
 if (debugDivergence) call fftw_destroy_plan(plan_divergence)
 if (debugFFTW) then
   call fftw_destroy_plan(plan_scalarField_forth)
   call fftw_destroy_plan(plan_scalarField_back)
 endif
 stop 0
end program DAMASK_spectral

!********************************************************************
! quit subroutine to satisfy IO_error
!
!********************************************************************
subroutine quit(id)
 use prec
 implicit none

 integer(pInt) id

 stop
end subroutine
