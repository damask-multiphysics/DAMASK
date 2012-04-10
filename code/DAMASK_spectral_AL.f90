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
program DAMASK_spectral_AL

 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran >4.6 at the moment)
 use DAMASK_interface
 use prec,             only: pInt, pReal, DAMASK_NaN
 use IO 
 use debug,            only: debug_spectral, &
                             debug_levelBasic, &
                             debug_spectralRestart, &
                             debug_spectralFFTW
 use math
 use CPFEM,            only: CPFEM_general, CPFEM_initAll
 use FEsolving,        only: restartWrite, restartInc
 use numerics,         only: err_div_tol, err_stress_tolrel, rotation_tol, itmax, itmin, &
                             memory_efficient, update_gamma, DAMASK_NumThreadsInt, &
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

!--------------------------------------------------------------------------------------------------
! variable storing information from load case file
 type bc_type
   real(pReal), dimension (3,3) :: deformation            = 0.0_pReal, &                            ! applied velocity gradient or time derivative of deformation gradient
                                   P                      = 0.0_pReal, &                            ! stress BC (if applicable)
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
 
 type(bc_type), allocatable, dimension(:) ::  bc

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
 real(pReal), dimension(3,3) ::           P_star_av = 0.0_pReal, &
                                          F_aim = math_I3, F_aim_lastInc = math_I3, lambda_av, &
                                          mask_stress, mask_defgrad, deltaF, F_star_av, &
                                          F_aim_lab                                                 ! quantities rotated to other coordinate system
 real(pReal), dimension(3,3,3,3) ::       C_inc0, C=0.0_pReal, S_lastInc, C_lastInc           ! stiffness and compliance
 real(pReal), dimension(6) ::             sigma                                                     ! cauchy stress
 real(pReal), dimension(6,6) ::           dsde                                                      ! small strain stiffness
 real(pReal), dimension(9,9) ::           s_prev99, c_prev99                                        ! compliance and stiffness in matrix notation 
 real(pReal), dimension(:,:), allocatable ::  s_reduced, c_reduced                                  ! reduced compliance and stiffness (only for stress BC)
 integer(pInt) ::                         size_reduced = 0_pInt                                     ! number of stress BCs

!--------------------------------------------------------------------------------------------------
! pointwise data 
 type(C_PTR) :: tensorField                                                                         ! fields in real an fourier space
 real(pReal),    dimension(:,:,:,:,:), pointer :: lambda_real, F_real                               ! fields in real space (pointer)
 complex(pReal), dimension(:,:,:,:,:), pointer :: lambda_fourier, F_fourier                         ! fields in fourier space (pointer)
 real(pReal),    dimension(:,:,:,:,:), allocatable ::  F_lastInc, F_star, lambda, P, F_star_lastIter
 real(pReal),    dimension(:,:,:,:,:,:,:), allocatable ::  dPdF
 real(pReal),    dimension(:,:,:,:),   allocatable ::  coordinates
 real(pReal),    dimension(:,:,:),     allocatable ::  temperature

!--------------------------------------------------------------------------------------------------
! variables storing information for spectral method and FFTW
 type(C_PTR) ::  plan_correction, plan_lambda                                                       ! plans for fftw
 real(pReal), dimension(3,3) ::                         xiDyad                                      ! product of wave vectors
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                   ! gamma operator (field) for spectral method
 real(pReal), dimension(:,:,:,:), allocatable ::        xi                                          ! wave vector field for divergence and for gamma operator
 integer(pInt), dimension(3) ::                         k_s

!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc = 1.0_pReal, timeinc_old = 0.0_pReal   ! elapsed time, begin of interval, time interval 
 real(pReal) :: guessmode, err_stress, err_stress_tol, err_f, err_p, err_crit, &
                err_f_point, err_p_point, pstress_av_L2, err_div_rms, err_div     
 real(pReal), dimension(3,3), parameter ::  ones = 1.0_pReal, zeroes = 0.0_pReal
 complex(pReal), dimension(3,3) ::          temp33_Complex
 real(pReal),    dimension(3,3) ::          temp33_Real
 integer(pInt) :: i, j, k, l, u, v, w, errorID = 0_pInt, ierr
 integer(pInt) :: N_Loadcases, loadcase, inc, iter, ielem, CPFEM_mode, guesses, guessmax=10_pInt,&
                  totalIncsCounter = 0_pInt,notConvergedCounter = 0_pInt, convergedCounter = 0_pInt
 logical          :: errmatinv, callCPFEM
 character(len=6) :: loadcase_string

!--------------------------------------------------------------------------------------------------
!variables controlling debugging
 logical :: debugGeneral, debugDivergence, debugRestart, debugFFTW

!##################################################################################################
! reading of information from load case file and geometry file
!##################################################################################################
 !$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                  ! set number of threads for parallel execution set by DAMASK_NUM_THREADS
 open (6, encoding='UTF-8')  
 call DAMASK_interface_init

 print '(a)', ''
 print '(a)', ' <<<+-  DAMASK_spectral_AL init  -+>>>'
 print '(a)', ' $Id$'
#include "compilation_info.f90"
 print '(a,a)', ' Working Directory:    ',trim(getSolverWorkingDirectoryName())
 print '(a,a)', ' Solver Job Name:      ',trim(getSolverJobName())
 print '(a)', ''
!--------------------------------------------------------------------------------------------------
! reading the load case file and allocate data structure containing load cases
 path = getLoadcaseName()
 call IO_open_file(myUnit,path)
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
   call IO_error(error_ID=837_pInt,ext_msg = trim(path))                                            ! error message for incomplete loadcase
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
                                                          IO_floatValue(line,positions,j+k)         ! assign values for the bc(loadcase)%P matrix
         enddo
         bc(loadcase)%maskStress = transpose(reshape(bc(loadcase)%maskStressVector,[ 3,3]))
         bc(loadcase)%P = math_plain9to33(temp_valueVector)
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
         u = 0_pInt                                                                                 ! assuming values given in radians
         v = 1_pInt                                                                                 ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,positions,j+1_pInt)))
           case('deg','degree')
             u = 1_pInt                                                                             ! for conversion from degree to radian           
           case('rad','radian') 
           case default               
             v = 0_pInt                                                                             ! immediately reading in angles, assuming radians
         end select
         forall(k = 1_pInt:3_pInt)  temp33_Real(k,1) = &
                                        IO_floatValue(line,positions,j+v+k) * real(u,pReal) * inRad
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
 if (update_gamma .and. .not. memory_efficient) call IO_error(error_ID = 847_pInt)

!--------------------------------------------------------------------------------------------------
! read header of geom file to get size information. complete geom file is intepretated by mesh.f90
 path = getModelName()
 call IO_open_file(myUnit,trim(path)//InputFileExtension)
 rewind(myUnit)
 read(myUnit,'(a1024)') line
 positions = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,positions,2_pInt))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,positions,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=842_pInt)
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
                                                                call IO_error(error_ID = 845_pInt)
 if (any(geomdim<=0.0_pReal)) call IO_error(error_ID = 802_pInt)
 if(mod(res(1),2_pInt)/=0_pInt .or.&
    mod(res(2),2_pInt)/=0_pInt .or.&
   (mod(res(3),2_pInt)/=0_pInt .and. res(3)/= 1_pInt)) call IO_error(error_ID = 803_pInt)

!--------------------------------------------------------------------------------------------------
! variables derived from resolution
 res1_red = res(1)/2_pInt + 1_pInt                                                                  ! size of complex array in first dimension (c2r, r2c)
 Npoints = res(1)*res(2)*res(3)
 wgt = 1.0_pReal/real(Npoints, pReal)

!--------------------------------------------------------------------------------------------------
! output of geometry
 print '(a)',  ''
 print '(a)',   '#############################################################'
 print '(a)',   'DAMASK spectral_AL:'
 print '(a)',   'The AL spectral method boundary value problem solver for'
 print '(a)',   'the Duesseldorf Advanced Material Simulation Kit'
 print '(a)',   '#############################################################'
 print '(a,a)', 'geometry file:        ',trim(path)//'.geom'
 print '(a)',   '============================================================='
 print '(a,3(i12  ))','resolution a b c:', res
 print '(a,3(f12.5))','dimension  x y z:', geomdim
 print '(a,i5)','homogenization:       ',homog
 print '(a)',   '#############################################################'
 print '(a,a)', 'loadcase file:        ',trim(getLoadcaseName())

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
 bc(1)%followFormerTrajectory = .false.                                                             ! cannot guess along trajectory for first inc of first loadcase

 do loadcase = 1_pInt, N_Loadcases
   write (loadcase_string, '(i6)' ) loadcase

   print '(a)', '============================================================='
   print '(a,i6)', 'loadcase:            ', loadcase

   if (.not. bc(loadcase)%followFormerTrajectory) print '(a)', 'drop guessing along trajectory'
   if (bc(loadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .true.) .and. &
           any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .false.)) errorID = 832_pInt               ! each row should be either fully or not at all defined
     enddo
     print '(a)','velocity gradient:'
   else
     print '(a)','deformation gradient rate:'
   endif
   write (*,'(3(3(f12.7,1x)/))',advance='no') merge(math_transpose33(bc(loadcase)%deformation),&
                  reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskDeformation))
   write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'stress / GPa:',&
        1e-9_pReal*merge(math_transpose33(bc(loadcase)%P),&
                         reshape(spread(DAMASK_NaN,1,9),[ 3,3]),transpose(bc(loadcase)%maskStress))
   if (any(bc(loadcase)%rotation /= math_I3)) &
     write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') ' rotation of loadframe:',&
                                                          math_transpose33(bc(loadcase)%rotation)
   print '(a,f12.6)','temperature:',bc(loadcase)%temperature
   print '(a,f12.6)','time:       ',bc(loadcase)%time
   print '(a,i5)'   ,'increments: ',bc(loadcase)%incs
   print '(a,i5)','output  frequency:  ',bc(loadcase)%outputfrequency
   print '(a,i5)','restart frequency:  ',bc(loadcase)%restartfrequency
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
 
!--------------------------------------------------------------------------------------------------
! debugging parameters
 debugRestart    = iand(debug_spectral,debug_spectralRestart)    > 0_pInt
 debugFFTW       = iand(debug_spectral,debug_spectralFFTW)       > 0_pInt
 debugGeneral = .true.
 
!##################################################################################################
! initialization 
!##################################################################################################

!--------------------------------------------------------------------------------------------------
! allocate more memory 
 allocate (P          (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (dPdF       (  res(1),  res(2),res(3),3,3,3,3),  source = 0.0_pReal)
 allocate (F_star     (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (F_star_lastIter     (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (F_lastInc  (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (lambda     (  res(1),  res(2),res(3),3,3),  source = 0.0_pReal)
 allocate (xi         (3,res1_red,res(2),res(3)),      source = 0.0_pReal)
 allocate (coordinates(  res(1),  res(2),res(3),3),    source = 0.0_pReal)
 allocate (temperature(  res(1),  res(2),res(3)),      source = bc(1)%temperature)                  ! start out isothermally
 tensorField = fftw_alloc_complex(int(res1_red*res(2)*res(3)*9_pInt,C_SIZE_T))                      ! allocate continous data using a C function, C_SIZE_T is of type integer(8)
 call c_f_pointer(tensorField, lambda_real,    [ res(1)+2_pInt,res(2),res(3),3,3])             ! place a pointer for the real representation
 call c_f_pointer(tensorField, F_real,         [ res(1)+2_pInt,res(2),res(3),3,3])             ! place a pointer for the real representation
 call c_f_pointer(tensorField, lambda_fourier, [ res1_red,     res(2),res(3),3,3])             ! place a pointer for the complex representation
 call c_f_pointer(tensorField, F_fourier,      [ res1_red,     res(2),res(3),3,3])             ! place a pointer for the complex representation

!--------------------------------------------------------------------------------------------------
! general initialization of fftw (see manual on fftw.org for more details)
 if (pReal /= C_DOUBLE .or. pInt /= C_INT) call IO_error(error_ID=808_pInt)                         ! check for correct precision in C
#ifdef _OPENMP
    if(DAMASK_NumThreadsInt > 0_pInt) then
      ierr = fftw_init_threads()
      if (ierr == 0_pInt) call IO_error(error_ID = 809_pInt)
      call fftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
    endif
#endif
 call fftw_set_timelimit(fftw_timelimit)                                                            ! set timelimit for plan creation

!--------------------------------------------------------------------------------------------------
! creating plans
 plan_lambda    =  fftw_plan_many_dft_r2c(3,[ res(3),res(2) ,res(1)],9,&                             ! dimensions , length in each dimension in reversed order
                                lambda_real,[ res(3),res(2) ,res(1)+2_pInt],&                        ! input data , physical length in each dimension in reversed order
                                        1,  res(3)*res(2)*(res(1)+2_pInt),&                        ! striding   , product of physical lenght in the 3 dimensions
                             lambda_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,fftw_planner_flag)   

 plan_correction  =  fftw_plan_many_dft_c2r(3,[ res(3),res(2) ,res(1)],9,&
                                    F_fourier,[ res(3),res(2) ,res1_red],&
                                         1,  res(3)*res(2)* res1_red,&
                                       F_real,[ res(3),res(2) ,res(1)+2_pInt],&
                                         1,  res(3)*res(2)*(res(1)+2_pInt),fftw_planner_flag)
 if (debugGeneral) print '(a)' , 'FFTW initialized'
 
!--------------------------------------------------------------------------------------------------
! calculation of discrete angular frequencies, ordered as in FFTW (wrap around) and remove the given highest frequencies
 do k = 1_pInt, res(3)
   k_s(3) = k - 1_pInt
   if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
     do j = 1_pInt, res(2)
       k_s(2) = j - 1_pInt
       if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
         do i = 1_pInt, res1_red
           k_s(1) = i - 1_pInt
           xi(1:3,i,j,k) = real(k_s, pReal)/geomdim
 enddo; enddo; enddo
 
!--------------------------------------------------------------------------------------------------
! calculate the gamma operator
 if(memory_efficient) then                                                                          ! allocate just single fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3), source = 0.0_pReal)
 else                                                                                               ! precalculation of gamma_hat field
   allocate (gamma_hat(res1_red ,res(2),res(3),3,3,3,3), source =0.0_pReal)
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res1_red
     if(any([i,j,k] /= 1_pInt)) then                                                                ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
       forall(l = 1_pInt:3_pInt, u = 1_pInt:3_pInt) &
         xiDyad(l,u) = xi(l, i,j,k)*xi(u, i,j,k)
       forall(l = 1_pInt:3_pInt, u = 1_pInt:3_pInt) &
         temp33_Real(l,u) = sum(C_inc0(l,1:3,u,1:3)*xiDyad)
       temp33_Real = math_inv33(temp33_Real)
       forall(l=1_pInt:3_pInt, u=1_pInt:3_pInt, v=1_pInt:3_pInt, w=1_pInt:3_pInt)&
         gamma_hat(i,j,k, l,u,v,w) =  temp33_Real(l,v)*xiDyad(u,w)
     endif  
   enddo; enddo; enddo
   gamma_hat(1,1,1, 1:3,1:3,1:3,1:3) = 0.0_pReal                                                    ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
 endif

!--------------------------------------------------------------------------------------------------
! init fields to no deformation
 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt 
   F_real(i,j,k,1:3,1:3) = math_I3; F_lastInc(i,j,k,1:3,1:3) = math_I3
   coordinates(i,j,k,1:3) = geomdim/real(res * [i,j,k], pReal) - geomdim/real(2_pInt*res,pReal)
   call CPFEM_general(3_pInt,coordinates(i,j,k,1:3),math_I3,math_I3,temperature(i,j,k),&
                      0.0_pReal,ielem,1_pInt,sigma,dsde,temp33_Real ,dPdF(i,j,k,1:3,1:3,1:3,1:3))
 enddo; enddo; enddo

 ielem = 0_pInt
 do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
   ielem = ielem + 1_pInt 
   call CPFEM_general(2_pInt,coordinates(i,j,k,1:3),math_I3,math_I3,temperature(i,j,k),&
                      0.0_pReal,ielem,1_pInt,sigma,dsde,temp33_Real ,dPdF(i,j,k,1:3,1:3,1:3,1:3))
   C = C + dPdF(i,j,k,1:3,1:3,1:3,1:3)
 enddo; enddo; enddo
 C_inc0 = C * wgt                                                                     ! linear reference material stiffness

!--------------------------------------------------------------------------------------------------
! possible restore deformation gradient from saved state
 if (restartInc > 1_pInt) then                                                                      ! using old values from file                                                      
   if (debugRestart) print '(a,i6,a)' , 'Reading values of increment ',&
                                             restartInc - 1_pInt,' from file' 
   call IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',&
                                                trim(getSolverJobName()),size(F_star))
   read (777,rec=1) F_star
   close (777)
   F_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = F_star
   F_lastInc = F_star
   F_aim = 0.0_pReal
   do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
     F_aim = F_aim + F_real(i,j,k,1:3,1:3)                                               ! calculating old average deformation
   enddo; enddo; enddo
   F_aim = F_aim * wgt
   F_aim_lastInc = F_aim
 endif

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
 write(538) 'frequencies', bc(1:N_Loadcases)%outputfrequency                                        ! one entry per loadcase
 write(538) 'times', bc(1:N_Loadcases)%time                                                         ! one entry per loadcase
 write(538) 'logscales',  bc(1:N_Loadcases)%logscale         
 write(538) 'increments', bc(1:N_Loadcases)%incs                                                    ! one entry per loadcase
 write(538) 'startingIncrement', restartInc - 1_pInt                                                ! start with writing out the previous inc
 write(538) 'eoh'                                                                                   ! end of header
 write(538) materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)                ! initial (non-deformed or read-in) results
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
   size_reduced = int(count(bc(loadcase)%maskStressVector), pInt)
   allocate (c_reduced(size_reduced,size_reduced), source =0.0_pReal)
   allocate (s_reduced(size_reduced,size_reduced), source =0.0_pReal)

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
         deltaF = timeinc * mask_defgrad * math_mul33x33(bc(loadcase)%deformation, F_aim)
       else                                                                                         ! deltaF = fDot *timeinc where applicable
         deltaF = timeinc * mask_defgrad * bc(loadcase)%deformation
       endif

!--------------------------------------------------------------------------------------------------
! coordinates at beginning of inc
       !call deformed_fft(res,geomdim,1.0_pReal,F_real(1:res(1),1:res(2),1:res(3),1:3,1:3),coordinates)! calculate current coordinates

!--------------------------------------------------------------------------------------------------
! winding forward of deformation aim in loadcase system
       temp33_Real = F_aim                                            
       F_aim = F_aim &                                                                         
                  + guessmode * mask_stress * (F_aim - F_aim_lastInc)*timeinc/timeinc_old &      
                  + deltaF
       F_aim_lastInc = temp33_Real
       F_star_av  = F_aim

!--------------------------------------------------------------------------------------------------
! update local deformation gradient
       deltaF = math_rotate_backward33(deltaF,bc(loadcase)%rotation)
       do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
         temp33_Real = F_real(i,j,k,1:3,1:3)
         F_real(i,j,k,1:3,1:3) = F_real(i,j,k,1:3,1:3) &                                          ! decide if guessing along former trajectory or apply homogeneous addon
                                + guessmode * (F_real(i,j,k,1:3,1:3) - F_lastInc(i,j,k,1:3,1:3))& ! guessing... 
                                            *timeinc/timeinc_old &
                                + (1.0_pReal-guessmode) * deltaF                                    ! if not guessing, use prescribed average deformation where applicable
         F_lastInc(i,j,k,1:3,1:3) = temp33_Real 
       enddo; enddo; enddo

!--------------------------------------------------------------------------------------------------
! Initialize / Update lambda to useful value
       P_star_av = P_star_av + math_mul3333xx33(C*wgt, F_aim-F_aim_lastInc)
       lambda_av = 0.0_pReal
       do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
          lambda(i,j,k,1:3,1:3) = P(i,j,k,1:3,1:3) +  math_mul3333xx33(dPdF(i,j,k,1:3,1:3,1:3,1:3), &
          F_real(i,j,k,1:3,1:3)-F_lastInc(i,j,k,1:3,1:3))
          lambda_av = lambda_av + lambda(i,j,k,1:3,1:3)
       enddo; enddo; enddo
       lambda_av=lambda_av*wgt

!--------------------------------------------------------------------------------------------------
!Initialize pointwise data for AL scheme: ToDo: good choice?
       F_star(1:res(1),1:res(2),1:res(3),1:3,1:3) = F_real(1:res(1),1:res(2),1:res(3),1:3,1:3)

!--------------------------------------------------------------------------------------------------
! calculate reduced compliance
       if(size_reduced > 0_pInt) then                                                               ! calculate compliance in case stress BC is applied
         C_lastInc = math_rotate_forward3333(C*wgt,bc(loadcase)%rotation)                      ! calculate stiffness from former inc
         c_prev99 = math_Plain3333to99(C_lastInc)
         k = 0_pInt                                                                                 ! build reduced stiffness
         do v = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(v)) then
             k = k + 1_pInt
             j = 0_pInt
             do u = 1_pInt,9_pInt
               if(bc(loadcase)%maskStressVector(u)) then
                 j = j + 1_pInt
                 c_reduced(k,j) = c_prev99(v,u)
         endif; enddo; endif; enddo
         call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)                         ! invert reduced stiffness
         if(errmatinv) call IO_error(error_ID=400_pInt)
         s_prev99 = 0.0_pReal                                                                       ! build full compliance
         k = 0_pInt
         do v = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(v)) then
             k = k + 1_pInt
             j = 0_pInt
             do u = 1_pInt,9_pInt
             if(bc(loadcase)%maskStressVector(u)) then
                   j = j + 1_pInt
                   s_prev99(v,u) = s_reduced(k,j)
         endif; enddo; endif; enddo
         S_lastInc = (math_Plain99to3333(s_prev99))
       endif

!--------------------------------------------------------------------------------------------------
! report begin of new increment
       print '(a)', '##################################################################'
       print '(A,I5.5,A,es12.5)', 'Increment ', totalIncsCounter, ' Time ',time
       
       guessmode = 1.0_pReal                                                                        ! keep guessing along former trajectory during same loadcase
       CPFEM_mode = 1_pInt                                                                          ! winding forward
       iter = 0_pInt
       err_crit = huge(1.0_pReal)                                                                   ! go into loop 
       callCPFEM=.true.
       guessmax = 2
       guesses = 0

!##################################################################################################
! convergence loop (looping over iterations)
!##################################################################################################
       do while((iter < itmax .and. (err_crit > 1.0_pReal)) .or. iter < itmin)
         iter = iter + 1_pInt

!--------------------------------------------------------------------------------------------------
! report begin of new iteration
         print '(a)', ''
         print '(a)', '=================================================================='
         print '(5(a,i6.6))', 'Loadcase ',loadcase,' Increment ',inc,'/',bc(loadcase)%incs,&
                                                                  ' @ Iteration ',iter,'/',itmax

!--------------------------------------------------------------------------------------------------
! stress BC handling
         if(size_reduced > 0_pInt) then                                                              ! calculate stress BC if applied
           err_stress = maxval(abs(mask_stress * (lambda_av - bc(loadcase)%P)))                      ! maximum deviaton (tensor norm not applicable)
           F_aim = F_aim  + math_mul3333xx33(S_lastInc,bc(loadcase)%P- lambda_av)
           err_stress_tol = maxval(abs(lambda_av)) * err_stress_tolrel                              ! don't use any tensor norm because the comparison should be coherent
         else
           err_stress_tol = + huge(1.0_pReal)
         endif
         F_aim_lab = math_rotate_backward33(F_aim,bc(loadcase)%rotation)
         write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'F aim  =',&
                                                math_transpose33(F_aim)

!--------------------------------------------------------------------------------------------------
! doing Fourier transform
         print '(a)', '... spectral method ...............................................'
         lambda_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = lambda(1:res(1),1:res(2),1:res(3),1:3,1:3)
         call fftw_execute_dft_r2c(plan_lambda,lambda_real,lambda_fourier)
         lambda_fourier(  res1_red,1:res(2) ,             1:res(3)              ,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         lambda_fourier(1:res1_red,  res(2)/2_pInt+1_pInt,1:res(3)              ,1:3,1:3)& 
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
         if(res(3)>1_pInt) &
          lambda_fourier(1:res1_red,1:res(2),                res(3)/2_pInt+1_pInt,1:3,1:3)&
                                                             = cmplx(0.0_pReal,0.0_pReal,pReal)
!--------------------------------------------------------------------------------------------------
! calculating RMS divergence criterion in Fourier space
         pstress_av_L2 = sqrt(maxval(math_eigenvalues33(math_mul33x33(lambda_av,&                    ! L_2 norm of average stress (http://mathworld.wolfram.com/SpectralNorm.html)
                                                     math_transpose33(lambda_av)))))
         err_div_RMS = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2)
           do i = 2_pInt, res1_red -1_pInt                                                          ! Has somewhere a conj. complex counterpart. Therefore count it twice.
             err_div_RMS = err_div_RMS &
                   + 2.0_pReal*(sum (real(math_mul33x3_complex(lambda_fourier(i,j,k,1:3,1:3),&           ! (sqrt(real(a)**2 + aimag(a)**2))**2 = real(a)**2 + aimag(a)**2. do not take square root and square again
                                                   xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal)&            ! --> sum squared L_2 norm of vector 
                               +sum(aimag(math_mul33x3_complex(lambda_fourier(i,j,k,1:3,1:3),& 
                                                                  xi(1:3,i,j,k))*TWOPIIMG)**2.0_pReal))
           enddo
           err_div_RMS = err_div_RMS &                                                              ! Those two layers (DC and Nyquist) do not have a conjugate complex counterpart
                         + sum( real(math_mul33x3_complex(lambda_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(lambda_fourier(1       ,j,k,1:3,1:3),&
                                                             xi(1:3,1       ,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum( real(math_mul33x3_complex(lambda_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)&
                         + sum(aimag(math_mul33x3_complex(lambda_fourier(res1_red,j,k,1:3,1:3),&
                                                             xi(1:3,res1_red,j,k))*TWOPIIMG)**2.0_pReal)
         enddo; enddo

         err_div_RMS = sqrt(err_div_RMS)*wgt
       !  if (err_div < err_div_RMS/pstress_av_L2 .and. guessmax<0) then
       !    print*, 'increasing div, stopping calc'
       !    iter = huge(1_pInt)
       !  endif
         err_div = err_div_RMS/pstress_av_L2
!--------------------------------------------------------------------------------------------------
! using gamma operator to update F 
         if(memory_efficient) then                                                                  ! memory saving version, on-the-fly calculation of gamma_hat
           do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res1_red
               if(any([i,j,k] /= 1_pInt)) then                                                      ! singular point at xi=(0.0,0.0,0.0) i.e. i=j=k=1       
                 forall(l = 1_pInt:3_pInt, u = 1_pInt:3_pInt) &
                   xiDyad(l,u) = xi(l, i,j,k)*xi(u, i,j,k)
                 forall(l = 1_pInt:3_pInt, u = 1_pInt:3_pInt) &
                   temp33_Real(l,u) = sum(C_inc0(l,1:3,u,1:3)*xiDyad)
                 temp33_Real = math_inv33(temp33_Real)
                 forall(l=1_pInt:3_pInt, u=1_pInt:3_pInt, v=1_pInt:3_pInt, w=1_pInt:3_pInt)&
                   gamma_hat(1,1,1, l,u,v,w) =  temp33_Real(l,v)*xiDyad(u,w)
                 forall(l = 1_pInt:3_pInt, u = 1_pInt:3_pInt) &
                   temp33_Complex(l,u) = sum(gamma_hat(1,1,1, l,u, 1:3,1:3) *&
                                                        lambda_fourier(i,j,k,1:3,1:3))
                 F_fourier(i,j,k,1:3,1:3) = - temp33_Complex 
             endif             
           enddo; enddo; enddo
         else                                                                                       ! use precalculated gamma-operator
           do k = 1_pInt, res(3);  do j = 1_pInt, res(2);  do i = 1_pInt,res1_red
             forall( u = 1_pInt:3_pInt, v = 1_pInt:3_pInt) &
               temp33_Complex(u,v) = sum(gamma_hat(i,j,k, u,v, 1:3,1:3) *&
                                                         lambda_fourier(i,j,k,1:3,1:3))
             F_fourier(i,j,k, 1:3,1:3) = - temp33_Complex
           enddo; enddo; enddo
         endif      
         F_fourier(1,1,1,1:3,1:3) = cmplx((F_aim_lab - F_star_av)*real(Npoints,pReal),0.0_pReal,pReal)

!--------------------------------------------------------------------------------------------------
! doing inverse Fourier transform
         call fftw_execute_dft_c2r(plan_correction,F_fourier,F_real)                                ! back transform of fluct deformation gradient
         F_real(1:res(1),1:res(2),1:res(3),1:3,1:3) = F_real(1:res(1),1:res(2),1:res(3),1:3,1:3) * wgt + &
                                                      F_star(1:res(1),1:res(2),1:res(3),1:3,1:3)

!--------------------------------------------------------------------------------------------------
!
         if(callCPFEM) then
           print '(a)', '... calling CPFEM to update P(F*) and F*.........................'
           F_star_lastIter = F_star
           ielem = 0_pInt
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             ielem = ielem + 1_pInt
             call CPFEM_general(3_pInt,&                                                              ! collect cycle
                                coordinates(i,j,k,1:3), F_lastInc(i,j,k,1:3,1:3),&
                                F_star(i,j,k,1:3,1:3),temperature(i,j,k),timeinc,ielem,1_pInt,&
                                sigma,dsde, P(i,j,k,1:3,1:3), dPdF(i,j,k,1:3,1:3,1:3,1:3))
           enddo; enddo; enddo
           ielem = 0_pInt
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             ielem = ielem + 1_pInt
             call CPFEM_general(CPFEM_mode,&
                                coordinates(i,j,k,1:3), F_lastInc(i,j,k,1:3,1:3),&
                                F_star(i,j,k,1:3,1:3),temperature(i,j,k),timeinc,ielem,1_pInt,&
                                sigma,dsde, P(i,j,k,1:3,1:3), dPdF(i,j,k,1:3,1:3,1:3,1:3))
             CPFEM_mode = 2_pInt                                                                  ! winding forward
             temp33_Real = lambda(i,j,k,1:3,1:3) - P(i,j,k,1:3,1:3) &
                         + math_mul3333xx33(C_inc0,F_real(i,j,k,1:3,1:3)- F_star(i,j,k,1:3,1:3))

             F_star(i,j,k,1:3,1:3) =  F_star(i,j,k,1:3,1:3) +  math_mul3333xx33(math_invSym3333(&
                                     C_inc0 + dPdF(i,j,k,1:3,1:3,1:3,1:3)), temp33_Real)
           enddo; enddo; enddo
         else
           guesses = guesses +1_pInt
           print*, '... linear approximation for P(F*) and F* ', guesses, ' of ', guessmax
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             temp33_Real = lambda(i,j,k,1:3,1:3) - (P(i,j,k,1:3,1:3)  + math_mul3333xx33(dPdF(i,j,k,1:3,1:3,1:3,1:3),&
                           F_star(i,j,k,1:3,1:3) -F_star_lastIter(i,j,k,1:3,1:3)))&
                         + math_mul3333xx33(C_inc0,F_real(i,j,k,1:3,1:3)- F_star(i,j,k,1:3,1:3))

             F_star(i,j,k,1:3,1:3) =  F_star(i,j,k,1:3,1:3) +  math_mul3333xx33(math_invSym3333(&
                                     C_inc0 + dPdF(i,j,k,1:3,1:3,1:3,1:3)), temp33_Real)
           enddo; enddo; enddo
         endif

         print '(a)', '... update  λ..........................'

         err_f = 0.0_pReal
         err_f_point = 0.0_pReal
         err_p = 0.0_pReal
         err_p_point = 0.0_pReal

         F_star_av = 0.0_pReal
         P_star_av = 0.0_pReal
         lambda_av = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           lambda(i,j,k,1:3,1:3) = lambda(i,j,k,1:3,1:3) + math_mul3333xx33(C_inc0,F_real(i,j,k,1:3,1:3) &
                                                                                 - F_star(i,j,k,1:3,1:3))
           F_star_av = F_star_av + F_star(i,j,k,1:3,1:3)
           lambda_av = lambda_av + lambda(i,j,k,1:3,1:3)
           P_star_av = P_star_av + P(i,j,k,1:3,1:3)

           temp33_real = F_star(i,j,k,1:3,1:3) - F_real(i,j,k,1:3,1:3)
           err_f_point = max(err_f_point, maxval(abs(temp33_real)))
           err_f = max(err_f, sqrt(math_mul33xx33(temp33_real,temp33_real)))
           
           temp33_real = lambda(i,j,k,1:3,1:3) - (P(i,j,k,1:3,1:3)  + math_mul3333xx33(dPdF(i,j,k,1:3,1:3,1:3,1:3),&
                           F_star(i,j,k,1:3,1:3) -F_star_lastIter(i,j,k,1:3,1:3)))
           err_p_point = max(err_p_point, maxval(abs(temp33_real)))
           err_p = max(err_p, sqrt(math_mul33xx33(temp33_real,temp33_real)))
         enddo; enddo; enddo

         F_star_av = F_star_av *wgt
         write (*,'(a,/,3(3(f12.7,1x)/))',advance='no') 'F* =',&
                                              math_transpose33(F_star_av)
         P_star_av = P_star_av *wgt
         write (*,'(a,/,3(3(es14.7,1x)/))',advance='no') 'P(F*) / GPa =',&
                                              math_transpose33(P_star_av) /1.e6_pReal
         lambda_av = lambda_av *wgt
         write (*,'(a,/,3(3(es14.7,1x)/))',advance='no') 'λ / GPa =',&
                                              math_transpose33(lambda_av) /1.e6_pReal
                                              
         err_f = err_f/sqrt(math_mul33xx33(F_star_av,F_star_av))
         err_p = err_p/sqrt(math_mul33xx33(P_star_av,P_star_av))

         write(6,'(a,es14.7,es14.7)') 'error F', err_f/1e-4, err_f
         write(6,'(a,es14.7,es14.7)') 'error P', err_p/1e-3, err_p
         write(6,'(a,es14.7,es14.7)') 'error stress     = ',err_stress/err_stress_tol, err_stress
         write(6,'(a,es14.7,es14.7)') 'error divergence = ',err_div/err_div_tol, err_div
         write(6,*) '  ' 
         write(6,'(a,es14.7)') 'error divergence  FT  RMS = ',err_div_RMS 
         write(6,'(a,es14.7)') 'max abs err F', err_f_point
         write(6,'(a,es14.7)') 'max abs err P', err_p_point
       err_crit = max(err_p/1e-3, err_f/1e-4,err_div/err_div_tol,err_stress/err_stress_tol)
       print*, 'critical error', err_crit

       if (.not. callCPFEM) then
         if(err_crit < 1.0_pReal .or. guesses >= guessmax) callCPFEM = .true.
         err_crit =huge(1.0_pReal)
       else
         if(iter >2 .and. iter< itmax-3) callCPFEM=.false.
         guesses = 0_pInt
       endif
       
       enddo    ! end looping when convergency is achieved 
       write(6,'(a)') '  ' 
       write(6,'(a)') '=================================================================='
       if(err_crit > 1.0_pReal) then
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' NOT converged'
         notConvergedCounter = notConvergedCounter + 1_pInt
       else
         convergedCounter = convergedCounter + 1_pInt
         write(6,'(A,I5.5,A)') 'increment ', totalIncsCounter, ' converged'
       endif

       if (mod(totalIncsCounter -1_pInt,bc(loadcase)%outputfrequency) == 0_pInt) then               ! at output frequency
         write(6,'(a)') '  ' 
         write(6,'(a)') '... writing results to file ......................................'
         write(538)  materialpoint_results(1_pInt:materialpoint_sizeResults,1,1_pInt:Npoints)       ! write result to file
       endif
       
       if( bc(loadcase)%restartFrequency > 0_pInt .and. &
                      mod(inc - 1_pInt,bc(loadcase)%restartFrequency) == 0_pInt) then               ! at frequency of writing restart information set restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
         restartWrite = .true.
         write(6,*) 'writing converged results for restart'
         call IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(F_star))                  ! writing deformation gradient field to file
         write (777,rec=1) F_star
         close (777)
         restartInc=totalIncsCounter
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
 call fftw_destroy_plan(plan_lambda); call fftw_destroy_plan(plan_correction)
 call quit(1_pInt)
end program DAMASK_spectral_AL

#include "DAMASK_quit.f90"