! Copyright 2011 Max-Planck-Institut fuer Eisenforschung GmbH
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
!##############################################################
!* $Id$
!********************************************************************
! Material subroutine for BVP solution using spectral method
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
!
!********************************************************************
!     Usage:
!             - start program with DAMASK_spectral
!                       -g (--geom, --geometry) PathToGeomFile/NameOfGeom.geom
!                       -l (--load, --loadcase) PathToLoadFile/NameOfLoadFile.load
!             - PathToGeomFile will be the working directory
!             - make sure the file "material.config" exists in the working
!               directory. For further configuration use "numerics.config" and
!               "numerics.config"
!********************************************************************
program DAMASK_spectral
!********************************************************************

 use DAMASK_interface
 use prec, only: pInt, pReal
 use IO
 use debug, only: spectral_debug_verbosity
 use math
 use mesh, only: mesh_ipCenterOfGravity
 use CPFEM, only: CPFEM_general, CPFEM_initAll
 use FEsolving, only: restartWrite, restartReadSpectral, restartReadStep
 use numerics, only: err_div_tol, err_stress_tolrel , rotation_tol,&
                     itmax, memory_efficient, DAMASK_NumThreadsInt, divergence_correction, &
                     fftw_planner_flag, fftw_timelimit
 use homogenization, only: materialpoint_sizeResults, materialpoint_results
!$ use OMP_LIB                                                                        ! the openMP function library

 implicit none
! variables to read from loadcase and geom file
 real(pReal), dimension(9) ::                      temp_valueVector                   ! stores information temporarily from loadcase file
 logical, dimension(9) ::                          temp_maskVector
 integer(pInt), parameter ::                       maxNchunksLoadcase = &
                                                   (1_pInt + 9_pInt)*3_pInt + &       ! deformation, rotation, and stress
                                                   (1_pInt + 1_pInt)*5_pInt + &       ! time, (log)incs, temp, restartfrequency, and outputfrequency
                                                   1_pInt                             ! dropguessing
 integer(pInt), dimension (1_pInt + maxNchunksLoadcase*2_pInt) :: posLoadcase
 integer(pInt), parameter :: maxNchunksGeom = 7_pInt                                  ! 4 identifiers, 3 values
 integer(pInt), dimension (1_pInt + maxNchunksGeom*2_pInt) ::  posGeom
 integer(pInt) :: headerLength, N_l=0_pInt, N_t=0_pInt, N_n=0_pInt, N_Fdot=0_pInt
 integer(pInt), parameter :: myUnit = 234_pInt
 character(len=1024) :: path, line, keyword
 logical ::  gotResolution =.false., gotDimension =.false., gotHomogenization = .false.
 type bc_type
   real(pReal), dimension (3,3) :: deformation, &           ! applied velocity gradient or time derivative of deformation gradient
                                   stress, &                ! stress BC (if applicable)
                                   rotation                 ! rotation of BC (if applicable)
   real(pReal) ::                  timeIncrement, &         ! length of increment
                                   temperature              ! isothermal starting conditions
   integer(pInt) ::                steps, &                 ! number of steps
                                   outputfrequency, &       ! frequency of result writes
                                   restartfrequency, &      ! frequency of result writes
                                   logscale                 ! linear/logaritmic time step flag
   logical ::                      followFormerTrajectory,& ! follow trajectory of former loadcase
                                   velGradApplied           ! decide wether velocity gradient or fdot is given 
   logical, dimension(3,3) ::      maskDeformation, &       ! mask of boundary conditions
                                   maskStress
   logical, dimension(9) ::        maskStressVector         ! linear mask of boundary conditions    
 end type
 type(bc_type), allocatable, dimension(:) :: bc
 type(bc_type) :: bc_init
 character(len=3) ::             loadcase_string

! variables storing information from geom file
 real(pReal) :: wgt
 real(pReal), dimension(3) ::  geomdimension = 0.0_pReal                    ! physical dimension of volume element in each direction
 integer(pInt) :: homog                                                     ! homogenization scheme used
 integer(pInt), dimension(3) :: res = 1_pInt                                ! resolution (number of Fourier points) in each direction
 logical :: spectralPictureMode = .false.                                   ! indicating 1 to 1 mapping of FP to microstructure

! stress, stiffness and compliance average etc.
 real(pReal), dimension(3,3) ::                pstress, pstress_av, defgrad_av, &
                                               defgradAim = math_I3, defgradAimOld= math_I3, defgradAimCorr= math_I3,&
                                               mask_stress, mask_defgrad, fDot, &
                                               pstress_av_load, defgradAim_lab                           ! quantities rotated to other coordinate system
 real(pReal), dimension(3,3,3,3) ::            dPdF, c0_reference, c_current = 0.0_pReal, s_prev, c_prev ! stiffness and compliance
 real(pReal), dimension(6) ::                  cstress                                                   ! cauchy stress
 real(pReal), dimension(6,6) ::                dsde                                                      ! small strain stiffness
 real(pReal), dimension(9,9) ::                s_prev99, c_prev99                                        ! compliance and stiffness in matrix notation 
 real(pReal), dimension(:,:), allocatable ::   s_reduced, c_reduced                                      ! reduced compliance and stiffness (only for stress BC)
 integer(pInt) ::                              size_reduced = 0.0_pReal                                  ! number of stress BCs

! pointwise data 
 real(pReal), dimension(:,:,:,:,:), allocatable ::  workfft, defgrad, defgradold
 real(pReal), dimension(:,:,:,:), allocatable ::    coordinates
 real(pReal), dimension(:,:,:), allocatable ::      temperature

! variables storing information for spectral method and FFTW
 real(pReal), dimension(3,3) ::                         xiDyad                                   ! product of wave vectors
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                ! gamma operator (field) for spectral method
 real(pReal), dimension(:,:,:,:), allocatable ::        xi                                       ! wave vector field 
 integer(pInt), dimension(3) ::                         k_s                                
 integer*8, dimension(3) ::                             fftw_plan                                ! plans for fftw (forward and backward)
 integer*8 ::                                           fftw_flag                                ! planner flag for fftw
 
! loop variables, convergence etc.
 real(pReal) :: time = 0.0_pReal, time0 = 0.0_pReal, timeinc                                     ! elapsed time, begin of interval, time interval 
 real(pReal) :: guessmode, err_div, err_stress, err_stress_tol, p_hat_avg
 complex(pReal), parameter ::               img = cmplx(0.0_pReal,1.0_pReal)
 real(pReal), dimension(3,3), parameter ::  ones = 1.0_pReal, zeroes = 0.0_pReal
 complex(pReal), dimension(3,3) ::          temp33_Complex
 real(pReal), dimension(3,3) ::             temp33_Real
 integer(pInt) :: i, j, k, l, m, n, p
 integer(pInt) :: N_Loadcases, loadcase, step, iter, ielem, CPFEM_mode, &
                  ierr, notConvergedCounter = 0_pInt, totalStepsCounter = 0_pInt
 logical :: errmatinv, regrid = .false.
 real(pReal) :: defgradDet, defgradDetMax, defgradDetMin
 real(pReal) :: correctionFactor

! debuging variables
 real(pReal), dimension(:,:,:,:), allocatable :: divergence
 real(pReal) :: p_real_avg, err_div_max, err_real_div_avg, err_real_div_max
 logical :: debugGeneral = .false., debugDivergence = .false., debugRestart = .false.

! initialize default value for loadcase
 bc_init%deformation = zeroes; bc_init%stress = zeroes; bc_init%rotation = zeroes
 bc_init%timeIncrement = 0.0_pReal; bc_init%temperature = 300.0_pReal
 bc_init%steps = 0_pInt; bc_init%logscale = 0_pInt
 bc_init%outputfrequency = 1_pInt; bc_init%restartfrequency = 1_pInt
 bc_init%maskDeformation = .false.;  bc_init%maskStress = .false.
 bc_init%maskStressVector = .false.; bc_init%velGradApplied = .false.
 bc_init%followFormerTrajectory = .true.
 bc_init%rotation = math_I3                                                               ! assume no rotation
 
! Initializing model size independed parameters
 !$ call omp_set_num_threads(DAMASK_NumThreadsInt)                          ! set number of threads for parallel execution set by DAMASK_NUM_THREADS
 if (.not.(command_argument_count()==4 .or. command_argument_count()==6)) &! check for correct number of given arguments
    call IO_error(error_ID=102_pInt)
 
 call DAMASK_interface_init()
 
 !$OMP CRITICAL (write2out)
 print '(a)', ''
 print '(a,a)', ' <<<+-  DAMASK_spectral init  -+>>>'
 print '(a,a)', ' $Id$'
 print '(a)', ''
 print '(a,a)', ' Working Directory:    ',trim(getSolverWorkingDirectoryName())
 print '(a,a)', ' Solver Job Name:      ',trim(getSolverJobName())
 print '(a)', ''
 !$OMP END CRITICAL (write2out)

! Reading the loadcase file and allocate variables for loadcases
 path = getLoadcaseName()
 if (.not. IO_open_file(myUnit,path)) call IO_error(error_ID=30_pInt,ext_msg = trim(path))
 rewind(myUnit)
 do
   read(myUnit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                                         ! skip empty lines
   posLoadcase = IO_stringPos(line,maxNchunksLoadcase)
   do i = 1_pInt, maxNchunksLoadcase, 1_pInt                           ! reading compulsory parameters for loadcase
       select case (IO_lc(IO_stringValue(line,posLoadcase,i)))
            case('l', 'velocitygrad', 'velgrad','velocitygradient')
                 N_l = N_l + 1_pInt
            case('fdot')
                 N_Fdot = N_Fdot + 1_pInt
            case('t', 'time', 'delta')
                 N_t = N_t + 1_pInt
            case('n', 'incs', 'increments', 'steps', 'logincs', 'logsteps')
                 N_n = N_n + 1_pInt
        end select
   enddo                                                               ! count all identifiers to allocate memory and do sanity check
 enddo

100 N_Loadcases = N_n
 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &                        ! sanity check
   call IO_error(error_ID=37_pInt,ext_msg = trim(path))                ! error message for incomplete loadcase

 allocate (bc(N_Loadcases))

! Reading the loadcase and assign values to the allocated data structure
 rewind(myUnit)
 loadcase = 0_pInt
 do
   read(myUnit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                                                ! skip empty lines
   loadcase = loadcase + 1_pInt
   bc(loadcase) = bc_init
   posLoadcase = IO_stringPos(line,maxNchunksLoadcase)
   do j = 1_pInt,maxNchunksLoadcase
     select case (IO_lc(IO_stringValue(line,posLoadcase,j)))
       case('fdot','l','velocitygrad','velgrad','velocitygradient')                                 ! assign values for the deformation BC matrix
         bc(loadcase)%velGradApplied = (IO_lc(IO_stringValue(line,posLoadcase,j)) == 'l' .or. &     ! in case of given L, set flag to true
                                        IO_lc(IO_stringValue(line,posLoadcase,j)) == 'velocitygrad' .or. &
                                        IO_lc(IO_stringValue(line,posLoadcase,j)) == 'velgrad' .or. &
                                        IO_lc(IO_stringValue(line,posLoadcase,j)) == 'velocitygradient')
         temp_valueVector = 0.0_pReal
         temp_maskVector = .false.
         forall (k = 1_pInt:9_pInt) temp_maskVector(k) = IO_stringValue(line,posLoadcase,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (temp_maskVector(k)) temp_valueVector(k) = IO_floatValue(line,posLoadcase,j+k)
         enddo
         bc(loadcase)%maskDeformation = transpose(reshape(temp_maskVector,(/3,3/)))
         bc(loadcase)%deformation = math_plain9to33(temp_valueVector)
       case('p', 'pk1', 'piolakirchhoff', 'stress')
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) bc(loadcase)%maskStressVector(k) = IO_stringValue(line,posLoadcase,j+k) /= '*'
         do k = 1_pInt,9_pInt
           if (bc(loadcase)%maskStressVector(k)) temp_valueVector(k) = IO_floatValue(line,posLoadcase,j+k)  ! assign values for the bc(loadcase)%stress matrix
         enddo
         bc(loadcase)%maskStress = transpose(reshape(bc(loadcase)%maskStressVector,(/3,3/)))
         bc(loadcase)%stress = math_plain9to33(temp_valueVector)
       case('t','time','delta')                                                  ! increment time
         bc(loadcase)%timeIncrement = IO_floatValue(line,posLoadcase,j+1_pInt)
       case('temp','temperature')                                                ! starting temperature
         bc(loadcase)%temperature = IO_floatValue(line,posLoadcase,j+1_pInt)
       case('n','incs','increments','steps')                                     ! steps
         bc(loadcase)%steps = IO_intValue(line,posLoadcase,j+1_pInt)
       case('logincs','logsteps')                                                ! true, if log scale
         bc(loadcase)%steps = IO_intValue(line,posLoadcase,j+1_pInt)
         bc(loadcase)%logscale = 1_pInt
       case('f','freq','frequency','outputfreq')                                 ! frequency of result writings
         bc(loadcase)%outputfrequency = IO_intValue(line,posLoadcase,j+1_pInt)                
       case('r','restart','restartwrite')                                        ! frequency of writing restart information
         bc(loadcase)%restartfrequency = IO_intValue(line,posLoadcase,j+1_pInt)                
       case('guessreset','dropguessing')
         bc(loadcase)%followFormerTrajectory = .false.                           ! do not continue to predict deformation along former trajectory
       case('euler')                                                             ! rotation of loadcase given in euler angles
         p = 0_pInt                                                              ! assuming values given in radians
         l = 1_pInt                                                              ! assuming keyword indicating degree/radians
         select case (IO_lc(IO_stringValue(line,posLoadcase,j+1_pInt)))
           case('deg','degree')
             p = 1_pInt                                                          ! for conversion from degree to radian           
           case('rad','radian') 
           case default               
             l = 0_pInt                                                          ! imediately reading in angles, assuming radians
         end select
         forall(k = 1_pInt:3_pInt)  temp33_Real(k,1) = IO_floatValue(line,posLoadcase,j+l+k) * real(p,pReal) * inRad
         bc(loadcase)%rotation = math_EulerToR(temp33_Real(:,1))
       case('rotation','rot')                                                    ! assign values for the rotation of loadcase matrix
         temp_valueVector = 0.0_pReal
         forall (k = 1_pInt:9_pInt) temp_valueVector(k) = IO_floatValue(line,posLoadcase,j+k)
         bc(loadcase)%rotation = math_plain9to33(temp_valueVector)
     end select
 enddo; enddo

101 close(myUnit)

!read header of geom file to get the information needed before the complete geom file is intepretated by mesh.f90
 path = getModelName()

 if (.not. IO_open_file(myUnit,trim(path)//InputFileExtension))&
        call IO_error(error_ID=101_pInt,ext_msg = trim(path)//InputFileExtension)
 rewind(myUnit)
 read(myUnit,'(a1024)') line
 posGeom = IO_stringPos(line,2_pInt)
 keyword = IO_lc(IO_StringValue(line,posGeom,2_pInt))
 if (keyword(1:4) == 'head') then
   headerLength = IO_intValue(line,posGeom,1_pInt) + 1_pInt
 else
   call IO_error(error_ID=42_pInt)
 endif
 
 rewind(myUnit)
 do i = 1_pInt, headerLength
   read(myUnit,'(a1024)') line
   posGeom = IO_stringPos(line,maxNchunksGeom)             
   select case ( IO_lc(IO_StringValue(line,posGeom,1)) )
     case ('dimension')
       gotDimension = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,posGeom,j)))
           case('x')
              geomdimension(1) = IO_floatValue(line,posGeom,j+1_pInt)
           case('y')
              geomdimension(2) = IO_floatValue(line,posGeom,j+1_pInt)
           case('z')
              geomdimension(3) = IO_floatValue(line,posGeom,j+1_pInt)
         end select
       enddo
     case ('homogenization')
       gotHomogenization = .true.
       homog = IO_intValue(line,posGeom,2_pInt)
     case ('resolution')
       gotResolution = .true.
       do j = 2_pInt,6_pInt,2_pInt
         select case (IO_lc(IO_stringValue(line,posGeom,j)))
           case('a')
             res(1) = IO_intValue(line,posGeom,j+1_pInt)
           case('b')
             res(2) = IO_intValue(line,posGeom,j+1_pInt)
           case('c')
             res(3) = IO_intValue(line,posGeom,j+1_pInt)
         end select
       enddo
     case ('picture')
       spectralPictureMode = .true.
   end select
 enddo
 close(myUnit)
 if (.not.(gotDimension .and. gotHomogenization .and. gotResolution)) call IO_error(error_ID=45_pInt)
 
 if(mod(res(1),2_pInt)/=0_pInt .or.&
    mod(res(2),2_pInt)/=0_pInt .or.&
   (mod(res(3),2_pInt)/=0_pInt .and. res(3)/= 1_pInt))  call IO_error(error_ID=103_pInt)

! Initialization of CPFEM_general (= constitutive law)
 call CPFEM_initAll(bc(1)%temperature,1_pInt,1_pInt)

! Get debugging parameters
 if (iand(spectral_debug_verbosity,1_pInt)==1_pInt) debugGeneral = .true.
 if (iand(spectral_debug_verbosity,2_pInt)==2_pInt) debugDivergence = .true.
 if (iand(spectral_debug_verbosity,4_pInt)==4_pInt) debugRestart = .true.

!Output of geometry
 !$OMP CRITICAL (write2out)
 print '(a)', ''
 print '(a)', '#############################################################'
 print '(a)', 'DAMASK spectral:'
 print '(a)', 'The spectral method boundary value problem solver for'
 print '(a)', 'the Duesseldorf Advanced Material Simulation Kit'
 print '(a)', '#############################################################'
 print '(a,a)', 'Geom File Name:       ',trim(path)//'.geom'
 print '(a)', '============================================================='
 print '(a,i12,i12,i12)','resolution a b c:', res
 print '(a,f12.5,f12.5,f12.5)','dimension x y z:', geomdimension
 print '(a,i5)','homogenization:     ',homog
 print '(a,L)','spectralPictureMode: ',spectralPictureMode
 print '(a)', '#############################################################'
 print '(a,a)','Loadcase File Name:   ',trim(getLoadcaseName())
 !$OMP END CRITICAL (write2out)
 if (bc(1)%followFormerTrajectory) then
   call IO_warning(warning_ID=33_pInt)                 ! cannot guess along trajectory for first step of first loadcase
   bc(1)%followFormerTrajectory = .false.
 endif
 
 ! consistency checks and output of loadcase
 do loadcase = 1_pInt, N_Loadcases
   !$OMP CRITICAL (write2out)
   print '(a)', '============================================================='
   print '(a,i5)', 'Loadcase:            ', loadcase
   write (loadcase_string, '(i3)' ) loadcase
   if (.not. bc(loadcase)%followFormerTrajectory) &
     print '(a)', 'Drop Guessing Along Trajectory'
   !$OMP END CRITICAL (write2out)
   if (any(bc(loadcase)%maskStress .eqv. bc(loadcase)%maskDeformation))&                ! exclusive or masking only
     call IO_error(error_ID=31_pInt,ext_msg=loadcase_string)
   if (any(bc(loadcase)%maskStress.and.transpose(bc(loadcase)%maskStress).and.&         !checking if no rotation is allowed by stress BC
     reshape((/.false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false./),(/3,3/))))&
     call IO_error(error_ID=38_pInt,ext_msg=loadcase_string)
   if (bc(loadcase)%velGradApplied) then
     do j = 1_pInt, 3_pInt
       if (any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .true.) .and.&
           any(bc(loadcase)%maskDeformation(j,1:3) .eqv. .false.)) call IO_error(error_ID=32_pInt,ext_msg=loadcase_string)     ! each line should be either fully or not at all defined
     enddo
     !$OMP CRITICAL (write2out)
     print '(a)','Velocity Gradient:'
     !$OMP END CRITICAL (write2out)
   else
     !$OMP CRITICAL (write2out)
     print '(a)','Change of Deformation Gradient:'
     !$OMP END CRITICAL (write2out)
   endif
   !$OMP CRITICAL (write2out)
   print '(3(3(f12.6,x)/)\)', merge(math_transpose3x3(bc(loadcase)%deformation),&
                                   reshape(spread(DAMASK_NaN,1,9),(/3,3/)),transpose(bc(loadcase)%maskDeformation))
   print '(a,/,3(3(f12.6,x)/)\)','Stress Boundary Condition/MPa:',merge(math_transpose3x3(bc(loadcase)%stress),&
                                              reshape(spread(DAMASK_NaN,1,9),(/3,3/)),&
                                              transpose(bc(loadcase)%maskStress))*1e-6
   !$OMP END CRITICAL (write2out)
   if (any(abs(math_mul33x33(bc(loadcase)%rotation,math_transpose3x3(bc(loadcase)%rotation))-math_I3)&                ! given rotation matrix contains strain
               >reshape(spread(rotation_tol,1,9),(/3,3/)))&
               .or. abs(math_det3x3(bc(loadcase)%rotation))>1.0_pReal + rotation_tol) call IO_error(error_ID=46_pInt,ext_msg=loadcase_string)
   !$OMP CRITICAL (write2out)
   if (any(bc(loadcase)%rotation/=math_I3)) &
               print '(a,3(3(f12.6,x)/)\)','Rotation of BCs:',math_transpose3x3(bc(loadcase)%rotation)
   !$OMP END CRITICAL (write2out)
   if (bc(loadcase)%timeIncrement < 0.0_pReal) call IO_error(error_ID=34_pInt,ext_msg=loadcase_string)                 ! negative time increment
   !$OMP CRITICAL (write2out)
   print '(a,f12.6)','Temperature:',bc(loadcase)%temperature
   print '(a,f12.6)','Time:       ',bc(loadcase)%timeIncrement
   !$OMP END CRITICAL (write2out)
   if (bc(loadcase)%steps < 1_pInt) call IO_error(error_ID=35_pInt,ext_msg=loadcase_string)                            ! non-positive increment count
   !$OMP CRITICAL (write2out)
   print '(a,i5)','Steps:       ',bc(loadcase)%steps
   !$OMP END CRITICAL (write2out)
   if (bc(loadcase)%outputfrequency < 1_pInt) call IO_error(error_ID=36_pInt,ext_msg=loadcase_string)                  ! non-positive result frequency
   !$OMP CRITICAL (write2out)
   print '(a,i5)','Freq. of Results Output: ',bc(loadcase)%outputfrequency
   !$OMP END CRITICAL (write2out)
   if (bc(loadcase)%restartfrequency < 1_pInt) call IO_error(error_ID=39_pInt,ext_msg=loadcase_string)                 ! non-positive restart frequency
   !$OMP CRITICAL (write2out)
   print '(a,i5)','Freq. of Restart Information Output: ',bc(loadcase)%restartfrequency
   !$OMP END CRITICAL (write2out)
 enddo

! Initialization of fftw (see manual on fftw.org for more details)
#ifdef _OPENMP
   if(DAMASK_NumThreadsInt>0_pInt) then
     call dfftw_init_threads(ierr)
     if(ierr == 0_pInt) call IO_error(error_ID=104_pInt)
     call dfftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
   endif
#endif
 !call dfftw_timelimit(fftw_timelimit)                          ! is not working, have to fix it in FFTW source file
 select case(IO_lc(fftw_planner_flag))                          ! setting parameters for the plan creation of FFTW. Basically a translation from fftw3.f
   case('estimate','fftw_estimate')                             ! ordered from slow execution (but fast plan creation) to fast execution
     fftw_flag = 64
   case('measure','fftw_measure')
     fftw_flag = 0
   case('patient','fftw_patient')
     fftw_flag= 32
   case('exhaustive','fftw_exhaustive')
     fftw_flag = 8 
   case default
     call IO_warning(warning_ID=47_pInt,ext_msg=trim(IO_lc(fftw_planner_flag)))
     fftw_flag = 32
 end select
!*************************************************************
! Loop over loadcases defined in the loadcase file
 do loadcase = 1_pInt,  N_Loadcases
!*************************************************************
   time0 = time                                                            ! loadcase start time                
   if (bc(loadcase)%followFormerTrajectory) then                           ! continue to guess along former trajectory where applicable
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                 ! change of load case, homogeneous guess for the first step
   endif
   
   mask_defgrad =  merge(ones,zeroes,bc(loadcase)%maskDeformation)
   mask_stress =  merge(ones,zeroes,bc(loadcase)%maskStress)
   size_reduced = count(bc(loadcase)%maskStressVector)
   allocate (c_reduced(size_reduced,size_reduced));          c_reduced = 0.0_pReal
   allocate (s_reduced(size_reduced,size_reduced));          s_reduced = 0.0_pReal

   timeinc = bc(loadcase)%timeIncrement/bc(loadcase)%steps                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
   fDot = bc(loadcase)%deformation                                        ! only valid for given fDot. will be overwritten later in case L is given

!*************************************************************
! loop oper steps defined in input file for current loadcase
   do step = 1_pInt,  bc(loadcase)%steps
!*************************************************************
! forwarding time
     if (bc(loadcase)%logscale == 1_pInt) then                                                 ! loglinear scale
       if (loadcase == 1_pInt) then                                                            ! 1st loadcase of loglinear scale            
         if (step == 1_pInt) then                                                              ! 1st step of 1st loadcase of loglinear scale
           timeinc = bc(1)%timeIncrement*(2.0_pReal**real(     1_pInt-bc(1)%steps ,pReal))     ! assume 1st step is equal to 2nd 
         else                                                                                  ! not-1st step of 1st loadcase of loglinear scale
           timeinc = bc(1)%timeIncrement*(2.0_pReal**real(step-1_pInt-bc(1)%steps ,pReal))
         endif
       else                                                                                      ! not-1st loadcase of loglinear scale
           timeinc = time0 *( (1.0_pReal + bc(loadcase)%timeIncrement/time0 )**real(          step/bc(loadcase)%steps ,pReal)  &
                             -(1.0_pReal + bc(loadcase)%timeIncrement/time0 )**real( (step-1_pInt)/bc(loadcase)%steps ,pReal) )
       endif
     endif
     time = time + timeinc
     totalStepsCounter = totalStepsCounter + 1_pInt

!*************************************************************
! Initialization Start
!*************************************************************
     if(totalStepsCounter >= restartReadStep) then                          ! Do calculations (otherwise just forwarding)

       if (regrid==.true. ) then                                            ! 'DeInitialize' the values changing in case of regridding
         regrid = .false.
         call dfftw_destroy_plan(fftw_plan(1)); call dfftw_destroy_plan(fftw_plan(2))
         if(debugDivergence) call dfftw_destroy_plan(fftw_plan(3))
         deallocate (defgrad)
         deallocate (defgradold)
         deallocate (coordinates)
         deallocate (temperature)
         deallocate (xi)
         deallocate (workfft)
        !ToDo:  here we have to create the new geometry and assign the values from the previous step
       endif 

       if(totalStepsCounter == restartReadStep) then                        ! Initialize values
         guessmode = 0.0_pReal                                              ! change of load case, homogeneous guess for the first step
         allocate (defgrad    (  res(1),res(2),res(3),3,3));  defgrad     = 0.0_pReal
         allocate (defgradold (  res(1),res(2),res(3),3,3));  defgradold  = 0.0_pReal
         allocate (coordinates(3,res(1),res(2),res(3)));      coordinates = 0.0_pReal
         allocate (temperature(  res(1),res(2),res(3)));      temperature = bc(1)%temperature  ! start out isothermally
         allocate (xi         (3,res(1)/2+1,res(2),res(3)));  xi          =0.0_pReal 
         allocate (workfft(res(1)+2,res(2),res(3),3,3)); workfft = 0.0_pReal
         if (debugDivergence) allocate (divergence(res(1)+2,res(2),res(3),3)); divergence = 0.0_pReal

         wgt = 1.0_pReal/real(res(1)*res(2)*res(3), pReal)
         call dfftw_plan_many_dft_r2c(fftw_plan(1),3,(/res(1),res(2),res(3)/),9,&
           workfft,(/res(1)       +2_pInt,res(2),res(3)/),1,(res(1)       +2_pInt)*res(2)*res(3),&
           workfft,(/res(1)/2_pInt+1_pInt,res(2),res(3)/),1,(res(1)/2_pInt+1_pInt)*res(2)*res(3),fftw_flag)   
         call dfftw_plan_many_dft_c2r(fftw_plan(2),3,(/res(1),res(2),res(3)/),9,&
           workfft,(/res(1)/2_pInt+1_pInt,res(2),res(3)/),1,(res(1)/2_pInt+1_pInt)*res(2)*res(3),&
           workfft,(/res(1)       +2_pInt,res(2),res(3)/),1,(res(1)       +2_pInt)*res(2)*res(3),fftw_flag)
         if (debugDivergence) &
           call dfftw_plan_many_dft_c2r(fftw_plan(3),3,(/res(1),res(2),res(3)/),3,&
             divergence,(/res(1)/2_pInt+1_pInt,res(2),res(3)/),1,(res(1)/2_pInt+1_pInt)*res(2)*res(3),&
             divergence,(/res(1)       +2_pInt,res(2),res(3)/),1,(res(1)       +2_pInt)*res(2)*res(3),fftw_flag)
         if (debugGeneral) then
           !$OMP CRITICAL (write2out)
           write (6,*) 'FFTW initialized'
           !$OMP END CRITICAL (write2out)
         endif

         if (restartReadStep==1_pInt) then                                            ! not restarting, no deformation at the beginning
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             defgrad(i,j,k,1:3,1:3) = math_I3            
             defgradold(i,j,k,1:3,1:3) = math_I3
           enddo; enddo; enddo
         else                                                                         ! using old values 
           if (IO_read_jobBinaryFile(777,'convergedSpectralDefgrad',trim(getSolverJobName()),size(defgrad))) then
             read (777,rec=1) defgrad
             close (777)
           endif
           defgradold = defgrad
           defgradAim = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             defgradAim = defgradAim + defgrad(i,j,k,1:3,1:3)                        ! calculating old average deformation
           enddo; enddo; enddo
           defgradAim = defgradAim * wgt
           defgradAimOld = defgradAim
           guessmode=0.0_pInt
         endif

         ielem = 0_pInt
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt 
           coordinates(1:3,i,j,k) = mesh_ipCenterOfGravity(1:3,1,ielem)  ! set to initial coordinates ToDo: SHOULD BE UPDATED TO CURRENT POSITION IN FUTURE REVISIONS!!! But do we know them? I don't think so. Otherwise we don't need geometry reconstruction
           call CPFEM_general(2_pInt,coordinates(1:3,i,j,k),math_I3,math_I3,temperature(i,j,k),0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
           c_current = c_current + dPdF
         enddo; enddo; enddo
         c0_reference = c_current * wgt                                                  ! linear reference material stiffness
         c_prev = math_rotate_forward3x3x3x3(c0_reference,bc(loadcase)%rotation) ! rotate_forward: lab -> load system 
   
         if (debugGeneral) then
           !$OMP CRITICAL (write2out)
           write (6,*) 'First Call to CPFEM_general finished'
           !$OMP END CRITICAL (write2out)
         endif
   
         do k = 1_pInt, res(3)                              ! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
           k_s(3) = k - 1_pInt
           if(k > res(3)/2_pInt + 1_pInt) k_s(3) = k_s(3) - res(3)
             do j = 1_pInt, res(2)
               k_s(2) = j - 1_pInt
               if(j > res(2)/2_pInt + 1_pInt) k_s(2) = k_s(2) - res(2) 
                 do i = 1, res(1)/2_pInt + 1_pInt
                   k_s(1) = i - 1_pInt
                   xi(3,i,j,k) = 0.0_pReal                                                  ! 2D case
                   if(res(3) > 1_pInt) xi(3,i,j,k) = real(k_s(3), pReal)/geomdimension(3) ! 3D case  
                                              xi(2,i,j,k) = real(k_s(2), pReal)/geomdimension(2)
                                              xi(1,i,j,k) = real(k_s(1), pReal)/geomdimension(1)
         enddo; enddo; enddo
  ! remove highest frequencies for calculation of divergence (CAREFULL, they will be used for pre calculatet gamma operator!) 
         do k = 1_pInt ,res(3); do j = 1_pInt ,res(2); do i = 1_pInt,res(1)/2_pInt + 1_pInt
           if(k==res(3)/2_pInt+1_pInt) xi(3,i,j,k)= 0.0_pReal
           if(j==res(2)/2_pInt+1_pInt) xi(2,i,j,k)= 0.0_pReal
           if(i==res(1)/2_pInt+1_pInt) xi(1,i,j,k)= 0.0_pReal
         enddo; enddo; enddo
   
         if(memory_efficient) then                            ! allocate just single fourth order tensor
           allocate (gamma_hat(1,1,1,3,3,3,3)); gamma_hat = 0.0_pReal
         else                                                 ! precalculation of gamma_hat field
           allocate (gamma_hat(res(1)/2_pInt + 1_pInt ,res(2),res(3),3,3,3,3)); gamma_hat = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)/2_pInt + 1_pInt
             if (any(xi(:,i,j,k) /= 0.0_pReal)) then     
               do l = 1_pInt ,3_pInt; do m = 1_pInt,3_pInt
                 xiDyad(l,m) = xi(l,i,j,k)*xi(m,i,j,k)
               enddo; enddo
               temp33_Real = math_inv3x3(math_mul3333xx33(c0_reference, xiDyad)) 
             else
               xiDyad  = 0.0_pReal
               temp33_Real = 0.0_pReal
             endif 
             do l=1_pInt,3_pInt; do m=1_pInt,3_pInt; do n=1_pInt,3_pInt; do p=1_pInt,3_pInt
               gamma_hat(i,j,k, l,m,n,p) = - 0.25*(temp33_Real(l,n)+temp33_Real(n,l)) *&
                                                   (xiDyad(m,p)+xiDyad(p,m))
             enddo; enddo; enddo; enddo         
           enddo; enddo; enddo
         endif

   ! write header of output file
         !$OMP CRITICAL (write2out)
         open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())&
                                                            //'.spectralOut',form='UNFORMATTED',status='REPLACE')
         write(538), 'load', trim(getLoadcaseName())
         write(538), 'workingdir', trim(getSolverWorkingDirectoryName())
         write(538), 'geometry', trim(getSolverJobName())//InputFileExtension
         write(538), 'resolution', res
         write(538), 'dimension', geomdimension
         write(538), 'materialpoint_sizeResults', materialpoint_sizeResults
         write(538), 'loadcases', N_Loadcases
         write(538), 'logscale', bc(1:N_Loadcases)%logscale                                   ! one entry per loadcase (0: linear, 1: log)
         write(538), 'frequencies', bc(1:N_Loadcases)%outputfrequency                         ! one entry per loadcase
         write(538), 'times', bc(1:N_Loadcases)%timeIncrement                                 ! one entry per loadcase
         bc(1)%steps= bc(1)%steps + 1_pInt                                                   
         write(538), 'increments', bc(1:N_Loadcases)%steps                                    ! one entry per loadcase ToDo: rename keyword to steps
         bc(1)%steps= bc(1)%steps - 1_pInt
         write(538), 'startingIncrement', restartReadStep -1_pInt                             ! start with writing out the previous step
         write(538), 'eoh'                                                                    ! end of header
         write(538),  materialpoint_results(materialpoint_sizeResults,1,res(1)*res(2)*res(3)) ! initial (non-deformed) results
        !$OMP END CRITICAL (write2out)
       endif
  !*************************************************************
  ! Initialization End
  !*************************************************************
          
     if (mod(step - 1_pInt,bc(loadcase)%restartFrequency)==0_pInt) then                    ! at frequency of writing restart information
       restartWrite = .true.                                                               ! setting restart parameter for FEsolving (first call to CPFEM_general will write ToDo: true?) 
     else
       restartWrite = .false.
     endif
       
     if (bc(loadcase)%velGradApplied) &                                                           ! calculate fDot from given L and current F
       fDot = math_mul33x33(bc(loadcase)%deformation, defgradAim)

  !winding forward of deformation aim in loadcase system
       temp33_Real = defgradAim                                            
       defgradAim = defgradAim &                                                                         
                  + guessmode * mask_stress * (defgradAim - defgradAimOld) &      
                  + mask_defgrad * fDot * timeinc 
       defgradAimOld = temp33_Real
   
  ! update local deformation gradient
       if (any(bc(loadcase)%rotation/=math_I3)) then                                               ! lab and loadcase coordinate system are NOT the same 
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           temp33_Real = defgrad(i,j,k,1:3,1:3)
           if (bc(loadcase)%velGradApplied) &                                                      ! use velocity gradient to calculate new deformation gradient (if not guessing)
                                  fDot = math_mul33x33(bc(loadcase)%deformation,&
                                  math_rotate_forward3x3(defgradold(i,j,k,1:3,1:3),bc(loadcase)%rotation))
             defgrad(i,j,k,1:3,1:3) = defgrad(i,j,k,1:3,1:3) &                                      ! decide if guessing along former trajectory or apply homogeneous addon
                                + guessmode * (defgrad(i,j,k,1:3,1:3) - defgradold(i,j,k,1:3,1:3))& ! guessing... 
                                + math_rotate_backward3x3((1.0_pReal-guessmode) * mask_defgrad * fDot,&
                                           bc(loadcase)%rotation) *timeinc                          ! apply the prescribed value where deformation is given if not guessing
             defgradold(i,j,k,1:3,1:3) = temp33_Real 
         enddo; enddo; enddo
       else                                                                                         ! one coordinate system for lab and loadcase, save some multiplications
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           temp33_Real = defgrad(i,j,k,1:3,1:3)
           if (bc(loadcase)%velGradApplied) &                                                       ! use velocity gradient to calculate new deformation gradient (if not guessing)
                                 fDot = math_mul33x33(bc(loadcase)%deformation,defgradold(i,j,k,1:3,1:3))
             defgrad(i,j,k,1:3,1:3) = defgrad(i,j,k,1:3,1:3) &                                      ! decide if guessing along former trajectory or apply homogeneous addon
                                + guessmode * (defgrad(i,j,k,1:3,1:3) - defgradold(i,j,k,1:3,1:3))& ! guessing... 
                                + (1.0_pReal-guessmode) * mask_defgrad * fDot * timeinc               ! apply the prescribed value where deformation is given if not guessing
             defgradold(i,j,k,1:3,1:3) = temp33_Real 
         enddo; enddo; enddo
       endif
       guessmode = 1.0_pReal                                                              ! keep guessing along former trajectory during same loadcase

       CPFEM_mode = 1_pInt                                                                ! winding forward
       iter = 0_pInt
       err_div = 2.0_pReal * err_div_tol                                                  ! go into loop 

       if(size_reduced > 0_pInt) then                                                     ! calculate compliance in case stress BC is applied
         c_prev99 = math_Plain3333to99(c_prev)
         k = 0_pInt                                                                       ! build reduced stiffness
         do n = 1_pInt,9_pInt
           if(bc(loadcase)%maskStressVector(n)) then
             k = k + 1_pInt
             j = 0_pInt
             do m = 1_pInt,9_pInt
               if(bc(loadcase)%maskStressVector(m)) then
                 j = j + 1_pInt
                 c_reduced(k,j) = c_prev99(n,m)
         endif; enddo; endif; enddo
         call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)               ! invert reduced stiffness
         if(errmatinv) call IO_error(error_ID=800)
         s_prev99 = 0.0_pReal                                                             ! build full compliance
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

       !$OMP CRITICAL (write2out)
       print '(a)', '#############################################################'
       print '(A,I5.5,A,es12.6)', 'Increment ', totalStepsCounter, ' Time ',time
       if (restartWrite ) then
         print '(A)', 'Writing converged Results of previous Step for Restart'
       if(IO_write_jobBinaryFile(777,'convergedSpectralDefgrad',size(defgrad))) then       ! and writing deformation gradient field to file
         write (777,rec=1) defgrad
         close (777)
       endif
     endif 
       !$OMP END CRITICAL (write2out)
  !*************************************************************
  ! convergence loop
       do while(iter < itmax .and. &
               (err_div     > err_div_tol    .or. &
                err_stress  > err_stress_tol)) 
         iter = iter + 1_pInt
  !*************************************************************
         print '(a)', '============================================================='
         print '(5(A,I5.5))', 'Loadcase ',loadcase,' Step ',step,'/',bc(loadcase)%steps,'@Iteration ',iter,'/',itmax
         do n = 1_pInt,3_pInt; do m = 1_pInt,3_pInt
           defgrad_av(m,n) = sum(defgrad(1:res(1),1:res(2),1:res(3),m,n)) * wgt
         enddo; enddo
         !$OMP CRITICAL (write2out)
         print '(a,/,3(3(f12.7,x)/)\)', 'Deformation Gradient:',math_transpose3x3(defgrad_av)
         print '(A)', '... Update Stress Field (Constitutive Evaluation P(F)) ......'
         !$OMP END CRITICAL (write2out)
         defgradDetMax = -huge(1.0_pReal)
         defgradDetMin = +huge(1.0_pReal)
         ielem = 0_pInt
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           defgradDet = math_det3x3(defgrad(i,j,k,1:3,1:3))
           defgradDetMax = max(defgradDetMax,defgradDet)
           defgradDetMin = min(defgradDetMin,defgradDet) 
           ielem = ielem + 1_pInt
           call CPFEM_general(3_pInt,&                                                       ! collect cycle
                              coordinates(1:3,i,j,k), defgradold(i,j,k,1:3,1:3), defgrad(i,j,k,1:3,1:3),&
                              temperature(i,j,k),timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo; enddo

         print '(a,x,es10.4)'       , 'Maximum Determinant of Deformation:', defgradDetMax
         print '(a,x,es10.4)'       , 'Minimum Determinant of Deformation:', defgradDetMin

         workfft = 0.0_pReal                                                            ! needed because of the padding for FFTW
         c_current = 0.0_pReal
         ielem = 0_pInt       
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(CPFEM_mode,&                                              ! first element in first iteration retains CPFEM_mode 1, 
                              coordinates(1:3,i,j,k),&
                              defgradold(i,j,k,1:3,1:3), defgrad(i,j,k,1:3,1:3),&       ! others get 2 (saves winding forward effort)
                              temperature(i,j,k),timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress,dPdF)
           CPFEM_mode = 2_pInt
           workfft(i,j,k,1:3,1:3) = pstress                                     ! build up average P-K stress 
           c_current = c_current + dPdF
         enddo; enddo; enddo
         restartWrite = .false. ! ToDo: don't know if we need it. Depends on how CPFEM_general is writing results 
         do n = 1_pInt,3_pInt; do m = 1_pInt,3_pInt
           pstress_av(m,n) = sum(workfft(1:res(1),1:res(2),1:res(3),m,n)) * wgt
         enddo; enddo
         
         !$OMP CRITICAL (write2out)
         print '(a,/,3(3(f12.7,x)/)\)', 'Piola-Kirchhoff Stress / MPa: ',math_transpose3x3(pstress_av)/1.e6
         
         err_stress_tol = 0.0_pReal
         pstress_av_load = math_rotate_forward3x3(pstress_av,bc(loadcase)%rotation)
         if(size_reduced > 0_pInt) then                                                                            ! calculate stress BC if applied
           err_stress = maxval(abs(mask_stress * (pstress_av_load - bc(loadcase)%stress)))                         ! maximum deviaton (tensor norm not applicable)
           err_stress_tol = maxval(abs(mask_defgrad * pstress_av_load)) * err_stress_tolrel                        ! don't use any tensor norm because the comparison should be coherent
           print '(A)', '... Correcting Deformation Gradient to Fullfill BCs .........'
           print '(2(a,es10.4))', 'Error Stress = ',err_stress, ', Tol. = ', err_stress_tol 
           defgradAimCorr = - math_mul3333xx33(s_prev, ((pstress_av_load - bc(loadcase)%stress)))                  ! residual on given stress components
           defgradAim = defgradAim + defgradAimCorr
           print '(a,/,3(3(f12.7,x)/)\)', 'New Deformation Aim:     ',math_transpose3x3(math_rotate_backward3x3(&
                                                                 defgradAim,bc(loadcase)%rotation))
           print '(a,x,es10.4)'         , 'Determinant of New Deformation Aim:', math_det3x3(defgradAim)
         endif
         print '(A)', '... Calculating Equilibrium Using Spectral Method ...........'
         !$OMP END CRITICAL (write2out)
         call dfftw_execute_dft_r2c(fftw_plan(1),workfft,workfft)                                                  ! FFT of pstress
   
         p_hat_avg = sqrt(maxval (math_eigenvalues3x3(math_mul33x33(workfft(1,1,1,1:3,1:3),&                       ! L_2 norm of average stress in fourier space,  
                                                  math_transpose3x3(workfft(1,1,1,1:3,1:3))))))                    ! ignore imaginary part as it is always zero for real only input))
         err_div = 0.0_pReal
         err_div_max = 0.0_pReal
         do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)/2_pInt+1_pInt
           err_div = err_div + sqrt(sum((&                                                                         ! avg of L_2 norm of div(stress) in fourier space (Suquet small strain)
                                         math_mul33x3_complex(workfft(i*2_pInt-1_pInt,j,k,1:3,1:3) + &
                                                              workfft(i*2_pInt       ,j,k,1:3,1:3)*img,&
                                                              xi(1:3,i,j,k))&
                                           )**2.0_pReal))                                                    
           if(debugDivergence) &
             err_div_max = max(err_div_max,abs(sqrt(sum((&                                                         !  maximum of L two norm of div(stress) in fourier space (Suquet large strain)
                                           math_mul33x3_complex(workfft(i*2_pInt-1_pInt,j,k,1:3,1:3)+&
                                                                workfft(i*2_pInt       ,j,k,1:3,1:3)*img,&
                                                                xi(1:3,i,j,k))&
                                                           )**2.0_pReal))))
         enddo; enddo; enddo
         correctionFactor = minval(geomdimension)*wgt**(-1.0_pReal/4.0_pReal)                                      ! multiplying by minimum dimension to get rid of dimension dependency and phenomenologigal factor wgt**(-1/4) to get rid of resolution dependency
         if (res(3)==1_pInt) correctionFactor = minval(geomdimension(1:2))*wgt**(-1.0_pReal/4.0_pReal)             ! 2D case, ToDo: correct?
         if (.not. divergence_correction) correctionFactor = 1.0_pReal   

         err_div = err_div*wgt/p_hat_avg*correctionFactor                                                          ! weighting by points and average stress and multiplying with correction factor
         err_div_max = err_div_max/p_hat_avg*correctionFactor                                                      ! weighting by average stress and multiplying with correction factor

         if(memory_efficient) then                                                                                 ! memory saving version, on-the-fly calculation of gamma_hat
           do k = 1_pInt, res(3); do j = 1_pInt, res(2) ;do i = 1_pInt, res(1)/2_pInt+1_pInt                         
             if (any(xi(:,i,j,k) /= 0.0_pReal)) then     
               do l = 1_pInt,3_pInt; do m = 1_pInt,3_pInt
                 xiDyad(l,m) = xi(l,i,j,k)*xi(m,i,j,k)
               enddo; enddo
               temp33_Real = math_inv3x3(math_mul3333xx33(c0_reference, xiDyad)) 
             else
               xiDyad = 0.0_pReal
               temp33_Real = 0.0_pReal
             endif 
             do l=1_pInt,3_pInt; do m=1_pInt,3_pInt; do n=1_pInt,3_pInt; do p=1_pInt,3_pInt
               gamma_hat(1,1,1, l,m,n,p) = - 0.25_pReal*(temp33_Real(l,n)+temp33_Real(n,l))*&
                                                         (xiDyad(m,p) +xiDyad(p,m))
             enddo; enddo; enddo; enddo   
             do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt
               temp33_Complex(m,n) = sum(gamma_hat(1,1,1,m,n,1:3,1:3) *(workfft(i*2_pInt-1_pInt,j,k,1:3,1:3)&                 
                                                                       +workfft(i*2_pInt       ,j,k,1:3,1:3)*img))
             enddo; enddo
             workfft(i*2_pInt-1_pInt,j,k,1:3,1:3) = real (temp33_Complex) 
             workfft(i*2_pInt       ,j,k,1:3,1:3) = aimag(temp33_Complex)
           enddo; enddo; enddo
         else                                                                                                       ! use precalculated gamma-operator
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)/2_pInt+1_pInt
             do m = 1_pInt,3_pInt; do n = 1_pInt,3_pInt
               temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n,1:3,1:3) *(workfft(i*2_pInt-1_pInt,j,k,1:3,1:3)&
                                                                       + workfft(i*2_pInt       ,j,k,1:3,1:3)*img))
             enddo; enddo
             workfft(i*2_pInt-1_pInt,j,k,1:3,1:3) = real (temp33_Complex)
             workfft(i*2_pInt       ,j,k,1:3,1:3) = aimag(temp33_Complex) 
           enddo; enddo; enddo
         endif
         if(debugDivergence) then
           divergence = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)/2_pInt+1_pInt
             ! real part at i*2-1, imaginary part at i*2 and multiply by i ==> switch and change sign
             divergence(i*2_pInt-1_pInt,j,k,1:3) = workfft(i*2_pInt       ,j,k,1:3,1)*xi(1:3,i,j,k)*pi*2.0_pReal&
                                                 + workfft(i*2_pInt       ,j,k,1:3,2)*xi(1:3,i,j,k)*pi*2.0_pReal&
                                                 + workfft(i*2_pInt       ,j,k,1:3,3)*xi(1:3,i,j,k)*pi*2.0_pReal
             divergence(i*2_pInt,j,k,1:3)      = - workfft(i*2_pInt-1_pInt,j,k,1:3,1)*xi(1:3,i,j,k)*pi*2.0_pReal&
                                                 - workfft(i*2_pInt-1_pInt,j,k,1:3,2)*xi(1:3,i,j,k)*pi*2.0_pReal&
                                                 - workfft(i*2_pInt-1_pInt,j,k,1:3,3)*xi(1:3,i,j,k)*pi*2.0_pReal
           enddo; enddo; enddo
           divergence = divergence*correctionFactor
           call dfftw_execute_dft_c2r(fftw_plan(3),divergence,divergence)
           divergence = divergence * wgt
           err_real_div_avg = 0.0_pReal
           err_real_div_max = 0.0_pReal
           do k = 1_pInt, res(3); do j = 1_pInt, res(2); do i = 1_pInt, res(1)
             err_real_div_avg = err_real_div_avg + sqrt(sum((divergence(i,j,k,1:3))**2.0_pReal))              ! avg of L_2 norm of div(stress) in fourier space (Suquet small strain)
             err_real_div_max = max(err_real_div_max,abs(sqrt(sum((divergence(i,j,k,1:3))**2.0_pReal))))      ! maximum of L two norm of div(stress) in fourier space (Suquet large strain)
           enddo; enddo; enddo
           p_real_avg = sqrt(maxval (math_eigenvalues3x3(math_mul33x33(pstress_av,&                          ! L_2 norm of average stress in fourier space,  
                                                      math_transpose3x3(pstress_av)))))                        ! ignore imaginary part as it is always zero for real only input))
           err_real_div_avg = err_real_div_avg*wgt/p_real_avg
           err_real_div_max = err_real_div_max/p_real_avg
         endif

  ! average strain 
         workfft(1,1,1,1:3,1:3) = defgrad_av - math_I3                                        ! zero frequency (real part)
         workfft(2,1,1,1:3,1:3) = 0.0_pReal                                                   ! zero frequency (imaginary part)
         
         call dfftw_execute_dft_c2r(fftw_plan(2),workfft,workfft)
         defgrad = defgrad + workfft(1:res(1),:,:,:,:)*wgt
         do m = 1,3; do n = 1,3 
           defgrad_av(m,n) = sum(defgrad(:,:,:,m,n))*wgt
         enddo; enddo
         defgradAim_lab = math_rotate_backward3x3(defgradAim,bc(loadcase)%rotation)
         do m = 1,3; do n = 1,3 
           defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + (defgradAim_lab(m,n) - defgrad_av(m,n))  ! anticipated target minus current state
         enddo; enddo
         !$OMP CRITICAL (write2out)
         if(.not. debugDivergence) then
           print '(2(a,es10.4))', 'Error Divergence = ',err_div,    ', Tol. = ', err_div_tol
         else
           print '(2(a,es10.4))', 'Error Divergence FT avg= ',err_div,    ', Tol. = ', err_div_tol
           print '(a,es10.4)', 'Error Divergence FT max= ',err_div_max
           print '(a,es10.4)', 'Error Divergence Real avg= ',err_real_div_avg
           print '(a,es10.4)', 'Error Divergence Real max= ',err_real_div_max
         endif
         !$OMP END CRITICAL (write2out)
         
       enddo    ! end looping when convergency is achieved 
       
       c_prev = math_rotate_forward3x3x3x3(c_current*wgt,bc(loadcase)%rotation)        ! calculate stiffness for next step
       !ToDo: Incfluence for next loadcase
       !$OMP CRITICAL (write2out)
       print '(a)', '============================================================='
       if(err_div<=err_div_tol .and. err_stress<=err_stress_tol) then
         print '(A,I5.5,A)', 'Increment ', totalStepsCounter, ' Converged'
       else
         print '(A,I5.5,A)', 'Increment ', totalStepsCounter, ' NOT Converged'
         notConvergedCounter = notConvergedCounter + 1
       endif
       if (mod(totalStepsCounter -1_pInt,bc(loadcase)%outputfrequency) == 0_pInt) then                 ! at output frequency
         print '(A)', '... Writing Results to File .................................'
         write(538),  materialpoint_results(materialpoint_sizeResults,1,res(1)*res(2)*res(3))          ! write result to file
       endif
       !$OMP END CRITICAL (write2out)
     endif
   enddo  ! end looping over steps in current loadcase
   deallocate(c_reduced)
   deallocate(s_reduced)
   enddo    ! end looping over loadcases
   !$OMP CRITICAL (write2out)
   print '(a)', '#############################################################'
   print '(a,i5.5,a,i5.5,a)', 'Of ', totalStepsCounter - restartReadStep + 1_pInt, ' Calculated Steps, ', notConvergedCounter, ' Steps did not Converge!'
   !$OMP END CRITICAL (write2out)
 close(538)
 call dfftw_destroy_plan(fftw_plan(1)); call dfftw_destroy_plan(fftw_plan(2))
 if(debugDivergence) call dfftw_destroy_plan(fftw_plan(3))
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
