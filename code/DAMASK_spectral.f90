! Copyright 2011 Max-Planck-Institut fuer Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Duesseldorf Advanced MAterial Simulation Kit.
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
!             - start program with DAMASK_spectral PathToGeomFile/NameOfGeom.geom
!               PathToLoadFile/NameOfLoadFile.load
!             - PathToGeomFile will be the working directory
!             - make sure the file "material.config" exists in the working
!               directory. For further configuration use "numerics.config"
!********************************************************************
program DAMASK_spectral
!********************************************************************

 use DAMASK_interface
 use prec, only: pInt, pReal
 use IO
 use math
 use mesh, only: mesh_ipCenterOfGravity
 use CPFEM, only: CPFEM_general, CPFEM_initAll 
 use numerics, only: err_div_tol, err_stress_tol, err_stress_tolrel,&
                     relevantStrain, itmax, memory_efficient, DAMASK_NumThreadsInt       
 use homogenization, only: materialpoint_sizeResults, materialpoint_results
!$ use OMP_LIB                                                           ! the openMP function library

 implicit none
 include 'include/fftw3.f' ! header file for fftw3 (declaring variables). Library files are also needed
                           ! compile FFTW 3.2.2 with ./configure --enable-threads
! variables to read from loadcase and geom file
 real(pReal), dimension(9) ::                      valuevector           ! stores information temporarily from loadcase file
 integer(pInt), parameter ::                       maxNchunksInput = 26  ! 5 identifiers, 18 values for the matrices and 3 scalars
 integer(pInt), dimension (1+maxNchunksInput*2) :: posInput
 integer(pInt), parameter ::                       maxNchunksGeom = 7    ! 4 identifiers, 3 values
 integer(pInt), dimension (1+2*maxNchunksGeom) ::  posGeom
 integer(pInt) unit, N_l, N_s, N_t, N_n, N_freq, N_Fdot, N_temperature   ! numbers of identifiers
 character(len=1024) path, line
 logical gotResolution,gotDimension,gotHomogenization

! variables storing information from loadcase file
 real(pReal)                                    time, time0, timeinc     ! elapsed time, begin of interval, time interval
 real(pReal), dimension (:,:,:), allocatable :: bc_deformation, &        ! applied velocity gradient or time derivative of deformation gradient
                                                bc_stress                ! stress BC (if applicable)
 real(pReal), dimension(:), allocatable ::      bc_timeIncrement, &      ! length of increment
                                                bc_temperature           ! isothermal starting conditions
 integer(pInt)                                  N_Loadcases, step        ! ToDo: rename?
 integer(pInt), dimension(:), allocatable ::    bc_steps, &              ! number of steps
                                                bc_frequency, &          ! frequency of result writes
                                                bc_logscale              ! linear/logaritmic time step flag
 logical, dimension(:), allocatable ::          followFormerTrajectory,& ! follow trajectory of former loadcase
                                                velGradApplied           ! decide wether velocity gradient or fdot is given 
 logical, dimension(:,:,:,:), allocatable ::    bc_mask                  ! mask of boundary conditions
 logical, dimension(:,:,:), allocatable ::      bc_maskvector            ! linear mask of boundary conditions

! variables storing information from geom file
 real(pReal) wgt
 real(pReal), dimension(3) ::  geomdimension                             ! physical dimension of volume element in each direction
 integer(pInt) homog                                                     ! homogenization scheme used
 integer(pInt), dimension(3) :: resolution                               ! resolution (number of Fourier points) in each direction

! stress etc.
 real(pReal), dimension(3,3) ::                         ones, zeroes, temp33_Real, &
                                                        pstress, pstress_av, defgrad_av,&
                                                        defgradAim, defgradAimOld,  defgradAimCorr,&
                                                        mask_stress, mask_defgrad, fDot                          
 real(pReal), dimension(3,3,3,3) ::                     dPdF, c0_reference, c_current, s_prev, c_prev ! stiffness and compliance
 real(pReal), dimension(6) ::                           cstress                                       ! cauchy stress
 real(pReal), dimension(6,6) ::                         dsde                                          ! small strain stiffness
 real(pReal), dimension(9,9) ::                         s_prev99, c_prev99                            ! compliance and stiffness in matrix notation         
 real(pReal), dimension(:,:,:,:,:), allocatable ::      workfft, defgrad, defgradold
 real(pReal), dimension(:,:,:,:), allocatable ::        coordinates
 real(pReal), dimension(:,:,:), allocatable ::          temperature
 real(pReal), dimension(:,:), allocatable ::            s_reduced, c_reduced                     ! reduced compliance and stiffness (only for stress BC)
 integer(pInt)                                          size_reduced                             ! number of stress BCs
 
! variables storing information for spectral method
 complex(pReal) ::                                      img
 complex(pReal), dimension(3,3) ::                      temp33_Complex
 real(pReal), dimension(3,3) ::                         xiDyad                                   ! product of wave vectors
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat                                ! gamma operator (field) for spectral method
 real(pReal), dimension(:,:,:,:), allocatable ::        xi                                       ! wave vector field 
 integer(pInt), dimension(3) ::                         k_s                                
 integer*8, dimension(2) ::                             plan_fft                                 ! plans for fftw (forward and backward)
 
! loop variables, convergence etc.
 real(pReal) guessmode, err_div, err_stress, p_hat_avg        
 integer(pInt)  i, j, k, l, m, n, p
 integer(pInt)  loadcase, ielem, iter, CPFEM_mode, ierr, not_converged_counter
 logical errmatinv

!Initializing
!$ call omp_set_num_threads(DAMASK_NumThreadsInt)         ! set number of threads for parallel execution set by DAMASK_NUM_THREADS

 print*, ''
 print*, '<<<+-  DAMASK_spectral init  -+>>>'
 print*, '$Id$'
 print*, ''
 
 unit = 234_pInt
 ones = 1.0_pReal; zeroes = 0.0_pReal
 img = cmplx(0.0,1.0)
 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_temperature = 0_pInt
 
 time = 0.0_pReal
 N_n = 0_pInt
 N_freq = 0_pInt
 N_Fdot = 0_pInt
 not_converged_counter = 0_pInt
 gotResolution =.false.; gotDimension =.false.; gotHomogenization = .false.
 resolution = 1_pInt
 geomdimension = 0.0_pReal
 

 if (IargC() /= 2) call IO_error(102)                     ! check for correct number of given arguments

! Reading the loadcase file and assign variables
 path = getLoadcaseName()
 print '(a,/,a)', 'Loadcase:      ',trim(path)
 print '(a,/,a)', 'Workingdir:    ',trim(getSolverWorkingDirectoryName())
 print '(a,/,a)', 'SolverJobName: ',trim(getSolverJobName())

 if (.not. IO_open_file(unit,path)) call IO_error(30,ext_msg = trim(path))
 
 rewind(unit)
 do
   read(unit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posInput = IO_stringPos(line,maxNchunksInput)
   do i = 1, maxNchunksInput, 1
       select case (IO_lc(IO_stringValue(line,posInput,i)))
            case('l', 'velocitygrad')
                 N_l = N_l+1
            case('fdot')
                 N_Fdot = N_Fdot+1
            case('s', 'stress', 'pk1', 'piolakirchhoff')
                 N_s = N_s+1
            case('t', 'time', 'delta')
                 N_t = N_t+1
            case('n', 'incs', 'increments', 'steps', 'logincs', 'logsteps')
                 N_n = N_n+1
            case('f', 'freq', 'frequency')
                 N_freq = N_freq+1
            case('temp','temperature')
                 N_temperature = N_temperature+1
        end select
   enddo                                                  ! count all identifiers to allocate memory and do sanity check
 enddo

101 N_Loadcases = N_n
 if ((N_l + N_Fdot /= N_n) .or. (N_n /= N_t)) &            ! sanity check
   call IO_error(31,ext_msg = trim(path))                  ! error message for incomplete loadcase

! allocate memory depending on lines in input file
 allocate (bc_deformation(3,3,N_Loadcases));        bc_deformation = 0.0_pReal
 allocate (bc_stress(3,3,N_Loadcases));             bc_stress = 0.0_pReal
 allocate (bc_mask(3,3,2,N_Loadcases));             bc_mask = .false.
 allocate (bc_maskvector(9,2,N_Loadcases));         bc_maskvector = .false.
 allocate (velGradApplied(N_Loadcases));            velGradApplied = .false.
 allocate (bc_timeIncrement(N_Loadcases));          bc_timeIncrement = 0.0_pReal
 allocate (bc_temperature(N_Loadcases));            bc_temperature = 300.0_pReal
 allocate (bc_steps(N_Loadcases));                  bc_steps = 0_pInt
 
 allocate (bc_logscale(N_Loadcases));               bc_logscale = 0_pInt
 allocate (bc_frequency(N_Loadcases));              bc_frequency = 1_pInt
 allocate (followFormerTrajectory(N_Loadcases));    followFormerTrajectory = .true.

 rewind(unit)
 loadcase = 0_pInt
 do
   read(unit,'(a1024)',END = 200) line
   if (IO_isBlank(line)) cycle                                                    ! skip empty lines
   loadcase = loadcase + 1
   posInput = IO_stringPos(line,maxNchunksInput)
   do j = 1,maxNchunksInput,2
     select case (IO_lc(IO_stringValue(line,posInput,j)))
       case('fdot','l','velocitygrad')                                            ! assign values for the deformation BC matrix
         velGradApplied(loadcase) = (IO_lc(IO_stringValue(line,posInput,j)) == 'l' .or. &
                                     IO_lc(IO_stringValue(line,posInput,j)) == 'velocitygrad')         ! in case of given L, set flag to true
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k,1,loadcase) = IO_stringValue(line,posInput,j+k) /= '*'
         do k = 1,9
           if (bc_maskvector(k,1,loadcase)) valuevector(k) = IO_floatValue(line,posInput,j+k)
         enddo
         bc_mask(:,:,1,loadcase) = transpose(reshape(bc_maskvector(1:9,1,loadcase),(/3,3/)))
         bc_deformation(:,:,loadcase) = math_plain9to33(valuevector)
       case('s', 'stress', 'pk1', 'piolakirchhoff')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k,2,loadcase) = IO_stringValue(line,posInput,j+k) /= '*'
         do k = 1,9
           if (bc_maskvector(k,2,loadcase)) valuevector(k) = IO_floatValue(line,posInput,j+k)  ! assign values for the bc_stress matrix
         enddo
         bc_mask(:,:,2,loadcase) = transpose(reshape(bc_maskvector(1:9,2,loadcase),(/3,3/)))
         bc_stress(:,:,loadcase) = math_plain9to33(valuevector)
       case('t','time','delta')                                            ! increment time
           bc_timeIncrement(loadcase) = IO_floatValue(line,posInput,j+1)
       case('temp','temperature')                                          ! starting temperature
           bc_temperature(loadcase) = IO_floatValue(line,posInput,j+1)
       case('n','incs','increments','steps')                               ! bc_steps
           bc_steps(loadcase) = IO_intValue(line,posInput,j+1)
       case('logincs','logsteps')                                          ! true, if log scale
           bc_steps(loadcase) = IO_intValue(line,posInput,j+1)
           bc_logscale(loadcase) = 1_pInt
       case('f','freq','frequency')                                        ! frequency of result writings
           bc_frequency(loadcase) = IO_intValue(line,posInput,j+1)                
       case('guessreset','dropguessing')
           followFormerTrajectory(loadcase) = .false.                       ! do not continue to predict deformation along former trajectory
     end select
 enddo; enddo

200 close(unit)

 if (followFormerTrajectory(1)) then
   call IO_warning(33)                 ! cannot guess along trajectory for first step of first loadcase
   followFormerTrajectory(1) = .false.
 endif
 
 do loadcase = 1, N_Loadcases                                                      ! consistency checks and output
   print *, '------------------------------------------------------'
   print '(a,i5)', 'Loadcase:', loadcase
   if (.not. followFormerTrajectory(loadcase)) &
     print '(a)', 'drop guessing along trajectory'
   if (any(bc_mask(:,:,1,loadcase) .eqv. bc_mask(:,:,2,loadcase)))&                ! exclusive or masking only
     call IO_error(31,loadcase)
   if (any(bc_mask(1:3,1:3,2,loadcase).and.transpose(bc_mask(1:3,1:3,2,loadcase)).and.&
           reshape((/.false.,.true.,.true.,.true.,.false.,.true.,.true.,.true.,.false./),(/3,3/))))&
     call IO_error(38,loadcase)
   if (velGradApplied(loadcase)) then
     do j = 1, 3
       if (any(bc_mask(j,:,1,loadcase) .eqv. .true.) .and.&
           any(bc_mask(j,:,1,loadcase) .eqv. .false.)) call IO_error(32,loadcase)     ! each line should be either fully or not at all defined
     enddo
     print '(a,/,3(3(f12.6,x)/))','L:'        ,math_transpose3x3(bc_deformation(:,:,loadcase))
     print '(a,/,3(3(l,x)/))',    'bc_mask for L:',transpose(bc_mask(:,:,1,loadcase))
   else
     print '(a,/,3(3(f12.6,x)/))','Fdot:'       ,math_transpose3x3(bc_deformation(:,:,loadcase))
     print '(a,/,3(3(l,x)/))',    'bc_mask for Fdot:',transpose(bc_mask(:,:,1,loadcase))
   endif
   print '(a,/,3(3(f12.6,x)/))','bc_stress/MPa:',math_transpose3x3(bc_stress(:,:,loadcase))*1e-6
   print '(a,/,3(3(l,x)/))',    'bc_mask for stress:'     ,transpose(bc_mask(:,:,2,loadcase))
   if (bc_timeIncrement(loadcase) < 0.0_pReal) call IO_error(34,loadcase)                 ! negative time increment
   print '(a,f12.6)','temperature: ',bc_temperature(loadcase)
   print '(a,f12.6)','time: ',bc_timeIncrement(loadcase)
   if (bc_steps(loadcase) < 1_pInt) call IO_error(35,loadcase)                            ! non-positive increment count
   print '(a,i6)','incs: ',bc_steps(loadcase)
   if (bc_frequency(loadcase) < 1_pInt) call IO_error(36,loadcase)                        ! non-positive result frequency
   print '(a,i6)','freq: ',bc_frequency(loadcase)
 enddo

!read header of geom file to get the information needed before the complete geom file is intepretated by mesh.f90
 path = getModelName()
 print *, '------------------------------------------------------'
 print '(a,a)', 'GeomName: ',trim(path)
 if (.not. IO_open_file(unit,trim(path)//InputFileExtension)) call IO_error(101,ext_msg = trim(path)//InputFileExtension)

 rewind(unit)
 do
   read(unit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posGeom = IO_stringPos(line,maxNchunksGeom)             

   select case ( IO_lc(IO_StringValue(line,posGeom,1)) )
     case ('dimension')
       gotDimension = .true.
       do i = 2,6,2
         select case (IO_lc(IO_stringValue(line,posGeom,i)))
           case('x')
              geomdimension(1) = IO_floatValue(line,posGeom,i+1)
           case('y')
              geomdimension(2) = IO_floatValue(line,posGeom,i+1)
           case('z')
              geomdimension(3) = IO_floatValue(line,posGeom,i+1)
         end select
       enddo
     case ('homogenization')
       gotHomogenization = .true.
       homog = IO_intValue(line,posGeom,2)
     case ('resolution')
       gotResolution = .true.
       do i = 2,6,2
         select case (IO_lc(IO_stringValue(line,posGeom,i)))
           case('a')
             resolution(1) = IO_intValue(line,posGeom,i+1)
           case('b')
             resolution(2) = IO_intValue(line,posGeom,i+1)
           case('c')
             resolution(3) = IO_intValue(line,posGeom,i+1)
         end select
       enddo
   end select
   if (gotDimension .and. gotHomogenization .and. gotResolution) exit
 enddo
 100 close(unit)
 
 if(mod(resolution(1),2_pInt)/=0_pInt .or.&
    mod(resolution(2),2_pInt)/=0_pInt .or.&
   (mod(resolution(3),2_pInt)/=0_pInt .and. resolution(3)/= 1_pInt))  call IO_error(103)

 print '(a,/,i5,i5,i5)','resolution a b c:', resolution
 print '(a,/,f12.5,f12.5,f12.5)','dimension x y z:', geomdimension
 print '(a,i4)','homogenization: ',homog
 
 allocate (defgrad    (  resolution(1),resolution(2),resolution(3),3,3));  defgrad     = 0.0_pReal
 allocate (defgradold (  resolution(1),resolution(2),resolution(3),3,3));  defgradold  = 0.0_pReal
 allocate (coordinates(3,resolution(1),resolution(2),resolution(3)));      coordinates = 0.0_pReal
 allocate (temperature(  resolution(1),resolution(2),resolution(3)));      temperature = bc_temperature(1)  ! start out isothermally
 allocate (xi         (3,resolution(1)/2+1,resolution(2),resolution(3)));  xi          =0.0_pReal 

 wgt = 1.0_pReal/real(resolution(1)*resolution(2)*resolution(3), pReal)
 defgradAim    = math_I3
 defgradAimOld = math_I3
 defgrad_av    = math_I3
 
! Initialization of CPFEM_general (= constitutive law) and of deformation gradient field
 call CPFEM_initAll(bc_temperature(1),1_pInt,1_pInt)
 ielem = 0_pInt
 c_current = 0.0_pReal 
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   defgradold(i,j,k,:,:) = math_I3                       ! no deformation at the beginning
   defgrad(i,j,k,:,:) = math_I3 
   ielem = ielem +1 
   coordinates(1:3,i,j,k) = mesh_ipCenterOfGravity(1:3,1,ielem)  ! set to initial coordinates ToDo: SHOULD BE UPDATED TO CURRENT POSITION IN FUTURE REVISIONS!!!
   call CPFEM_general(2,coordinates(1:3,i,j,k),math_I3,math_I3,temperature(i,j,k),0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
   c_current = c_current + dPdF
 enddo; enddo; enddo
 c0_reference = c_current * wgt  ! linear reference material stiffness
 c_prev = c0_reference
 
 do k = 1, resolution(3)                              ! calculation of discrete angular frequencies, ordered as in FFTW (wrap around)
    k_s(3) = k-1
    if(k > resolution(3)/2+1) k_s(3) = k_s(3)-resolution(3)
    do j = 1, resolution(2)
      k_s(2) = j-1
      if(j > resolution(2)/2+1) k_s(2) = k_s(2)-resolution(2) 
      do i = 1, resolution(1)/2+1
        k_s(1) = i-1
        xi(3,i,j,k) = 0.0_pReal                                                  ! 2D case
        if(resolution(3) > 1) xi(3,i,j,k) = real(k_s(3), pReal)/geomdimension(3) ! 3D case  
                              xi(2,i,j,k) = real(k_s(2), pReal)/geomdimension(2)
                              xi(1,i,j,k) = real(k_s(1), pReal)/geomdimension(1)
 enddo; enddo; enddo
! remove highest frequencies for calculation of divergence (CAREFULL, they will be used for pre calculatet gamma operator!) 
 do k = 1,resolution(3); do j = 1,resolution(2); do i = 1,resolution(1)/2+1
   if(k==resolution(3)/2+1) xi(3,i,j,k)= 0.0_pReal
   if(j==resolution(2)/2+1) xi(2,i,j,k)= 0.0_pReal
   if(i==resolution(1)/2+1) xi(1,i,j,k)= 0.0_pReal
 enddo; enddo; enddo
 
 if(memory_efficient) then                            ! allocate just single fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3)); gamma_hat = 0.0_pReal
 else                                                 ! precalculation of gamma_hat field
   allocate (gamma_hat(resolution(1)/2+1,resolution(2),resolution(3),3,3,3,3)); gamma_hat = 0.0_pReal
   do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
     if (any(xi(:,i,j,k) /= 0.0_pReal)) then     
       do l = 1,3; do m = 1,3
          xiDyad(l,m) = xi(l,i,j,k)*xi(m,i,j,k)
       enddo; enddo
       temp33_Real = math_inv3x3(math_mul3333xx33(c0_reference, xiDyad)) 
     else
        xiDyad  = 0.0_pReal
        temp33_Real = 0.0_pReal
     endif 
     do l=1,3; do m=1,3; do n=1,3; do p=1,3
       gamma_hat(i,j,k, l,m,n,p) = - 0.25*(temp33_Real(l,n)+temp33_Real(n,l)) *&
                                          (xiDyad(m,p)+xiDyad(p,m))
     enddo; enddo; enddo; enddo         
   enddo; enddo; enddo
 endif

 allocate (workfft(resolution(1)+2,resolution(2),resolution(3),3,3)); workfft = 0.0_pReal

! Initialization of fftw (see manual on fftw.org for more details)

 if(DAMASK_NumThreadsInt>0_pInt) then
   call dfftw_init_threads(ierr)
   if(ierr == 0_pInt) call IO_error(104,ierr)
   call dfftw_plan_with_nthreads(DAMASK_NumThreadsInt) 
 endif

 call dfftw_plan_many_dft_r2c(plan_fft(1),3,(/resolution(1),resolution(2),resolution(3)/),9,&
   workfft,(/resolution(1)  +2,resolution(2),resolution(3)/),1,(resolution(1)  +2)*resolution(2)*resolution(3),&
   workfft,(/resolution(1)/2+1,resolution(2),resolution(3)/),1,(resolution(1)/2+1)*resolution(2)*resolution(3),FFTW_PATIENT)   
 call dfftw_plan_many_dft_c2r(plan_fft(2),3,(/resolution(1),resolution(2),resolution(3)/),9,&
   workfft,(/resolution(1)/2+1,resolution(2),resolution(3)/),1,(resolution(1)/2+1)*resolution(2)*resolution(3),&
   workfft,(/resolution(1)  +2,resolution(2),resolution(3)/),1,(resolution(1)  +2)*resolution(2)*resolution(3),FFTW_PATIENT)

! write header of output file
 open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())&
                                                    //'.spectralOut',form='UNFORMATTED')       
 write(538), 'load', trim(getLoadcaseName())
 write(538), 'workingdir', trim(getSolverWorkingDirectoryName())
 write(538), 'geometry', trim(getSolverJobName())//InputFileExtension
 write(538), 'resolution', resolution
 write(538), 'dimension', geomdimension
 write(538), 'materialpoint_sizeResults', materialpoint_sizeResults
 write(538), 'loadcases', N_Loadcases
 write(538), 'logscale', bc_logscale                                       ! one entry per loadcase (0: linear, 1: log)
 write(538), 'frequencies', bc_frequency                                   ! one entry per loadcase
 write(538), 'times', bc_timeIncrement                                     ! one entry per loadcase
 bc_steps(1) = bc_steps(1)+1                                               ! +1 to store initial situation
 write(538), 'increments', bc_steps                                        ! one entry per loadcase
 bc_steps(1) = bc_steps(1)-1                                               ! re-adjust for correct looping
 write(538), 'eoh'                                                         ! end of header

 write(538)  materialpoint_results(:,1,:)                                  ! initial (non-deformed) results
! Initialization done

!*************************************************************
! Loop over loadcases defined in the loadcase file
 do loadcase = 1, N_Loadcases
!*************************************************************
   time0 = time                                                            ! loadcase start time                

   if (followFormerTrajectory(loadcase)) then                              ! continue to guess along former trajectory where applicable
     guessmode = 1.0_pReal
   else
     guessmode = 0.0_pReal                                                 ! change of load case, homogeneous guess for the first step
   endif
   
   mask_defgrad =  merge(ones,zeroes,bc_mask(:,:,1,loadcase))
   mask_stress =  merge(ones,zeroes,bc_mask(:,:,2,loadcase))
   size_reduced = count(bc_maskvector(1:9,2,loadcase))
   allocate (c_reduced(size_reduced,size_reduced));          c_reduced = 0.0_pReal
   allocate (s_reduced(size_reduced,size_reduced));          s_reduced = 0.0_pReal

   timeinc = bc_timeIncrement(loadcase)/bc_steps(loadcase)                ! only valid for given linear time scale. will be overwritten later in case loglinear scale is used
   fDot = bc_deformation(:,:,loadcase)                                  ! only valid for given fDot. will be overwritten later in case L is given
!*************************************************************
! loop oper steps defined in input file for current loadcase
   do step = 1, bc_steps(loadcase)
!*************************************************************
     if (bc_logscale(loadcase) == 1_pInt) then                                                    ! loglinear scale
        if (loadcase == 1_pInt) then                                                              ! 1st loadcase of loglinear scale            
            if (step == 1_pInt) then                                                              ! 1st step of 1st loadcase of loglinear scale
                timeinc = bc_timeIncrement(1)*(2.0**(1 - bc_steps(1)))                            ! assume 1st step is equal to 2nd 
            else                                                                                  ! not-1st step of 1st loadcase of loglinear scale
                timeinc = bc_timeIncrement(1)*(2.0**(step - (1 + bc_steps(1))))
            endif
        else                                                                                      ! not-1st loadcase of loglinear scale
            timeinc = time0 * (  ((1.0_pReal+bc_timeIncrement(loadcase)/time0)**(float( step   )/(bc_steps(loadcase))))  &
                               - ((1.0_pReal+bc_timeIncrement(loadcase)/time0)**(float((step-1))/(bc_steps(loadcase)))) )
        endif
     endif
     time = time + timeinc
     
     if (velGradApplied(loadcase)) &                                                 ! calculate fDot from given L and current F
       fDot = math_mul33x33(bc_deformation(1:3,1:3,loadcase), defgradAim)

!winding forward of deformation aim
     temp33_Real = defgradAim                                            
     defgradAim = defgradAim &                                                                         
                + guessmode * mask_stress * (defgradAim - defgradAimOld) &      
                + mask_defgrad * fDot * timeinc 
     defgradAimOld = temp33_Real
 
! update local deformation gradient
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       temp33_Real = defgrad(i,j,k,:,:)
       if (velGradApplied(loadcase)) &                                                  ! use velocity gradient to calculate new deformation gradient (if not guessing)
         fDot = math_mul33x33(bc_deformation(1:3,1:3,loadcase),defgradold(i,j,k,1:3,1:3))
         defgrad(i,j,k,1:3,1:3) = defgrad(i,j,k,1:3,1:3) &                                      ! decide if guessing along former trajectory or apply homogeneous addon
                            + guessmode * (defgrad(i,j,k,1:3,1:3) - defgradold(i,j,k,1:3,1:3))& ! guessing... 
                            + (1.0_pReal-guessmode) * mask_defgrad * fDot *timeinc    ! apply the prescribed value where deformation is given if not guessing
         defgradold(i,j,k,1:3,1:3) = temp33_Real   
     enddo; enddo; enddo
     guessmode = 1.0_pReal                                                              ! keep guessing along former trajectory during same loadcase

     CPFEM_mode = 1_pInt                                                                ! winding forward
     iter = 0_pInt
     err_div = 2.0_pReal * err_div_tol                                                  ! go into loop 

     if(size_reduced > 0_pInt) then                                                     ! calculate compliance in case stress BC is applied
       c_prev99 = math_Plain3333to99(c_prev)
       k = 0_pInt                                                                       ! build reduced stiffness
       do n = 1,9
         if(bc_maskvector(n,2,loadcase)) then
           k = k + 1_pInt
           j = 0_pInt
           do m = 1,9
             if(bc_maskvector(m,2,loadcase)) then
               j = j + 1_pInt
               c_reduced(k,j) = c_prev99(n,m)
       endif; enddo; endif; enddo
       call math_invert(size_reduced, c_reduced, s_reduced, i, errmatinv)               ! invert reduced stiffness
       if(errmatinv) call IO_error(800)
       s_prev99 = 0.0_pReal                                                             ! build full compliance
       k = 0_pInt
       do n = 1,9
         if(bc_maskvector(n,2,loadcase)) then
           k = k + 1_pInt
           j = 0_pInt
           do m = 1,9
           if(bc_maskvector(m,2,loadcase)) then
                 j = j + 1_pInt
                 s_prev99(n,m) = s_reduced(k,j)
       endif; enddo; endif; enddo
       s_prev = (math_Plain99to3333(s_prev99))
     endif
!*************************************************************
! convergence loop
     do while(iter < itmax .and. &
             (err_div     > err_div_tol    .or. &
              err_stress  > err_stress_tol)) 
       iter = iter + 1_pInt
       print '(A)', '************************************************************'
       print '(3(A,I5.5,tr2)A)', '**** Loadcase = ',loadcase, 'Step = ',step, 'Iteration = ',iter,'****'
       print '(A)', '************************************************************'
       workfft = 0.0_pReal                                                            ! needed because of the padding for FFTW
!*************************************************************
       do n = 1,3; do m = 1,3
         defgrad_av(m,n) = sum(defgrad(:,:,:,m,n)) * wgt
       enddo; enddo
       print '(a,/,3(3(f12.7,x)/))', 'Deformation Gradient:',math_transpose3x3(defgrad_av)
       
       print '(A,/)', '== Update Stress Field (Constitutive Evaluation P(F)) ======'
       ielem = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1
         call CPFEM_general(3,&                                                       ! collect cycle
                            coordinates(1:3,i,j,k), defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                            temperature(i,j,k),timeinc,ielem,1_pInt,&
                            cstress,dsde, pstress, dPdF)
       enddo; enddo; enddo
       
       c_current = 0.0_pReal
       ielem = 0_pInt       
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1_pInt
         call CPFEM_general(CPFEM_mode,&                                              ! first element in first iteration retains CPFEM_mode 1, 
                            coordinates(1:3,i,j,k),&
                            defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&               ! others get 2 (saves winding forward effort)
                            temperature(i,j,k),timeinc,ielem,1_pInt,&
                            cstress,dsde, pstress,dPdF)
         CPFEM_mode = 2_pInt
         workfft(i,j,k,:,:) = pstress                                     ! build up average P-K stress 
         c_current = c_current + dPdF
       enddo; enddo; enddo
       
       do n = 1,3; do m = 1,3
         pstress_av(m,n) = sum(workfft(1:resolution(1),:,:,m,n)) * wgt
       enddo; enddo
       print '(a,/,3(3(f12.7,x)/))', 'Piola-Kirchhoff Stress / MPa: ',math_transpose3x3(pstress_av)/1.e6

       err_stress_tol = 0.0_pReal
       if(size_reduced > 0_pInt) then                                                                            ! calculate stress BC if applied
         err_stress = maxval(abs(mask_stress * (pstress_av - bc_stress(1:3,1:3,loadcase))))                      ! maximum deviaton (tensor norm not applicable)
         err_stress_tol = maxval(abs(mask_defgrad * pstress_av)) * err_stress_tolrel                             ! don't use any tensor norm because the comparison should be coherent
         err_stress_tol = err_stress_tol * err_stress_tolrel                                                     ! weighting by relative criterion 
         print '(A,/)', '== Correcting Deformation Gradient to Fullfill BCs ========='
         print '(2(a,E10.5)/)', 'Error Stress = ',err_stress, ', Tol. = ', err_stress_tol 
         defgradAimCorr = - math_mul3333xx33(s_prev, ((pstress_av - bc_stress(1:3,1:3,loadcase))))               ! residual on given stress components
         defgradAim = defgradAim + defgradAimCorr
         print '(a,/,3(3(f12.7,x)/))', 'Deformation Aim:     ',math_transpose3x3(defgradAim)
         print '(a,x,f12.7,/)'       , 'Determinant of Deformation Aim: ', math_det3x3(defgradAim)
       endif
       print '(A,/)', '== Calculating Equilibrium Using Spectral Method ===========' 
       call dfftw_execute_dft_r2c(plan_fft(1),workfft,workfft)                                                   ! FFT of pstress
 
       p_hat_avg = sqrt(maxval (math_eigenvalues3x3(math_mul33x33(workfft(1,1,1,1:3,1:3),&                       ! L_2 norm of average stress in fourier space,  
                                                math_transpose3x3(workfft(1,1,1,1:3,1:3))))))                    ! ignore imaginary part as it is always zero for real only input))
       err_div = 0.0_pReal
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
         err_div = err_div + sqrt(sum((math_mul33x3_complex(workfft(i*2-1,j,k,1:3,1:3)+&                         ! avg of L_2 norm of div(stress) in fourier space (Suquet small strain)
                                                            workfft(i*2,  j,k,1:3,1:3)*img,xi(1:3,i,j,k)))**2.0))                                                    
       enddo; enddo; enddo

       err_div = err_div*wgt/p_hat_avg*(minval(geomdimension)*wgt**(-1/4))                                       ! weigthting,  multiplying by minimum dimension to get rid of dimension dependency and phenomenologigal factor wgt**(-1/4) to get rid of resolution dependency

       if(memory_efficient) then                                                                                 ! memory saving version, on-the-fly calculation of gamma_hat
         do k = 1, resolution(3); do j = 1, resolution(2) ;do i = 1, resolution(1)/2+1                         
           if (any(xi(:,i,j,k) /= 0.0_pReal)) then     
             do l = 1,3; do m = 1,3
               xiDyad(l,m) = xi(l,i,j,k)*xi(m,i,j,k)
             enddo; enddo
             temp33_Real = math_inv3x3(math_mul3333xx33(c0_reference, xiDyad)) 
           else
             xiDyad = 0.0_pReal
             temp33_Real = 0.0_pReal
           endif 
           do l=1,3; do m=1,3; do n=1,3; do p=1,3
             gamma_hat(1,1,1, l,m,n,p) = - 0.25_pReal*(temp33_Real(l,n)+temp33_Real(n,l))*&
                                                       (xiDyad(m,p) +xiDyad(p,m))
           enddo; enddo; enddo; enddo   
           do m = 1,3; do n = 1,3
             temp33_Complex(m,n) = sum(gamma_hat(1,1,1,m,n,:,:) *(workfft(i*2-1,j,k,:,:)&                 
                                                                 +workfft(i*2  ,j,k,:,:)*img))
           enddo; enddo
           workfft(i*2-1,j,k,:,:) = real (temp33_Complex) 
           workfft(i*2  ,j,k,:,:) = aimag(temp33_Complex)
         enddo; enddo; enddo
       else                                                                                                       ! use precalculated gamma-operator
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
           do m = 1,3; do n = 1,3
             temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n,:,:) *(workfft(i*2-1,j,k,:,:)&
                                                                 + workfft(i*2  ,j,k,:,:)*img))
           enddo; enddo
           workfft(i*2-1,j,k,:,:) = real (temp33_Complex)
           workfft(i*2  ,j,k,:,:) = aimag(temp33_Complex) 
         enddo; enddo; enddo
       endif

! average strain 
       workfft(1,1,1,:,:) = defgrad_av - math_I3                                        ! zero frequency (real part)
       workfft(2,1,1,:,:) = 0.0_pReal                                                   ! zero frequency (imaginary part)
       
       call dfftw_execute_dft_c2r(plan_fft(2),workfft,workfft)
       defgrad = defgrad + workfft(1:resolution(1),:,:,:,:)*wgt
       do m = 1,3; do n = 1,3
         defgrad_av(m,n) = sum(defgrad(:,:,:,m,n))*wgt
         defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + (defgradAim(m,n) - defgrad_av(m,n))  ! anticipated target minus current state
       enddo; enddo

       print '(2(a,E10.5)/)', 'Error Divergence = ',err_div,    ', Tol. = ', err_div_tol
       
     enddo    ! end looping when convergency is achieved 
     
     c_prev = c_current*wgt        ! calculate stiffness for next step
     
     if (mod(step,bc_frequency(loadcase)) == 0_pInt) &                          ! at output frequency
       write(538)  materialpoint_results(:,1,:)                                 ! write result to file
     if(err_div<err_div_tol .and. err_stress<err_stress_tol) then
       print '(2(A,I5.5),A,/)', '== Step = ',step, ' of Loadcase = ',loadcase, ' Converged =============='
     else
       print '(2(A,I5.5),A,/)', '== Step = ',step, ' of Loadcase = ',loadcase, ' NOT Converged =========='
       not_converged_counter = not_converged_counter + 1
     endif
   enddo  ! end looping over steps in current loadcase
   deallocate(c_reduced)
   deallocate(s_reduced)
   enddo    ! end looping over loadcases
   print '(A,/)', '############################################################'
   print '(a,i5.5,a)', 'A Total of ', not_converged_counter, ' Steps did not Converge!'
 close(538)
 call dfftw_destroy_plan(plan_fft(1)); call dfftw_destroy_plan(plan_fft(2))
    
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
