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
!             - start program with mpie_spectral PathToGeomFile/NameOfGeom.geom
!               PathToLoadFile/NameOfLoadFile.load
!             - PathToGeomFile will be the working directory
!             - make sure the file "material.config" exists in the working
!               directory. For further configuration use "numerics.config"
!********************************************************************
program mpie_spectral
!********************************************************************

 use mpie_interface
 use prec, only: pInt, pReal
 use IO
 use math
 use CPFEM, only: CPFEM_general, CPFEM_initAll 
 use numerics, only: err_div_tol, err_stress_tol, err_stress_tolrel, err_defgrad_tol,&
                     itmax, fast_execution, mpieNumThreadsInt       
 use homogenization, only: materialpoint_sizeResults, materialpoint_results
!$ use OMP_LIB                                                          ! the openMP function library

 implicit none
 include 'fftw3.f' !header file for fftw3 (declaring variables). Library files are also needed
! compile FFTW 3.2.2 with ./configure --enable-threads
! variables to read from loadcase and geom file
 real(pReal), dimension(9) ::                      valuevector           ! stores information temporarily from loadcase file
 integer(pInt), parameter ::                       maxNchunksInput = 24  ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunksInput*2) :: posInput
 integer(pInt), parameter ::                       maxNchunksGeom = 7    ! 4 identifiers, 3 values
 integer(pInt), dimension (1+2*maxNchunksGeom) ::  posGeom
 integer(pInt) unit, N_l, N_s, N_t, N_n                                  ! numbers of identifiers
 character(len=1024) path, line
 logical gotResolution,gotDimension,gotHomogenization
 logical, dimension(9) :: bc_maskvector

! variables storing information from loadcase file
 real(pReal)                                    timeinc
 real(pReal), dimension (:,:,:), allocatable :: bc_velocityGrad, &
                                                bc_stress             ! velocity gradient and stress BC
 real(pReal), dimension(:), allocatable ::      bc_timeIncrement      ! length of increment
 integer(pInt)                                  N_Loadcases, steps
 integer(pInt), dimension(:), allocatable ::    bc_steps              ! number of steps
 logical, dimension(:,:,:,:), allocatable ::    bc_mask               ! mask of boundary conditions

! variables storing information from geom file
 real(pReal) wgt
 real(pReal), dimension(3) ::  geomdimension
 integer(pInt) homog
 integer(pInt), dimension(3) :: resolution

! stress etc.
 real(pReal), dimension(3,3) ::                         ones, zeroes, temp33_Real, damper,&
                                                        pstress, pstress_av, cstress_av, defgrad_av,&
                                                        defgradAim, defgradAimOld, defgradAimCorr, defgradAimCorrPrev,&
                                                        mask_stress, mask_defgrad                          
 real(pReal), dimension(3,3,3,3) ::                     dPdF, c0, s0
 real(pReal), dimension(6) ::                           cstress                                ! cauchy stress in Mandel notation
 real(pReal), dimension(6,6) ::                         dsde, c066, s066                       ! Mandel notation of 4th order tensors
 real(pReal), dimension(:,:,:,:,:), allocatable ::      workfft, defgrad, defgradold
 
! variables storing information for spectral method
 complex(pReal) ::                                      img
 complex(pReal), dimension(3,3) ::                      temp33_Complex
 real(pReal), dimension(3,3) ::                         xinormdyad
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat
 real(pReal), dimension(3) ::                           xi, xi_middle
 integer(pInt), dimension(3) ::                         k_s
 integer*8, dimension(2) ::                             plan_fft
 
! loop variables, convergence etc.
 real(pReal) guessmode, err_div, err_stress, err_defgrad, sigma0        
 integer(pInt)  i, j, k, l, m, n, p
 integer(pInt)  loadcase, ielem, iter, calcmode, CPFEM_mode, ierr
 logical errmatinv
 
 real(pReal) temperature                               ! not used, but needed for call to CPFEM_general

!Initializing
!$ call omp_set_num_threads(mpieNumThreadsInt)         ! set number of threads for parallel execution set by MPIE_NUM_THREADS

 bc_maskvector = ''
 unit = 234_pInt
 ones   = 1.0_pReal; zeroes = 0.0_pReal
 img = cmplx(0.0,1.0)
 
 N_l = 0_pInt; N_s = 0_pInt
 N_t = 0_pInt; N_n = 0_pInt
 gotResolution =.false.; gotDimension =.false.; gotHomogenization = .false.
 resolution = 1_pInt; geomdimension = 0.0_pReal
 
 temperature = 300.0_pReal

 if (IargC() /= 2) call IO_error(102)                     ! check for correct number of given arguments

! Reading the loadcase file and assign variables
 path = getLoadcaseName()
 print*,'Loadcase: ',trim(path)
 print*,'Workingdir: ',trim(getSolverWorkingDirectoryName())

 if (.not. IO_open_file(unit,path)) call IO_error(45,ext_msg = path)
 
 rewind(unit)
 do
   read(unit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posInput = IO_stringPos(line,maxNchunksInput)
   do i = 1, maxNchunksInput, 1
       select case (IO_lc(IO_stringValue(line,posInput,i)))
            case('l','velocitygrad')
                 N_l = N_l+1
            case('s','stress')
                 N_s = N_s+1
            case('t','time','delta')
                 N_t = N_t+1
            case('n','incs','increments','steps')
                 N_n = N_n+1
        end select
   enddo                                                  ! count all identifiers to allocate memory and do sanity check
   if ((N_l /= N_s).or.(N_s /= N_t).or.(N_t /= N_n)) &    ! sanity check
     call IO_error(46,ext_msg = path)                     ! error message for incomplete input file
 enddo

101 N_Loadcases = N_l

! allocate memory depending on lines in input file
 allocate (bc_velocityGrad(3,3,N_Loadcases));       bc_velocityGrad = 0.0_pReal
 allocate (bc_stress(3,3,N_Loadcases));             bc_stress = 0.0_pReal
 allocate (bc_mask(3,3,2,N_Loadcases));             bc_mask = .false.
 allocate (bc_timeIncrement(N_Loadcases));          bc_timeIncrement = 0.0_pReal
 allocate (bc_steps(N_Loadcases));                  bc_steps = 0_pInt

 rewind(unit)
 i = 0_pInt
 do
   read(unit,'(a1024)',END = 200) line
   if (IO_isBlank(line)) cycle                           ! skip empty lines
   i = i + 1
   posInput = IO_stringPos(line,maxNchunksInput)         ! ToDo: Add error message for case that information is not complete
   do j = 1,maxNchunksInput,2
     select case (IO_lc(IO_stringValue(line,posInput,j)))
       case('l','velocitygrad')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,j+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,j+k) ! assign values for the velocity gradient matrix
         enddo
         bc_mask(:,:,1,i) = transpose(reshape(bc_maskvector,(/3,3/)))
         bc_velocityGrad(:,:,i) = math_transpose3x3(reshape(valuevector,(/3,3/)))
       case('s','stress')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,j+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,j+k)  ! assign values for the bc_stress matrix
         enddo
         bc_mask(:,:,2,i) = transpose(reshape(bc_maskvector,(/3,3/)))
         bc_stress(:,:,i) = math_transpose3x3(reshape(valuevector,(/3,3/)))
       case('t','time','delta')                                            ! increment time
           bc_timeIncrement(i) = IO_floatValue(line,posInput,j+1)
       case('n','incs','increments','steps')                               ! bc_steps
           bc_steps(i) = IO_intValue(line,posInput,j+1)
     end select
 enddo; enddo

200 close(unit)

 do i = 1, N_Loadcases
   if (any(bc_mask(:,:,1,i) == bc_mask(:,:,2,i))) call IO_error(47,i)     ! bc_mask consistency
   print '(a,/,3(3(f12.6,x)/))','L'        ,math_transpose3x3(bc_velocityGrad(:,:,i))
   print '(a,/,3(3(f12.6,x)/))','bc_stress',math_transpose3x3(bc_stress(:,:,i))
   print '(a,/,3(3(l,x)/))',    'bc_mask for velocitygrad',transpose(bc_mask(:,:,1,i))
   print '(a,/,3(3(l,x)/))',    'bc_mask for stress'      ,transpose(bc_mask(:,:,2,i))
   print *,'time',bc_timeIncrement(i)
   print *,'incs',bc_steps(i)
   print *, ''
 enddo

!read header of geom file to get the information needed before the complete geom file is intepretated by mesh.f90
 path = getSolverJobName()
 print*,'JobName: ',trim(path)
 if (.not. IO_open_file(unit,trim(path)//InputFileExtension)) call IO_error(101,ext_msg = path)

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
 
 if(mod(resolution(1),2)/=0 .or. mod(resolution(2),2)/=0 .or. mod(resolution(3),2)/=0)  call IO_error(102)  !!ToDo: add correct error to IO
 
 print '(a,/,i4,i4,i4)','resolution a b c', resolution
 print '(a,/,f6.1,f6.1,f6.1)','dimension x y z', geomdimension
 print *,'homogenization',homog
 
 allocate (defgrad   (resolution(1),resolution(2),resolution(3),3,3)); defgrad    = 0.0_pReal
 allocate (defgradold(resolution(1),resolution(2),resolution(3),3,3)); defgradold = 0.0_pReal
 
 wgt = 1.0_pReal/real(resolution(1)*resolution(2)*resolution(3), pReal)
 defgradAim = math_I3
 defgradAimOld = math_I3
 defgrad_av = math_I3
 
! Initialization of CPFEM_general (= constitutive law) and of deformation gradient field
 call CPFEM_initAll(temperature,1_pInt,1_pInt)
 ielem = 0_pInt
 c066 = 0.0_pReal 
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   defgradold(i,j,k,:,:) = math_I3                    !no deformation at the beginning
   defgrad(i,j,k,:,:) = math_I3 
   ielem = ielem +1 
   call CPFEM_general(2,math_I3,math_I3,temperature,0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
   c066 = c066 + dsde
 enddo; enddo; enddo
 c066 = c066 * wgt
 c0 = math_mandel66to3333(c066)
 call math_invert(6, c066, s066,i, errmatinv)
 s0 = math_mandel66to3333(s066)
 
!calculation of calculate gamma_hat field in the case of fast execution (needs a lot of memory) 
 if(fast_execution) then
   allocate (gamma_hat(resolution(1)/2+1,resolution(2),resolution(3),3,3,3,3)); gamma_hat = 0.0_pReal
   do k = 1, resolution(3)
     k_s(3) = k-1
     if(k > resolution(3)/2+1) k_s(3) = k_s(3)-resolution(3)
     do j = 1, resolution(2)
       k_s(2) = j-1
       if(j > resolution(2)/2+1) k_s(2) = k_s(2)-resolution(2)     
       do i = 1, resolution(1)/2+1
         k_s(1) = i-1
         xi(3) = 0.0_pReal     !for the 2D case
         if(resolution(3) > 1) xi(3) = real(k_s(3), pReal)/geomdimension(3) !3D case
                               xi(2) = real(k_s(2), pReal)/geomdimension(2)
                               xi(1) = real(k_s(1), pReal)/geomdimension(1)
         if (any(xi /= 0.0_pReal)) then     
           do l = 1,3; do m = 1,3
              xinormdyad(l,m) = xi(l)*xi(m)/sum(xi**2)
           enddo; enddo
         else
           xinormdyad = 0.0_pReal
         endif 
         temp33_Real = math_mul3333xx33(c0, xinormdyad)            
         temp33_Real = math_inv3x3(temp33_Real) 
         do l=1,3; do m=1,3; do n=1,3; do p=1,3
           gamma_hat(i,j,k, l,m,n,p) = - (0.5*temp33_Real(l,n)+0.5*temp33_Real(n,l)) *&
                                         (0.5*xinormdyad(m,p)+0.5*xinormdyad(p,m))
         enddo; enddo; enddo; enddo         
   enddo; enddo; enddo
 else ! or allocate just one fourth order tensor
   allocate (gamma_hat(1,1,1,3,3,3,3)); gamma_hat = 0.0_pReal
 endif
 
! calculate xi for the calculation of divergence in Fourier space (middle frequency)
 xi_middle(3) = 0.0_pReal !2D case
 if(resolution(3) > 1) xi_middle(3) = real(resolution(3)/2, pReal)/geomdimension(3) !3D case
                       xi_middle(2) = real(resolution(2)/2, pReal)/geomdimension(2)
                       xi_middle(1) = real(resolution(1)/2, pReal)/geomdimension(1)
 
 allocate (workfft(resolution(1)+2,resolution(2),resolution(3),3,3)); workfft = 0.0_pReal

! Initialization of fftw (see manual on fftw.org for more details) 
 call dfftw_init_threads(ierr) !toDo: add error code
 call dfftw_plan_with_nthreads(mpieNumThreadsInt) 
 call dfftw_plan_many_dft_r2c(plan_fft(1),3,(/resolution(1),resolution(2),resolution(3)/),9,&
   workfft,(/resolution(1)  +2,resolution(2),resolution(3)/),1,(resolution(1)  +2)*resolution(2)*resolution(3),&
   workfft,(/resolution(1)/2+1,resolution(2),resolution(3)/),1,(resolution(1)/2+1)*resolution(2)*resolution(3),FFTW_PATIENT)   
 call dfftw_plan_many_dft_c2r(plan_fft(2),3,(/resolution(1),resolution(2),resolution(3)/),9,&
   workfft,(/resolution(1)/2+1,resolution(2),resolution(3)/),1,(resolution(1)/2+1)*resolution(2)*resolution(3),&
   workfft,(/resolution(1)  +2,resolution(2),resolution(3)/),1,(resolution(1)  +2)*resolution(2)*resolution(3),FFTW_PATIENT)

! write header of output file
 open(538,file=trim(getSolverWorkingDirectoryName())//trim(getSolverJobName())&
                                               //'_'//trim(getLoadcase())&
                                                    //'.spectralOut',form='UNFORMATTED')       
 write(538), 'load',trim(getLoadcaseName())
 write(538), 'workingdir',trim(getSolverWorkingDirectoryName())
 write(538), 'geometry',trim(getSolverJobName())//InputFileExtension
 write(538), 'resolution',resolution
 write(538), 'dimension',geomdimension
 write(538), 'materialpoint_sizeResults', materialpoint_sizeResults
 write(538), 'increments', sum(bc_steps)
 write(538), 'eoh'
 write(538)  materialpoint_results(:,1,:) 
 write(538)  materialpoint_results(:,1,:) !to be conform with t16 Marc format
! Initialization done

!*************************************************************
!Loop over loadcases defined in the loadcase file
 do loadcase = 1, N_Loadcases
!*************************************************************

   timeinc = bc_timeIncrement(loadcase)/bc_steps(loadcase)
   guessmode = 0.0_pReal                                   ! change of load case, homogeneous guess for the first step
  
   mask_defgrad =  merge(ones,zeroes,bc_mask(:,:,1,loadcase))
   mask_stress =  merge(ones,zeroes,bc_mask(:,:,2,loadcase))
   damper = ones/10
!*************************************************************
! loop oper steps defined in input file for current loadcase
   do steps = 1, bc_steps(loadcase)
!*************************************************************
     temp33_Real = defgradAim
     defgradAim = defgradAim &                        ! update macroscopic displacement gradient (defgrad BC)
                  + guessmode * mask_stress * (defgradAim - defgradAimOld) &
                  + math_mul33x33(bc_velocityGrad(:,:,loadcase), defgradAim)*timeinc
     defgradAimOld = temp33_Real
     
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       temp33_Real = defgrad(i,j,k,:,:)
       defgrad(i,j,k,:,:) = defgrad(i,j,k,:,:)&         ! old fluctuations as guess for new step, no fluctuations for new loadcase
                          + guessmode * (defgrad(i,j,k,:,:) - defgradold(i,j,k,:,:))&                                           
                          + (1.0_pReal-guessmode) * math_mul33x33(bc_velocityGrad(:,:,loadcase),defgradold(i,j,k,:,:))*timeinc
       defgradold(i,j,k,:,:) = temp33_Real                                                                  
     enddo; enddo; enddo

     guessmode = 1.0_pReal                              ! keep guessing along former trajectory during same loadcase
     if(all(bc_mask(:,:,1,loadcase))) then
       calcmode = 1_pInt                                ! if no stress BC is given (calmode 0 is not needed)
     else
       calcmode = 0_pInt                                ! start calculation of BC fulfillment
     endif
     CPFEM_mode = 1_pInt                                ! winding forward
     iter = 0_pInt
     err_div= 2_pReal * err_div_tol                     ! go into loop 
     defgradAimCorr = 0.0_pReal                         ! reset damping calculation
     damper = damper * 0.9_pReal
    
!*************************************************************
! convergence loop
     do while(iter < itmax .and. &     
             (err_div > err_div_tol .or. &
              err_stress > err_stress_tol .or. &
              err_defgrad > err_defgrad_tol))    
       iter = iter + 1_pInt
       print*, ' '
       print '(3(A,I5.5,tr2))', ' Loadcase = ',loadcase, ' Step = ',steps,'Iteration = ',iter
       cstress_av = 0.0_pReal
       workfft = 0.0_pReal !needed because of the padding for FFTW
!*************************************************************
       
! adjust defgrad to fulfill BCs 
       select case (calcmode)       
       case (0)                                                                      
         print *, 'Update Stress Field (constitutive evaluation P(F))'  
         ielem = 0_pInt
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1
           call CPFEM_general(3, defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo; enddo
         
         ielem = 0_pInt       
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(CPFEM_mode,&                                  ! first element in first iteration retains CPFEM_mode 1, 
                              defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&   ! others get 2 (saves winding forward effort)
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
           CPFEM_mode = 2_pInt
           workfft(i,j,k,:,:) = pstress
           cstress_av = cstress_av + math_mandel6to33(cstress)          
         enddo; enddo; enddo
         
         cstress_av = cstress_av * wgt
         do m = 1,3; do n = 1,3
           pstress_av(m,n) = sum(workfft(1:resolution(1),:,:,m,n)) * wgt   
           defgrad_av(m,n) = sum(defgrad(:,:,:,m,n)) * wgt
         enddo; enddo

         err_stress = maxval(abs(mask_stress * (pstress_av - bc_stress(:,:,loadcase))))
         err_stress_tol = maxval(abs(pstress_av))*err_stress_tolrel                               
         
         print*, 'Correcting deformation gradient to fullfill BCs'
         defgradAimCorrPrev = defgradAimCorr
         defgradAimCorr     = -mask_stress * math_mul3333xx33(s0, (mask_stress*(pstress_av - bc_stress(:,:,loadcase))))

         do m=1,3; do n =1,3                                        ! calculate damper (correction is far to strong)
           if ( sign(1.0_pReal,defgradAimCorr(m,n))/=sign(1.0_pReal,defgradAimCorrPrev(m,n))) then
             damper(m,n) = max(0.01_pReal,damper(m,n)*0.8)
           else
             damper(m,n) = min(1.0_pReal,damper(m,n) *1.2)
           endif
         enddo; enddo
         defgradAimCorr = mask_Stress*(damper * defgradAimCorr)
         defgradAim = defgradAim + defgradAimCorr                                                                  

         do m = 1,3; do n = 1,3                        
           defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + (defgradAim(m,n) - defgrad_av(m,n)) !anticipated target minus current state
         enddo; enddo
         err_div = 2 * err_div_tol 
         err_defgrad = maxval(abs(mask_defgrad * (defgrad_av - defgradAim)))         
         print '(a,/,3(3(f12.7,x)/))', ' Deformation Gradient:   ',math_transpose3x3(defgrad_av)                           
         print '(a,/,3(3(f10.4,x)/))', ' Cauchy Stress [MPa]: '   ,math_transpose3x3(cstress_av)/1.e6        
         print '(2(a,E8.2))', ' error stress               ',err_stress, '  Tol. = ', err_stress_tol 
         print '(2(a,E8.2))', ' error deformation gradient ',err_defgrad,'  Tol. = ', err_defgrad_tol*0.8  
         if(err_stress < err_stress_tol*0.8) then 
           calcmode = 1_pInt
         endif
             
! Using the spectral method to calculate the change of deformation gradient, check divergence of stress field in fourier space        
       case (1)
         print *, 'Update Stress Field (constitutive evaluation P(F))'
         ielem = 0_pInt
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(3, defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo; enddo
         ielem = 0_pInt                 
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1_pInt
           call CPFEM_general(2,&    
                              defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
           workfft(i,j,k,:,:) = pstress 
           cstress_av = cstress_av + math_mandel6to33(cstress)          
         enddo; enddo; enddo       
         cstress_av = cstress_av * wgt
         do m = 1,3; do n = 1,3
           pstress_av(m,n) = sum(workfft(1:resolution(1),:,:,m,n))*wgt
         enddo; enddo
         
         print *, 'Calculating equilibrium using spectral method'
         err_div = 0.0_pReal; sigma0 = 0.0_pReal
         call dfftw_execute_dft_r2c(plan_fft(1),workfft,workfft)       ! FFT of pstress
         do m = 1,3 ! L infinity Norm of stress tensor 
           sigma0 = max(sigma0, sum(abs(workfft(1,1,1,m,:) + (workfft(2,1,1,m,:))*img)))     
         enddo
         err_div = (maxval(abs(math_mul33x3_complex(workfft(resolution(1)+1,resolution(2)/2+1,resolution(3)/2+1,:,:)+& ! L infinity Norm of div(stress)
                                                    workfft(resolution(1)+2,resolution(2)/2+1,resolution(3)/2+1,:,:)*img,xi_middle))))
         err_div = err_div/sigma0     !weighting of error

         if(fast_execution) then     ! fast execution with stored gamma_hat
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
           do m = 1,3; do n = 1,3
             temp33_Complex(m,n) = sum(gamma_hat(i,j,k, m,n,:,:) *(workfft(i*2-1,j,k,:,:)&                 
                                                                 + workfft(i*2  ,j,k,:,:)*img))
           enddo; enddo
           workfft(i*2-1,j,k,:,:) = real (temp33_Complex) 
           workfft(i*2  ,j,k,:,:) = aimag(temp33_Complex) 
          enddo; enddo; enddo
          
         else ! memory saving version, in-time calculation of gamma_hat
         do k = 1, resolution(3)
           k_s(3) = k-1
           if(k > resolution(3)/2+1) k_s(3) = k_s(3)-resolution(3)
           do j = 1, resolution(2)
             k_s(2) = j-1
             if(j > resolution(2)/2+1) k_s(2) = k_s(2)-resolution(2)
             do i = 1, resolution(1)/2+1
               k_s(1) = i-1
               xi(3) = 0.0_pReal     !for the 2D case
               if(resolution(3) > 1) xi(3) = real(k_s(3), pReal)/geomdimension(3) !3D case
                                     xi(2) = real(k_s(2), pReal)/geomdimension(2)
                                     xi(1) = real(k_s(1), pReal)/geomdimension(1)
               if (any(xi(:) /= 0.0_pReal)) then     
                 do l = 1,3; do m = 1,3
                    xinormdyad(l,m) = xi(l)*xi( m)/sum(xi**2)
                 enddo; enddo
               else
                 xinormdyad = 0.0_pReal
               endif 
               temp33_Real = math_mul3333xx33(c0, xinormdyad)            
               temp33_Real = math_inv3x3(temp33_Real) 
               do l=1,3; do m=1,3; do n=1,3; do p=1,3
                gamma_hat(1,1,1, l,m,n,p) = - (0.5*temp33_Real(l,n)+0.5*temp33_Real(n,l))*&
                                              (0.5*xinormdyad(m,p) +0.5*xinormdyad(p,m))
               enddo; enddo; enddo; enddo   
               do m = 1,3; do n = 1,3
                 temp33_Complex(m,n) = sum(gamma_hat(1,1,1,m,n,:,:) *(workfft(i*2-1,j,k,:,:)&                 
                                                                     +workfft(i*2  ,j,k,:,:)*img))
               enddo; enddo
               workfft(i*2-1,j,k,:,:) = real (temp33_Complex) 
               workfft(i*2  ,j,k,:,:) = aimag(temp33_Complex)
           enddo; enddo; enddo
         endif
         
         workfft(1,1,1,:,:) = defgrad_av - math_I3 !zero frequency (real part)
         workfft(2,1,1,:,:) = 0.0_pReal            !zero frequency (imaginary part)
         
         call dfftw_execute_dft_c2r(plan_fft(2),workfft,workfft)
         defgrad = defgrad + workfft(1:resolution(1),:,:,:,:)*wgt
         do m = 1,3; do n = 1,3
           defgrad_av(m,n) = sum(defgrad(:,:,:,m,n))*wgt
           defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + (defgradAim(m,n) - defgrad_av(m,n)) !anticipated target minus current state
         enddo; enddo
         
         err_stress = maxval(abs(mask_stress * (pstress_av - bc_stress(:,:,loadcase))))
         err_stress_tol = maxval(abs(pstress_av))*err_stress_tolrel                             !accecpt relativ error specified
         err_defgrad = maxval(abs(mask_defgrad * (defgrad_av - defgradAim)))  
         
         print '(2(a,E8.2))', ' error divergence           ',err_div,    '  Tol. = ', err_div_tol
         print '(2(a,E8.2))', ' error stress               ',err_stress, '  Tol. = ', err_stress_tol 
         print '(2(a,E8.2))', ' error deformation gradient ',err_defgrad,'  Tol. = ', err_defgrad_tol 

         if((err_stress > err_stress_tol .or. err_defgrad > err_defgrad_tol) .and. err_div < err_div_tol) then  ! change to calculation of BCs, reset damper etc.
           calcmode = 0_pInt
           defgradAimCorr = 0.0_pReal
           damper = damper * 0.9_pReal
         endif     
       end select 
     enddo    ! end looping when convergency is achieved 
     
     write(538)  materialpoint_results(:,1,:) !write to output file
     
     print '(a,x,f12.7)'         , ' Determinant of Deformation Aim:', math_det3x3(defgradAim)
     print '(a,/,3(3(f12.7,x)/))', ' Deformation Aim:        ',math_transpose3x3(defgradAim)
     print '(a,/,3(3(f12.7,x)/))', ' Deformation Gradient:   ',math_transpose3x3(defgrad_av) 
     print '(a,/,3(3(f10.4,x)/))', ' Cauchy Stress [MPa]:    ',math_transpose3x3(cstress_av)/1.e6
     print '(A)', '************************************************************'
   enddo  ! end looping over steps in current loadcase
 enddo    ! end looping over loadcases
 close(538)
 call dfftw_destroy_plan(plan_fft(1)); call dfftw_destroy_plan(plan_fft(2))
    
end program mpie_spectral

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