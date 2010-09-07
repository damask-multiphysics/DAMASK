!* $Id$
!********************************************************************
! Material subroutine for BVP solution using spectral method
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts
!            D.D. Tjahjanto
!            C. Kords
!            M. Diehl
!            R. Lebensohn
!
! MPI fuer Eisenforschung, Duesseldorf
!
!********************************************************************
!     Usage:
!             - start program with mpie_spectral PathToMeshFile/NameOfMesh.mesh
!               PathToLoadFile/NameOfLoadFile.load
!             - PathToLoadFile will be the working directory
!             - make sure the file "material.config" exists in the working
!               directory
!********************************************************************
program mpie_spectral
!********************************************************************

 use mpie_interface
 use prec, only: pInt, pReal
 use IO
 use math
 use CPFEM, only: CPFEM_general

 implicit none
 include 'fftw3.f' !header file for fftw3 (declaring variables). Library file is also needed

 !variables to read in from loadcase and mesh file
 real(pReal), dimension(9) ::                      valuevector           ! stores information temporarily from loadcase file
 integer(pInt), parameter ::                       maxNchunksInput = 24  ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunksInput*2) :: posInput
 integer(pInt), parameter ::                       maxNchunksMesh = 7    ! 4 identifiers, 3 values
 integer(pInt), dimension (1+2*maxNchunksMesh) ::  posMesh
 integer(pInt) unit, N_l, N_s, N_t, N_n                                  ! numbers of identifiers
 logical gotResolution,gotDimension,gotHomogenization
 logical, dimension(9) :: bc_maskvector
 character(len=1024) path, line

! variables storing information from loadcase file
 integer(pInt)                                  N_Loadcases, steps
 integer(pInt), dimension(:), allocatable ::    bc_steps              ! number of steps
 real(pReal)                                    timeinc
 real(pReal), dimension (:,:,:), allocatable :: bc_velocityGrad, &
                                                bc_stress             ! velocity gradient and stress BC
 real(pReal), dimension(:), allocatable ::      bc_timeIncrement      ! length of increment
 logical, dimension(:,:,:,:), allocatable ::    bc_mask               ! mask of boundary conditions

! variables storing information from mesh file
 integer(pInt) homog, prodnn
 real(pReal) wgt
 integer(pInt), dimension(3) :: resolution
 real(pReal), dimension(3) ::  meshdimension

! stress etc.
 real(pReal), dimension(6) ::                           cstress                   ! cauchy stress in Mandel notation (not needed)
 real(pReal), dimension(3,3) ::                         pstress                   ! Piola-Kirchhoff stress in Matrix notation
 real(pReal), dimension(3,3,3,3) ::                     dPdF, c0, s0              ! ??, reference stiffnes, compliance
 real(pReal), dimension(6,6) ::                         dsde, s066            
 real(pReal), dimension(3,3) ::                         defgradmacro
 real(pReal), dimension(3,3) ::                         pstress_av, defgrad_av, temp33_Real
 real(pReal), dimension(:,:,:,:,:), allocatable ::      pstress_field, defgrad, defgradold, ddefgrad

! variables storing information for spectral method
 complex(pReal), dimension(:,:,:,:,:), allocatable ::   workfft
 complex(pReal), dimension(3,3) ::                      temp33_Complex
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable ::  gamma_hat
 real(pReal), dimension(:,:,:,:,:), allocatable ::      xknormdyad
 real(pReal), dimension(3) ::                           xk
 integer(pInt), dimension(3) ::                         k_s
 integer*8, dimension(2,3,3) ::                         plan_fft

! convergency etc.
 logical errmatinv
 integer(pInt) itmax, ierr
 real(pReal) error, err_stress_av, err_stress_max, err_strain_av, err_strain_max
 real(pReal), dimension(3,3) :: strain_err, pstress_err
 
! loop variables etc.
 integer(pInt)  i, j, k, l, m, n, p
 integer(pInt)  loadcase, ielem, iter, calcmode
 real(pReal) guessmode                                             ! flip-flop to guess defgrad fluctuation field evolution
 
 real(pReal) temperature                                           ! not used, but needed

!gmsh output
 character(len=1024) :: nriter
 character(len=1024) :: nrstep
!gmsh output

!Initializing
 bc_maskvector = ''
 unit = 234_pInt

 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_n = 0_pInt

 pstress_err = .0_pReal; strain_err = .0_pReal
 resolution = 1_pInt; meshdimension = .0_pReal

 error = 0.001_pReal
 itmax = 50_pInt

 temperature = 300.0_pReal

 gotResolution =.false.; gotDimension =.false.; gotHomogenization = .false.

 if (IargC() < 2) call IO_error(102)                     ! check for correct number of arguments given

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

! allocate memory depending on lines in input file
101 N_Loadcases = N_l

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
   posInput = IO_stringPos(line,maxNchunksInput)
   do j = 1,maxNchunksInput,2
     select case (IO_lc(IO_stringValue(line,posInput,j)))
       case('l','velocitygrad')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,j+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,j+k)  ! assign values for the velocity gradient matrix
         enddo
         bc_mask(:,:,1,i) = reshape(bc_maskvector,(/3,3/))
         bc_velocityGrad(:,:,i) = reshape(valuevector,(/3,3/))
       case('s','stress')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,j+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,j+k)  ! assign values for the bc_stress matrix
         enddo
         bc_mask(:,:,2,i) = reshape(bc_maskvector,(/3,3/))
         bc_stress(:,:,i) = reshape(valuevector,(/3,3/))
       case('t','time','delta')                                            ! increment time
           bc_timeIncrement(i) = IO_floatValue(line,posInput,j+1)
       case('n','incs','increments','steps')                               ! bc_steps
           bc_steps(i) = IO_intValue(line,posInput,j+1)
     end select
 enddo; enddo

200 close(unit)

 do i = 1, N_Loadcases
   if (any(bc_mask(:,:,1,i) == bc_mask(:,:,2,i))) call IO_error(47,i)     ! bc_mask consistency
   print '(a,/,3(3(f12.6,x)/))','L',bc_velocityGrad(:,:,i)
   print '(a,/,3(3(f12.6,x)/))','bc_stress',bc_stress(:,:,i)
   print '(a,/,3(3(l,x)/))','bc_mask for velocitygrad',bc_mask(:,:,1,i)
   print '(a,/,3(3(l,x)/))','bc_mask for stress',bc_mask(:,:,2,i)
   print *,'time',bc_timeIncrement(i)
   print *,'incs',bc_steps(i)
   print *, ''
 enddo

!read header of mesh file to get the information needed before the complete mesh file is intepretated by mesh.f90
 path = getSolverJobName()
 print*,'JobName: ',trim(path)
 if (.not. IO_open_file(unit,trim(path)//InputFileExtension)) call IO_error(101,ext_msg = path)

 rewind(unit)
 do
   read(unit,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posMesh = IO_stringPos(line,maxNchunksMesh)

   select case ( IO_lc(IO_StringValue(line,posMesh,1)) )
     case ('dimension')
         gotDimension = .true.
         do i = 2,6,2
           select case (IO_lc(IO_stringValue(line,posMesh,i)))
             case('x')
                 meshdimension(1) = IO_floatValue(line,posMesh,i+1)
             case('y')
                 meshdimension(2) = IO_floatValue(line,posMesh,i+1)
             case('z')
                 meshdimension(3) = IO_floatValue(line,posMesh,i+1)
           end select
         enddo
     case ('homogenization')
         gotHomogenization = .true.
         homog = IO_intValue(line,posMesh,2)
     case ('resolution')
         gotResolution = .true.
         do i = 2,6,2
           select case (IO_lc(IO_stringValue(line,posMesh,i)))
             case('a')
                 resolution(1) = 2**IO_intValue(line,posMesh,i+1)
             case('b')
                 resolution(2) = 2**IO_intValue(line,posMesh,i+1)
             case('c')
                 resolution(3) = 2**IO_intValue(line,posMesh,i+1)
           end select
         enddo
   end select
   if (gotDimension .and. gotHomogenization .and. gotResolution) exit
 enddo
 100 close(unit)

 print '(a,/,i3,i3,i3)','resolution a b c', resolution
 print '(a,/,f6.2,f6.2,f6.2)','dimension x y z', meshdimension
 print *,'homogenization',homog
 print *, ''

 allocate (workfft(resolution(1)/2+1,resolution(2),resolution(3),3,3));          workfft              = .0_pReal
 allocate (gamma_hat(resolution(1)/2+1,resolution(2),resolution(3),3,3,3,3));    gamma_hat            = .0_pReal
 allocate (xknormdyad(resolution(1)/2+1,resolution(2),resolution(3),3,3));       xknormdyad           = .0_pReal
 allocate (pstress_field(resolution(1),resolution(2),resolution(3),3,3));        pstress_field        = .0_pReal
 allocate (defgrad(resolution(1),resolution(2),resolution(3),3,3));              defgrad              = .0_pReal
 allocate (defgradold(resolution(1),resolution(2),resolution(3),3,3));           defgradold           = .0_pReal
 allocate (ddefgrad(resolution(1),resolution(2),resolution(3),3,3));            ddefgrad              = .0_pReal
 
 call dfftw_init_threads(ierr)
 call dfftw_plan_with_nthreads(4)
 do m = 1,3; do n = 1,3
   call dfftw_plan_dft_r2c_3d(plan_fft(1,m,n),resolution(1),resolution(2),resolution(3),& 
                    pstress_field(:,:,:,m,n), workfft(:,:,:,m,n), FFTW_PATIENT, FFTW_DESTROY_INPUT)
   call dfftw_plan_dft_c2r_3d(plan_fft(2,m,n),resolution(1),resolution(2),resolution(3),& 
                    workfft(:,:,:,m,n), ddefgrad(:,:,:,m,n), FFTW_PATIENT)
 enddo; enddo
 
 prodnn = resolution(1)*resolution(2)*resolution(3)
 wgt = 1._pReal/real(prodnn, pReal)
 defgradmacro = math_I3
 
 c0 = .0_pReal
 ielem = 0_pInt
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   defgradold(i,j,k,:,:) = math_I3                                             !no deformation at the beginning
   defgrad(i,j,k,:,:) = math_I3 
   ielem = ielem +1                                                             
   call CPFEM_general(2,math_I3,math_I3,temperature,0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
   c0 = c0 + dPdF
 enddo; enddo; enddo

 call math_invert(6,math_Mandel3333to66(c0),s066,i,errmatinv) !i is just a dummy variable
 if(errmatinv) call IO_error(45,ext_msg = "problem in c0 inversion")  ! todo: change number and add message to io.f90 (and remove No. 48)
 s0 = math_Mandel66to3333(s066)*real(prodnn, pReal)
   
!calculation of xknormdyad (needed to calculate gamma_hat)
 do k = 1, resolution(3)
   k_s(3) = k-1
   if(k > resolution(3)/2+1) k_s(3) = k_s(3)-resolution(3)
     xk(3) = .0_pReal
     if(resolution(3) > 1) xk(3) = real(k_s(3), pReal)/meshdimension(3)
     do j = 1, resolution(2)
       k_s(2) = j-1
       if(j > resolution(2)/2+1) k_s(2) = k_s(2)-resolution(2)
       xk(2) = real(k_s(2), pReal)/meshdimension(2)
       do i = 1, resolution(1)/2+1
         k_s(1) = i-1
         xk(1) = real(k_s(1), pReal)/meshdimension(1)

         if (any(xk /= .0_pReal)) then
           do l = 1,3; do m = 1,3
             xknormdyad(i,j,k, l,m) = xk(l)*xk(m)/(xk(1)**2+xk(2)**2+xk(3)**2)
           enddo; enddo
         endif

 enddo; enddo; enddo
 
! Initialization done
 open(539,file='stress-strain.out')

!*************************************************************
!Loop over loadcases defined in the loadcase file
 do loadcase = 1, N_Loadcases
!*************************************************************

   timeinc = bc_timeIncrement(loadcase)/bc_steps(loadcase)
   guessmode = 0.0_pReal                                                            ! change of load case

!*************************************************************
! loop oper steps defined in input file for current loadcase
   do steps = 1, bc_steps(loadcase)
!*************************************************************
     write(*,*) '***************************************************'
     write(*,*) 'STEP = ',steps

     defgradmacro = defgradmacro& 
                  + math_mul33x33(bc_velocityGrad(:,:,loadcase), defgradmacro)*timeinc   !update macroscopic displacement gradient (stores the desired BCs of defgrad) 
  
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       temp33_Real = defgrad(i,j,k,:,:)
       defgrad(i,j,k,:,:) = defgrad(i,j,k,:,:)& 
                          + guessmode  * (defgrad(i,j,k,:,:) - defgradold(i,j,k,:,:))&                                           ! old fluctuations as guess for new step
                          + (1.0_pReal-guessmode) * math_mul33x33(bc_velocityGrad(:,:,loadcase),defgradold(i,j,k,:,:))*timeinc   ! no fluctuations for new loadcase  
       defgradold(i,j,k,:,:) = temp33_Real                                                                  
     enddo; enddo; enddo

     guessmode = 1.0_pReal                                                          ! keep guessing along former trajectory
     calcmode = 1_pInt
     iter = 0_pInt
     err_stress_av = 2.*error; err_strain_av = 2.*error

!*************************************************************
! convergency loop
     do while((iter <= itmax).and.((err_stress_av > error).or.(err_strain_av > error)))
       iter = iter + 1
       write(*,*) 'ITER = ',iter
!*************************************************************
       err_strain_av = .0_pReal; err_stress_av = .0_pReal
       err_strain_max = .0_pReal; err_stress_max = .0_pReal
       pstress_av = .0_pReal; defgrad_av=.0_pReal

       write(*,*) 'UPDATE STRESS FIELD'
       ielem = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1
         call CPFEM_general(3, defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                            temperature,timeinc,ielem,1_pInt,&
                            cstress,dsde, pstress, dPdF)
       enddo; enddo; enddo
   
       c0 = .0_pReal
       l = 0_pInt
       ielem = 0_pInt      
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1
         call CPFEM_general(calcmode,&    ! first element in first iteration retains calcMode 1, others get 2 (saves winding forward effort)
                            defgradold(i,j,k,:,:), defgrad(i,j,k,:,:),&
                            temperature,timeinc,ielem,1_pInt,&
                            cstress,dsde, pstress, dPdF)
         calcmode = 2
         c0 = c0 + dPdF
         temp33_Real = pstress

         do m = 1,3; do n = 1,3                                                       ! calculate stress error
           if(abs(temp33_Real(m,n)) > 0.1_pReal * abs(pstress_err(m,n))) then         ! only stress components larger than 10% are taking under consideration
             err_stress_av  = err_stress_av + abs((pstress_field(i,j,k,m,n)-temp33_Real(m,n))/temp33_Real(m,n))
             err_stress_max = max(err_stress_max, abs((pstress_field(i,j,k,m,n)-temp33_Real(m,n))/temp33_Real(m,n)))
             l=l+1
           endif
         enddo; enddo
         pstress_field(i,j,k,:,:) = temp33_Real
         pstress_av = pstress_av + temp33_Real          ! average stress
       enddo; enddo; enddo
  
       err_stress_av = err_stress_av/l        ! do the weighting of the error
       pstress_av = pstress_av*wgt            ! do the weighting of average stress
       pstress_err = pstress_av 

       if(iter==1) then                        !update gamma_hat with new reference stiffness
         c0 = c0 *wgt
         do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
           temp33_Real = .0_pReal
           do l = 1,3; do m = 1,3; do n = 1,3; do p = 1,3
             temp33_Real(l,m) = temp33_Real(l,m)+c0(l,n,m,p)*xknormdyad(i,j,k, n,p)
           enddo; enddo; enddo; enddo
           temp33_Real = math_inv3x3(temp33_Real)
           do l=1,3; do m=1,3; do n=1,3; do p=1,3
             gamma_hat(i,j,k, l,m,n,p) = -temp33_Real(l,n)*xknormdyad(i,j,k, m,p)
           enddo; enddo; enddo; enddo
         enddo; enddo; enddo
       endif  

       write(*,*) 'SPECTRAL METHOD TO GET CHANGE OF DEFORMATION GRADIENT FIELD'
       do m = 1,3; do n = 1,3
         call dfftw_execute_dft_r2c(plan_fft(1,m,n), pstress_field(:,:,:,m,n),workfft(:,:,:,m,n))
       enddo; enddo 
       
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)/2+1
         temp33_Complex = .0_pReal
         do m = 1,3; do n = 1,3
           temp33_Complex(m,n) = sum(gamma_hat(i,j,k,m,n,:,:) * workfft(i,j,k,:,:))
         enddo; enddo
         workfft(i,j,k,:,:) = temp33_Complex(:,:) 
       enddo; enddo; enddo
  
       do m = 1,3; do n = 1,3
         call dfftw_execute_dft_c2r(plan_fft(2,m,n), workfft(:,:,:,m,n),ddefgrad(:,:,:,m,n))
       enddo; enddo
       
       ddefgrad = ddefgrad * wgt
       defgrad = defgrad + ddefgrad

       l = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         defgrad_av(:,:) = defgrad_av(:,:) + defgrad(i,j,k,:,:)       ! calculate average strain
         do m = 1,3; do n = 1,3                                       ! calculate strain error
           if(abs(defgrad(i,j,k,m,n)) > 0.1 * abs(strain_err(m,n))) then
             err_strain_av = err_strain_av + abs(real(ddefgrad(i,j,k,m,n), pReal)/defgrad(i,j,k,m,n))
             err_strain_max = max(err_strain_max, abs(real(ddefgrad(i,j,k,m,n), pReal)/defgrad(i,j,k,m,n)))
             l=l+1
           endif
         enddo; enddo
       enddo; enddo; enddo

       err_strain_av = err_strain_av/l                 ! weight by number of non-zero strain components
       defgrad_av = defgrad_av * wgt                   ! weight by number of points
       strain_err = defgrad_av
       
       do m = 1,3; do n = 1,3
         if(bc_mask(m,n,1,loadcase)) then              ! adjust defgrad to fulfill displacement BC (defgradmacro)
           defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + (defgradmacro(m,n)-defgrad_av(m,n))
         else                                          ! adjust defgrad to fulfill stress BC
           defgrad(:,:,:,m,n) = defgrad(:,:,:,m,n) + sum( s0(m,n,:,:)*(bc_stress(:,:,loadcase)-pstress_av(:,:)), &
                                                          mask = bc_mask(:,:,2,loadcase) ) !works at the moment only for 0 Stress as BC
         endif
       enddo; enddo

       write(*,*) 'STRESS FIELD ERROR AV  = ',err_stress_av
       write(*,*) 'STRAIN FIELD ERROR AV  = ',err_strain_av
       write(*,*) 'STRESS FIELD ERROR MAX = ',err_stress_max
       write(*,*) 'STRAIN FIELD ERROR MAX = ',err_strain_max

     enddo    ! end looping when convergency is achieved
     write(539,'(f12.6,a,f12.6)'),defgrad_av(3,3)-1,'   ',pstress_av(3,3)
     write(*,*) 'U11 U22 U33'
     write(*,*) defgrad_av(1,1)-1,defgrad_av(2,2)-1,defgrad_av(3,3)-1
     write(*,*) 'U11/U33'
     write(*,*) (defgrad_av(1,1)-1)/(defgrad_av(3,3)-1)
     write(*,*) 'S11 S22 S33'
     write(*,*) pstress_av(1,1),pstress_av(2,2),pstress_av(3,3)

!gsmh output
     write(nriter, *) iter
     write(nrstep, *) steps
     nrstep = 'stress'//trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh'
     open(589,file = nrstep)
     write(589, '(A, /, A, /, A, /, A, /, I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn
     ielem = 0_pInt
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       ielem = ielem + 1
       write(589, '(I10,f16.8,tr2,f16.8,tr2,f16.8)'), ielem,math_mul33x3(defgrad(i,j,k,:,:),real((/i, j, k/), pReal))
     enddo; enddo; enddo
     write(589, '(A, /, A, /, I10)'), '$EndNodes', '$Elements', prodnn
     do i = 1, prodnn
       write(589, '(I10, A, I10)'), i, ' 15 2    1     2', i
     enddo
     write(589, '(A)'), '$EndElements'
     write(589, '(A, /, A, /, A, /, A, /, A, /, A, /, A, /, A, /, I10)'), '$NodeData', '1',&
                             '"'//trim(adjustl(nrstep))//'"', '1','0.0', '3', '0', '9', prodnn
     ielem = 0_pInt
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
            ielem = ielem + 1
            write(589, '(i10,f16.8,tr2,f16.8,tr2,f16.8,tr2,f16.8,tr2,f16.8,tr2,f16.8,&
                                  tr2,f16.8,tr2,f16.8,tr2,f16.8,tr2)'), ielem, pstress_field(i,j,k,:,:)
     enddo; enddo; enddo
     write(589, *), '$EndNodeData'
     close(589)
!end gmsh

   enddo  ! end looping over steps in current loadcase
 enddo    ! end looping over loadcases
close(539)

do i=1,2; do m = 1,3; do n = 1,3
  call dfftw_destroy_plan(plan_fft(i,m,n))
enddo; enddo; enddo

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