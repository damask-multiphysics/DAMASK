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
 include 'fftw3.f' !header file for fftw3 (declaring variables) library is also needed

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
 integer(pInt), dimension (3,3) ::              bc_stress_i           ! conversion from bc_mask
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
 real(pReal), dimension(6) ::           cstress                   ! cauchy stress in Mandel notation
 real(pReal), dimension (3,3) ::        pstress                   ! Piola-Kirchhoff stress in Matrix notation
 real(pReal), dimension (3,3,3,3) ::    dPdF,c0,s0                ! ??, reference stiffnes, (reference stiffness)^-1
 real(pReal), dimension(6,6) ::         dsde,c066,s066            ! mandel notation
 real(pReal), dimension(3,3) ::         disgradmacro
 real(pReal), dimension(3,3) ::         cstress_av, defgrad_av, aux33
 real(pReal), dimension(:,:,:,:,:), allocatable :: cstress_field, defgrad, defgradold, start

! variables storing information for spectral method
 complex(pReal), dimension(:,:,:,:,:), allocatable :: workfft
 complex(pReal), dimension(3,3) :: ddefgrad
 real(pReal), dimension(:,:,:,:,:,:,:), allocatable :: gamma_hat
 real(pReal), dimension(3) :: xk
 real(pReal), dimension(3,3) :: xknormdyad
 integer(pInt), dimension(3) :: k_s
 integer*8, dimension(2) :: plan

! convergency etc.
 logical errmatinv
 integer(pInt) itmax, ierr
 real(pReal) error, err_stress_av, err_stress_max, err_strain_av, err_strain_max
 real(pReal), dimension(3,3) :: strain_err, cstress_err

! loop variables etc.
 integer(pInt)  i, j, k, l, m, n, p
 integer(pInt)  loadcase, ielem, ial, iter, calcmode

 real(pReal) temperature                                           ! not used, but needed

!gmsh
 character(len=1024) :: nriter
 character(len=1024) :: nrstep
!gmsh

!Initializing
 bc_maskvector = ''
 unit = 234_pInt

 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_n = 0_pInt

 disgradmacro = .0_pReal
 c0 = .0_pReal; c066 = .0_pReal
 s0 = .0_pReal; s066 = .0_pReal
 cstress_err = .0_pReal; strain_err = .0_pReal
 
 cstress = .0_pReal
 dsde = .0_pReal

 resolution = 1_pInt; meshdimension = .0_pReal

 error = 0.001_pReal
 itmax = 50_pInt

 temperature = 300.0_pReal

 gotResolution =.false.; gotDimension =.false.; gotHomogenization = .false.

 if (IargC() < 2) call IO_error(102)                     ! check for correct Nr. of arguments given

! Reading the loadcase file and assingn variables
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

! consistency checks
 do i = 1, N_Loadcases
   if (any(bc_mask(:,:,1,i) == bc_mask(:,:,2,i))) &
     call IO_error(47,i)                                                  ! bc_mask consistency
   if (any(math_transpose3x3(bc_stress(:,:,i)) + bc_stress(:,:,i) /= 2.0_pReal * bc_stress(:,:,i))) &
     call IO_error(48,i)                                                  ! bc_stress symmetry

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

 allocate (workfft(resolution(1),resolution(2),resolution(3),3,3));                 workfft              = .0_pReal
 allocate (gamma_hat(3,3,3,3,resolution(1),resolution(2),resolution(3)));           gamma_hat            = .0_pReal
 allocate (cstress_field(resolution(1),resolution(2),resolution(3),3,3));           cstress_field        = .0_pReal
 allocate (defgrad(3,3,resolution(1),resolution(2),resolution(3)));                 defgrad              = .0_pReal
 allocate (defgradold(3,3,resolution(1),resolution(2),resolution(3)));              defgradold           = .0_pReal
 allocate (start(3,3,resolution(1),resolution(2),resolution(3)));                   start                = .0_pReal
 
 call dfftw_plan_dft_r2c_3d(plan(1),resolution(1),resolution(2),resolution(3), cstress_field(:,:,:,:,:),workfft(1:resolution(1)/2+1,:,:,:,:), FFTW_PATIENT)
 call dfftw_plan_dft_3d(plan(2), resolution(1),resolution(2),resolution(3), workfft(:,:,:,:,:),workfft(:,:,:,:,:), FFTW_FORWARD, FFTW_PATIENT)
 
 prodnn = resolution(1)*resolution(2)*resolution(3)
 wgt = 1._pReal/real(prodnn, pReal)

 ielem = 0_pInt

 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   defgradold(:,:,i,j,k) = math_I3                                             !to fit calculation of first step to calculation of following steps
   defgrad(:,:,i,j,k) = math_I3                                                !to fit calculation of first step to calculation of following steps
   ielem = ielem +1                                                             !loop over FPs and determine elastic constants of reference material
   call CPFEM_general(2,math_I3,math_I3,temperature,0.0_pReal,ielem,1_pInt,cstress,dsde, pstress, dPdF)
   c066 = c066+dsde
 enddo; enddo; enddo

 c066 = c066*wgt

 call math_invert(6,c066,s066,i,errmatinv) !i is just a dummy variable
 if(errmatinv) call IO_error(45,ext_msg = "problem in c0 inversion")  ! todo: change number and add message to io.f90

 s0 = math_Mandel66to3333(s066)
 c0 = math_Mandel66to3333(c066)

 do k = 1, resolution(3)
   k_s(3) = k-1
   if(k > resolution(3)/2) k_s(3) = k_s(3)-resolution(3)
     xk(3) = .0_pReal
     if(resolution(3) > 1) xk(3) = real(k_s(3), pReal)/meshdimension(3)
     do j = 1, resolution(2)
       k_s(2) = j-1
       if(j > resolution(2)/2) k_s(2) = k_s(2)-resolution(2)
       xk(2) = real(k_s(2), pReal)/meshdimension(2)
       do i = 1, resolution(1)
         k_s(1) = i-1
         if(i >  resolution(1)/2) k_s(1) = k_s(1) -resolution(1)
         xk(1) = real(k_s(1), pReal)/meshdimension(1)

         xknormdyad=.0_pReal

         if (any(xk /= .0_pReal)) then
           do l = 1,3; do m = 1,3
             xknormdyad(l,m) = xk(l)*xk(m)/(xk(1)**2+xk(2)**2+xk(3)**2)
           enddo; enddo
         endif

!forall loops don't work for the next 2 loop constructs!!!
         aux33 = .0_pReal
         do l = 1,3; do m = 1,3; do n = 1,3; do p = 1,3
           aux33(l,m) = aux33(l,m)+c0(l,n,m,p)*xknormdyad(n,p)
         enddo; enddo; enddo; enddo

         aux33 = math_inv3x3(aux33)

         do l=1,3; do m=1,3; do n=1,3; do p=1,3
           gamma_hat(l,m,n,p,i,j,k) = -aux33(l,n)*xknormdyad(m,p)
         enddo; enddo; enddo; enddo

 enddo; enddo; enddo

! Initialization done
 open(539,file='stress-strain.out')
!*************************************************************
!Loop over loadcases defined in the loadcase file
 do loadcase = 1, N_Loadcases
!*************************************************************
   bc_stress_i = 0_pInt  !convert information about stress BC's from bc_mask in an integer-array
   do m = 1,3; do n = 1,3
       if(bc_mask(m,n,2,loadcase)) bc_stress_i(m,n) = 1
   enddo; enddo

   timeinc = bc_timeIncrement(loadcase)/bc_steps(loadcase)

   do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)    !no fluctuation as guess for new loadcase (last two summands will disappear)
     start(:,:,i,j,k) = bc_velocityGrad(:,:,loadcase)*timeinc -defgrad(:,:,i,j,k) + defgradold(:,:,i,j,k)
   enddo; enddo; enddo

!*************************************************************
! loop oper steps defined in input file for current loadcase
   do steps = 1, bc_steps(loadcase)
!*************************************************************
     write(*,*) '***************************************************'
     write(*,*) 'STEP = ',steps

     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       aux33 = defgrad(:,:,i,j,k)
       defgrad(:,:,i,j,k)    = 2 * defgrad(:,:,i,j,k) - defgradold(:,:,i,j,k) + start(:,:,i,j,k)   ! old fluctuations as guess
       defgradold(:,:,i,j,k) = aux33                                                               ! wind forward
     enddo; enddo; enddo

     disgradmacro = disgradmacro + bc_velocityGrad(:,:,loadcase)*timeinc  !update macroscopic displacementgradient (stores the desired BCs of defgrad)
     start = .0_pReal
     calcmode = 1_pInt
     iter = 0_pInt
     err_stress_av = 2.*error; err_strain_av = 2.*error

!*************************************************************
! convergency loop
     do while((iter <= itmax).and.((err_stress_av > error).or.(err_strain_av > error)))
       iter = iter+1
       write(*,*) 'ITER = ',iter
!*************************************************************
       err_strain_av = .0_pReal; err_stress_av = .0_pReal
       err_strain_max = .0_pReal; err_stress_max = .0_pReal
       cstress_av = .0_pReal; defgrad_av=.0_pReal

       write(*,*) 'UPDATE STRESS FIELD'
       ielem = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1
         call CPFEM_general(3, defgradold(:,:,i,j,k), defgrad(:,:,i,j,k),&
                          temperature,timeinc,ielem,1_pInt,&
                          cstress,dsde, pstress, dPdF)
       enddo; enddo; enddo

       l = 0_pInt
       ielem = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ielem = ielem + 1
         call CPFEM_general(calcmode,&    ! first element in first iteration retains calcMode 1, others get 2 (saves winding forward effort)
                             defgradold(:,:,i,j,k), defgrad(:,:,i,j,k),&
                             temperature,timeinc,ielem,1_pInt,&
                             cstress,dsde, pstress, dPdF)
         calcmode = 2

         aux33 = math_Mandel6to33(cstress)

         do m = 1,3; do n = 1,3                                           ! calculate stress error
           if(abs(aux33(m,n)) > 0.1 * abs(cstress_err(m,n))) then         ! only stress components larger than 10% are taking under consideration
             err_stress_av  = err_stress_av + abs((cstress_field(i,j,k,m,n)-aux33(m,n))/aux33(m,n)) !any find maxval> leave loop
             err_stress_max = max(err_stress_max, abs((cstress_field(i,j,k,m,n)-aux33(m,n))/aux33(m,n)))
             l=l+1
           endif
         enddo; enddo
         cstress_field(i,j,k,:,:) = aux33
         cstress_av = cstress_av + aux33          ! average stress
       enddo; enddo; enddo
   
       err_stress_av = err_stress_av/l        ! do the weighting of the error
       cstress_av = cstress_av*wgt            ! do the weighting of average stress
       cstress_err = cstress_av 
    
       write(*,*) 'SPECTRAL METHOD TO GET CHANGE OF DEFORMATION GRADIENT FIELD'
       do m = 1,3; do n = 1,3
         call dfftw_execute_dft_r2c(plan(1), cstress_field(:,:,:,m,n),workfft(1:resolution(1)/2+1,:,:,m,n))
       enddo; enddo 
       workfft=conjg(workfft)
       do i = 0, resolution(1)/2-2 !unpack fft data for conj complex symmetric part. can be directly used in calculation of cstress_field
         m = 1
         do k = 1, resolution(3)
           n = 1
           do j = 1, resolution(2)
             workfft(resolution(1)-i,j,k,:,:) = conjg(workfft(2+i,n,m,:,:))
             if(n == 1) n = resolution(2) +1
             n = n-1
          enddo
          if(m == 1) m = resolution(3) +1
            m = m -1
       enddo; enddo

       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         ddefgrad = .0_pReal
         do m = 1,3; do n = 1,3
           ddefgrad(m,n) = ddefgrad(m,n) +sum(gamma_hat(m,n,:,:,i,j,k)*workfft(i,j,k,:,:))
         enddo; enddo
         workfft(i,j,k,:,:) = ddefgrad(:,:) 
       enddo; enddo; enddo
  
       do m = 1,3; do n = 1,3
          call dfftw_execute_dft(plan(2), workfft(:,:,:,m,n), workfft(:,:,:,m,n))
          defgrad(m,n,:,:,:) = defgrad(m,n,:,:,:) + real(workfft(:,:,:,m,n), pReal)*wgt
       enddo; enddo
   
       l = 0_pInt
       do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
         defgrad_av(:,:) = defgrad_av(:,:) + defgrad(:,:,i,j,k)       ! calculate average strain
         do m = 1,3; do n = 1,3                                       ! calculate strain error
         if(abs(defgrad(m,n,i,j,k)) > 0.1 * abs(strain_err(m,n))) then
             err_strain_av = err_strain_av + abs((real(workfft(i,j,k,m,n), pReal)*wgt)/defgrad(m,n,i,j,k))
             err_strain_max = max(err_strain_max, abs((real(workfft(i,j,k,m,n), pReal)*wgt)/defgrad(m,n,i,j,k)))
             l=l+1
           endif
         enddo; enddo
       enddo; enddo; enddo

       err_strain_av = err_strain_av/l                 ! weight by number of non-zero strain components
       defgrad_av = defgrad_av * wgt                   ! weight by number of points
       strain_err = defgrad_av
       
       do m = 1,3; do n = 1,3
         if(bc_mask(m,n,1,loadcase)) then !adjust defgrad to achieve displacement BC (disgradmacro)
           defgrad(m,n,:,:,:) = defgrad(m,n,:,:,:) + (disgradmacro(m,n)+math_I3(m,n)-defgrad_av(m,n))
         endif
         if(bc_mask(m,n,2,loadcase)) then !adjust defgrad to achieve convergency in stress
           defgrad(m,n,:,:,:) = defgrad(m,n,:,:,:) + sum(s0(m,n,:,:)*bc_stress_i(:,:)*(bc_stress(:,:,loadcase)-cstress_av(:,:)))
         endif
       enddo; enddo

       write(*,*) 'STRESS FIELD ERROR AV  = ',err_strain_av
       write(*,*) 'STRAIN FIELD ERROR AV  = ',err_stress_av
       write(*,*) 'STRESS FIELD ERROR MAX = ',err_strain_max
       write(*,*) 'STRAIN FIELD ERROR MAX = ',err_stress_max

     enddo    ! end looping when convergency is achieved
     write(539,'(f12.6,a,f12.6)'),defgrad_av(3,3)-1,'   ',cstress_av(3,3)
     write(*,*) 'U11 U22 U33'
     write(*,*) defgrad_av(1,1)-1,defgrad_av(2,2)-1,defgrad_av(3,3)-1
     write(*,*) 'U11/U33'
     write(*,*) (defgrad_av(1,1)-1)/(defgrad_av(3,3)-1)
     write(*,*) 'S11 S22 S33'
     write(*,*) cstress_av(1,1),cstress_av(2,2),cstress_av(3,3)

!gsmh output
     write(nriter, *) iter
     write(nrstep, *) steps
     nrstep='defgrad'//trim(adjustl(nrstep))//trim(adjustl(nriter))//'_cpfem.msh'
     open(589,file=nrstep)
     write(589, '(A, /, A, /, A, /, A, /, I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn
     do i = 1, prodnn
       write(589, '(I10, I10, I10, I10)'), i, mod((i-1), resolution(1)) +1, mod(((i-1)/resolution(1)),&
                       resolution(2)) +1, mod(((i-1)/(resolution(1)*resolution(2))), resolution(3)) +1
     enddo
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
                                 tr2,f16.8,tr2,f16.8,tr2,f16.8,tr2)'), ielem, defgrad(:,:,i,j,k)
     enddo; enddo; enddo
     write(589, *), '$EndNodeData'
     close(589)

     write(nriter, *) iter
     write(nrstep, *) steps
     nrstep = 'stress'//trim(adjustl(nrstep))//trim(adjustl(nriter))//'_cpfem.msh'
     open(589,file = nrstep)
     write(589, '(A, /, A, /, A, /, A, /, I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn
     do i = 1, prodnn
       write(589, '(I10, I10, I10, I10)'), i, mod((i-1), resolution(1)) +1, mod(((i-1)/resolution(1)),&
                       resolution(2)) +1, mod(((i-1)/(resolution(1)*resolution(2))), resolution(3)) +1
     enddo
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
                                  tr2,f16.8,tr2,f16.8,tr2,f16.8,tr2)'), ielem, cstress_field(i,j,k,:,:)
     enddo; enddo; enddo
     write(589, *), '$EndNodeData'
     close(589)
!end gmsh

   enddo  ! end loping over steps in current loadcase
 enddo  ! end looping over loadcases
close(539)
call dfftw_destroy_plan(plan(2))
call dfftw_destroy_plan(plan(1))
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