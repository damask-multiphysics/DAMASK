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
!##############################################################
!* $Id: mpie_spectral2d.f90 665 2010-10-13 16:04:44Z MPIE\m.diehl $
!********************************************************************
! Material subroutine for BVP solution using spectral method
! Version for 2D Problems, suitable for testing and implementation of new features
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Not working by now due to changes on other routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 use numerics, only: relevantStrain, rTol_crystalliteStress

 implicit none
 include 'fftw3.f' !header file for fftw3 (declaring variables). Library file is also needed

! variables to read from loadcase and mesh file
 real(pReal), dimension(9) ::                      valuevector           ! stores information temporarily from loadcase file
 integer(pInt), parameter ::                       maxNchunksInput = 24  ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunksInput*2) :: posInput
 integer(pInt), parameter ::                       maxNchunksMesh = 7    ! 4 identifiers, 3 values
 integer(pInt), dimension (1+2*maxNchunksMesh) ::  posMesh
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

! variables storing information from mesh file
 real(pReal) wgt
 real(pReal), dimension(3) ::  meshdimension
 integer(pInt) homog, prodnn
 integer(pInt), dimension(3) :: resolution

! stress etc.
 real(pReal), dimension(3,3) ::                         ones, zeroes, temp33_Real, damper,&
                                                        pstress, pstress_av, cstress_av, defgrad_av,&
                                                        defgradAim, defgradAimOld, defgradAimCorr, defgradAimCorrPrev,&
                                                        mask_stress, mask_defgrad
 real(pReal), dimension(3,3,3,3) ::                     dPdF, c0, s0
 real(pReal), dimension(6) ::                           cstress                                ! cauchy stress in Mandel notation
 real(pReal), dimension(6,6) ::                         dsde, c066, s066                         
 real(pReal), dimension(:,:,:,:), allocatable ::        defgradold, cstress_field, defgrad 
 complex(pReal), dimension(:,:,:,:), allocatable ::     pstress_field, tau
 complex(pReal), dimension(:,:), allocatable ::         ddefgrad       
 
! variables storing information for spectral method
 complex(pReal), dimension(:,:,:,:), allocatable ::     workfft
 complex(pReal), dimension(3,3) ::                      temp33_Complex
 real(pReal), dimension(3,3) ::                         xinormdyad
 real(pReal), dimension(:,:,:,:,:,:), allocatable ::    gamma_hat
 real(pReal), dimension(:,:,:), allocatable ::          xi
 integer(pInt), dimension(2) ::                         k_s
 integer*8, dimension(3,3,3) ::                         plan_fft

! convergence etc.
 real(pReal) err_div, err_stress, err_defgrad
 real(pReal) err_div_tol, err_stress_tol, err_stress_tolrel, sigma0
 integer(pInt) itmax, ierr
 logical errmatinv

! loop variables etc.
 real(pReal) guessmode                                 ! flip-flop to guess defgrad fluctuation field evolution
 integer(pInt)  i, j, k, l, m, n, p
 integer(pInt)  loadcase, ielem, iter, calcmode, CPFEM_mode
 
 real(pReal) temperature                               ! not used, but needed for call to CPFEM_general

!gmsh output
 character(len=1024) :: nriter
 character(len=1024) :: nrstep
 character(len=1024) :: nrloadcase
 real(pReal), dimension(:,:,:), allocatable ::    displacement
!gmsh output

!Initializing
 bc_maskvector = ''
 unit = 234_pInt

 ones   = 1.0_pReal
 zeroes = 0.0_pReal
 
 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_n = 0_pInt

 resolution = 1_pInt; meshdimension = 0.0_pReal
 
 err_div_tol = 1.0e-4
 itmax = 250_pInt
 err_stress_tolrel=0.01
 
 temperature = 300.0_pReal

 gotResolution =.false.; gotDimension =.false.; gotHomogenization = .false.

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
     call IO_error(46,ext_msg = path)                     !error message for incomplete input file

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
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,j+k) ! assign values for the velocity gradient matrix
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
             resolution(1) = IO_intValue(line,posMesh,i+1)
           case('b')
             resolution(2) = IO_intValue(line,posMesh,i+1)
           case('c')
             resolution(3) = IO_intValue(line,posMesh,i+1)
         end select
       enddo
   end select
   if (gotDimension .and. gotHomogenization .and. gotResolution) exit
   if (resolution(3) /=1) exit
 enddo
 100 close(unit)

 print '(a,/,i4,i4,i4)','resolution a b c', resolution
 print '(a,/,f6.1,f6.1,f6.1)','dimension x y z', meshdimension
 print *,'homogenization',homog
 print *, ''
 
 allocate (workfft(resolution(1),resolution(2),3,3));              workfft              = 0.0_pReal
 allocate (gamma_hat(resolution(1),resolution(2),3,3,3,3));        gamma_hat            = 0.0_pReal
 allocate (xi(resolution(1),resolution(2),3));                     xi                   = 0.0_pReal
 allocate (pstress_field(resolution(1),resolution(2),3,3));        pstress_field        = 0.0_pReal
 allocate (cstress_field(resolution(1),resolution(2),3,3));        cstress_field        = 0.0_pReal
 allocate (tau(resolution(1),resolution(2),3,3));                  tau                  = 0.0_pReal
 allocate (displacement(resolution(1),resolution(2),3));         displacement         = 0.0_pReal
 allocate (defgrad(resolution(1),resolution(2),3,3));              defgrad              = 0.0_pReal
 allocate (defgradold(resolution(1),resolution(2),3,3));           defgradold           = 0.0_pReal
 allocate (ddefgrad(resolution(1),resolution(2)));                ddefgrad              = 0.0_pReal
 
! Initialization of fftw (see manual on fftw.org for more details) 
 call dfftw_init_threads(ierr)
 call dfftw_plan_with_nthreads(4)
 do m = 1,3; do n = 1,3
   call dfftw_plan_dft_2d(plan_fft(1,m,n),resolution(1),resolution(2),& 
                    cstress_field(:,:,m,n), workfft(:,:,m,n), FFTW_PATIENT, FFTW_FORWARD) !only for calculation of div (P)
   call dfftw_plan_dft_2d(plan_fft(2,m,n),resolution(1),resolution(2),& 
                    tau(:,:,m,n), workfft(:,:,m,n), FFTW_PATIENT, FFTW_FORWARD)                    
   call dfftw_plan_dft_2d(plan_fft(3,m,n),resolution(1),resolution(2),& 
                    workfft(:,:,m,n), ddefgrad(:,:), FFTW_PATIENT, FFTW_BACKWARD)
 enddo; enddo
 
 prodnn = resolution(1)*resolution(2)
 wgt = 1_pReal/real(prodnn, pReal)
 defgradAim = math_I3
 defgradAimOld = math_I3
 defgrad_av = math_I3
! Initialization of CPFEM_general (= constitutive law) and of deformation gradient field
 ielem = 0_pInt
 c066 = 0.0_pReal
 do j = 1, resolution(3); do i = 1, resolution(2)
   defgradold(i,j,:,:) = math_I3                    !no deformation at the beginning
   defgrad(i,j,:,:) = math_I3 
   ielem = ielem +1                                                             
   call CPFEM_general(2,math_I3,math_I3,temperature,0.0_pReal,ielem,1_pInt,cstress,dsde,pstress,dPdF)
   c066 = c066 + dsde
 enddo; enddo
   c066 = c066 * wgt
   c0 = math_mandel66to3333(c066)
   call math_invert(6, c066, s066,i, errmatinv)
   s0 = math_mandel66to3333(s066)
 
!calculation of xinormdyad (to calculate gamma_hat) and xi (waves, for proof of equilibrium)
 
 do j = 1, resolution(2)
   k_s(2) = j-1
   if(j > resolution(2)/2+1) k_s(2) = k_s(2)-resolution(2)     
   do i = 1, resolution(1)
     k_s(1) = i-1
     if(i > resolution(1)/2+1) k_s(1) = k_s(1)-resolution(1)
                           xi(i,j,3) = 0.0_pReal
                           xi(i,j,2) = real(k_s(2), pReal)/meshdimension(2)
                           xi(i,j,1) = real(k_s(1), pReal)/meshdimension(1)
       if (any(xi(i,j,:) /= 0.0_pReal)) then     
         do l = 1,2; do m = 1,2
            xinormdyad(l,m) = xi(i,j, l)*xi(i,j, m)/sum(xi(i,j,:)**2)
         enddo; enddo
       else
         xinormdyad = 0.0_pReal
       endif
       temp33_Real = math_mul3333xx33(c0, xinormdyad)            
       temp33_Real = math_inv3x3(temp33_Real) 
       do l=1,3; do m=1,3; do n=1,3; do p=1,3
         gamma_hat(i,j, l,m,n,p) = - temp33_Real(l,n) * xinormdyad(m,p)
       ! gamma_hat(i,j,k, l,m,n,p) = - 0.5_pReal * temp33_Real(l,n) * xinormdyad(m,p)& ! symmetrization????????
                          !            - 0.5_pReal * temp33_Real(m,n) * xinormdyad(l,p)
       enddo; enddo; enddo; enddo     
 enddo; enddo
 
 open(539,file='stress-strain.out')
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
     
     do j = 1, resolution(2); do i = 1, resolution(1)
       temp33_Real = defgrad(i,j,:,:)
       defgrad(i,j,:,:) = defgrad(i,j,:,:)&         ! old fluctuations as guess for new step, no fluctuations for new loadcase
                          + guessmode * (defgrad(i,j,:,:) - defgradold(i,j,:,:))&                                           
                          + (1.0_pReal-guessmode) * math_mul33x33(bc_velocityGrad(:,:,loadcase),defgradold(i,j,:,:))*timeinc
       defgradold(i,j,:,:) = temp33_Real                                                                  
     enddo; enddo

     guessmode = 1.0_pReal                             ! keep guessing along former trajectory during same loadcase
     calcmode = 0_pInt                                 ! start calculation of BC fullfillment
     CPFEM_mode = 1_pInt                               ! winding forward
     iter = 0_pInt
     err_div= 2_pReal * err_div_tol                    ! go into loop 
     defgradAimCorr = 0.0_pReal                        ! reset damping calculation
     damper = damper * 0.9_pReal

!*************************************************************
! convergence loop
     do while( iter <= itmax .and. &     
             (err_div > err_div_tol .or. &
              err_stress > err_stress_tol))    
       iter = iter + 1
       print '(A,I5.5,tr2,A,I5.5)', ' Step = ',steps,'Iteration = ',iter
!*************************************************************
       
! adjust defgrad to fulfill BCs 
       select case (calcmode)         
       case (0)                                                                      
         print *, 'Update Stress Field (constitutive evaluation P(F))'  
         ielem = 0_pInt
         do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1
           call CPFEM_general(3, defgradold(i,j,:,:), defgrad(i,j,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo

         ielem = 0_pInt         
         do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1
           call CPFEM_general(CPFEM_mode,&                                  ! first element in first iteration retains CPFEM_mode 1, 
                              defgradold(i,j,:,:), defgrad(i,j,:,:),&   ! others get 2 (saves winding forward effort)
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
           CPFEM_mode = 2_pInt
           pstress_field(i,j,:,:) = pstress 
           cstress_field(i,j,:,:) = math_mandel6to33(cstress)          
         enddo; enddo
         
         do m = 1,3; do n = 1,3
           pstress_av(m,n) = sum(pstress_field(:,:,m,n)) * wgt
           cstress_av(m,n) = sum(cstress_field(:,:,m,n)) * wgt
           defgrad_av(m,n) = sum(defgrad(:,:,m,n)) * wgt
         enddo; enddo
         
         err_stress = maxval(abs(mask_stress * (cstress_av - bc_stress(:,:,loadcase))))
         err_stress_tol = maxval(abs(cstress_av))*err_stress_tolrel                               
         
         print*, 'Correcting deformation gradient to fullfill BCs'
         defgradAimCorrPrev = defgradAimCorr
         defgradAimCorr     = -mask_stress * math_mul3333xx33(s0, (mask_stress*(cstress_av - bc_stress(:,:,loadcase))))

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
           defgrad(:,:,m,n) = defgrad(:,:,m,n) + (defgradAim(m,n) - defgrad_av(m,n)) !anticipated target minus current state
         enddo; enddo
         err_div = 2 * err_div_tol 
         err_defgrad = maxval(abs(mask_defgrad * (defgrad_av - defgradAim)))         
         print '(a,/,3(3(f12.7,x)/))', ' Deformation Gradient:   ',defgrad_av(1:3,:)               
         print '(a,/,3(3(f10.4,x)/))', ' Cauchy Stress [MPa]: ',cstress_av(1:3,:)/1.e6        
         print '(a,E8.2)', ' error defgrad           ',err_defgrad
         print '(2(a,E8.2))', ' error stress            ',err_stress,'  Tol. = ', err_stress_tol*0.8   
         if(err_stress < err_stress_tol*0.8) then 
           calcmode = 1
         endif
             
! Using the spectral method to calculate the change of deformation gradient, check divergence of stress field in fourier space        
       case (1)   
         print *, 'Update Stress Field (constitutive evaluation P(F))'
         ielem = 0_pInt
         do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1
           call CPFEM_general(3, defgradold(i,j,:,:), defgrad(i,j,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
         enddo; enddo

         ielem = 0_pInt                 
         do j = 1, resolution(2); do i = 1, resolution(1)
           ielem = ielem + 1
           call CPFEM_general(2,&    
                              defgradold(i,j,:,:), defgrad(i,j,:,:),&
                              temperature,timeinc,ielem,1_pInt,&
                              cstress,dsde, pstress, dPdF)
           pstress_field(i,j,:,:) = pstress 
           cstress_field(i,j,:,:) = math_mandel6to33(cstress)          
         enddo; enddo
         
         do j = 1, resolution(2); do i = 1, resolution(1)
           tau(i,j,:,:) = cstress_field(i,j,:,:) - math_mul3333xx33(c0, defgrad(i,j,:,:)-math_I3)  
         enddo; enddo           
         
         print *, 'Calculating equilibrium using spectral method'
         err_div = 0.0_pReal; sigma0 = 0.0_pReal
         do m = 1,3; do n = 1,3
           call dfftw_execute_dft(plan_fft(1,m,n), cstress_field(:,:,m,n),workfft(:,:,m,n))
           if(n==3) sigma0 = max(sigma0, sum(abs(workfft(1,1,m,:))))                     ! L infinity Norm of stress tensor 
         enddo; enddo

        do j = 1, resolution(2); do i = 1, resolution(1)
           err_div = err_div + (maxval(abs(math_mul33x3_complex(workfft(i,j,:,:),xi(i,j,:))))) ! L infinity Norm of div(stress)
         enddo; enddo
         err_div = err_div/real(prodnn, pReal)/sigma0       !weighting of error
         
         do m = 1,3; do n = 1,3
           call dfftw_execute_dft(plan_fft(2,m,n), tau(:,:,m,n), workfft(:,:,m,n))
         enddo; enddo
         
         do j = 1, resolution(2); do i = 1, resolution(1)
           temp33_Complex = 0.0_pReal
           do m = 1,3; do n = 1,3
             temp33_Complex(m,n) = sum(gamma_hat(i,j,m,n,:,:) * workfft(i,j,:,:))
           enddo; enddo
           workfft(i,j,:,:) = temp33_Complex(:,:) 
         enddo; enddo
         workfft(1,1,:,:) = zeroes!????        
         
         do m = 1,3; do n = 1,3
            call dfftw_execute_dft(plan_fft(3,m,n), workfft(:,:,m,n),ddefgrad(:,:))
            defgrad(:,:,m,n) = defgrad_av(m,n) + real(ddefgrad, pReal) * wgt
            pstress_av(m,n) = sum(pstress_field(:,:,m,n))*wgt
            cstress_av(m,n) = sum(cstress_field(:,:,m,n))*wgt
            defgrad_av(m,n) = sum(defgrad(:,:,m,n))*wgt
            defgrad(:,:,m,n) = defgrad(:,:,m,n) + (defgradAim(m,n) - defgrad_av(m,n)) !anticipated target minus current state
         enddo; enddo
         
         err_stress = maxval(abs(mask_stress * (cstress_av - bc_stress(:,:,loadcase))))
         err_stress_tol = maxval(abs(cstress_av))*err_stress_tolrel                             !accecpt relativ error specified
            
         print '(2(a,E8.2))', ' error divergence        ',err_div,'  Tol. = ', err_div_tol
         print '(2(a,E8.2))', ' error stress            ',err_stress,'  Tol. = ', err_stress_tol           
         if(err_stress > err_stress_tol .and. err_div < err_div_tol) then         ! change to calculation of BCs, reset damper etc.
           calcmode = 0
           defgradAimCorr = 0.0_pReal
           damper = damper * 0.9_pReal
         endif     
       end select 
     enddo    ! end looping when convergency is achieved 

     write(539,'(E12.6,a,E12.6)'),defgrad_av(3,3)-1,'   ',  cstress_av(3,3)
     print '(a,/,3(3(f12.7,x)/))', ' Deformation Gradient:   ',defgrad_av(1:3,:)
     print *, ''
     print '(a,/,3(3(f10.4,x)/))', ' Cauchy Stress [MPa]: ',cstress_av(1:3,:)/1.e6
     print '(A)', '************************************************************'

! Postprocessing (gsmh output)

     temp33_Real(1,:) = 0.0_pReal; temp33_Real(1,3) = -(real(resolution(3))/meshdimension(3)) ! start just below origin
     k=1
     do j = 1, resolution(2); do i = 1, resolution(1)
       if((j==1).and.(i==1)) then
         temp33_Real(1,:) = temp33_Real(1,:) + math_mul33x3(defgrad(i,j,:,:),&
                          (/0.0_pReal,0.0_pReal,(real(resolution(3))/meshdimension(3))/))
         temp33_Real(2,:) = temp33_Real(1,:)
         temp33_Real(3,:) = temp33_Real(1,:)
         displacement(i,j,:) = temp33_Real(1,:)
       else 
         if(i==1) then
           temp33_Real(2,:) = temp33_Real(2,:) + math_mul33x3(defgrad(i,j,:,:),&
                     (/0.0_pReal,(real(resolution(2))/meshdimension(2)),0.0_pReal/))
           temp33_Real(3,:) = temp33_Real(2,:)
           displacement(i,j,:) = temp33_Real(2,:)
         else   
           temp33_Real(3,:) = temp33_Real(3,:) + math_mul33x3(defgrad(i,j,:,:),&
                     (/(real(resolution(1))/meshdimension(1)),0.0_pReal,0.0_pReal/))
           displacement(i,j,:) = temp33_Real(3,:)   
         endif
       endif
     enddo; enddo

     write(nrloadcase, *) loadcase; write(nriter, *) iter; write(nrstep, *) steps
     open(589,file = 'stress' //trim(adjustl(nrloadcase))//'-'//trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh')
     open(588,file = 'disgrad'//trim(adjustl(nrloadcase))//'-'//trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh')
     write(589, '(4(A, /), I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn
     write(588, '(4(A, /), I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn

     ielem = 0_pInt
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       ielem = ielem + 1
       write(589, '(I10, 3(tr2, E12.6))'), ielem, displacement(i,j,:) !for deformed configuration
       write(588, '(I10, 3(tr2, E12.6))'), ielem, displacement(i,j,:)        
   !    write(589, '(4(I10,tr2))'), ielem, i-1,j-1,k-1 !for undeformed configuration
   !    write(588, '(4(I10,tr2))'), ielem, i-1,j-1,k-1 
     enddo; enddo; enddo

     write(589, '(2(A, /), I10)'), '$EndNodes', '$Elements', prodnn
     write(588, '(2(A, /), I10)'), '$EndNodes', '$Elements', prodnn

     do i = 1, prodnn
       write(589, '(I10, A, I10)'), i, ' 15 2    1     2', i
       write(588, '(I10, A, I10)'), i, ' 15 2    1     2', i
     enddo

     write(589, '(A)'), '$EndElements'
     write(588, '(A)'), '$EndElements'
     write(589, '(8(A, /), I10)'), '$NodeData', '1','"'//trim(adjustl('stress'//trim(adjustl(nrloadcase))//'-'//&
                            trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh'))//'"','1','0.0', '3', '0', '9', prodnn
     write(588, '(8(A, /), I10)'), '$NodeData', '1','"'//trim(adjustl('disgrad'//trim(adjustl(nrloadcase))//'-'//&
                            trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh'))//'"','1','0.0', '3', '0', '9', prodnn
     ielem = 0_pInt
     do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
            ielem = ielem + 1
             write(589, '(i10, 9(tr2, E14.8))'), ielem, cstress_field(i,j,:,:)
             write(588, '(i10, 9(tr2, E14.8))'), ielem, defgrad(i,j,:,:) - math_I3
     enddo; enddo; enddo

     write(589, *), '$EndNodeData'
     write(588, *), '$EndNodeData'
     close(589); close(588) 
   enddo  ! end looping over steps in current loadcase
 enddo    ! end looping over loadcases
close(539)

do i=1,3; do m = 1,3; do n = 1,3
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