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
!
include "prec.f90"             ! uses nothing else


MODULE mpie_interface
 use prec, only: pInt, pReal
 character(len=64), parameter :: FEsolver = 'Spectral'
 character(len=5),  parameter :: InputFileExtension = '.mesh'

CONTAINS

!********************************************************************
! initialize interface module
!
!********************************************************************
subroutine mpie_interface_init()

 write(6,*)
 write(6,*) '<<<+-  mpie_spectral init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)

 return
endsubroutine

!********************************************************************
! extract working directory from loadcase file
! possibly based on current working dir
!********************************************************************
function getSolverWorkingDirectoryName()

 implicit none

 character(len=1024) cwd,outname,getSolverWorkingDirectoryName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \

 call getarg(2,outname)                                ! path to loadFile
 
 if (scan(outname,pathSep) == 1) then                  ! absolute path given as command line argument
   getSolverWorkingDirectoryName = outname(1:scan(outname,pathSep,back=.true.))
 else
   call getcwd(cwd)
   getSolverWorkingDirectoryName = trim(cwd)//'/'//outname(1:scan(outname,pathSep,back=.true.))
 endif
 
 getSolverWorkingDirectoryName = rectifyPath(getSolverWorkingDirectoryName)
 
 return
 
endfunction

!********************************************************************
! basename of meshfile from command line arguments
!
!********************************************************************
function getSolverJobName()

 use prec, only: pInt

 implicit none

 character(1024) getSolverJobName, outName, cwd
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 integer(pInt) posExt,posSep

 getSolverJobName = ''

 call getarg(1,outName)
 posExt = scan(outName,'.',back=.true.)
 posSep = scan(outName,pathSep,back=.true.)
 
 if (posExt <= posSep) posExt = len_trim(outName)+1       ! no extension present
 getSolverJobName = outName(1:posExt-1)                   ! path to mesh file (excl. extension)
 
 if (scan(getSolverJobName,pathSep) /= 1) then            ! relative path given as command line argument
   call getcwd(cwd)
   getSolverJobName = rectifyPath(trim(cwd)//'/'//getSolverJobName)
 else
   getSolverJobName = rectifyPath(getSolverJobName)
 endif
 
 getSolverJobName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                    getSolverJobName)
 return
endfunction


!********************************************************************
! relative path of loadcase from command line arguments
!
!********************************************************************
function getLoadcaseName()

 use prec, only: pInt

 implicit none

 character(len=1024) getLoadcaseName, outName, cwd
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 integer(pInt) posExt,posSep
 posExt = 0                                                 !not sure if its needed

 call getarg(2,getLoadcaseName)
 posExt = scan(getLoadcaseName,'.',back=.true.)
 posSep = scan(getLoadcaseName,pathSep,back=.true.)
 
 if (posExt <= posSep) getLoadcaseName = trim(getLoadcaseName)//('.load')   ! no extension present
 if (scan(getLoadcaseName,pathSep) /= 1) then          ! relative path given as command line argument
   call getcwd(cwd)
   getLoadcaseName = rectifyPath(trim(cwd)//'/'//getLoadcaseName)
 else
   getLoadcaseName = rectifyPath(getLoadcaseName)
 endif
 
 getLoadcaseName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                    getLoadcaseName)
 return
endfunction


!********************************************************************
! remove ../ and ./ from path
!
!********************************************************************
function rectifyPath(path)

 use prec, only: pInt

 implicit none
 
 character(len=*) path
 character(len=len_trim(path)) rectifyPath
 integer(pInt) i,j,k,l

 !remove ./ from path
 l = len_trim(path)
 rectifyPath = path
 do i = l,2,-1
    if ( rectifyPath(i-1:i) == './' .and. rectifyPath(i-2:i-2) /= '.' ) &
      rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 enddo

 !remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 
 i = index(rectifyPath(i:l),'../')
 j = 0_pInt

 do while (i > j)
    j = scan(rectifyPath(:i-2),'/',back=.true.)
    rectifyPath(j+1:l) = rectifyPath(i+3:l)//repeat(' ',2+i-j)
    i = j+index(rectifyPath(j+1:l),'../')
 enddo
 if(len_trim(rectifyPath) == 0) rectifyPath = '/'
 return
 
 endfunction rectifyPath


!********************************************************************
! relative path from absolute a to absolute b
!
!********************************************************************
function makeRelativePath(a,b)

 use prec, only: pInt

 implicit none
 
 character (len=*) :: a,b
 character (len=1024) :: makeRelativePath
 integer(pInt) i,posLastCommonSlash,remainingSlashes

 posLastCommonSlash = 0
 remainingSlashes = 0
 do i = 1,min(1024,len_trim(a),len_trim(b))
   if (a(i:i) /= b(i:i)) exit
   if (a(i:i) == '/') posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1,len_trim(a)
   if (a(i:i) == '/') remainingSlashes = remainingSlashes + 1
 enddo
 makeRelativePath = repeat('../',remainingSlashes)//b(posLastCommonSlash+1:len_trim(b))
 return
endfunction makeRelativePath

END MODULE



include "IO.f90"               ! uses prec
include "numerics.f90"         ! uses prec, IO
include "math.f90"             ! uses prec, numerics
include "debug.f90"            ! uses prec, numerics
include "FEsolving.f90"        ! uses prec, IO
include "mesh.f90"             ! uses prec, math, IO, FEsolving
include "material.f90"         ! uses prec, math, IO, mesh
include "lattice.f90"          ! uses prec, math, IO, material
include "constitutive_phenopowerlaw.f90" ! uses prec, math, IO, lattice, material, debug
include "constitutive_j2.f90"            ! uses prec, math, IO, lattice, material, debug
include "constitutive_dislotwin.f90"    ! uses prec, math, IO, lattice, material, debug
include "constitutive_nonlocal.f90"      ! uses prec, math, IO, lattice, material, debug
include "constitutive.f90"     ! uses prec, IO, math, lattice, mesh, debug
include "crystallite.f90"      ! uses prec, math, IO, numerics 
include "homogenization_isostrain.f90"   ! uses prec, math, IO, 
include "homogenization_RGC.f90"         ! uses prec, math, IO, numerics, mesh: added <<<updated 31.07.2009>>>
include "homogenization.f90"   ! uses prec, math, IO, numerics
include "CPFEM.f90"            ! uses prec, math, IO, numerics, debug, FEsolving, mesh, lattice, constitutive, crystallite



!********************************************************************
program mpie_spectral
!********************************************************************

 use mpie_interface
 use prec, only: pInt, pReal
 use IO
! use math, only: math_I3,math_transpose3x3,math_Mandel66to3333
  use math
 use CPFEM, only: CPFEM_general

 implicit none

 real(pReal), dimension (:,:,:), allocatable :: bc_velocityGrad, &
                                                bc_stress          ! velocity gradient and stress BC
 real(pReal), dimension(:), allocatable :: bc_timeIncrement        ! length of increment
 integer(pInt), dimension(:), allocatable :: bc_steps              ! number of steps
 logical, dimension(:,:,:,:), allocatable :: bc_mask               ! mask
 
 real(pReal) temperature
 real(pReal), dimension(6) :: stress
 real(pReal), dimension(6,6) :: dsde
 
 character(len=1024) path,line
 logical, dimension(9) :: bc_maskvector
 logical gotResolution,gotDimension,gotHomogenization
 
 integer(pInt), parameter :: maxNchunksInput = 24                 ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunksInput*2) :: posInput
 integer(pInt), parameter :: maxNchunksMesh = 7                   ! 4 identifiers, 3 values  
 integer(pInt), dimension (1+2*maxNchunksMesh) :: posMesh
 
 real(pReal), dimension(9) :: valuevector
 
 integer(pInt) unit, N_l, N_s, N_t, N_n, N, i, j, k, l            ! numbers of identifiers, loop variables
 integer(pInt) a, b, c, e, homog
 real(pReal) x, y, z

!-------------------------
!begin RL
!-------------------------
 
 real(pReal), dimension(:), allocatable :: datafft
 real(pReal), dimension(:,:,:,:,:), allocatable :: workfft,workfftim,sg,disgrad,defgradold

 integer(pInt), dimension (3,3) :: iudot,iscau
  
 real(pReal), dimension(3,3) :: disgradmacro,disgradmacroactual
 real(pReal), dimension(3,3) :: ddisgradmacro,ddisgradmacroacum,ddisgrad,ddisgradim
 real(pReal), dimension(3,3) :: defgrad0,defgrad
 real(pReal), dimension(3,3) :: udot,scauchy,scauav,aux33

 integer(pInt), dimension(3) :: nn
 integer(pInt), dimension(2) :: nn2

 real(pReal), dimension(3) :: delt,xk,xk2
 real(pReal), dimension(6) :: aux6
 real(pReal), dimension(3,3,3,3) :: c0,s0,g1
 real(pReal), dimension(6,6) :: c066,s066

 integer(pInt) itmax, jload, ielem, ii, jj, k1, kxx, kyy, kzz, kx, ky, kz, idum, iter, imicro, m1, n1, p, q

 real(pReal) prodnn,wgt,error,tdot,erre,errs,evm,svm,det,xknorm

 logical errmatinv

!-------------------------
!end RL
!-------------------------
 
 if (IargC() < 2) call IO_error(102)

 path = getLoadcaseName()
 bc_maskvector = ''
 unit = 234_pInt
 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_n = 0_pInt
 
 print*,'Loadcase: ',trim(path)
 print*,'Workingdir: ',trim(getSolverWorkingDirectoryName())

 if (.not. IO_open_file(unit,path)) call IO_error(45,ext_msg = path)

 rewind(unit)
 do
   read(unit,'(a1024)',END = 101) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posInput = IO_stringPos(line,maxNchunksInput)
   do i = 1,maxNchunksInput,1
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
     call IO_error(46,ext_msg = path)                       ! error message for incomplete input file

 enddo

! allocate memory depending on lines in input file
101 N = N_l
 allocate (bc_velocityGrad(3,3,N));       bc_velocityGrad = 0.0_pReal
 allocate (bc_stress(3,3,N));             bc_stress = 0.0_pReal
 allocate (bc_mask(3,3,2,N));             bc_mask = .false.
 allocate (bc_timeIncrement(N));          bc_timeIncrement = 0.0_pReal
 allocate (bc_steps(N));                  bc_steps = 0_pInt

 rewind(unit)
 j = 0_pInt
 do
   read(unit,'(a1024)',END = 200) line
   if (IO_isBlank(line)) cycle                           ! skip empty lines
   j = j+1
   posInput = IO_stringPos(line,maxNchunksInput)
   do i = 1,maxNchunksInput,2
     select case (IO_lc(IO_stringValue(line,posInput,i)))
       case('l','velocitygrad')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,i+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,i+k)  ! assign values for the velocity gradient matrix
         enddo
         bc_mask(:,:,1,j) = reshape(bc_maskvector,(/3,3/))
         bc_velocityGrad(:,:,j) = reshape(valuevector,(/3,3/))
       case('s','stress')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,posInput,i+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,posInput,i+k)  ! assign values for the bc_stress matrix
         enddo
         bc_mask(:,:,2,j) = reshape(bc_maskvector,(/3,3/))
         bc_stress(:,:,j) = reshape(valuevector,(/3,3/))
       case('t','time','delta')                                            ! increment time
           bc_timeIncrement(j) = IO_floatValue(line,posInput,i+1)
       case('n','incs','increments','steps')                               ! bc_steps
           bc_steps(j) = IO_intValue(line,posInput,i+1)
     end select
   enddo
 enddo
200 close(unit)

 ! consistency checks
 do j = 1,N
   if (any(bc_mask(:,:,1,j) == bc_mask(:,:,2,j))) &
     call IO_error(47,j)                                                  ! bc_mask consistency
   if (any(math_transpose3x3(bc_stress(:,:,j)) + bc_stress(:,:,j) /= 2.0_pReal * bc_stress(:,:,j))) &
     call IO_error(48,j)                                                  ! bc_stress symmetry
 
   print '(a,/,3(3(f12.6,x)/))','L',bc_velocityGrad(:,:,j)
   print '(a,/,3(3(f12.6,x)/))','bc_stress',bc_stress(:,:,j)
   print '(a,/,3(3(l,x)/))','bc_mask',bc_mask(:,:,1,j)
   print *,'time',bc_timeIncrement(j)
   print *,'incs',bc_steps(j)
   print *, ''
 enddo
 
!read header of mesh file
 a = 1_pInt
 b = 1_pInt
 c = 1_pInt
 x = 1_pReal
 y = 1_pReal
 z = 1_pReal
 
 gotResolution =     .false.
 gotDimension =      .false.
 gotHomogenization = .false.
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
                 x = IO_floatValue(line,posMesh,i+1)
             case('y')
                 y = IO_floatValue(line,posMesh,i+1)
             case('z')
                 z = IO_floatValue(line,posMesh,i+1)
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
                 a = 2**IO_intValue(line,posMesh,i+1)
             case('b')
                 b = 2**IO_intValue(line,posMesh,i+1)
             case('c')
                 c = 2**IO_intValue(line,posMesh,i+1)
           end select
         enddo
   end select
   if (gotDimension .and. gotHomogenization .and. gotResolution) exit
 enddo
 100 close(unit)
 
 print '(a,/,i3,i3,i3)','resolution a b c',a,b,c
 print '(a,/,f6.2,f6.2,f6.2)','dimension x y z',x, y, z
 print *,'homogenization',homog
 print *, ''

 temperature = 300.0_pReal

!-------------------------
!begin RL
!-------------------------
 allocate (datafft(2*a*b*c))
 allocate (workfft(3,3,a,b,c))
 allocate (workfftim(3,3,a,b,c))
 allocate (sg(3,3,a,b,c))
 allocate (disgrad(3,3,a,b,c))
 allocate (defgradold(3,3,a,b,c))

 error = 0.00001
 itmax = 100

 delt(1) = 1.
 delt(2) = 1.
 delt(3) = 1.

 nn(1) = a
 nn(2) = b
 nn(3) = c
 
 nn2(1) = a
 nn2(2) = b

 prodnn = nn(1)*nn(2)*nn(3)
 wgt = 1./prodnn

! C_0 and S_0 CALCULATION
!!
!! PHILIP: FE_exec_elem?
!!
 stress = 0. !should be initialized somewhere
 dsde = 0. !should be initialized somewhere
 c0 = 0.  !stiffness of reference material
 c066 = 0. ! other way of notating c0

 do ielem = 1, prodnn !call each element with identity (math_i3) to initialize with high stress
   call CPFEM_general(3,math_i3,math_i3,temperature,0.0_pReal,ielem,1_pInt,stress,dsde)
   c066 = c066+dsde   
   c0 = c0+math_Mandel66to3333(dsde)
 enddo

 c066 = c066*wgt   
 c0 = c0*wgt

  
 call math_invert(6,c066,s066,idum,errmatinv)
 if(errmatinv) then
   write(*,*) 'ERROR IN C0 INVERSION'
   stop
 endif

 s0 = math_Mandel66to3333(s066)

! INITIALIZATION BEFORE STARTING WITH LOADINGS

 disgradmacro = 0.

 do k = 1, c
   do j = 1, b
     do i = 1, a
       defgradold(:,:,i,j,k) = math_i3(:,:)
     enddo
   enddo
 enddo

 do jload = 1,N !Loop over loadcases defined in the loadcase file
   udot(:,:) = bc_velocityGrad(:,:,jload)
   scauchy(:,:) = bc_stress(:,:,jload)
   iudot = 0
   iscau = 0
  
   do i = 1, 3 !convert information about rb's from bc_mask in corresponding arrays
     do j = 1, 3
       if(bc_mask(i,j,1,jload)) iudot(i,j) = 1
       if(bc_mask(i,j,2,jload)) iscau(i,j) = 1
     enddo
   enddo
  
   tdot = bc_timeIncrement(jload)/bc_steps(jload)
   ! evm = dvm*tdot ?
   do imicro = 1, bc_steps(jload)  ! loop oper steps defined in input file for current loadcase
     write(*,*) '***************************************************'
     write(*,*) 'STEP = ',imicro
   
! INITIALIZATION BEFORE NEW TIME STEP
     disgradmacro = disgradmacro+udot*tdot  !update macroscopic displacementgradient
     ddisgradmacro = 0.
     ielem=0  

     do k = 1, c  !loop over FPs
       do j = 1, b
         do i = 1, a
           ielem = ielem+1
           disgrad(:,:,i,j,k)=disgradmacro(:,:) 

           defgrad0(:,:)=defgradold(:,:,i,j,k)
           defgrad(:,:)=math_i3(:,:)+disgrad(:,:,i,j,k) ! calculate new defgrad 
           call CPFEM_general(3,defgrad0,defgrad,temperature,0.0_pReal,ielem,1_pInt,stress,dsde)
!    sg(:,:,i,j,k)=math_Mandel6to33(stress) ?
         enddo
       enddo
     enddo
	 
	 ielem = 0
	
     do k = 1, c !loop over FPs
       do j = 1, b
         do i = 1, a
           ielem = ielem+1
           disgrad(:,:,i,j,k) = disgradmacro(:,:)
		   
           defgrad0(:,:) = defgradold(:,:,i,j,k)
           defgrad(:,:) = math_i3(:,:)+disgrad(:,:,i,j,k)
           call CPFEM_general(1,defgrad0,defgrad,temperature,0.0_pReal,ielem,1_pInt,stress,dsde)
           sg(:,:,i,j,k) = math_Mandel6to33(stress)
         enddo
       enddo
     enddo
  
     ddisgradmacroacum = 0.

     iter = 0
     erre = 2.*error
     errs = 2.*error

     do while(iter <= itmax.and.(errs.gt.error.or.erre.gt.error))
       iter = iter+1
       write(*,*)'ITER = ',iter
       write(*,*) 'DIRECT FFT OF STRESS FIELD'
       do ii = 1,3
         do jj = 1,3
           k1 = 0
           do k = 1,c
!      write(*,'(1H+,a,i2,2(a,i4))') 
!     #   'STRESS - COMPONENT',ii,jj,'  -   Z = ',k,'   OUT OF ',npts
             do j = 1,b
               do i = 1,a
                 k1 = k1+1
                 datafft(k1) = sg(ii,jj,i,j,k)
                 k1 = k1+1
                 datafft(k1) = 0.
               enddo
             enddo
           enddo
           if(c.gt.1) then
             call fourn(datafft,nn,3,1)
           else
             call fourn(datafft,nn2,2,1)
           endif
           k1 = 0
           do k = 1, c
             do j = 1, b
               do i = 1, a
                 k1 = k1+1
                 workfft(ii,jj,i,j,k) = datafft(k1)
                 k1 = k1+1  
                 workfftim(ii,jj,i,j,k) = datafft(k1)
               enddo
             enddo
           enddo
         enddo
       enddo
       write(*,*) 'CALCULATING G^pq,ij : SG^ij ...'
       do kzz = 1,c
         do kyy = 1,b
           do kxx = 1,a
             if(kxx.le.a/2) kx = kxx-1
             if(kxx.gt.a/2) kx = kxx-a-1
 
             if(kyy.le.b/2) ky = kyy-1
             if(kyy.gt.b/2) ky = kyy-b-1
 
             if(kzz.le.c/2) kz = kzz-1
             if(kzz.gt.c/2) kz = kzz-c-1

             xk(1) = kx/(delt(1)*nn(1))
             xk(2) = ky/(delt(2)*nn(2))

             if(c.gt.1) then
               xk(3) = kz/(delt(3)*nn(3))
             else
               xk(3) = 0.
             endif
		   
             xknorm = sqrt(xk(1)**2+xk(2)**2+xk(3)**2)

             if (xknorm.ne.0.) then
               do i = 1,3
                 xk2(i) = xk(i)/(xknorm*xknorm*2.*pi)
                 xk(i) = xk(i)/xknorm
               enddo
             endif
	 
             do i = 1,3
               do k = 1,3
                 aux33(i,k) = 0.
                 do j = 1,3
                   do l = 1,3
                     aux33(i,k) = aux33(i,k)+c0(i,j,k,l)*xk(j)*xk(l)
                   enddo
                 enddo
               enddo
             enddo
             aux33 = math_inv3x3(aux33)
	  
!     call minv(aux33,3,det,minv1,minv2)
!     if(det.eq.0) then
!      write(*,*) kx,ky,kz,'  --> SINGULAR SYSTEM'
!      stop
!      pause
!     endif

             do p = 1,3
               do q = 1,3
                 do i = 1,3
                   do j = 1,3
                      g1(p,q,i,j) = -aux33(p,i)*xk(q)*xk(j)
                   enddo
                 enddo
               enddo
             enddo

             do i = 1,3
               do j = 1,3
                 ddisgrad(i,j) = 0.
                 ddisgradim(i,j) = 0.
!      if(kx.eq.0.and.ky.eq.0.and.kz.eq.0.) goto 4
                 if(.not.(kx.eq.0.and.ky.eq.0.and.kz.eq.0.)) then
                   do k = 1,3
                     do l = 1,3
                       ddisgrad(i,j) = ddisgrad(i,j)+g1(i,j,k,l)*workfft(k,l,kxx,kyy,kzz)
                       ddisgradim(i,j) = ddisgradim(i,j)+g1(i,j,k,l)*workfftim(k,l,kxx,kyy,kzz)
                     enddo
                   enddo
                 endif
               enddo
             enddo

             workfft(:,:,kxx,kyy,kzz) = ddisgrad(:,:)
             workfftim(:,:,kxx,kyy,kzz) = ddisgradim(:,:)
	
           enddo
         enddo
       enddo

       write(*,*) 'INVERSE FFT TO GET DISPLACEMENT GRADIENT FIELD'

       do m1 = 1,3
         do n1 = 1,3

           k1 = 0

           do k = 1,c
             do j = 1,b
               do i = 1,a
			   
                 k1 = k1+1
                 datafft(k1) = workfft(m1,n1,i,j,k)
				 
                 k1 = k1+1
                 datafft(k1) = workfftim(m1,n1,i,j,k)
               enddo
             enddo
           enddo

           if(c.gt.1) then
             call fourn(datafft,nn,3,-1)
           else
             call fourn(datafft,nn2,2,-1)
           endif
           datafft = datafft*wgt

            k1 = 0
            do k = 1,c
              do j = 1,b
                do i = 1,a

                  k1 = k1+1
                  disgrad(m1,n1,i,j,k) = disgrad(m1,n1,i,j,k)+ddisgradmacro(m1,n1)+datafft(k1)
                  k1 = k1+1
                enddo
              enddo
            enddo

          enddo
       enddo

       write(*,*) 'UPDATE STRESS FIELD'
!!!    call evpal
       ielem = 0  
       do k = 1,c
         do j = 1,b
           do i = 1,a

             ielem = ielem+1

             defgrad0(:,:) = defgradold(:,:,i,j,k)
             defgrad(:,:) = math_i3(:,:)+disgrad(:,:,i,j,k)

             call CPFEM_general(3,defgrad0,defgrad,temperature,0.0_pReal,ielem,1_pInt,stress,dsde)

             sg(:,:,i,j,k) = math_Mandel6to33(stress)
           enddo
         enddo
       enddo

       ielem = 0  

       do k = 1,c
         do j = 1,b
           do i = 1,a

             ielem = ielem+1

             defgrad0(:,:) = defgradold(:,:,i,j,k)
             defgrad(:,:) = math_i3(:,:)+disgrad(:,:,i,j,k)
			 
             call CPFEM_general(2,defgrad0,defgrad,temperature,0.0_pReal,ielem,1_pInt,stress,dsde)

              sg(:,:,i,j,k) = math_Mandel6to33(stress)
           enddo
         enddo
       enddo

!!!    call get_smacro
!  AVERAGE STRESS

       scauav = 0.
   
       do k = 1,c
         do j = 1,b
           do i = 1,a
             scauav(:,:) = scauav(:,:)+sg(:,:,i,j,k)*wgt
           enddo
         enddo
       enddo
!  MIXED BC

       do i = 1,3
         do j = 1,3
           ddisgradmacro(i,j) = 0.
           if(iudot(i,j).eq.0) then
             do k = 1,3
               do l = 1,3
                 ddisgradmacro(i,j) = ddisgradmacro(i,j)+s0(i,j,k,l)*iscau(k,l)*(scauchy(k,l)-scauav(k,l))
               enddo
             enddo
           endif
         enddo
       enddo
       ddisgradmacroacum = ddisgradmacroacum+ddisgradmacro
!      write(*,*) 'DDEFGRADMACRO(1,1),(2,2) = ',ddefgradmacro(1,1),ddefgradmacro(2,2)
!!  svm = ?
!    erre = erre/evm
!    errs = errs/svm
       write(*,*) 'STRESS FIELD ERROR  = ',errs
       write(*,*) 'STRAIN FIELD ERROR  = ',erre
!      write(21,101) iter,erre,errs,svm
!101   format(i3,4(1x,e10.4),10(1x,F7.4))

     enddo    ! WHILE  ENDDO

     disgradmacroactual = disgradmacro+ddisgradmacroacum

!!      write(*,*) 'defgradmacro(1,1),defgradmacro(2,2),defgradmacro(3,3)'
!!      write(*,*) defgradmacroactual(1,1),defgradmacroactual(2,2),defgradmacroactual(3,3)
!!      write(*,*) 'defgradmacro(1,1)/defgradmacro(3,3)'
!!      write(*,*) defgradmacroactual(1,1)/defgradmacroactual(3,3)
!!      write(*,*) 'scauav(1,1),scauav(2,2),scauav(3,3)'
!!      write(*,*) scauav(1,1),scauav(2,2),scauav(3,3)

     do k = 1,c
       do j = 1,b
         do i = 1,a
           defgradold(:,:,i,j,k) = math_i3(:,:)+disgrad(:,:,i,j,k)
         enddo
       enddo
     enddo
   enddo  ! IMICRO ENDDO
 enddo  ! JLOAD ENDDO Ende looping over loadcases

!-------------------------
!end RL
!-------------------------
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



!********************************************************************
! fourn subroutine (fourier transform)
!     FROM NUMERICAL RECIPES IN F77 (FIXED FORMAT), 
!     CONVERTED INTO FREE FORMAT (RL @ MPIE, JUNE 2010)
!********************************************************************
subroutine fourn(data,nn,ndim,isign)

      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
!      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
!     *k2,n,nprev,nrem,ntot
 INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot = 1
!      do 11 idim = 1,ndim
 do idim = 1,ndim  ! 11
        ntot = ntot*nn(idim)
!11    continue
 enddo ! 11
!
      nprev = 1
!      do 18 idim = 1,ndim
 do idim = 1,ndim  ! 18
        n = nn(idim)
        nrem = ntot/(n*nprev)
        ip1 = 2*nprev
        ip2 = ip1*n
        ip3 = ip2*nrem
        i2rev = 1
!        do 14 i2 = 1,ip2,ip1
   do i2 = 1,ip2,ip1  !  14
          if(i2.lt.i2rev)then

!            do 13 i1 = i2,i2+ip1-2,2
!              do 12 i3 = i1,ip3,ip2
     do i1 = i2,i2+ip1-2,2  ! 13
       do i3 = i1,ip3,ip2  ! 12
                i3rev = i2rev+i3-i2
                tempr = data(i3)
                tempi = data(i3+1)
                data(i3) = data(i3rev)
                data(i3+1) = data(i3rev+1)
                data(i3rev) = tempr
                data(i3rev+1) = tempi

!13          continue
     enddo  !  12
       enddo  !  13
          endif
          ibit = ip2/2

!1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
     do while ((ibit.ge.ip1).and.(i2rev.gt.ibit))  ! if 1
            i2rev = i2rev-ibit
            ibit = ibit/2
!          goto 1
!          endif
     enddo  ! do while (if 1)

          i2rev = i2rev+ibit
!14      continue
   enddo  ! 14
        ifp1 = ip1

!2       if(ifp1.lt.ip2)then
    do while (ifp1.lt.ip2)  !  if 2
          ifp2 = 2*ifp1
          theta = isign*6.28318530717959d0/(ifp2/ip1)
          wpr = -2.d0*sin(0.5d0*theta)**2
          wpi = sin(theta)
          wr = 1.d0
          wi = 0.d0
!          do 17 i3 = 1,ifp1,ip1
!            do 16 i1 = i3,i3+ip1-2,2
!              do 15 i2 = i1,ip3,ifp2
     do i3 = 1,ifp1,ip1  !  17
       do i1 = i3,i3+ip1-2,2  !  16
         do i2 = i1,ip3,ifp2  !  15
                k1 = i2
                k2 = k1+ifp1
                tempr = sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi = sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2) = data(k1)-tempr

                data(k2+1) = data(k1+1)-tempi
                data(k1) = data(k1)+tempr
                data(k1+1) = data(k1+1)+tempi
!15            continue
!16          continue
         enddo  !  15
       enddo  !  16
            wtemp = wr
            wr = wr*wpr-wi*wpi+wr
            wi = wi*wpr+wtemp*wpi+wi
!17        continue
     enddo  !  17
          ifp1 = ifp2
!        goto 2
!        endif
    enddo  !  do while (if 2)

        nprev = n*nprev
!18    continue
 enddo  ! 18
      return
      END subroutine
