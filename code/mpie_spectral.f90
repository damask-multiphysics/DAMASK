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
end subroutine

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
 
end function

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
 
 if (posExt <= posSep) posExt = len_trim(outName)       ! no extension present

 getSolverJobName = outName(posSep+1:posExt-1)          ! path to mesh file (excl. extension)

 if (scan(getSolverJobName,pathSep) /= 1) then          ! relative path given as command line argument
   call getcwd(cwd)
   getSolverJobName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                       rectifyPath(trim(cwd)//'/'//getSolverJobName))
 endif
  
 return
end function


!********************************************************************
! realive path of loadcase from command line arguments
!
!********************************************************************
function getLoadcaseName()

 use prec, only: pInt

 implicit none

 character(len=1024) getLoadcaseName, outName, cwd
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 integer(pInt) posExt,posSep

 call getarg(2,getLoadcaseName)

 if (scan(getLoadcaseName,pathSep) /= 1) then          ! relative path given as command line argument
   call getcwd(cwd)
   getLoadcaseName = rectifyPath(trim(cwd)//'/'//getLoadcaseName)
 endif
 
 getLoadcaseName = makeRelativePath(getSolverWorkingDirectoryName(),&
                                    getLoadcaseName)
  
 return
end function


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
 do i=l,2,-1
    if ( rectifyPath(i-1:i)=='./' .and. rectifyPath(i-2:i-2) /= '.' ) &
      rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
 end do

 !remove ../ and corresponding directory from rectifyPath
 l = len_trim(rectifyPath)
 
 i = index(rectifyPath(i:l),'../')
 j = 0_pInt

 do while (i > j)
    j = scan(rectifyPath(:i-2),'/',back=.true.)
    rectifyPath(j+1:l) = rectifyPath(i+3:l)//repeat(' ',2+i-j)
    i = j+index(rectifyPath(j+1:l),'../')
 end do
 return
 
 end function rectifyPath


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
end function makeRelativePath

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
 use math, only: math_I3,math_transpose3x3
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
 integer(pInt), parameter :: maxNchunks = 24                 ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunks*2) :: pos
 real(pReal), dimension(9) :: valuevector
 integer(pInt) unit, N_l, N_s, N_t, N_n, N, i,j,k,l          ! numbers of identifiers, loop variables

 if (IargC() < 2) call IO_error(102)

 path = getLoadcaseName()
 bc_maskvector = ''
 unit = 234_pInt
 N_l = 0_pInt
 N_s = 0_pInt
 N_t = 0_pInt
 N_n = 0_pInt

 if (.not. IO_open_file(unit,path)) call IO_error(45,ext_msg=path)

 rewind(unit)
 do
   read(unit,'(a1024)',END=101) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   pos = IO_stringPos(line,maxNchunks)
   do i = 1,maxNchunks,1
       select case (IO_lc(IO_stringValue(line,pos,i)))
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
     call IO_error(46,ext_msg=path)                       ! error message for incomplete input file

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
   read(unit,'(a1024)',END=200) line
   if (IO_isBlank(line)) cycle                           ! skip empty lines
   j=j+1
   pos = IO_stringPos(line,maxNchunks)
   do i = 1,maxNchunks,2
     select case (IO_lc(IO_stringValue(line,pos,i)))
       case('l','velocitygrad')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,pos,i+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,pos,i+k)  ! assign values for the velocity gradient matrix
         enddo
         bc_mask(:,:,1,j) = reshape(bc_maskvector,(/3,3/))
         bc_velocityGrad(:,:,j) = reshape(valuevector,(/3,3/))
       case('s','stress')
         valuevector = 0.0_pReal
         forall (k = 1:9) bc_maskvector(k) = IO_stringValue(line,pos,i+k) /= '#'
         do k = 1,9
           if (bc_maskvector(k)) valuevector(k) = IO_floatValue(line,pos,i+k)  ! assign values for the bc_stress matrix
         enddo
         bc_mask(:,:,2,j) = reshape(bc_maskvector,(/3,3/))
         bc_stress(:,:,j) = reshape(valuevector,(/3,3/))
       case('t','time','delta')                                            ! increment time
           bc_timeIncrement(j) = IO_floatValue(line,pos,i+1)
       case('n','incs','increments','steps')                               ! bc_steps
           bc_steps(j) = IO_intValue(line,pos,i+1)
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

 temperature = 300.0_pReal
 call CPFEM_general(2,math_i3,math_i3,temperature,0.0_pReal,1_pInt,1_pInt,stress,dsde)

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
