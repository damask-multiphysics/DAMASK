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

character(len=64), parameter :: FEsolver = 'Spectral'
character(len=5),  parameter :: InputFileExtension = '.mesh'

CONTAINS

subroutine mpie_interface_init
 write(6,*)
 write(6,*) '<<<+-  mpie_spectral init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 return
end subroutine

function getSolverWorkingDirectoryName()
 implicit none
 character(1024) cwd,outname,getSolverWorkingDirectoryName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \

 print *, 'start of func'

 call getarg(2,outname)                                ! path to loadFile
 
 if (scan(outname,pathSep) == 1) then                  ! absolute path given as command line argument
   getSolverWorkingDirectoryName = outname(1:scan(outname,pathSep,back=.true.))
   print *, 'abs',scan(outname,pathSep,back=.true.)
 else
   call getcwd(cwd)
   getSolverWorkingDirectoryName = trim(cwd)//'/'//outname(1:scan(outname,pathSep,back=.true.))
   print *, 'rel',scan(outname,pathSep,back=.true.),getSolverWorkingDirectoryName
 endif
 
 return
! getSolverWorkingDirectoryName = rectifyPath(getSolverWorkingdirectory) 
end function

function getSolverJobName()
 use prec
 implicit none

 character(1024) getSolverJobName, outName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 integer(pInt) extPos

 getSolverJobName=''
 call getarg(1,outName)
 extPos = len_trim(outName)-4
 getSolverJobName = outName(scan(outName,pathSep,back=.true.)+1:extPos)
! write(6,*) 'getSolverJobName', getSolverJobName
end function

!function removes ../ and ./ in Path
function rectifyPath(Path)
implicit none
character(len=1024) Path, rectifyPath
integer i,j,k,l

!remove ./ from path
l = min(1024,len_trim(Path))
rectifyPath = path
do i=l,2,-1
    if ( rectifyPath(i-1:i)=='./' .and. rectifyPath(i-2:i-2) /= '.' ) &
      rectifyPath(i-1:l) = rectifyPath(i+1:l)//'  '
end do

!remove ../ from rectifyPath and change directory
i=0
do while (i<=len_trim(rectifyPath))
    k=0
    j=0
    if(rectifyPath(len_trim(rectifyPath)-2-i:len_trim(rectifyPath)-i)=='../') then           !search for ../ (directory down)
         j=i
         k=i
         if(len_trim(rectifyPath)-5-i<=0) print*, 'enter error message here'   !invalid rectifyPath (below root)
         do while(rectifyPath(len_trim(rectifyPath)-j-5:len_trim(rectifyPath)-j-3)=='../')   !search for more ../ in front of first hit
              j=j+3
         end do
         j=((j-i)/3+1)*2                                                !calculate number of ../
         do while(j/=0)                                                 !find position to break rectifyPath
              i=i+1
              if(rectifyPath(len_trim(rectifyPath)-i:len_trim(rectifyPath)-i)=='/') j=j-1
         end do
         if(i>len_trim(rectifyPath)) print*, 'enter error message here'        !invalid path (below root)
         rectifyPath(len_trim(rectifyPath)-i:len_trim(rectifyPath))=rectifyPath(len_trim(rectifyPath)-k:len_trim(rectifyPath)+i-k)
         i=k-1
     end if
     i=i+1
end do
end function rectifyPath




! make out of two absolute Paths (a,b) relative Path from a to b
function makeRelativePath(a,b)
 implicit none
 character (len=1024) :: makeRelativePath,a,b
 integer i,posLastCommonSlash,remainingSlashes

 posLastCommonSlash = 0
 remainingSlashes = 0
 do i = 1,min(len_trim(a),len_trim(b))
   if (a(i:i) /= b(i:i)) exit
   if (a(i:i) == '/') posLastCommonSlash = i
 enddo
 do i = posLastCommonSlash+1,len_trim(a)
   if (a(i:i) == '/') remainingSlashes = remainingSlashes + 1
 enddo
 makeRelativePath=repeat('../',remainingSlashes)//b(posLastCommonSlash+1:len_trim(b))
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



program mpie_spectral

 use prec
 use mpie_interface
 
 implicit none
 character(len=1024) path, line
 integer(pInt), parameter :: maxNchunks = 24                ! 4 identifiers, 18 values for the matrices and 2 scalars
 integer(pInt), dimension (1+maxNchunks*2) :: pos
 real(pReal), dimension (:,:), allocatable :: l, s          ! velocity gradient and stress BC
 real(pReal), dimension(:), allocatable :: t, n             ! length of time and step number
 character, dimension(:,:), allocatable :: mask             ! BC mask
 integer(pInt) unit, N_l, N_s, N_t, N_n, i, j, k            ! numbers of identifiers, loop variables
 
 call mpie_interface_init()
 if (IargC() < 2) then
   print *,'buh'
 else
   path = rectifyPath(getSolverWorkingDirectoryName()) 
 endif
 
! initialize variables
unit=2
N_l=0
N_s=0
N_t=0
N_n=0
if (IO_open_file(unit,path)) rewind(unit)                  ! open file
do
    read(unit,'(a1024)',END=101) line
    if (IO_isBlank(line)) cycle                            ! skip empty lines
    pos = IO_stringPos(line,maxNchunks)
    do i = 1,maxNchunks,1
         select case (IO_lc(IO_stringValue(line,pos,i)))
              case('l')
                   N_l=N_l+1
              case('s')
                   N_s=N_s+1
              case('t')
                   N_t=N_t+1
              case('n')
                   N_n=N_n+1
          end select
    enddo                                                  ! count all identifiers to allocate memory and do sanity check
   if ((N_l /= N_s).or.(N_s /= N_t).or.(N_t /= N_n)) then  ! sanity check
   print*, 'insert error message code here'                ! error message for incomplete input file
   end if
enddo

! allocate memory depending on lines in input file
101 allocate (l(9,N_l))
allocate (s(9,N_s))
allocate (mask(9,N_s))
allocate (t(N_t))
allocate (n(N_n))

! initialize variables
do i=1,9
    do j=1,N_l
         mask(i,j)='x'
    end do
end do
do i=1,9
    do j=1,N_l
         l(i,j)=0
    end do
end do
do i=1,9
    do j=1,N_s
         s(i,j)=0
    end do
end do
i=0
j=0

rewind(unit)
do
    read(unit,'(a1024)',END=100) line
    if (IO_isBlank(line)) cycle                                            ! build BC mask from input file
    j=j+1
    pos = IO_stringPos(line,maxNchunks)
         do i = 1,maxNchunks,2
              select case (IO_lc(IO_stringValue(line,pos,i)))
                   case('l')
                   do k=1,9
                        if((IO_lc(IO_stringValue(line,pos,i+k)))/='-') mask(k,j)='l'
                   end do
                   if(((mask(2,j)/='x')).and.((mask(4,j)=='x'))&                ! if one non-diagonal element is defined, the
                       &.or.((mask(4,j)/='x')).and.((mask(2,j)=='x'))&          ! correspondig one should not be empty
                       &.or.((mask(3,j)/='x')).and.((mask(7,j)=='x'))&
                       &.or.((mask(7,j)/='x')).and.((mask(3,j)=='x'))&
                       &.or.((mask(6,j)/='x')).and.((mask(8,j)=='x'))&
                       &.or.((mask(8,j)/='x')).and.((mask(6,j)=='x'))) print*, 'enter error message here'
                   case('s')
                   do k=1,9
                        if((IO_lc(IO_stringValue(line,pos,i+k)))/='-') then
                             if(mask(k,j)=='l') then
                                  print*, 'enter error message here'           ! stress and velocity gradient bc at the same place
                             else
                                  mask(k,j)='s'
                             end if
                        end if
                   end do
              end select
         enddo
    enddo
100 rewind(unit)

do i=1,9
    do j=1,N_l
         if(mask(i,j)=='x') print*,'enter error message here' !check if sufficient Nr. of BCs are found
    end do
end do

j=0
i=0
do
   read(unit,'(a1024)',END=200) line
   if (IO_isBlank(line)) cycle                           ! skip empty lines
   j=j+1
   pos = IO_stringPos(line,maxNchunks)
     do i = 1,maxNchunks,2
       select case (IO_lc(IO_stringValue(line,pos,i)))
         case('l')
              do k=1,9
                   if(mask(k,j)=='l')  L(k,j) = IO_floatValue(line,pos,i+k)  ! assign values for the velocity gradient matrix
              end do
         case('s')
              do k=1,9
                   if(mask(k,j)=='s') then
                        s(k,j) = IO_floatValue(line,pos,i+k)  ! assign values for the stress BC
                        select case(k)
                             case(4)
                                  if(s(4,j)/=s(2,j)) print*, 'enter error code here' !non-symmetric stress BC
                        end select
                   else
                   end if
              end do
         case('t')                                       ! assign the scalars
             t(j) = IO_floatValue(line,pos,i+1)
         case('n')
             n(j) = IO_floatValue(line,pos,i+1)
       end select
     enddo
   enddo
200 end program mpie_spectral

subroutine quit(id)
 use prec

 implicit none

 integer(pInt) id

 stop
end subroutine
