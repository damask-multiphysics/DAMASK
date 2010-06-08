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
 character(len=1024) path
 
 call mpie_interface_init()
 if (IargC() < 2) then
   print *,'buh'
 else
   path = getSolverWorkingDirectoryName() 
   print *, path 
 endif
end program mpie_spectral

subroutine quit(id)
 use prec

 implicit none

 integer(pInt) id

 stop
end subroutine
