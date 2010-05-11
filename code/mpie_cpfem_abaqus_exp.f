!* $Id$
!********************************************************************
! Material subroutine for Abaqus
!
! written by P. Eisenlohr,
!            F. Roters,
!            K. Janssens
!
! MPI fuer Eisenforschung, Duesseldorf
! PSI, Switzerland
!
! REMARK: change compile command to include the switch
!         "-free" in "abaqus_v6.env"
!
!********************************************************************

include "prec.f90"             ! uses nothing else


MODULE mpie_interface

character(len=64), parameter :: FEsolver = 'Abaqus'
character(len=4),  parameter :: InputFileExtension = '.inp'

CONTAINS

subroutine mpie_cpfem_init
  write(6,*)
  write(6,*) '<<<+-  mpie_cpfem_abaqus init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
 return
end subroutine

function getSolverWorkingDirectoryName
 use prec
 implicit none
 character(1024) getSolverWorkingDirectoryName
 integer(pInt) LENOUTDIR

 getSolverWorkingDirectoryName=''
 CALL VGETOUTDIR( getSolverWorkingDirectoryName, LENOUTDIR )
! write(6,*) 'getSolverWorkingDirectoryName', getSolverWorkingDirectoryName
end function

function getSolverJobName
 use prec
 implicit none

 character(1024) getSolverJobName, JOBNAME
 integer(pInt) LENJOBNAME

 getSolverJobName=''
 CALL VGETJOBNAME(getSolverJobName , LENJOBNAME )
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
 include "constitutive_phenopowerlaw.f90" ! uses prec, math, IO, latt ice, material, debug
 include "constitutive_j2.f90"            ! uses prec, math, IO, latt ice, material, debug
 include "constitutive_dislotwin.f90"    ! uses prec, math, IO, latt ice, material, debug
 include "constitutive_nonlocal.f90"      ! uses prec, math, IO, latt ice, material, debug
 include "constitutive.f90"     ! uses prec, IO, math, lattice, mesh, debug
 include "crystallite.f90"      ! uses prec, math, IO, numerics
 include "homogenization_isostrain.f90"   ! uses prec, math, IO, 
 include "homogenization_RGC.f90"         ! uses prec, math, IO, numerics, mesh: added <<<updated 31.07.2009>>>
 include "homogenization.f90"   ! uses prec, math, IO, numerics
 include "CPFEM.f90"            ! uses prec, math, IO, numerics, debug, FEsolving, mesh, lattice, constitutive, crystallite

subroutine vumat (jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                  stepTime, totalTime, dt, cmname, coordMp, charLength, &
                  props, density, strainInc, relSpinInc, &
                  tempOld, stretchOld, defgradOld, fieldOld, &
                  stressOld, stateOld, enerInternOld, enerInelasOld, &
                  tempNew, stretchNew, defgradNew, fieldNew, &
                  stressNew, stateNew, enerInternNew, enerInelasNew )

 include 'vaba_param.inc'

 dimension jblock(*), props(nprops), density(*), coordMp(*), &
           charLength(*), strainInc(*), &
           relSpinInc(*), tempOld(*), &
           stretchOld(*), &
           defgradOld(*), &
           fieldOld(*), stressOld(*), &
           stateOld(*), enerInternOld(*), &
           enerInelasOld(*), tempNew(*), &
           stretchNew(*), &
           defgradNew(*), &
           fieldNew(*), &
           stressNew(*), stateNew(*), &
           enerInternNew(*), enerInelasNew(*)

 character*80 cmname


 call vumatXtrArg ( jblock(1), &
                    ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                    stepTime, totalTime, dt, cmname, coordMp, charLength, &
                    props, density, strainInc, relSpinInc, &
                    tempOld, stretchOld, defgradOld, fieldOld, &
                    stressOld, stateOld, enerInternOld, enerInelasOld, &
                    tempNew, stretchNew, defgradNew, fieldNew, &
                    stressNew, stateNew, enerInternNew, enerInelasNew, &
                    jblock(5), jblock(2))

 return
 end subroutine


 subroutine vumatXtrArg (nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
                         stepTime, totalTime, dt, cmname, coordMp, charLength, &
                         props, density, strainInc, relSpinInc, &
                         tempOld, stretchOld, defgradOld, fieldOld, &
                         stressOld, stateOld, enerInternOld, enerInelasOld, &
                         tempNew, stretchNew, defgradNew, fieldNew, &
                         stressNew, stateNew, enerInternNew, enerInelasNew, &
                         nElement, nMatPoint)

 use prec, only:      pReal, &
                      pInt
 use FEsolving, only: cycleCounter, &
                      theTime, &
                      outdatedByNewInc, &
                      outdatedFFN1, &
                      terminallyIll, &
                      symmetricSolver
 use math, only:      invnrmMandel
 use debug, only:     debug_info, &
                      debug_reset, &
                      verboseDebugger
 use mesh, only:      mesh_FEasCP
 use CPFEM, only:     CPFEM_general,CPFEM_init_done
 use homogenization, only: materialpoint_sizeResults, materialpoint_results

 include 'vaba_param.inc'      ! Abaqus exp initializes a first step in single prec. for this a two-step compilation is used.
                               ! symbolic linking switches between .._sp.inc and .._dp.inc for both consecutive compilations...


 dimension props(nprops), density(nblock), &
           strainInc(nblock,ndir+nshr), &
           relSpinInc(nblock,nshr), defgradOld(nblock,ndir+nshr+nshr), &
           stressOld(nblock,ndir+nshr), &
           stateOld(nblock,nstatev), enerInternOld(nblock), &
           enerInelasOld(nblock), tempNew(nblock), tempOld(nblock), &
           stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), &
           stressNew(nblock,ndir+nshr) 

 dimension enerInelasNew(nblock),stateNew(nblock,nstatev),enerInternNew(nblock)
 dimension nElement(nblock),nMatPoint(nblock)

 character*80 cmname

! local variables
 real(pReal), dimension(3,3) :: defgrd0,defgrd1
 real(pReal), dimension(6) ::   stress
 real(pReal), dimension(6,6) :: ddsdde
 real(pReal) temp, timeInc
 integer(pInt) computationMode, n, i

 do n = 1,nblock                                                       ! loop over vector of IPs

   if ( .not. CPFEM_init_done ) then
     outdatedByNewInc = .false.

     if ( verboseDebugger ) then
       !$OMP CRITICAL (write2out)
         write(6,'(i6,x,i2,x,a)') nElement(n),nMatPoint(n),'first call special case..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   else if (theTime < totalTime) then                                  ! reached convergence
     outdatedByNewInc = .true.

     if ( verboseDebugger ) then
       !$OMP CRITICAL (write2out)
         write (6,'(i6,x,i2,x,a)') nElement(n),nMatPoint(n),'lastIncConverged + outdated'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   endif

   outdatedFFN1 = .false.
   terminallyIll = .false.
   cycleCounter = 1
   if ( outdatedByNewInc ) then
     outdatedByNewInc = .false.
     call debug_info()                                             ! first after new inc reports debugging
     call debug_reset()                                            ! resets debugging
     computationMode = 8                                           ! calc and age results with implicit collection
   else
     computationMode = 9                                           ! plain calc with implicit collection
   endif

   theTime  = totalTime                                            ! record current starting time

   if ( verboseDebugger ) then
     !$OMP CRITICAL (write2out)
       write(6,'(a16,x,i2,x,a,i5,x,i5,a)') 'computationMode',computationMode,'(',nElement(n),nMatPoint(n),')'; call flush(6)
     !$OMP END CRITICAL (write2out)
   endif
  
   defgrd0 = 0.0_pReal
   defgrd1 = 0.0_pReal
   temp    = tempOld(n)
   timeInc = dt

  !     ABAQUS explicit:     deformation gradient as vector 11, 22, 33, 12, 23, 31, 21, 32, 13
  !     ABAQUS explicit:     deformation gradient as vector 11, 22, 33, 12, 21
  
   forall (i=1:ndir)
     defgrd0(i,i) = defgradOld(n,i)
     defgrd1(i,i) = defgradNew(n,i)
   end forall
   if (nshr == 1) then
     defgrd0(1,2) = defgradOld(n,4)
     defgrd1(1,2) = defgradNew(n,4)
     defgrd0(2,1) = defgradOld(n,5)
     defgrd1(2,1) = defgradNew(n,5)
   else
     defgrd0(1,2) = defgradOld(n,4)
     defgrd1(1,2) = defgradNew(n,4)
     defgrd0(1,3) = defgradOld(n,9)
     defgrd1(1,3) = defgradNew(n,9)
     defgrd0(2,1) = defgradOld(n,7)
     defgrd1(2,1) = defgradNew(n,7)
     defgrd0(2,3) = defgradOld(n,5)
     defgrd1(2,3) = defgradNew(n,5)
     defgrd0(3,1) = defgradOld(n,6)
     defgrd1(3,1) = defgradNew(n,6)
     defgrd0(3,2) = defgradOld(n,8)
     defgrd1(3,2) = defgradNew(n,8)
   endif

   call CPFEM_general(computationMode,defgrd0,defgrd1,temp,timeInc,nElement(n),nMatPoint(n),stress,ddsdde)
  
  !     Mandel:     11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
  !     straight:   11, 22, 33, 12, 23, 13
  !     ABAQUS implicit:     11, 22, 33, 12, 13, 23
  !     ABAQUS explicit:     11, 22, 33, 12, 23, 13
  !     ABAQUS explicit:     11, 22, 33, 12

   stressNew(n,1:ndir+nshr) = stress(1:ndir+nshr)*invnrmMandel(1:ndir+nshr)
  
   stateNew(n,:) = materialpoint_results(1:min(nstatev,materialpoint_sizeResults),nMatPoint(n),mesh_FEasCP('elem', nElement(n)))
   tempNew(n) = temp
  
 enddo

 return
 end subroutine

!********************************************************************
!     This subroutine replaces the corresponding Marc subroutine
!********************************************************************
 subroutine quit(mpie_error)

 use prec, only:      pReal, &
                      pInt
 implicit none
 integer(pInt) mpie_error

 call xit
 end subroutine
