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

subroutine mpie_interface_init()
  write(6,*)
  write(6,*) '<<<+-  mpie_cpfem_abaqus init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
 return
end subroutine

function getSolverWorkingDirectoryName()
 use prec
 implicit none
 character(1024) getSolverWorkingDirectoryName
 integer(pInt) LENOUTDIR

 getSolverWorkingDirectoryName=''
 CALL GETOUTDIR( getSolverWorkingDirectoryName, LENOUTDIR )
! write(6,*) 'getSolverWorkingDirectoryName', getSolverWorkingDirectoryName
end function

function getSolverJobName()
 use prec
 implicit none

 character(1024) getSolverJobName, JOBNAME
 integer(pInt) LENJOBNAME

 getSolverJobName=''
 CALL GETJOBNAME(getSolverJobName , LENJOBNAME )
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
 
subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
       RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,&
       TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,&
       NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,&
       DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
 
 use prec, only:      pReal, &
                      pInt
 use FEsolving, only: cycleCounter, &
                      theInc, &
                      cutBack, &
                      calcMode, &
                      lastMode, &
                      theTime, &
                      theDelta, &
                      lastIncConverged, &
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


 implicit none
 
 CHARACTER*80 CMNAME
 integer(pInt) ndi, nshr, ntens, nstatv, nprops, noel, npt,&
               kslay, kspt, kstep, kinc
 real(pReal) STRESS(NTENS),STATEV(NSTATV),&
       DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
       STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
       PROPS(NPROPS),COORDS(3),DROT(3,3),&
       DFGRD0(3,3),DFGRD1(3,3)
 real(pReal) SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP,&
        DTEMP, PNEWDT, CELENT
 real(pReal), dimension (3,3) :: pstress             ! not used, but needed for call of cpfem_general
 real(pReal), dimension (3,3,3,3) :: dPdF            ! not used, but needed for call of cpfem_general

! local variables
 real(pReal), dimension(6) ::   stress_h
 real(pReal), dimension(6,6) :: ddsdde_h
 integer(pInt) computationMode, i, cp_en

 if (verboseDebugger .and. noel == 1 .and. npt == 1) then
   !$OMP CRITICAL (write2out)
     write(6,*) 'el',noel,'ip',npt
     write(6,*) 'got kinc as',kinc
     write(6,*) 'got dStran',dstran
     call flush(6)
   !$OMP END CRITICAL (write2out)
 endif

 if ( .not. CPFEM_init_done ) then

   computationMode = 2                                                    ! calc + init
   if ( verboseDebugger ) then
     !$OMP CRITICAL (write2out)
       write(6,'(i6,x,i2,x,a)') noel,npt,'first call special case..!'; call flush(6)
     !$OMP END CRITICAL (write2out)
   endif

 else
   cp_en = mesh_FEasCP('elem',noel)
   if (theTime < time(2) .or. theInc /= kinc) then                        ! reached convergence
     lastIncConverged = .true.
     outdatedByNewInc = .true.
     terminallyIll = .false.
     cycleCounter = 0

     if ( verboseDebugger ) then
       !$OMP CRITICAL (write2out)
         write (6,'(i6,x,i2,x,a)') noel,npt,'lastIncConverged + outdated'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   else if ( dtime < theDelta ) then                                      ! or check for cutBack
     calcMode = .true.                                                    ! pretend last step was calculation
     cutBack = .true.
     terminallyIll = .false.
     cycleCounter = 0

     if ( verboseDebugger ) then
       !$OMP CRITICAL (write2out)
         write(6,'(i6,x,i2,x,a)') noel,npt,'cutback detected..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   endif

   calcMode(npt,cp_en) = .not. calcMode(npt,cp_en)                        ! ping pong (calc <--> collect)

   if ( calcMode(npt,cp_en) ) then                                        ! now calc
     if ( lastMode .ne. calcMode(npt,cp_en) ) then                        ! first after ping pong
       call debug_reset()                                                 ! resets debugging
       outdatedFFN1 = .false.
       cycleCounter = cycleCounter + 1
     endif
     if ( outdatedByNewInc ) then
       outdatedByNewInc = .false.
       computationMode = 1                                                ! calc and age results
     else
       computationMode = 2                                                ! plain calc
     endif
   else                                                                   ! now collect
     if ( lastMode .ne. calcMode(npt,cp_en) ) call debug_info()           ! first after ping pong reports debugging
     if ( lastIncConverged ) then
       lastIncConverged = .false.
       computationMode = 4                                                ! collect and backup Jacobian after convergence
     elseif ( cutBack ) then
       cutBack = .false.
       computationMode = 5                                                ! collect and restore Jacobian after cutback
     else
       computationMode = 3                                                ! plain collect
     endif
   endif

 endif
 
 theTime  = time(2)                                                       ! record current starting time
 theDelta = dtime                                                         ! record current time increment
 theInc   = kinc                                                          ! record current increment number
 if (CPFEM_init_done) lastMode = calcMode(npt,cp_en)                      ! record calculationMode

 if ( verboseDebugger ) then
   !$OMP CRITICAL (write2out)
     write(6,'(a16,x,i2,x,a,i5,a,i5,x,i5,a)') 'computationMode',computationMode,'(',cp_en,':',noel,npt,')'; call flush(6)
   !$OMP END CRITICAL (write2out)
 endif
   
 call CPFEM_general(computationMode,dfgrd0,dfgrd1,temp,dtime,noel,npt,stress_h,ddsdde_h, pstress, dPdF)

!     Mandel:              11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     straight:            11, 22, 33, 12, 23, 13
!     ABAQUS explicit:     11, 22, 33, 12, 23, 13
!     ABAQUS implicit:     11, 22, 33, 12, 13, 23
!     ABAQUS implicit:     11, 22, 33, 12

 forall(i=1:ntens) ddsdde(1:ntens,i) = invnrmMandel(i)*ddsdde_h(1:ntens,i)*invnrmMandel(1:ntens)
 stress(1:ntens) = stress_h(1:ntens)*invnrmMandel(1:ntens)
 if(symmetricSolver) ddsdde(1:ntens,1:ntens) = 0.5_pReal*(ddsdde(1:ntens,1:ntens) + transpose(ddsdde(1:ntens,1:ntens)))
 if(ntens == 6) then
   stress_h = stress
   stress(5) = stress_h(6)
   stress(6) = stress_h(5)
   ddsdde_h = ddsdde
   ddsdde(:,5) = ddsdde_h(:,6)
   ddsdde(:,6) = ddsdde_h(:,5)
   ddsdde_h = ddsdde
   ddsdde(5,:) = ddsdde_h(6,:)
   ddsdde(6,:) = ddsdde_h(5,:)
 end if

 statev = materialpoint_results(1:min(nstatv,materialpoint_sizeResults),npt,mesh_FEasCP('elem', noel))

 if ( terminallyIll ) pnewdt = 0.5_pReal                                  ! force cutback directly ?

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
