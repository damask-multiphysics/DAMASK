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

MODULE cpfem_interface

character(len=64), parameter :: FEsolver = 'Abaqus'

CONTAINS

subroutine mpie_cpfem_init ()

!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  mpie_cpfem_abaqus init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
  call flush(6)
!$OMP END CRITICAL (write2out)
  return
end subroutine

END MODULE

 include "prec.f90"             ! uses nothing else
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
 include "constitutive_dislobased.f90"    ! uses prec, math, IO, latt ice, material, debug
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
                      debug_reset
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

! local variables
 real(pReal), dimension(6) ::   stress_h
 real(pReal), dimension(6,6) :: ddsdde_h
 integer(pInt) computationMode, i, cp_en

 if (noel == 1 .and. npt == 1) then
!$OMP CRITICAL (write2out)
   write(6,*) 'el',noel,'ip',npt
   write(6,*) 'got kinc as',kinc
   write(6,*) 'got dStran',dstran
   call flush(6)
!$OMP END CRITICAL (write2out)
 endif

 if ( .not. CPFEM_init_done ) then

   computationMode = 2                                                    ! calc + init
!$OMP CRITICAL (write2out)
   write(6,'(i6,x,i2,x,a)') noel,npt,'first call special case..!'; call flush(6)
!$OMP END CRITICAL (write2out)

 else
   cp_en = mesh_FEasCP('elem',noel)
   if (theTime < time(2) .or. theInc /= kinc) then                        ! reached convergence
     lastIncConverged = .true.
     outdatedByNewInc = .true.
     terminallyIll = .false.
     cycleCounter = 0
!$OMP CRITICAL (write2out)
     write (6,'(i6,x,i2,x,a)') noel,npt,'lastIncConverged + outdated'; call flush(6)
!$OMP END CRITICAL (write2out)
   endif

   if ( dtime < theDelta ) then                                           ! cutBack
     calcMode = .true.                                                    ! pretend last step was calculation
     cutBack = .true.
     terminallyIll = .false.
     cycleCounter = 0
!$OMP CRITICAL (write2out)
     write(6,'(i6,x,i2,x,a)') noel,npt,'cutback detected..!'; call flush(6)
!$OMP END CRITICAL (write2out)
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

!$OMP CRITICAL (write2out)
 write(6,'(a16,x,i2,x,a,i5,a,i5,x,i5,a)') 'computationMode',computationMode,'(',cp_en,':',noel,npt,')'; call flush(6)
!$OMP END CRITICAL (write2out)
   
 call CPFEM_general(computationMode,dfgrd0,dfgrd1,temp,dtime,noel,npt,stress,ddsdde,ntens)

!     Mandel:     11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     straight:   11, 22, 33, 12, 23, 13
 forall(i=1:ntens) ddsdde(1:ntens,i) = invnrmMandel(i)*ddsdde(1:ntens,i)*invnrmMandel(1:ntens)
 stress(1:ntens) = stress(1:ntens)*invnrmMandel(1:ntens)
 if(symmetricSolver) ddsdde(1:ntens,1:ntens) = 0.5_pReal*(ddsdde(1:ntens,1:ntens) + transpose(ddsdde(1:ntens,1:ntens)))
!     ABAQUS:     11, 22, 33, 12, 13, 23
 if(ntens == 6) then
   stress_h=stress
   stress(5)=stress_h(6)
   stress(6)=stress_h(5)
   ddsdde_h=ddsdde
   ddsdde(:,5)=ddsdde_h(:,6)
   ddsdde(:,6)=ddsdde_h(:,5)
   ddsdde_h=ddsdde
   ddsdde(5,:)=ddsdde_h(6,:)
   ddsdde(6,:)=ddsdde_h(5,:)
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

!############################################################################
!

!      include "KJ_Disp.f"
      subroutine disp(u,kstep,kinc,time,node,noel,jdof,coords)

!     hardwired aba_param.inc
      implicit real*8(a-h,o-z)
      parameter (nprecd=2)
!
      dimension u(3),time(2), coords(3)

      real ktime,ktcl,ktmax,ktmin,ktdeltaup,ktdeltadown
      real klmin,klmax,kldelta,klbegincycle
      real kdeltau, kru, kret
! When using stdb_abqerr for debugging
!      dimension intv(2),realv(4)
!      character*8 charv(1)

! ratchet displacement per cycle
      kru = 0.1
! ratcheting ends at time kret
      kret = 400.
! displacement amplitude
      kdeltau = 0.5
! time cycle length:
      ktcl = 4.
      ktmax = ktcl/4
      ktmin = 3.*ktmax
      ktdeltadown = ktmin - ktmax
      ktdeltaup = ktcl - ktdeltadown
! load minimum & maximum:
      klmin = -kdeltau
      klmax = kdeltau
      kldelta = klmax - klmin
      klbegincycle = klmin + kldelta * (ktcl-ktmin) / ktdeltaup
! load as a function of (total time); trianglar loading cycle
      ktime = time(2)
      kru = kru * MIN(kret, ktime) / ktcl
      if ( ktime .lt. ktmax ) then
! special case for path to first maximum
          u(1) = kru + (ktime/ktmax) * klmax
      else
          do while ( ktime .ge. ktcl )
             ktime = ktime - ktcl
          end do
          if ( ktime .le. ktmax ) then
             u(1) = kru + klbegincycle + ktime * (klmax-klbegincycle) / ktmax
          else if ( ktime .lt. ktmin ) then
             u(1) = kru + klmax - (ktime-ktmax) * kldelta / (ktmin-ktmax)
          else
             u(1) = kru + klmin + (ktime-ktmin) * (klbegincycle-klmin)/(ktcl-ktmin)
          end if
      endif
      return
      end
!
