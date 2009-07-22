!********************************************************************
! Material subroutine for MSC.Marc
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts
!            D.D. Tjahjanto
!
! MPI fuer Eisenforschung, Duesseldorf
!
!********************************************************************
!     Usage:
!             - choose material as hypela2
!             - set statevariable 2 to index of homogenization
!             - set statevariable 3 to index of microstructure
!             - make sure the file "material.config" exists in the working
!               directory
!             - make sure the file "numerics.config" exists in the working 
!               directory
!             - use nonsymmetric option for solver (e.g. direct 
!               profile or multifrontal sparse, the latter seems
!               to be faster!)
!             - in case of ddm (domain decomposition)a SYMMETRIC
!               solver has to be used, i.e uncheck "non-symmetric"
!********************************************************************
!     Marc subroutines used:
!             - hypela2
!             - plotv
!             - quit
!********************************************************************
!     Marc common blocks included:
!             - concom: lovl, ncycle, inc, incsub
!             - creeps: timinc
!********************************************************************
!
 include "prec.f90"             ! uses nothing else
 include "math.f90"             ! uses prec
 include "IO.f90"               ! uses prec, math
 include "numerics.f90"         ! uses prec, IO
 include "debug.f90"            ! uses prec, numerics
 include "FEsolving.f90"        ! uses prec, IO
 include "mesh.f90"             ! uses prec, math, IO, FEsolving
 include "material.f90"         ! uses prec, math, IO, mesh
 include "lattice.f90"          ! uses prec, math, IO, material
 include "constitutive_phenopowerlaw.f90"     ! uses prec, math, IO, lattice, material, debug
 include "constitutive_j2.f90"                   ! uses prec, math, IO, lattice, material, debug
 include "constitutive_dislobased.f90"           ! uses prec, math, IO, lattice, material, debug
 include "constitutive.f90"     ! uses prec, IO, math, lattice, mesh, debug
 include "crystallite.f90"      ! uses prec, math, IO, numerics 
 include "homogenization_isostrain.f90"          ! uses prec, math, IO, 
 include "homogenization.f90"   ! uses prec, math, IO, numerics
 include "CPFEM.f90"            ! uses prec, math, IO, numerics, debug, FEsolving, mesh, lattice, constitutive, crystallite


!********************************************************************
! This is the Marc material routine
!********************************************************************
!
! *************   user subroutine for defining material behavior  **************
!
!
! CAUTION : Due to calculation of the Deformation gradients, Stretch Tensors and
!         Rotation tensors at previous and current states, the analysis can be
!         computationally expensive. Please use the user subroutine ->  hypela
!         if these kinematic quantities are not needed in the constitutive model
!
!
! IMPORTANT NOTES :
!
! (1) F,R,U are only available for continuum and membrane elements (not for
!     shells and beams).
!
! (2) For total Lagrangian formulation use the -> 'Elasticity,1' card(=
!     total Lagrange with large disp) in the parameter section of input deck.
!     For updated Lagrangian formulation use the -> 'Plasticity,3' card(=
!     update+finite+large disp+constant d) in the parameter section of
!     input deck.
!
!     The following operation obtains U (stretch tensor) at t=n+1 :
!
!     call scla(un1,0.d0,itel,itel,1)
!     do 3 k=1,3
!      do 2 i=1,3
!       do 1 j=1,3
!        un1(i,j)=un1(i,j)+dsqrt(strechn1(k))*eigvn1(i,k)*eigvn1(j,k)
!1      continue
!2     continue
!3    continue
!
!********************************************************************
subroutine hypela2(&
     d,&          ! stress strain law to be formed
     g,&          ! change in stress due to temperature effects
     e,&          ! total elastic strain
     de,&         ! increment of strain
     s,&          ! stress - should be updated by user
     t,&          ! state variables (comes in at t=n, must be updated to have state variables at t=n+1)
     dt,&         ! increment of state variables
     ngens,&      ! size of stress - strain law
     n,&          ! element number
     nn,&         ! integration point number
     kcus,&       ! (1) layer number, (2) internal layer number
     matus,&      ! (1) user material identification number, (2) internal material identification number
     ndi,&        ! number of direct components
     nshear,&     ! number of shear components
     disp,&       ! incremental displacements
     dispt,&      ! displacements at t=n (at assembly, lovl=4) and displacements at t=n+1 (at stress recovery, lovl=6)
     coord,&      ! coordinates
     ffn,&        ! deformation gradient
     frotn,&      ! rotation tensor
     strechn,&    ! square of principal stretch ratios, lambda(i)
     eigvn,&      ! i principal direction components for j eigenvalues
     ffn1,&       ! deformation gradient
     frotn1,&     ! rotation tensor
     strechn1,&   ! square of principal stretch ratios, lambda(i)
     eigvn1,&     ! i principal direction components for j eigenvalues
     ncrd,&       ! number of coordinates
     itel,&       ! dimension of F and R, either 2 or 3
     ndeg,&       ! number of degrees of freedom  ==> is this at correct list position ?!?
     ndm,&        !
     nnode,&      ! number of nodes per element
     jtype,&      ! element type
     lclass,&     ! element class
     ifr,&        ! set to 1 if R has been calculated
     ifu &        ! set to 1 if stretch has been calculated
   )

 use prec, only:      pReal, &
                      pInt
 use FEsolving, only: cycleCounter, &
                      theInc, &
                      theCycle, &
                      theLovl, &
                      theTime, &
                      lastIncConverged, &
                      outdatedByNewInc, &
                      outdatedFFN1, &
                      symmetricSolver
 use CPFEM, only:     CPFEM_general
 use math, only:      invnrmMandel
 use debug, only:     debug_info, &
                      debug_reset
 implicit none
 
!     ** Start of generated type statements **
 real(pReal) coord, d, de, disp, dispt, dt, e, eigvn, eigvn1, ffn, ffn1
 real(pReal) frotn, frotn1, g
 integer(pInt) ifr, ifu, itel, jtype, kcus, lclass, matus, n, ncrd, ndeg
 integer(pInt) ndi, ndm, ngens, nn, nnode, nshear
 real(pReal) s, strechn, strechn1, t
 !     ** End of generated type statements **

 dimension e(*),de(*),t(*),dt(*),g(*),d(ngens,*),s(*), n(2),coord(ncrd,*),disp(ndeg,*),matus(2),dispt(ndeg,*),ffn(itel,*),&
           frotn(itel,*),strechn(itel),eigvn(itel,*),ffn1(itel,*),frotn1(itel,*),strechn1(itel),eigvn1(itel,*),kcus(2), lclass(2)

! Marc common blocks are in fixed format so they have to be reformated to free format (f90)
! Beware of changes in newer Marc versions

 include "concom%%MARCVERSION%%"     ! concom is needed for inc, subinc, ncycle, lovl
 include "creeps%%MARCVERSION%%"     ! creeps is needed for timinc (time increment)

 integer(pInt)        computationMode, i

 if (inc == 0) then
   cycleCounter = 4
 else
   if (theCycle > ncycle .or. theInc /= inc) cycleCounter = 0                  ! reset counter for each cutback or new inc
   if (theCycle /= ncycle .or. theLovl /= lovl) then
     cycleCounter = cycleCounter+1   ! ping pong
     outdatedFFN1 = .false.
     write (6,*) n(1),nn,'cycleCounter',cycleCounter
     call debug_info()                          ! output of debugging/performance statistics of former
     call debug_reset()
   endif
 endif
 if (cptim > theTime .or. theInc /= inc) then                                   ! reached convergence
   lastIncConverged = .true.
   outdatedByNewInc = .true.
   write (6,*) n(1),nn,'lastIncConverged + outdated'
 endif

 if (mod(cycleCounter,2) /= 0) computationMode = 4   ! recycle in odd cycles
 if (mod(cycleCounter,4) == 2) computationMode = 3   ! collect in 2,6,10,...
 if (mod(cycleCounter,4) == 0) computationMode = 2   ! compute in 0,4,8,...
 if (computationMode == 4 .and. ncycle == 0 .and. .not. lastIncConverged) &
   computationMode = 6    ! recycle but restore known good consistent tangent
 if (computationMode == 4 .and. lastIncConverged) then
   computationMode  = 5   ! recycle and record former consistent tangent
   lastIncConverged = .false.
 endif
 if (computationMode == 2 .and. outdatedByNewInc) then
   computationMode  = 1   ! compute and age former results
   outdatedByNewInc = .false.
 endif

 theTime  = cptim                                   ! record current starting time
 theInc   = inc                                     ! record current increment number
 theCycle = ncycle                                  ! record current cycle count
 theLovl  = lovl                                    ! record current lovl

 call CPFEM_general(computationMode,ffn,ffn1,t(1),timinc,n(1),nn,s,d,ngens)

!     Mandel: 11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     Marc:   11, 22, 33, 12, 23, 13
 forall(i=1:ngens) d(1:ngens,i) = invnrmMandel(i)*d(1:ngens,i)*invnrmMandel(1:ngens)
 s(1:ngens) = s(1:ngens)*invnrmMandel(1:ngens)
 if(symmetricSolver) d(1:ngens,1:ngens) = 0.5_pReal*(d(1:ngens,1:ngens)+transpose(d(1:ngens,1:ngens)))

 return
 
end subroutine


!********************************************************************
!     This routine sets user defined output variables for Marc
!********************************************************************
!
!     select a variable contour plotting (user subroutine).
!
!********************************************************************
subroutine plotv(&
     v,&          ! variable
     s,&          ! stress array
     sp,&         ! stresses in preferred direction
     etot,&       ! total strain (generalized)
     eplas,&      ! total plastic strain
     ecreep,&     ! total creep strain
     t,&          ! current temperature
     m,&          ! element number
     nn,&         ! integration point number
     layer,&      ! layer number
     ndi,&        ! number of direct stress components
     nshear,&     ! number of shear stress components
     jpltcd &     ! user variable index
    )
 use prec,  only: pReal,pInt
 use mesh,  only: mesh_FEasCP
 use homogenization, only: materialpoint_results
 implicit none

 real(pReal) s(*),etot(*),eplas(*),ecreep(*),sp(*)
 real(pReal) v, t(*)
 integer(pInt) m, nn, layer, ndi, nshear, jpltcd

 v = materialpoint_results(jpltcd,nn,mesh_FEasCP('elem', m))
 return

end subroutine



! subroutine utimestep(timestep,timestepold,icall,time,timeloadcase)
!********************************************************************
!     This routine modifies the addaptive time step of Marc
!********************************************************************
! use prec, only: pReal,pInt
! use CPFEM, only : CPFEM_timefactor_max
! implicit none
!
! real(pReal) timestep, timestepold, time,timeloadcase 
! integer(pInt) icall
!
! user subroutine for modifying the time step in auto step
!
!   timestep    :  the current time step as suggested by marc
!                  to be modified in this routine
!   timestepold :  the current time step before it was modified by marc
!   icall       :  =1 for setting the initial time step
!                  =2 if this routine is called during an increment
!                  =3 if this routine is called at the beginning
!                     of the increment
!   time        :  time at the start of the current increment
!   timeloadcase:  time period of the current load case
!
!   it is in general not recommended to increase the time step
!   during the increment.
!   this routine is called right after the time step has (possibly)
!   been updated by marc.
!
!  user coding
!     reduce timestep during increment in case mpie_timefactor is too large
! if(icall==2_pInt) then
!    if(mpie_timefactor_max>1.25_pReal) then
!        timestep=min(timestep,timestepold*0.8_pReal)
!    end if
! return
!     modify timestep at beginning of new increment
! else if(icall==3_pInt) then
!    if(mpie_timefactor_max<=0.8_pReal) then
!        timestep=min(timestep,timestepold*1.25_pReal)
!    else if (mpie_timefactor_max<=1.0_pReal) then
!        timestep=min(timestep,timestepold/mpie_timefactor_max)
!    else if (mpie_timefactor_max<=1.25_pReal) then
!        timestep=min(timestep,timestepold*1.01_pReal)
!    else
!        timestep=min(timestep,timestepold*0.8_pReal)
!    end if
! end if
! return
! end