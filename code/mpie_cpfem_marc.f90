!* $Id$
!********************************************************************
! Material subroutine for MSC.Marc
!
! written by P. Eisenlohr,
!            F. Roters,
!            L. Hantcherli,
!            W.A. Counts
!            D.D. Tjahjanto
!            C. Kords
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


MODULE mpie_interface

character(len=64), parameter :: FEsolver = 'Marc'
character(len=4),  parameter :: InputFileExtension = '.dat'
character(len=4),  parameter :: LogFileExtension = '.log'

CONTAINS

subroutine mpie_interface_init()


!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  mpie_cpfem_marc init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
!$OMP END CRITICAL (write2out)
 return
end subroutine

function getSolverWorkingDirectoryName()
 implicit none
 character(1024) getSolverWorkingDirectoryName, outName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forward and backward slash

 getSolverWorkingDirectoryName=''
 outName=''
 inquire(6, name=outName) ! determine outputfile
 getSolverWorkingDirectoryName=outName(1:scan(outName,pathSep,back=.true.))
! write(6,*) 'getSolverWorkingDirectoryName', getSolverWorkingDirectoryName
end function

function getSolverJobName()
 use prec
 implicit none

 character(1024) getSolverJobName, outName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forward and backward slash
 integer(pInt) extPos

 getSolverJobName=''
 outName=''
 inquire(6, name=outName) ! determine outputfile
 extPos = len_trim(outName)-4
 getSolverJobName=outName(scan(outName,pathSep,back=.true.)+1:extPos)
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
 include "constitutive_j2.f90"            ! uses prec, math, IO, lattice, material, debug
 include "constitutive_phenopowerlaw.f90" ! uses prec, math, IO, lattice, material, debug
 include "constitutive_titanmod.f90"      ! uses prec, math, IO, lattice, material, debug
 include "constitutive_dislotwin.f90"     ! uses prec, math, IO, lattice, material, debug
 include "constitutive_nonlocal.f90"      ! uses prec, math, IO, lattice, material, debug
 include "constitutive.f90"     ! uses prec, IO, math, lattice, mesh, debug
 include "crystallite.f90"      ! uses prec, math, IO, numerics, Fesolving, material, mesh, constitutive
 include "homogenization_isostrain.f90"   ! uses prec, math, IO
 include "homogenization_RGC.f90"         ! uses prec, math, IO, numerics, mesh: added <<<updated 31.07.2009>>>
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
 use CPFEM, only:     CPFEM_initAll,CPFEM_general,CPFEM_init_done
 use OMP_LIB
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

 real(pReal), dimension(6) ::   stress
 real(pReal), dimension(6,6) :: ddsdde
 
 real(pReal), dimension (3,3) :: pstress                                  ! dummy argument for call of cpfem_general (used by mpie_spectral)
 real(pReal), dimension (3,3,3,3) :: dPdF                                 ! dummy argument for call of cpfem_general (used by mpie_spectral)

 integer(pInt) computationMode, i, cp_en, &
               defaultNumThreadsInt, &                                    ! default value stored in environment variable OMP_NUM_THREADS
               mpieNumThreadsInt                                          ! value stored in environment variable MPIE_NUM_THREADS

 character mpieNumThreadsString
 
 
 defaultNumThreadsInt = omp_get_num_threads()                             ! remember default number of available threads defined by environment variable OMP_NUM_THREADS
 call GetEnv('MPIE_NUM_THREADS',mpieNumThreadsString)                     ! get environment variable MPIE_NUM_THREADS...
 read(mpieNumThreadsString,'(i2)') mpieNumThreadsInt                      ! ...convert it to integer...
 if (mpieNumThreadsInt < 1) mpieNumThreadsInt = 1                         ! ...ensure that its at least one...
 call omp_set_num_threads(mpieNumThreadsInt)                              ! ...and use it as number of threads for parallel execution
 
 
 if (.not. CPFEM_init_done) call CPFEM_initAll(t(1),n(1),nn)

 if (lovl == 4) then                                                      ! Marc requires stiffness in separate call
   if ( timinc < theDelta .and. theInc == inc ) then                      ! first after cutback
     computationMode = 7                                                  !  --> restore tangent and return it
   else
     computationMode = 6                                                  !  --> just return known tangent
   endif
 else                                                                     ! stress requested (lovl == 6)
   cp_en = mesh_FEasCP('elem',n(1))

   if (cptim > theTime .or. inc /= theInc) then                           ! reached "convergence"

     terminallyIll = .false.
     cycleCounter = -1                                                    ! first calc step increments this to cycle = 0
     if (inc == 0) then                                                   ! >> start of analysis <<
       lastIncConverged = .false.                                         ! no Jacobian backup
       outdatedByNewInc = .false.                                         ! no aging of state
       lastMode = .false.                                                 ! pretend last step was collection
       calcMode = .false.                                                 ! pretend last step was collection
       !$OMP CRITICAL (write2out)
       write (6,'(i6,x,i2,x,a)') n(1),nn,'<< hypela2 >> start of analysis..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     else if (inc - theInc > 1) then                                      ! >> restart of broken analysis <<
       lastIncConverged = .false.                                         ! no Jacobian backup
       outdatedByNewInc = .false.                                         ! no aging of state
       lastMode = .true.                                                  ! pretend last step was calculation
       calcMode = .true.                                                  ! pretend last step was calculation
       !$OMP CRITICAL (write2out)
       write (6,'(i6,x,i2,x,a)') n(1),nn,'<< hypela2 >> restart of analysis..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     else                                                                 ! >> just the next inc <<
       lastIncConverged = .true.                                          ! request Jacobian backup
       outdatedByNewInc = .true.                                          ! request aging of state
       lastMode = .true.                                                  ! assure last step was calculation
       calcMode = .true.                                                  ! assure last step was calculation
       !$OMP CRITICAL (write2out)
       write (6,'(i6,x,i2,x,a)') n(1),nn,'<< hypela2 >> new increment..!'; call flush(6)
       !$OMP END CRITICAL (write2out)
     endif

   else if ( timinc < theDelta ) then                                     ! >> cutBack <<

     terminallyIll = .false.
     cycleCounter = -1                                                    ! first calc step increments this to cycle = 0
     calcMode = .true.                                                    ! pretend last step was calculation
     !$OMP CRITICAL (write2out)
     write(6,'(i6,x,i2,x,a)') n(1),nn,'<< hypela2 >> cutback detected..!'; call flush(6)
     !$OMP END CRITICAL (write2out)

   endif                                                                  ! convergence treatment end

   calcMode(nn,cp_en) = .not. calcMode(nn,cp_en)                          ! ping pong (calc <--> collect)

   if ( calcMode(nn,cp_en) ) then                                         ! now --- CALC ---
     if ( lastMode /= calcMode(nn,cp_en) ) then                           ! first after ping pong
       call debug_reset()                                                 ! resets debugging
       outdatedFFN1  = .false.
       cycleCounter  = cycleCounter + 1_pInt
     endif
     if ( outdatedByNewInc ) then
       outdatedByNewInc = .false.                                         ! reset flag
       computationMode = 1                                                ! calc and age results
     else
       computationMode = 2                                                ! plain calc
     endif
   else                                                                   ! now --- COLLECT ---
     if ( lastMode /= calcMode(nn,cp_en) .and. &
          .not. terminallyIll ) call debug_info()                         ! first after ping pong reports (meaningful) debugging
     if ( lastIncConverged ) then
       lastIncConverged = .false.                                         ! reset flag
       computationMode = 4                                                ! collect and backup Jacobian after convergence
     else
       computationMode = 3                                                ! plain collect
     endif
   endif

   theTime  = cptim                                                       ! record current starting time
   theDelta = timinc                                                      ! record current time increment
   theInc   = inc                                                         ! record current increment number
   lastMode = calcMode(nn,cp_en)                                          ! record calculationMode
 endif

 call CPFEM_general(computationMode,ffn,ffn1,t(1),timinc,n(1),nn,stress,ddsdde, pstress, dPdF)

!     Mandel: 11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     Marc:   11, 22, 33, 12, 23, 13
!     Marc:   11, 22, 33, 12

 forall(i=1:ngens) d(1:ngens,i) = invnrmMandel(i)*ddsdde(1:ngens,i)*invnrmMandel(1:ngens)
 s(1:ngens) = stress(1:ngens)*invnrmMandel(1:ngens)
 if(symmetricSolver) d(1:ngens,1:ngens) = 0.5_pReal*(d(1:ngens,1:ngens)+transpose(d(1:ngens,1:ngens)))

 call omp_set_num_threads(defaultNumThreadsInt)                               ! reset number of threads to stored default value

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
