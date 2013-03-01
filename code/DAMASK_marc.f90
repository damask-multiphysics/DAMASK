! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced Material Simulation Kit.
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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luc Hantcherli, Max-Planck-Institut für Eisenforschung GmbH
!> @author W.A. Counts
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Material subroutine for MSC.Marc
!> @details     Usage:
!> @details             - choose material as hypela2
!> @details             - set statevariable 2 to index of homogenization
!> @details             - set statevariable 3 to index of microstructure
!> @details             - make sure the file "material.config" exists in the working
!> @details               directory
!> @details             - make sure the file "numerics.config" exists in the working 
!> @details               directory
!> @details             - use nonsymmetric option for solver (e.g. direct 
!> @details               profile or multifrontal sparse, the latter seems
!> @details               to be faster!)
!> @details             - in case of ddm (domain decomposition)a SYMMETRIC
!> @details               solver has to be used, i.e uncheck "non-symmetric"
!> @details     Marc subroutines used:
!> @details             - hypela2
!> @details             - plotv
!> @details             - quit
!> @details     Marc common blocks included:
!> @details             - concom: lovl, ncycle, inc, incsub
!> @details             - creeps: timinc
!--------------------------------------------------------------------------------------------------

#ifndef INT
#define INT 4
#endif

#ifndef FLOAT
#define FLOAT 8
#endif

#define Marc

#include "prec.f90"

module DAMASK_interface
 
 implicit none
 character(len=4), parameter :: InputFileExtension = '.dat'
 character(len=4), parameter :: LogFileExtension = '.log'

contains


!--------------------------------------------------------------------------------------------------
!> @brief only output of current version
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init

implicit none

!$OMP CRITICAL (write2out)
 write(6,'(/,a)') ' <<<+-  DAMASK_marc init  -+>>>'
 write(6,'(a)')   ' $Id$'
#include "compilation_info.f90"
!$OMP END CRITICAL (write2out)

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the current workingDir 
!--------------------------------------------------------------------------------------------------
function getSolverWorkingDirectoryName()

 implicit none
 character(1024) getSolverWorkingDirectoryName, inputName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forward and backward slash

 getSolverWorkingDirectoryName=''
 inputName=''
 inquire(5, name=inputName) ! determine inputputfile
 getSolverWorkingDirectoryName=inputName(1:scan(inputName,pathSep,back=.true.))
! write(6,*) 'getSolverWorkingDirectoryName', getSolverWorkingDirectoryName

end function getSolverWorkingDirectoryName


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
function getSolverJobName()
 use prec, only: &
   pReal, &
   pInt

 implicit none
 character(1024) :: getSolverJobName, inputName
 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! forward and backward slash
 integer(pInt) :: extPos

 getSolverJobName=''
 inputName=''
 inquire(5, name=inputName) ! determine inputfile
 extPos = len_trim(inputName)-4
 getSolverJobName=inputName(scan(inputName,pathSep,back=.true.)+1:extPos)

end function getSolverJobName


end module DAMASK_interface

#include "IO.f90"
#include "numerics.f90"
#include "debug.f90"
#include "math.f90"
#include "FEsolving.f90"
#include "mesh.f90"
#include "material.f90"
#include "lattice.f90"
#include "constitutive_none.f90"
#include "constitutive_j2.f90"
#include "constitutive_phenopowerlaw.f90"
#include "constitutive_titanmod.f90"
#include "constitutive_dislotwin.f90"
#include "constitutive_nonlocal.f90"
#include "constitutive.f90"
#include "crystallite.f90"
#include "homogenization_isostrain.f90"
#include "homogenization_RGC.f90"
#include "homogenization.f90"
#include "CPFEM.f90"


!--------------------------------------------------------------------------------------------------
!> @brief This is the MSC.Marc user subroutine for defining material behavior
!> @details CAUTION : Due to calculation of the Deformation gradients, Stretch Tensors and
!> @details         Rotation tensors at previous and current states, the analysis can be
!> @details         computationally expensive. Please use the user subroutine ->  hypela
!> @details         if these kinematic quantities are not needed in the constitutive model
!> @details
!> @details IMPORTANT NOTES :
!> @details
!> @details (1) F,R,U are only available for continuum and membrane elements (not for
!> @details     shells and beams).
!> @details
!> @details (2) For total Lagrangian formulation use the -> 'Elasticity,1' card(=
!> @details     total Lagrange with large disp) in the parameter section of input deck.
!> @details     For updated Lagrangian formulation use the -> 'Plasticity,3' card(=
!> @details     update+finite+large disp+constant d) in the parameter section of
!> @details     input deck.
!> @details
!> @details     The following operation obtains U (stretch tensor) at t=n+1 :
!> @details
!> @details     call scla(un1,0.d0,itel,itel,1)
!> @details     do k=1,3
!> @details      do i=1,3
!> @details       do  j=1,3
!> @details        un1(i,j)=un1(i,j)+dsqrt(strechn1(k))*eigvn1(i,k)*eigvn1(j,k)
!> @details       enddo
!> @details      enddo
!> @details     enddo
!--------------------------------------------------------------------------------------------------
subroutine hypela2(&
     d,&          !< stress strain law to be formed
     g,&          !< change in stress due to temperature effects
     e,&          !< total elastic strain
     de,&         !< increment of strain
     s,&          !< stress - should be updated by user
     t,&          !< state variables (comes in at t=n, must be updated to have state variables at t=n+1)
     dt,&         !< increment of state variables
     ngens,&      !< size of stress - strain law
     n,&          !< element number
     nn,&         !< integration point number
     kcus,&       !< (1) layer number, (2) internal layer number
     matus,&      !< (1) user material identification number, (2) internal material identification number
     ndi,&        !< number of direct components
     nshear,&     !< number of shear components
     disp,&       !< incremental displacements
     dispt,&      !< displacements at t=n (at assembly, lovl=4) and displacements at t=n+1 (at stress recovery, lovl=6)
     coord,&      !< coordinates
     ffn,&        !< deformation gradient
     frotn,&      !< rotation tensor
     strechn,&    !< square of principal stretch ratios, lambda(i)
     eigvn,&      !< i principal direction components for j eigenvalues
     ffn1,&       !< deformation gradient
     frotn1,&     !< rotation tensor
     strechn1,&   !< square of principal stretch ratios, lambda(i)
     eigvn1,&     !< i principal direction components for j eigenvalues
     ncrd,&       !< number of coordinates
     itel,&       !< dimension of F and R, either 2 or 3
     ndeg,&       !< number of degrees of freedom  ==> is this at correct list position ?!?
     ndm,&        !< 
     nnode,&      !< number of nodes per element
     jtype,&      !< element type
     lclass,&     !< element class
     ifr,&        !< set to 1 if R has been calculated
     ifu &        !< set to 1 if stretch has been calculated
   )

 use prec, only:      pReal, &
                      pInt
 use numerics, only:  numerics_unitlength
 use FEsolving, only: cycleCounter, &
                      theInc, &
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
 use mesh, only:      mesh_FEasCP, &
                      mesh_element, &
                      mesh_node0, &
                      mesh_node, &
                      mesh_build_subNodeCoords, &
                      mesh_build_ipCoordinates, &
                      FE_Nnodes, &
                      FE_geomtype
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_init_done, &
   CPFEM_initAll, &
   CPFEM_CALCRESULTS, &
   CPFEM_AGERESULTS, &
   CPFEM_COLLECT, &
   CPFEM_RESTOREJACOBIAN, &
   CPFEM_BACKUPJACOBIAN
   
!$ use numerics, only: DAMASK_NumThreadsInt                                   ! number of threads set by DAMASK_NUM_THREADS
 
 implicit none
!$ include "omp_lib.h"                                                          ! the openMP function library
!     ** Start of generated type statements **
 real(pReal) coord, d, de, disp, dispt, dt, e, eigvn, eigvn1, ffn, ffn1
 real(pReal) frotn, frotn1, g
 integer(pInt) ifr, ifu, itel, jtype, kcus, lclass, matus, n, ncrd, ndeg
 integer(pInt) ndi, ndm, ngens, nn, nnode, nshear
 real(pReal) s, strechn, strechn1, t
 logical :: cutBack
 !     ** End of generated type statements **

 dimension e(*),de(*),t(*),dt(*),g(*),d(ngens,*),s(*), n(2),coord(ncrd,*),disp(ndeg,*),matus(2),dispt(ndeg,*),ffn(itel,*),&
           frotn(itel,*),strechn(itel),eigvn(itel,*),ffn1(itel,*),frotn1(itel,*),strechn1(itel),eigvn1(itel,*),kcus(2), lclass(2)

! Marc common blocks are in fixed format so they have to be reformated to free format (f90)
! Beware of changes in newer Marc versions

 include "include/concom%%MARCVERSION%%"     ! concom is needed for inc, subinc, ncycle, lovl
 include "include/creeps%%MARCVERSION%%"     ! creeps is needed for timinc (time increment)

 real(pReal), dimension(6) ::   stress
 real(pReal), dimension(6,6) :: ddsdde
 
 real(pReal), dimension (3,3) :: pstress                                  ! dummy argument for call of cpfem_general (used by mpie_spectral)
 real(pReal), dimension (3,3,3,3) :: dPdF                                 ! dummy argument for call of cpfem_general (used by mpie_spectral)

 integer(pInt) computationMode, i, cp_en
 integer(pInt) node, FEnodeID

! OpenMP variable
!$ integer(pInt) defaultNumThreadsInt                                     ! default value set by Marc
 
 
!$ defaultNumThreadsInt = omp_get_num_threads()                           ! remember number of threads set by Marc

 if (.not. CPFEM_init_done) call CPFEM_initAll(t(1),n(1),nn)

!$ call omp_set_num_threads(DAMASK_NumThreadsInt)                         ! set number of threads for parallel execution set by DAMASK_NUM_THREADS

 if (lovl == 4 ) then
   if(timinc < theDelta .and. theInc == inc ) &                           ! first after cutback
     computationMode = CPFEM_RESTOREJACOBIAN
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
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> start of analysis..! ',n(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     else if (inc - theInc > 1) then                                      ! >> restart of broken analysis <<
       lastIncConverged = .false.                                         ! no Jacobian backup
       outdatedByNewInc = .false.                                         ! no aging of state
       lastMode = .true.                                                  ! pretend last step was calculation
       calcMode = .true.                                                  ! pretend last step was calculation
       !$OMP CRITICAL (write2out)
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> restart of analysis..! ',n(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     else                                                                 ! >> just the next inc <<
       lastIncConverged = .true.                                          ! request Jacobian backup
       outdatedByNewInc = .true.                                          ! request aging of state
       lastMode = .true.                                                  ! assure last step was calculation
       calcMode = .true.                                                  ! assure last step was calculation
       !$OMP CRITICAL (write2out)
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> new increment..! ',n(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     endif
   else if ( timinc < theDelta ) then                                     ! >> cutBack <<
     terminallyIll = .false.
     cycleCounter = -1                                                    ! first calc step increments this to cycle = 0
     calcMode = .true.                                                    ! pretend last step was calculation
     !$OMP CRITICAL (write2out)
       write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> cutback detected..! ',n(1),nn
       flush(6)
     !$OMP END CRITICAL (write2out)
   endif                                                                  ! convergence treatment end

   calcMode(nn,cp_en) = .not. calcMode(nn,cp_en)                          ! ping pong (calc <--> collect)

   if ( calcMode(nn,cp_en) ) then                                         ! now --- CALC ---
     if ( lastMode /= calcMode(nn,cp_en) ) then                           ! first after ping pong
       call debug_reset()                                                 ! resets debugging
       outdatedFFN1  = .false.
       cycleCounter  = cycleCounter + 1_pInt
       call mesh_build_subNodeCoords()                                    ! update subnodal coordinates
       call mesh_build_ipCoordinates()                                    ! update ip coordinates
     endif
     if ( outdatedByNewInc ) then
       computationMode = ior(CPFEM_CALCRESULTS,CPFEM_AGERESULTS)
       outdatedByNewInc = .false.                                         ! reset flag
     else
       computationMode = CPFEM_CALCRESULTS
     endif
   else                                                                   ! now --- COLLECT ---
     if ( lastMode /= calcMode(nn,cp_en) .and. &
          .not. terminallyIll ) then
       call debug_info()                                                  ! first after ping pong reports (meaningful) debugging
     endif
     if ( lastIncConverged ) then
       computationMode = ior(CPFEM_COLLECT,CPFEM_BACKUPJACOBIAN)          ! collect and backup Jacobian after convergence
       lastIncConverged = .false.                                         ! reset flag
     else
       computationMode = CPFEM_COLLECT                                    ! plain collect
     endif
     do node = 1,FE_Nnodes(FE_geomtype(mesh_element(2,cp_en)))
       FEnodeID = mesh_FEasCP('node',mesh_element(4+node,cp_en))
       mesh_node(1:3,FEnodeID) = mesh_node0(1:3,FEnodeID) + numerics_unitlength * dispt(1:3,node)
     enddo
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
 
!$ call omp_set_num_threads(defaultNumThreadsInt)                               ! reset number of threads to stored default value

end subroutine hypela2


!--------------------------------------------------------------------------------------------------
!> @brief sets user defined output variables for Marc
!> @details select a variable contour plotting (user subroutine).
!--------------------------------------------------------------------------------------------------
subroutine plotv(&
     v,&          !< variable
     s,&          !< stress array
     sp,&         !< stresses in preferred direction
     etot,&       !< total strain (generalized)
     eplas,&      !< total plastic strain
     ecreep,&     !< total creep strain
     t,&          !< current temperature
     m,&          !< element number
     nn,&         !< integration point number
     layer,&      !< layer number
     ndi,&        !< number of direct stress components
     nshear,&     !< number of shear stress components
     jpltcd &     !< user variable index
    )
 use prec,  only: pReal,pInt
 use mesh,  only: mesh_FEasCP
 use IO,    only: IO_error
 use homogenization, only: materialpoint_results,materialpoint_sizeResults
 
 implicit none
 real(pReal) s(*),etot(*),eplas(*),ecreep(*),sp(*)
 real(pReal) v, t(*)
 integer(pInt) m, nn, layer, ndi, nshear, jpltcd

 if (jpltcd > materialpoint_sizeResults) call IO_error(700_pInt,jpltcd)                             ! complain about out of bounds error

 v = materialpoint_results(jpltcd,nn,mesh_FEasCP('elem', m))

end subroutine plotv
